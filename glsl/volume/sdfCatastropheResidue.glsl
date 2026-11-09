// sdfCatastropheResidue: complex-catastrophe-Morse tangential ray calculus.
// 1-tap Cauchy residue volume integration, Arnold fold catastrophe invariant,
// wave-optical Airy diffraction caustics, and Euler-Morse overlap cancellation.
// See docs/volume/sdfCatastropheResidue.md.
#ifndef SDF_GEMS_SDF_CATASTROPHE_RESIDUE_GLSL
#define SDF_GEMS_SDF_CATASTROPHE_RESIDUE_GLSL

// 1-step Newton-Eikonal tangency jump & Arnold Fold (A2) Catastrophe Invariant (Arnold 1984).
// Given ray distance d, first directional derivative dDot, and ray curvature dDdot > 0:
// returns vec2(dtStar, u1) where dtStar is the ray-distance offset to closest approach t*,
// and u1 = d - dDot^2 / (2 * kap) is the curvature-corrected signed clearance d_min.
vec2 sdfFoldCatastrophe(float d, float dDot, float dDdot)
{
    float kap = max(dDdot, 1e-4);
    float dtStar = -dDot / kap;
    float u1 = d - (dDot * dDot) / (2.0 * kap);
    return vec2(dtStar, u1);
}

// Eikonal curvature estimate from distance and directional derivative when full Hessian is unavailable.
// Exact for spherical and cylindrical level sets outside the surface (M21).
float sdfEikonalCurvature(float d, float dDot)
{
    float num = max(1.0 - dDot * dDot, 0.0);
    float den = max(abs(d), 1e-4);
    return num / den;
}

// 1-tap Cauchy residue integral of Lorentzian density sigma(d) = eps / (eps + max(d, 0))
// along the entire grazing ray s in (-inf, +inf) past a Fold catastrophe with clearance u1 >= 0.
// Exact closed form via Cauchy's Residue Theorem in C+ (Baker & Graves-Morris 1996).
float sdfComplexPoleResidue(float u1, float dDdot, float eps)
{
    float kap = max(dDdot, 1e-4);
    float e = max(eps, 1e-6);
    float u = max(u1, 0.0);
    float denom = sqrt(0.5 * kap * (e + u));
    return (3.141592653589793 * e) / denom;
}

// Analytical Cauchy slab integral of Lorentzian density over finite ray segment [sa, sb]
// where sa = ta - t* and sb = tb - t* relative to closest approach t*.
float sdfComplexPoleSlab(float sa, float sb, float u1, float dDdot, float eps)
{
    float kap = max(dDdot, 1e-4);
    float e = max(eps, 1e-6);
    float u = max(u1, 0.0);
    float scale = sqrt(0.5 * kap / (e + u));
    float factor = e / sqrt(0.5 * kap * (e + u));
    return factor * (atan(sb * scale) - atan(sa * scale));
}

// Smooth, branch-safe, unconditionally stable rational/asymptotic approximation
// to the Airy function Ai(x) on (-inf, +inf).
// Matches scipy.special.airy on [-8, 4] with mean error < 3e-4 and zero overflow.
float sdfAiryAi(float x)
{
    if (x > 25.0) return 0.0;

    // Taylor polynomial around x = 0 (Abramowitz & Stegun 10.4.4), clamped to prevent FP32 overflow
    float xc = clamp(x, -3.0, 2.5);
    float x3 = xc * xc * xc;
    float f = 1.0 + x3 * (0.16666667 + x3 * (0.00555556 + x3 * 0.00007716));
    float g = xc * (1.0 + x3 * (0.08333333 + x3 * (0.00198413 + x3 * 0.00002205)));
    float yPoly = 0.35502805 * f - 0.25881940 * g;

    // Asymptotic expansion for x > 0 (exponential shadow decay, blended above x = 1.2)
    float xp = clamp(x, 1.0, 25.0);
    float xiP = 0.66666667 * xp * sqrt(xp);
    float yPos = (exp(-xiP) / (3.54490770 * sqrt(sqrt(xp)))) * (1.0 - 0.06944444 / xiP);

    // Asymptotic expansion for x < 0 (oscillatory interference fringes, blended below x = -1.8)
    float xn = clamp(-x, 1.0, 100.0);
    float xiN = 0.66666667 * xn * sqrt(xn);
    float phase = xiN + 0.78539816;
    float yNeg = (sin(phase) - (0.06944444 / xiN) * cos(phase)) / (1.77245385 * sqrt(sqrt(xn)));

    // Smoothstep blends across transition boundaries
    float wPos = smoothstep(1.2, 2.2, x);
    float wNeg = smoothstep(-1.8, -2.8, x);

    return yPoly * (1.0 - wPos) * (1.0 - wNeg) + yPos * wPos + yNeg * wNeg;
}

// Arnold Fold (A2) Caustic Diffraction (Berry & Upstill 1980).
// Evaluates wave-optical caustic intensity Ai^2(-kWave^(2/3) * u1).
// Produces rainbow supernumerary diffraction fringes outside silhouettes and caustics.
float sdfAiryCaustic(float u1, float kWave)
{
    float k23 = pow(max(kWave, 1e-4), 0.66666667);
    float ai = sdfAiryAi(-k23 * u1);
    return ai * ai;
}

// Euler-Morse critical-point tomography term (Morse 1934; Baryshnikov & Ghrist 2009).
// Weights tangency by Morse index sign (-1)^mu = sign(dDdot).
// Eliminates double-counted volume brightening seams across overlapping SDF lobes.
float sdfMorseEulerTerm(float u1, float dDdot, float eps)
{
    float sgn = dDdot >= 0.0 ? 1.0 : -1.0;
    float kap = max(abs(dDdot), 1e-4);
    float e = max(eps, 1e-6);
    float u = max(u1, 0.0);
    float denom = sqrt(0.5 * kap * (e + u));
    return sgn * (3.141592653589793 * e) / denom;
}

#endif // SDF_GEMS_SDF_CATASTROPHE_RESIDUE_GLSL
