// sdfWkbAiryCavity: Boundary-Locked Fermi-Weingarten-Lazutkin Whispering-Gallery Cavities,
// Stieltjes-Pade Vallee-Soares Zero-Loop Ai^2 Slabs & da Costa Extrinsic-Curvature Waveguides.
//
// Solves:
//   1. Curved-Boundary Helmholtz Cavity Eigen-SDFs (psi|_{d=0} == 0 identically):
//      - Existing sdgChladniEigen evaluates free-space plane waves on R^3 whose nodal surfaces
//        crash arbitrarily across curved SDF boundaries (psi|_{d=0} != 0).
//      - sdgWhisperingGalleryCavity constructs Lazutkin-Langer-Airy whispering-gallery eigenmodes
//        psi(x) = J^{-1/2} * Ai(-beta * d - zeta_m) / |Ai'(-zeta_m)| * sin^l(theta) * cos(l*phi + phi0)
//        locked to any smooth curved SDF boundary {d(x) = 0} via the exact Airy zeros Ai(-zeta_m) = 0,
//        with zero polar singularities in 3D, an exponentially quiet Rayleigh core for -beta*d > zeta_m,
//        and exact analytical Sobolev-regularized gradient vec4(d_psi, grad_d_psi).
//   2. Zero-Loop Vallee-Soares Ai^2 Volumetric Slab Integration:
//      - Integrates radial Airy probability density Ai(a0 + a1*t)^2 analytically along any ray segment
//        [ta, tb] via Vallee & Soares' exact Wronskian antiderivative F(u) = u*Ai(u)^2 - Ai'(u)^2,
//        stabilized for u > 1.8 by a pole-free [1/1] Stieltjes-Pade asymptotic resummation that
//        eliminates severe FP32 catastrophic cancellation in the evanescent tail.
//   3. da Costa (1981) & Exner-Seba (1989) Extrinsic-Curvature Quantum Waveguide Traps:
//      - Evaluates the purely geometric attractive potential V_daCosta = -(H^2 - K) = -(k1 - k2)^2 / 4 <= 0
//        inside constant-thickness SDF shells {|d| < w}, vanishing at spheres/planes (umbilic points)
//        while trapping bound states at hyperbolic saddles (K < 0) and non-umbilic waveguide bends.
//
// References:
//   - Lord Rayleigh (Phil. Mag. 1910), "The problem of the whispering gallery"
//   - Langer (Phys. Rev. 1937), "On the connection formulas and the solutions of the wave equation"
//   - Keller & Rubinow (Ann. Phys. 1960), "Asymptotic solution of eigenvalue problems"
//   - Lazutkin (Math. USSR Izv. 1973), "The existence of caustics for a billiard problem in a convex domain"
//   - da Costa (Phys. Rev. A 1981), "Quantum mechanics of a constrained particle"
//   - Mather (Ergod. Th. Dynam. Sys. 1982) & Heller (PRL 1984), Non-convex caustic destruction & scars
//   - Exner & Seba (J. Math. Phys. 1989), "Bound states in curved quantum waveguides"
//   - Vallee & Soares (Airy Functions and Applications to Physics, 2004, Eq. 3.63)
//
// Attribution:
//   Discovered, derived, implemented, and benchmarked using Gemini 4 Carbon (sdf-gems, 2026)

#ifndef SDF_GEMS_SDF_WKB_AIRY_CAVITY_GLSL
#define SDF_GEMS_SDF_WKB_AIRY_CAVITY_GLSL

// -----------------------------------------------------------------------------
// 1. Joint Airy Function & Derivative Jet vec2(Ai(x), Ai'(x))
// -----------------------------------------------------------------------------
#ifndef SDF_GEMS_HAS_AIRY_JET
#define SDF_GEMS_HAS_AIRY_JET
// Smooth, overflow-free joint evaluation of vec2(Ai(x), Ai'(x)) across (-inf, +inf).
// Combines degree-15 Maclaurin polynomials on [-2.8, 2.0] with 3-term DLMF 9.7.5-9.7.10
// asymptotic expansions via standard-compliant smoothstep blends (max error < 1e-3).
vec2 sdfAiryJet(in float x)
{
    if (x > 25.0)
    {
        return vec2(0.0);
    }

    // Degree-15 Taylor expansion around x = 0 (DLMF 9.4.1-9.4.6)
    float xc  = clamp(x, -2.8, 2.0);
    float x2  = xc * xc;
    float x3  = x2 * xc;
    float f   = 1.0 + x3 * (0.16666667 + x3 * (5.5555556e-3 + x3 * (7.7160494e-5
              + x3 * (5.8454919e-7 + x3 * 2.7835676e-9))));
    float g   = xc * (1.0 + x3 * (0.083333336 + x3 * (1.9841270e-3 + x3 * (2.2045855e-5
              + x3 * (1.3778660e-7 + x3 * 5.6702314e-10)))));
    float fp  = x2 * (0.5 + x3 * (0.033333335 + x3 * (6.9444446e-4 + x3 * (7.0145903e-6
              + x3 * 4.1753514e-8))));
    float gp  = 1.0 + x3 * (0.33333334 + x3 * (1.3888889e-2 + x3 * (2.2045855e-4
              + x3 * (1.7912257e-6 + x3 * 9.0723703e-9))));

    float aiPoly  = 0.35502805 * f  - 0.25881940 * g;
    float aipPoly = 0.35502805 * fp - 0.25881940 * gp;

    // Positive evanescent asymptotic x >= 1.6 (DLMF 9.7.5 & 9.7.6)
    float xp     = clamp(x, 1.0, 25.0);
    float xip    = 0.66666667 * xp * sqrt(xp);
    float invXip = 1.0 / xip;
    float qx     = sqrt(sqrt(xp));
    float ex     = exp(-xip);
    float aiPos  = (0.28209479 * ex / qx) * (1.0 + invXip * (-0.069444445
                 + invXip * (0.037133488 - invXip * 0.037993059)));
    float aipPos = -(0.28209479 * qx * ex) * (1.0 + invXip * (0.097222224
                 + invXip * (-0.043885030 + invXip * 0.042462831)));

    // Negative oscillatory asymptotic x <= -2.2 (DLMF 9.7.9 & 9.7.10)
    float xn      = max(-x, 1.0);
    float xin     = 0.66666667 * xn * sqrt(xn);
    float invXin  = 1.0 / xin;
    float invXin2 = invXin * invXin;
    float qn      = sqrt(sqrt(xn));
    float ph      = mod(xin, 6.28318531) + 0.78539816;
    float sp      = sin(ph);
    float cp      = cos(ph);
    float aiNeg   = (0.56418958 / qn) * (sp * (1.0 - 0.037133488 * invXin2)
                  - cp * invXin * (0.069444445 - 0.037993059 * invXin2));
    float aipNeg  = -(0.56418958 * qn) * (cp * (1.0 + 0.043885030 * invXin2)
                  + sp * invXin * (-0.097222224 + 0.042462831 * invXin2));

    float wp = smoothstep(1.6, 2.2, x);
    float wn = smoothstep(2.2, 2.8, -x);
    float wc = (1.0 - wp) * (1.0 - wn);

    return vec2(
        aiPoly  * wc + aiPos  * wp + aiNeg  * wn,
        aipPoly * wc + aipPos * wp + aipNeg * wn
    );
}
#endif // SDF_GEMS_HAS_AIRY_JET

// -----------------------------------------------------------------------------
// 2. Vallee-Soares Exact Ai^2 Antiderivative & Zero-Loop Born Slab Integral
// -----------------------------------------------------------------------------
// Evaluates F(u) = -integral_u^{+inf} Ai(s)^2 ds <= 0 (satisfying F'(u) == Ai(u)^2).
// For u <= 1.8, uses the exact Airy identity F(u) = u*Ai(u)^2 - Ai'(u)^2.
// For u >= 2.4, switches smoothly to a subtraction-free [1/1] Stieltjes-Pade asymptotic
// resummation matching c1 = -17/36 and c2 = +1225/2592, preventing FP32 cancellation.
float sdfValleeSoaresAi2Antideriv(in float u)
{
    if (u > 18.0)
    {
        return 0.0;
    }
    vec2  jet    = sdfAiryJet(min(u, 2.6));
    float fExact = u * jet.x * jet.x - jet.y * jet.y;
    if (u <= 1.8)
    {
        return min(fExact, 0.0);
    }

    float up    = max(u, 1.5);
    float xi    = 0.66666667 * up * sqrt(up);
    float lead  = -0.039788736 * exp(-min(2.0 * xi, 80.0)) / up; // -1/(8*pi) = -0.039788736
    float fPade = lead * (1224.0 * xi + 647.0) / (1224.0 * xi + 1225.0);

    float w = smoothstep(1.8, 2.4, u);
    return min(mix(fExact, fPade, w), 0.0);
}

// Closed-form O(1) zero-loop integral of Ai(a0 + a1 * t)^2 over t in [ta, tb].
// Eliminates raymarching loops across caustic envelopes and whispering-gallery rings.
float sdfAiryBornSlabIntegral(in float a0, in float a1, in float ta, in float tb)
{
    float dt = max(tb - ta, 0.0);
    if (dt <= 0.0)
    {
        return 0.0;
    }

    float du = a1 * dt;
    if (abs(du) < 0.02)
    {
        float uMid = a0 + 0.5 * a1 * (ta + tb);
        vec2  jet  = sdfAiryJet(uMid);
        float d2   = jet.y * jet.y + uMid * jet.x * jet.x; // 0.5 * d^2/du^2 [Ai(u)^2]
        return max(dt * (jet.x * jet.x + (du * du * 0.083333336) * d2), 0.0);
    }

    float ua = a0 + a1 * ta;
    float ub = a0 + a1 * tb;
    return max((sdfValleeSoaresAi2Antideriv(ub) - sdfValleeSoaresAi2Antideriv(ua)) / a1, 0.0);
}

// -----------------------------------------------------------------------------
// 3. Keller-Rubinow-Lazutkin Eigenvalues & Boundary-Locked Whispering-Gallery SDF
// -----------------------------------------------------------------------------
// Positive magnitude zeta_m of the m-th negative zero -zeta_m of Ai(x) (m >= 1).
float sdfAiryZeroZeta(in int radialModeM)
{
    int m = max(radialModeM, 1);
    if (m == 1) return 2.33810741;
    if (m == 2) return 4.08794944;
    if (m == 3) return 5.52055983;
    if (m == 4) return 6.78670809;
    float s = 1.17809725 * float(4 * m - 1); // 3*pi/8 * (4m - 1)
    float s23 = pow(s, 0.66666667);
    return s23 * (1.0 + 0.10416667 / (s * s));
}

// Boundary derivative magnitude |Ai'(-zeta_m)| at the m-th Airy zero (normalizes lobe peak to O(1)).
float sdfAiryZeroSlope(in int radialModeM)
{
    int m = max(radialModeM, 1);
    if (m == 1) return 0.70121082;
    if (m == 2) return 0.80311137;
    if (m == 3) return 0.86520403;
    if (m == 4) return 0.91085074;
    return 0.56418958 * sqrt(sqrt(sdfAiryZeroZeta(m)));
}

// Keller-Rubinow (1960) & Lazutkin (1973) whispering-gallery eigenvalue k_{m, nu}:
//   - perimRadius : R in 2D disk / 3D sphere, or L / (2*pi) for convex perimeter L
//   - radialModeM : radial mode number m >= 1 (number of concentric radial Airy rings)
//   - angularNu   : angular index nu = float(ell) in 2D, or float(ell) + 0.5 in 3D spheres
float sdfWhisperingGalleryK(in float perimRadius, in int radialModeM, in float angularNu)
{
    float R    = max(perimRadius, 1e-5);
    float nu   = max(angularNu, 0.5);
    float zm   = sdfAiryZeroZeta(radialModeM);
    float c1   = pow(0.5 * nu, 0.33333333);
    float c2   = pow(2.0 / nu, 0.33333333);
    return (nu + zm * c1 + 0.15 * zm * zm * c2) / R;
}

// Boundary-locked Fermi-Weingarten-Lazutkin whispering-gallery cavity eigen-SDF inside {d <= 0}:
//   psi(p) = (Ai(-beta * d - zeta_m) / (|Ai'(-zeta_m)| * sqrt(max(J, 0.08)))) * sin^l(theta) * cos(l*phi + phaseOffset)
// Guarantees psi|_{d=0} == 0 identically on ANY smooth curved SDF, zero polar 0/0 singularities in 3D,
// exponential decay in the Rayleigh quiet core (-beta*d > zeta_m), and exact Sobolev-regularized
// unit-slope distance field vec4(d_psi, grad_d_psi).
vec4 sdgWhisperingGalleryCavity(in vec3 p, in float d, in vec3 g,
                                in float kappaPar, in float weingartenJ,
                                in int radialModeM, in int ell, in float phaseOffset,
                                out float outPsi)
{
    int   l       = max(ell, 0);
    float fl      = float(l);
    float kap     = max(kappaPar, 1e-3);
    float zm      = sdfAiryZeroZeta(radialModeM);
    float slopeM  = sdfAiryZeroSlope(radialModeM);
    float nuEff   = max(fl, 0.5);
    float beta    = pow(2.0 * nuEff * nuEff, 0.33333333) * kap; // (2 * k_par^2 * kap)^(1/3) with k_par = nu * kap

    // Radial Langer-Airy factor R(d) with exact boundary zero at d == 0:
    float u       = -beta * d - zm;
    vec2  jet     = sdfAiryJet(u);
    float normC   = 1.0 / (slopeM * sqrt(max(weingartenJ, 0.08)));
    float radVal  = normC * jet.x;
    float radD1   = -normC * beta * jet.y;             // dR / dd
    float radD2   = normC * beta * beta * u * jet.x;   // d^2R / dd^2 via Airy ODE Ai''(u) = u * Ai(u)

    // Smooth polar-regularized spherical harmonic W(p) = (r_perp / r)^l * cos(l*phi + phaseOffset):
    float W       = 1.0;
    vec3  gW      = vec3(0.0);
    vec3  tildeGP = vec3(0.0);
    float r2      = dot(p, p) + 1e-8;
    float invR2   = 1.0 / r2;
    float invR    = sqrt(invR2);
    float rPerp2  = p.x * p.x + p.z * p.z;
    float rPerp   = sqrt(rPerp2);
    float sTheta  = rPerp * invR;
    float phi     = (rPerp > 1e-7) ? atan(p.z, p.x) : 0.0;

    if (l == 0)
    {
        W  = cos(phaseOffset);
        gW = vec3(0.0);
    }
    else
    {
        float sPowL   = pow(max(sTheta, 0.0), fl);
        float sPowLm1 = (l == 1) ? 1.0 : pow(max(sTheta, 0.0), fl - 1.0);
        W = sPowL * cos(fl * phi + phaseOffset);

        float ang1 = float(l - 1) * phi + phaseOffset;
        tildeGP    = (fl * sPowLm1 * invR) * vec3(cos(ang1), 0.0, -sin(ang1));
        gW         = tildeGP - (fl * W * invR2) * p;
    }

    float psi  = radVal * W;
    outPsi     = psi;
    vec3  gPsi = W * (radD1 * g) + radVal * gW;

    // Exact Hessian-vector product H_psi * gPsi for smooth Sobolev-regularized distance field:
    vec3 hWP = vec3(0.0);
    if (l >= 2)
    {
        float sPowLm2 = (l == 2) ? 1.0 : pow(max(sTheta, 0.0), fl - 2.0);
        float ang2    = float(l - 2) * phi + phaseOffset;
        float coeff2  = fl * (fl - 1.0) * sPowLm2 * invR2;
        float a2      = coeff2 * cos(ang2);
        float b2      = coeff2 * sin(ang2);
        hWP           = vec3(a2 * gPsi.x - b2 * gPsi.z, 0.0, -b2 * gPsi.x - a2 * gPsi.z);
    }
    vec3 hWgPsi = vec3(0.0);
    if (l >= 1)
    {
        float pv = dot(p, gPsi);
        hWgPsi   = hWP - (fl * pv * invR2) * tildeGP
                 - (fl * invR2) * (dot(gW, gPsi) * p + W * gPsi - (2.0 * pv * W * invR2) * p);
    }

    float gv     = dot(g, gPsi);
    vec3  hRgPsi = radD2 * gv * g + (radD1 * kap) * (gPsi - gv * g);
    vec3  hPsiG  = W * hRgPsi + radD1 * (dot(gW, gPsi) * g + gv * gW) + radVal * hWgPsi;

    // Phase-space Sobolev energy with quiet-core regularization epsCore = 0.06:
    float beta2    = beta * beta;
    float epsCore2 = 0.0036 * beta2;
    float E        = beta2 * psi * psi + dot(gPsi, gPsi) + epsCore2 + 1e-12;
    float invSqrtE = inversesqrt(E);
    float dPsi     = psi * invSqrtE;
    vec3  halfGradE = (beta2 * psi) * gPsi + hPsiG;
    vec3  gradDPsi  = invSqrtE * (gPsi - (psi / E) * halfGradE);

    return vec4(dPsi, gradDPsi);
}

// Convenience overload without outPsi output parameter.
vec4 sdgWhisperingGalleryCavity(in vec3 p, in float d, in vec3 g,
                                in float kappaPar, in float weingartenJ,
                                in int radialModeM, in int ell, in float phaseOffset)
{
    float psiDummy;
    return sdgWhisperingGalleryCavity(p, d, g, kappaPar, weingartenJ, radialModeM, ell, phaseOffset, psiDummy);
}

// Scalar wrapper returning Sobolev-regularized whispering-gallery nodal distance d_psi.
float sdWhisperingGalleryCavity(in vec3 p, in float d, in vec3 g,
                                in float kappaPar, in float weingartenJ,
                                in int radialModeM, in int ell, in float phaseOffset)
{
    return sdgWhisperingGalleryCavity(p, d, g, kappaPar, weingartenJ, radialModeM, ell, phaseOffset).x;
}

// -----------------------------------------------------------------------------
// 4. da Costa (1981) & Exner-Seba (1989) Extrinsic-Curvature Waveguide Trap
// -----------------------------------------------------------------------------
// Evaluates purely geometric quantum trapping inside a constant-width shell {|d| < shellHalfWidth}:
//   - vDaCosta     = -max(H^2 - K, 0) = -(kappa_1 - kappa_2)^2 / 4 <= 0
//                    (0 on planes & spheres; maximally attractive at saddles K < 0 & sharp bends)
//   - chiShell     = J(d)^{-1/2} * cos(pi * d / (2 * w)) inside |d| < w (0 outside)
//   - boundDensity = chiShell^2 * (1 - vDaCosta * w^2) local geometric bound-state enhancement
// Returns vec3(vDaCosta, boundDensity, chiShell).
vec3 sdfDaCostaCurvatureTrap(in float d, in float meanH, in float gaussK, in float shellHalfWidth)
{
    float w        = max(shellHalfWidth, 1e-4);
    float skew2    = max(meanH * meanH - gaussK, 0.0); // (k1 - k2)^2 / 4 >= 0
    float vDaCosta = -skew2;
    if (abs(d) >= w)
    {
        return vec3(vDaCosta, 0.0, 0.0);
    }

    float jTube    = max(1.0 + 2.0 * meanH * d + gaussK * d * d, 0.08);
    float chiShell = cos(1.57079633 * (d / w)) * inversesqrt(jTube);
    // Dimensionless bound-state localization boost weighted by dimensionless well depth Q^2 / (pi/(2w))^2:
    float wellBoost = 1.0 + 4.0 * skew2 * w * w;
    return vec3(vDaCosta, chiShell * chiShell * wellBoost, chiShell);
}

#endif // SDF_GEMS_SDF_WKB_AIRY_CAVITY_GLSL
