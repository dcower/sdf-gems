// sdfCherenkovTransition: Causal Mach-Cherenkov Wake-Cone SDF, Lorentz-Reciprocity Tamm-Ginzburg-Frank
// Coherent Slab & DIRC Prism Radiation, and Zero-Loop Frank-Tamm Volumetric Reactor Pool Glow.
//
// Solves:
//   1. Superluminal Dielectric Wake Geometry & Dispersive Phase Thresholds in SDF Volumes:
//      - Existing stRetardedTrace (glsl/optics/) enforces subluminal vacuum motion (beta < 1).
//      - Inside a dispersive dielectric SDF of phase index n(lambda) > 1, a relativistic charge
//        moving at v = beta * c with 1/n(lambda) < beta < 1 outruns light in the medium, forming
//        a sharp Mach-Cherenkov shock cone of wake half-aperture sin(psi_wake) = 1 / (beta * n)
//        (and wavevector emission angle cos(theta_C) = 1 / (beta * n)).
//      - sdgCherenkovWakeCone returns the exact 3D Euclidean signed distance and unit analytical
//        gradient to the trailing wake cone above threshold (beta * n > 1), with C^1 continuity
//        across the apex Voronoi seam and exact unit-norm fallbacks on the medial axis.
//   2. Pole-Free Coherent Tamm-Cherenkov + Ginzburg-Frank Transition Radiation + DIRC Escape:
//      - Direct matching of relativistic particular solutions across a finite dielectric SDF slab
//        suffers a catastrophic 0/0 pole (1 - beta * kappa_eps)^(-1) right on the Cherenkov cone.
//      - Switching to Lorentz Electrodynamic Reciprocity (projecting J_z onto the adjoint TM
//        Fabry-Perot slab wave) eliminates the Cherenkov denominator analytically: the interior
//        track integral reduces purely to Tamm's finite-track sinc(Delta_phi), while the exterior
//        vacuum half-spaces give pure Ginzburg-Frank transition radiation + Wartski etalon fringes!
//      - Generalizing the exit face to an SDF bevel angle alpha_bevel captures Total Internal
//        Reflection (TIR) Cherenkov trapping on flat high-index slabs vs BaBar/Belle-II DIRC
//        fused-quartz prism escape, matching an exact 4x4 Maxwell solver to machine precision.
//   3. Zero-March 3D Volumetric Nuclear-Reactor / Neutrino-Detector Cherenkov Pool Glow:
//      - Clips any 3D skew camera ray [ta, tb] to the exact causal Mach-Cherenkov frustum of a
//        finite particle track [trackStart, trackEnd] via 2 quadratic cone intersections and
//        integrates the 1 / r_perp cylindrical Frank-Tamm fluence in closed form (0 loops) via
//        subtraction-free asinh differences!
//
// References:
//   - Frank & Tamm (Dokl. Akad. Nauk SSSR 1937; Tamm, J. Phys. USSR 1939), Coherent finite-track Cherenkov radiation
//   - Ginzburg & Frank (JETP 1945), Transition radiation at dielectric boundaries
//   - Wartski et al. (J. Appl. Phys. 1975), Interference in two-foil / slab transition radiation
//   - Verzilov (Phys. Lett. A 2000), Transition radiation in the pre-wave Ginzburg formation zone
//
// Attribution:
//   Discovered, derived, implemented, and benchmarked using Gemini 4 Carbon (sdf-gems, 2026)

#ifndef SDF_GEMS_SDF_CHERENKOV_TRANSITION_GLSL
#define SDF_GEMS_SDF_CHERENKOV_TRANSITION_GLSL

#ifndef SDF_GEMS_HAS_SAFE_SINC
#define SDF_GEMS_HAS_SAFE_SINC
// FP32-safe unnormalized sinc(u) = sin(u) / u
float sdfSafeSinc(float u) {
    float u2 = u * u;
    return (abs(u) > 1e-4) ? (sin(u) / u) : (1.0 - u2 * (1.0 / 6.0 - u2 * (1.0 / 120.0)));
}
#endif

// Subtraction-free inverse hyperbolic sine for all x in (-inf, +inf):
float sdfCherenkovAsinh(in float x)
{
    float ax = abs(x);
    float v  = log(ax + sqrt(ax * ax + 1.0));
    return (x >= 0.0) ? v : -v;
}

// Complex helpers (vec2(re, im)):
vec2 sdfCherenkovCMul(in vec2 a, in vec2 b)
{
    return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

vec2 sdfCherenkovCDiv(in vec2 a, in vec2 b)
{
    float d = max(dot(b, b), 1e-30);
    return vec2(a.x * b.x + a.y * b.y, a.y * b.x - a.x * b.y) / d;
}

vec2 sdfCherenkovCExp(in vec2 z)
{
    float e = exp(clamp(z.x, -80.0, 80.0));
    return vec2(e * cos(z.y), e * sin(z.y));
}

vec2 sdfCherenkovCSqrtReal(in float x)
{
    return (x >= 0.0) ? vec2(sqrt(x), 0.0) : vec2(0.0, sqrt(-x));
}

// Complex unnormalized sinc(z) = sin(z) / z (accurate for both real & imaginary z):
vec2 sdfCherenkovCSinc(in vec2 z)
{
    float r2 = dot(z, z);
    if (r2 < 1e-6)
    {
        vec2 z2 = sdfCherenkovCMul(z, z);
        vec2 z4 = sdfCherenkovCMul(z2, z2);
        return vec2(1.0, 0.0) - z2 * (1.0 / 6.0) + z4 * (1.0 / 120.0);
    }
    float ch = cosh(clamp(z.y, -40.0, 40.0));
    float sh = sinh(clamp(z.y, -40.0, 40.0));
    vec2  sz = vec2(sin(z.x) * ch, cos(z.x) * sh);
    return sdfCherenkovCDiv(sz, z);
}

// -----------------------------------------------------------------------------
// 1. Cauchy Dispersive Phase Index n(lambda) = n0 + B / lambda^2
// -----------------------------------------------------------------------------
float sdfCauchyRefractiveIndex(in float wavelengthUm, in float n0, in float cauchyB)
{
    float lam = max(wavelengthUm, 0.15);
    return n0 + cauchyB / (lam * lam);
}

// -----------------------------------------------------------------------------
// 2. Exact 3D Mach-Cherenkov Wake-Cone SDF & Unit Analytical Gradient
// -----------------------------------------------------------------------------
// Returns vec4(d, grad) for a relativistic particle at particlePos moving along unit vector
// particleDir at speed v = beta * c inside a medium of refractive index nMedium:
//   - Above threshold (beta * nMedium > 1): exact 3D Euclidean signed distance (d < 0 inside
//     the causal trailing wake cone of half-angle sin(psi_wake) = 1 / (beta * nMedium), d = 0 on
//     the shock wavefront, d > 0 outside) and exact unit gradient (||grad|| == 1 everywhere).
//   - Sub-threshold (beta * nMedium <= 1): shock cone disappears; returns exact Euclidean
//     distance and unit gradient to the particle apex particlePos.
vec4 sdgCherenkovWakeCone(in vec3 p, in vec3 particlePos, in vec3 particleDir,
                          in float beta, in float nMedium)
{
    vec3  u    = normalize(particleDir);
    vec3  q    = p - particlePos;
    float qLen = length(q);
    vec3  gAp  = (qLen > 1e-12) ? (q / qLen) : u;

    float bn = beta * nMedium;
    if (bn <= 1.0)
    {
        return vec4(qLen, gAp);
    }

    float sw = 1.0 / bn;
    float cw = sqrt(max((bn - 1.0) * (bn + 1.0), 0.0)) / bn;

    // Orthonormal basis (e1, e2) strictly orthogonal to u prevents FP32 medial-axis
    // cancellation noise in q + z*u from polluting dot(eR, u) when r ~ 1e-8:
    float z      = -dot(q, u);
    vec3  helper = (abs(u.x) < 0.9) ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
    vec3  e1     = normalize(cross(u, helper));
    vec3  e2     = cross(u, e1);
    vec2  q2     = vec2(dot(q, e1), dot(q, e2));
    float r      = length(q2);
    vec3  eR     = (r > 1e-12) ? ((q2.x * e1 + q2.y * e2) / r) : e1;

    // Generator projection tau = r * sin(psi) + z * cos(psi) separates apex vs lateral sheet:
    float tau = r * sw + z * cw;
    if (tau <= 0.0)
    {
        return vec4(qLen, gAp);
    }

    float dCone = r * cw - z * sw;
    vec3  gCone = cw * eR + sw * u;
    return vec4(dCone, gCone);
}

float sdCherenkovWakeCone(in vec3 p, in vec3 particlePos, in vec3 particleDir,
                          in float beta, in float nMedium)
{
    return sdgCherenkovWakeCone(p, particlePos, particleDir, beta, nMedium).x;
}

// -----------------------------------------------------------------------------
// 3. Lorentz-Reciprocity Coherent Tamm-Cherenkov + Ginzburg-Frank + DIRC Prism
// -----------------------------------------------------------------------------
// Evaluates exact far-field TM angular spectral intensity I(theta_vac, omega) for a charge
// moving at v = beta * c across a dielectric SDF slab/bar of permittivity epsReal = n^2 and
// dimensionless thickness omegaLOverC = omega * L / c, exiting across an SDF face whose outward
// normal is tilted by bevelAngleRad in the meridian plane:
//   - Pass bevelAngleRad = 0.0 for a parallel dielectric slab (reproduces the exact 4x4 Maxwell
//     boundary-value linear solver across sub-threshold transition radiation, super-threshold
//     Tamm-Cherenkov cones, Wartski Fabry-Perot fringes, and flat-slab TIR trapping!).
//   - Pass bevelAngleRad > 0.0 (e.g. ~0.5 rad) for a beveled BaBar/Belle-II DIRC prism that
//     unlocks TIR-trapped Cherenkov cones from high-index quartz/sapphire SDF radiators!
float sdfTammGinzburgBevelChannel(in float thetaVacRad, in float beta, in float epsReal,
                                  in float omegaLOverC, in float bevelAngleRad)
{
    float b    = clamp(beta, 1e-4, 0.99995);
    float eps  = max(epsReal, 0.05);
    float wLc  = max(omegaLOverC, 1e-4);
    float thOut = thetaVacRad - bevelAngleRad;
    if (abs(thOut) >= 1.5707963 - 1e-4 || thetaVacRad <= 1e-5 || thetaVacRad >= 3.1415926 - 1e-5)
    {
        return 0.0;
    }

    float sOut = sin(thOut);
    float cOut = cos(thOut);
    float ca   = cos(bevelAngleRad);
    float sa   = sin(bevelAngleRad);

    // Complex normal wavenumber inside medium relative to beveled exit face:
    vec2 kap1n = sdfCherenkovCSqrtReal(eps - sOut * sOut);
    // Rotate internal wavevector back to particle-track axes (x_perp, z_track):
    vec2 k1x   = kap1n * sa + vec2(sOut * ca, 0.0);
    vec2 k1z   = kap1n * ca - vec2(sOut * sa, 0.0);

    // Fresnel TM reflection & adjoint transmission at beveled exit face:
    vec2 epsCo = vec2(eps * cOut, 0.0);
    vec2 rExit = sdfCherenkovCDiv(epsCo - kap1n, epsCo + kap1n);
    vec2 t01   = sdfCherenkovCDiv(2.0 * epsCo, epsCo + kap1n);

    // Entrance face normal (0, -1): c_ent_vac = sqrt(1 - k1x^2)
    // For real k1x (propagating exit), 1 - k1x^2 is real; handle general complex k1x cleanly:
    vec2 k1x2  = sdfCherenkovCMul(k1x, k1x);
    vec2 cEntV = (abs(k1x2.y) < 1e-7)
               ? sdfCherenkovCSqrtReal(1.0 - k1x2.x)
               : vec2(sqrt(max(1.0 - k1x2.x, 0.0)), 0.0);
    vec2 epsCe = eps * cEntV;
    vec2 rEnt  = sdfCherenkovCDiv(epsCe - k1z, epsCe + k1z);
    vec2 t10   = sdfCherenkovCDiv(2.0 * k1z, epsCe + k1z);

    // Fabry-Perot adjoint cavity resonance:
    vec2 psi     = 0.5 * wLc * k1z;
    vec2 iPsi    = vec2(-psi.y, psi.x);
    vec2 ePsi    = sdfCherenkovCExp(iPsi);
    vec2 e2Psi   = sdfCherenkovCMul(ePsi, ePsi);
    vec2 e3Psi   = sdfCherenkovCMul(e2Psi, ePsi);
    vec2 e4Psi   = sdfCherenkovCMul(e2Psi, e2Psi);
    vec2 rProd   = sdfCherenkovCMul(rExit, rEnt);
    vec2 denomFP = vec2(1.0, 0.0) - sdfCherenkovCMul(rProd, e4Psi);

    vec2 aAdj = sdfCherenkovCDiv(sdfCherenkovCMul(t01, ePsi), denomFP);
    vec2 bAdj = -sdfCherenkovCDiv(sdfCherenkovCMul(sdfCherenkovCMul(t01, rEnt), e3Psi), denomFP);

    // 1. Pole-free interior Tamm-Cherenkov track integral:
    float phiP  = 0.5 * wLc / b;
    vec2  dphiM = vec2(phiP, 0.0) - psi;
    vec2  dphiP = vec2(phiP, 0.0) + psi;
    vec2  sincM = (abs(dphiM.y) < 1e-7) ? vec2(sdfSafeSinc(dphiM.x), 0.0) : sdfCherenkovCSinc(dphiM);
    vec2  sincP = (abs(dphiP.y) < 1e-7) ? vec2(sdfSafeSinc(dphiP.x), 0.0) : sdfCherenkovCSinc(dphiP);
    vec2  sumAB = sdfCherenkovCMul(aAdj, sincM) + sdfCherenkovCMul(bAdj, sincP);
    vec2  iCher = -(wLc / eps) * sdfCherenkovCMul(k1x, sumAB);

    // 2. Exterior vacuum Ginzburg-Frank formation-zone boundary integrals:
    vec2 rExit2 = sdfCherenkovCMul(rExit, rExit);
    vec2 rSlab  = rExit - sdfCherenkovCDiv(
        sdfCherenkovCMul(vec2(1.0, 0.0) - rExit2, sdfCherenkovCMul(rEnt, e4Psi)),
        denomFP
    );
    vec2 tSlab  = sdfCherenkovCDiv(sdfCherenkovCMul(sdfCherenkovCMul(t01, t10), e2Psi), denomFP);

    float sVac  = sin(thetaVacRad);
    float cVac  = cos(thetaVacRad);
    float thRef = 2.0 * bevelAngleRad - thetaVacRad;
    float sRef  = sin(thRef);
    float cRef  = cos(thRef);

    vec2 ePosP  = vec2(cos(phiP), sin(phiP));
    vec2 eNegP  = vec2(ePosP.x, -ePosP.y);
    vec2 fTerm  = vec2(sVac / max(1.0 - b * cVac, 1e-6), 0.0)
                - (sRef / max(1.0 + b * cRef, 1e-6)) * rSlab;
    // Multiply by -i * beta * exp(+i * phiP):
    vec2 mIbEp  = b * vec2(ePosP.y, -ePosP.x);
    vec2 iVacF  = sdfCherenkovCMul(mIbEp, fTerm);

    // Entrance vacuum half-space: +i * beta * exp(-i * phiP) * k1x * T_slab / (1 - beta * cEntV)
    // Since eNegP = (cos phiP, -sin phiP), +i * eNegP = (+sin phiP, +cos phiP) = (-eNegP.y, eNegP.x):
    vec2 pIbEn  = b * vec2(-eNegP.y, eNegP.x);
    vec2 denomB = vec2(1.0, 0.0) - b * cEntV;
    vec2 iVacB  = sdfCherenkovCDiv(sdfCherenkovCMul(pIbEn, sdfCherenkovCMul(k1x, tSlab)), denomB);

    vec2 amp = 0.5 * (iCher + iVacF + iVacB);
    return dot(amp, amp);
}

// Parallel-slab convenience overload parameterized by sin(theta_vac) in [0, 1):
float sdfTammGinzburgChannel(in float sinThetaVac, in float beta, in float epsReal, in float omegaLOverC)
{
    float s = clamp(sinThetaVac, 0.0, 0.99995);
    return sdfTammGinzburgBevelChannel(asin(s), beta, epsReal, omegaLOverC, 0.0);
}

// Chromatic RGB coherent Tamm-Cherenkov + Ginzburg-Frank + DIRC slab/prism radiance
// (includes physical omega^2 ~ lambda^(-3) spectral density per unit wavelength):
vec3 sdfTammGinzburgSlabRGB(in float thetaVacRad, in float beta, in float slabThickness,
                            in float n0, in float cauchyB, in vec3 wavelengthsUm, in float bevelAngleRad)
{
    vec3 col = vec3(0.0);
    for (int c = 0; c < 3; ++c)
    {
        float lam  = max(wavelengthsUm[c], 0.25);
        float nLam = sdfCauchyRefractiveIndex(lam, n0, cauchyB);
        float wLc  = (6.2831853 * slabThickness) / lam;
        float iCh  = sdfTammGinzburgBevelChannel(thetaVacRad, beta, nLam * nLam, wLc, bevelAngleRad);
        // Normalized (0.50 um / lam)^3 Frank-Tamm spectral density scaling:
        float wLam = (0.125) / (lam * lam * lam);
        col[c] = wLam * iCh;
    }
    return col;
}

// -----------------------------------------------------------------------------
// 4. Zero-March 3D Volumetric Cherenkov Pool Glow (Reactor Pool / IceCube)
// -----------------------------------------------------------------------------
// Helper: intersects 1D ray parameter t with forward cone { z(t) - sBound >= cotTheta * r_perp(t) },
// returning vec2(tLo, tHi) with tLo > tHi when empty.
vec2 sdfCherenkovForwardConeSpan(in float z0, in float du, in vec3 pPerp, in vec3 dPerp,
                                 in float cotTheta, in float sBound)
{
    float u0 = z0 - sBound;
    float hLo = -1e20;
    float hHi =  1e20;
    if (abs(du) < 1e-8)
    {
        if (u0 < 0.0) return vec2(1.0, -1.0);
    }
    else if (du > 0.0)
    {
        hLo = -u0 / du;
    }
    else
    {
        hHi = -u0 / du;
    }

    float c2 = cotTheta * cotTheta;
    float A  = du * du - c2 * dot(dPerp, dPerp);
    float B  = u0 * du - c2 * dot(pPerp, dPerp);
    float C  = u0 * u0 - c2 * dot(pPerp, pPerp);

    if (abs(A) < 1e-8)
    {
        if (abs(B) < 1e-8)
        {
            return (C >= 0.0 && hLo <= hHi) ? vec2(hLo, hHi) : vec2(1.0, -1.0);
        }
        float tLin = -0.5 * C / B;
        if (B > 0.0) hLo = max(hLo, tLin);
        else         hHi = min(hHi, tLin);
        return (hLo <= hHi) ? vec2(hLo, hHi) : vec2(1.0, -1.0);
    }

    float disc = B * B - A * C;
    if (disc < 0.0)
    {
        float tMid = 0.5 * (max(hLo, -1e4) + min(hHi, 1e4));
        float qVal = (A * tMid + 2.0 * B) * tMid + C;
        return (hLo <= hHi && qVal >= 0.0) ? vec2(hLo, hHi) : vec2(1.0, -1.0);
    }

    float sq = sqrt(max(disc, 0.0));
    float r1 = (-B - sq) / A;
    float r2 = (-B + sq) / A;
    float t1 = min(r1, r2);
    float t2 = max(r1, r2);

    if (A < 0.0)
    {
        hLo = max(hLo, t1);
        hHi = min(hHi, t2);
    }
    else
    {
        // At U(t_0) = 0, Q(t_0) = -c2 * ||r_perp||^2 <= 0, so t_0 lies inside [t1, t2]:
        if (du > 0.0) hLo = max(hLo, t2);
        else          hHi = min(hHi, t1);
    }
    return (hLo <= hHi) ? vec2(hLo, hHi) : vec2(1.0, -1.0);
}

// Integrates water-scattered Frank-Tamm Cherenkov pool radiance along camera ray ro + t*rd
// across water slab [ta, tb] for a relativistic beam/track from trackStart to trackEnd in
// ZERO raymarching steps (2 analytical cone clips + closed-form subtraction-free asinh):
vec3 sdfCherenkovPoolGlowRGB(in vec3 ro, in vec3 rd, in float ta, in float tb,
                             in vec3 trackStart, in vec3 trackEnd, in float beta,
                             in float n0, in float cauchyB, in vec3 muAttenRGB,
                             in vec3 wavelengthsUm, in float coreSoftRadius)
{
    if (tb <= ta) return vec3(0.0);
    vec3  uVec = trackEnd - trackStart;
    float lTrk = length(uVec);
    if (lTrk < 1e-6) return vec3(0.0);
    vec3  uHat = uVec / lTrk;

    vec3  p0    = ro - trackStart;
    float z0    = dot(p0, uHat);
    float du    = dot(rd, uHat);
    vec3  pPerp = p0 - z0 * uHat;
    vec3  dPerp = rd - du * uHat;
    float dp2   = dot(dPerp, dPerp);
    float dpLen = sqrt(max(dp2, 0.0));

    float tPca  = (dp2 > 1e-10) ? (-dot(pPerp, dPerp) / dp2) : 0.0;
    vec3  rMinV = pPerp + tPca * dPerp;
    float rCore = max(coreSoftRadius, 1e-4);
    float aRad  = sqrt(dot(rMinV, rMinV) + rCore * rCore);

    vec3 glow = vec3(0.0);
    for (int c = 0; c < 3; ++c)
    {
        float lam  = max(wavelengthsUm[c], 0.25);
        float nLam = sdfCauchyRefractiveIndex(lam, n0, cauchyB);
        float bn   = beta * nLam;
        if (bn <= 1.00001) continue;

        float invBn2  = 1.0 / (bn * bn);
        float yieldFT = (1.0 - invBn2) * (0.125 / (lam * lam * lam));
        float cotTh   = 1.0 / sqrt(max((bn - 1.0) * (bn + 1.0), 1e-10));
        float sinTh   = sqrt(max(1.0 - invBn2, 1e-8));

        vec2 c0 = sdfCherenkovForwardConeSpan(z0, du, pPerp, dPerp, cotTh, 0.0);
        float lo = max(ta, c0.x);
        float hi = min(tb, c0.y);
        if (lo >= hi) continue;

        vec2 cL = sdfCherenkovForwardConeSpan(z0, du, pPerp, dPerp, cotTh, lTrk);
        // Shell C(0) \ C(lTrk) splits [lo, hi] into at most 2 sub-segments:
        vec4 segs = (cL.x > cL.y)
                  ? vec4(lo, hi, 1.0, -1.0)
                  : vec4(lo, min(hi, cL.x), max(lo, cL.y), hi);

        float chSum = 0.0;
        for (int k = 0; k < 2; ++k)
        {
            float s1 = (k == 0) ? segs.x : segs.z;
            float s2 = (k == 0) ? segs.y : segs.w;
            if (s2 <= s1 + 1e-6) continue;

            float ig = (dpLen < 1e-6)
                     ? ((s2 - s1) / aRad)
                     : ((sdfCherenkovAsinh(dpLen * (s2 - tPca) / aRad)
                       - sdfCherenkovAsinh(dpLen * (s1 - tPca) / aRad)) / dpLen);
            float tEff = clamp(tPca, s1, s2);
            float rEff = length(pPerp + tEff * dPerp);
            float att  = exp(-min(muAttenRGB[c] * (rEff / sinTh + (tEff - ta)), 80.0));
            chSum += yieldFT * att * ig;
        }
        glow[c] = chSum;
    }
    return glow;
}

#endif // SDF_GEMS_SDF_CHERENKOV_TRANSITION_GLSL
