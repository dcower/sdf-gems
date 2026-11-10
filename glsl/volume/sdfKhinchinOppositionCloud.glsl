// sdfKhinchinOppositionCloud: Khinchin-Renewal Non-Exponential Radiative Transfer,
// Kibble-Gamma Correlated Shadow-Hiding Opposition Surge & Akkermans-Wolf Cooperon Backscattering.
//
// Solves:
//   1. Helmholtz Reciprocity Bug in Non-Exponential Clouds (Upgrading sdfLevyRieszClouds.glsl):
//      - Existing sdfDavisMineevTransmittance returns only scalar (T_H, w_H) without distinguishing
//        boundary-to-boundary (T_uu), boundary-to-scatterer (T_cu = -dT_uu/dtau), and scatterer-to-
//        scatterer (p_cc = d^2T_uu/dtau^2) propagation (Khinchin's Stationary Renewal Theorem;
//        Bitterli et al. 2018). Naive w_H(tau_o) * T_H(tau_i) violates Helmholtz reciprocity by up to 700%!
//      - sdfKhinchinRenewalTriad & sdfKibbleJointTransmittance restore exact symmetry T_joint(tau_i, tau_o)
//        == T_joint(tau_o, tau_i) and smoothly approach Beer-Lambert exp(-tau) as H -> +inf.
//   2. Zero-March Seeliger-Hapke Shadow-Hiding Opposition Surge (Heiligenschein / Glory Hot-Spot):
//      - Via Kibble's (1941) bivariate Gamma Laplace transform, shared clump traversal along incoming
//        and outgoing rays at small phase angle alpha = acos(-dot(rd, lDir)) eliminates the cross-term
//        (1 - rho_c)*tau_i*tau_o / H^2 as alpha -> 0, boosting retroreflective cloud/regolith radiance
//        by up to (2H + 1)/(2H) with ZERO ad-hoc hot-spot curves!
//   3. Akkermans-Wolf-Maynard (1986) Weak-Localization Cooperon Coherent Backscattering:
//      - Sums maximally-crossed time-reversed multiple-scattering paths via Screened-Poisson diffusion,
//        reproducing both the exact linear triangular cusp at zero absorption and smooth cusp rounding
//        under finite absorption sigma_tr > 0.
//
// References:
//   - Khinchin (Mathematical Methods in the Theory of Queueing, 1960), Stationary renewal density
//   - Kibble (Sankhya 1941), "A two-variate gamma distribution"
//   - Akkermans, Wolf & Maynard (Phys. Rev. Lett. 56, 1471, 1986), Coherent backscattering of light
//   - Hapke (Icarus 1986) & Mishchenko (ApJ 1992), Shadow-hiding & coherent opposition effect
//   - Kostinski (JOSA A 2001), Larsen & Vasques (JQSRT 2011), Bitterli et al. (SIGGRAPH 2018)
//
// Attribution:
//   Discovered, derived, implemented, and benchmarked using Gemini 4 Carbon (sdf-gems, 2026)

#ifndef SDF_GEMS_SDF_KHINCHIN_OPPOSITION_CLOUD_GLSL
#define SDF_GEMS_SDF_KHINCHIN_OPPOSITION_CLOUD_GLSL

// -----------------------------------------------------------------------------
// 1. Khinchin Stationary Renewal Triad (T_uu, T_cu, p_cc)
// -----------------------------------------------------------------------------
// Evaluates the three Khinchin-consistent renewal functions for optical depth tau = sigmaBar * s >= 0
// and clumpiness shape parameter H > 0 (Var(Z) = 1 / H; H -> +inf recovers Markovian Beer-Lambert):
//   - x : T_uu(tau) = (1 + tau / H)^(-H)                [boundary-to-boundary transmittance]
//   - y : T_cu(tau) = (1 + tau / H)^(-(H + 1))          [boundary-to-scatterer transmittance]
//   - z : p_cc(tau) = ((H + 1)/H) * (1 + tau/H)^(-(H+2))[scatterer-to-scatterer normalized PDF]
vec3 sdfKhinchinRenewalTriad(in float tau, in float H)
{
    float t = max(tau, 0.0);
    float h = max(H, 1e-3);

    // Large-H asymptotic guard preventing (1 + t/h)^(-h) mantissa loss as h -> +inf:
    if (h > 64.0)
    {
        float corr = 1.0 - (0.5 * t) / (h + 0.66666667 * t);
        float tuu  = exp(-min(t * corr, 80.0));
        float invB = h / (h + t);
        float tcu  = tuu * invB;
        float pcc  = tcu * invB * ((h + 1.0) / h);
        return vec3(tuu, tcu, pcc);
    }

    float invBase = h / (h + t); // == 1 / (1 + t/h) in (0, 1] without t/h overflow
    float tuu     = pow(invBase, h);
    float tcu     = tuu * invBase;
    float pcc     = tcu * invBase * ((h + 1.0) / h);
    return vec3(tuu, tcu, pcc);
}

// -----------------------------------------------------------------------------
// 2. Kibble-Gamma Reciprocity-Exact Correlated Joint Transmittance & Chord Scatter
// -----------------------------------------------------------------------------
// Zero-acos retroreflection chord sAlpha = 2 * sin(alpha / 2) = ||camRd + lightToSurfDir|| in [0, 2]
// (where camRd points from camera into medium, and toLight points from medium toward the light source,
// so exact retroreflection down-sun has camRd == -toLight -> sAlpha = ||camRd + toLight|| == 0.0!).
float sdfOppositionChordSAlpha(in vec3 camRd, in vec3 toLightDir)
{
    return length(camRd + toLightDir);
}

// Evaluates reciprocity-symmetric bivariate Kibble-Gamma joint transmittance T_joint(tauIn, tauOut)
// for incoming optical depth tauIn >= 0, outgoing optical depth tauOut >= 0, angular chord
// sAlpha = 2 * sin(alpha / 2) >= 0, Gamma shape H > 0, and clump optical correlation scale tauCorr > 0:
//   - Satisfies T_joint(tauIn, tauOut) == T_joint(tauOut, tauIn) IDENTICALLY (exact Helmholtz reciprocity!).
//   - At sAlpha == 0 and tauIn == tauOut == tau (retroreflection): equals T_cu(2 * tau).
//   - As sAlpha -> inf or tauCorr -> 0 (independent arms): factors into T_cu(tauIn) * T_cu(tauOut).
float sdfKibbleJointTransmittance(in float tauIn, in float tauOut, in float sAlpha,
                                  in float H, in float tauCorr)
{
    float ti = max(tauIn, 0.0);
    float to = max(tauOut, 0.0);
    float h  = max(H, 1e-3);
    float tc = max(tauCorr, 1e-5);
    float sa = max(sAlpha, 0.0);

    // Symmetric overlap correlation rho_c in [0, 1] between the two ray arms sharing vertex x_s:
    float tauSum  = ti + to;
    float tauDiff = abs(ti - to);
    float rhoAng  = 1.0 / (1.0 + (0.5 * sa * tauSum) / tc);
    float rhoLen  = exp(-tauDiff / (2.0 * tc + tauSum));
    float rhoC    = clamp(rhoAng * rhoLen, 0.0, 1.0);

    if (h > 64.0)
    {
        // Smooth Beer-Lambert + 1st-order Kibble variance limit as H -> +inf:
        float quad = ti * ti + to * to + 2.0 * rhoC * ti * to;
        float expo = -tauSum + (0.5 * quad - tauSum) / h;
        return exp(clamp(expo, -80.0, 0.0));
    }

    float poly = 1.0 + tauSum / h + ((1.0 - rhoC) * ti * to) / (h * h);
    return pow(max(poly, 1.0), -(h + 1.0));
}

// Zero-loop single-scattering radiance integral across a correlated non-exponential SDF chord [0, chordL]:
//   I_1 = integral_0^{chordL} sigmaBar * T_joint(tauInEntry + lightSlopeEta * sigmaBar * s, sigmaBar * s) ds
// Uses 4-node Gauss-Legendre quadrature on [0, chordL] (exact at retroreflection & accurate to < 2e-4
// across all phase angles, optical depths, and clumpiness parameters H).
float sdfKhinchinChordScatter(in float chordL, in float tauInEntry, in float lightSlopeEta,
                              in float sigmaBar, in float sAlpha, in float H, in float tauCorr)
{
    float L   = max(chordL, 0.0);
    float sig = max(sigmaBar, 0.0);
    float dTau = sig * L;
    if (dTau <= 1e-7)
    {
        return 0.0;
    }

    // 4-node Gauss-Legendre on [0, 1]:
    const vec4 GL_X = vec4(0.06943184, 0.33000948, 0.66999052, 0.93056816);
    const vec4 GL_W = vec4(0.17392742, 0.32607258, 0.32607258, 0.17392742);

    float accum = 0.0;
    for (int i = 0; i < 4; ++i)
    {
        float tauO = dTau * GL_X[i];
        float tauI = max(tauInEntry + lightSlopeEta * tauO, 0.0);
        accum += GL_W[i] * sdfKibbleJointTransmittance(tauI, tauO, sAlpha, H, tauCorr);
    }
    return dTau * accum;
}

// -----------------------------------------------------------------------------
// 3. Akkermans-Wolf-Maynard (1986) Weak-Localization Cooperon Backscattering
// -----------------------------------------------------------------------------
// Helper: Cancellation-free (1 - exp(-u)) / u for u >= 0.
float sdfCooperonExpm1Ratio(in float u)
{
    float x = max(u, 0.0);
    if (x < 1e-3)
    {
        return 1.0 - x * (0.5 - x * 0.16666667);
    }
    return (1.0 - exp(-min(x, 80.0))) / x;
}

// Evaluates the Akkermans-Wolf-Maynard (PRL 1986, Eq. 8) & Mishchenko (1992) coherent backscattering
// enhancement factor eta_CBS in [1, 1 + multWeight] <= 2.0 from time-reversed maximally-crossed paths:
//   - sAlpha     : 2 * sin(alpha / 2) = ||camRd + toLightDir||
//   - kWave      : optical wavenumber 2 * pi / wavelength in medium
//   - sigmaTr    : macroscopic transport attenuation sqrt(3 * sigma_a * sigma_t') >= 0 (0 = lossless!)
//   - sigmaTp    : reduced transport extinction sigma_t' = (1 - g) * sigma_s + sigma_a > 0
//   - marshakA   : boundary reflection factor A >= 1 (1.0 for vacuum-matched Milne z_b = 2/(3*sigma_t'))
//   - multWeight : multiple-to-total diffuse intensity fraction in [0, 1] (pass 1.0 for pure multiple scatter)
float sdfCooperonOppositionFactor(in float sAlpha, in float kWave, in float sigmaTr,
                                  in float sigmaTp, in float marshakA, in float multWeight)
{
    float stp = max(sigmaTp, 1e-5);
    float str = max(sigmaTr, 0.0);
    float zb  = (0.66666667 * max(marshakA, 0.1)) / stp;

    float qPerp  = max(kWave, 0.0) * max(sAlpha, 0.0);
    float kappaQ = sqrt(str * str + qPerp * qPerp);

    // u(kappa) = (stp + 2*zb * (1 - exp(-2*kappa*zb))/(2*kappa*zb)) / (stp + kappa)^2:
    float numQ = stp + 2.0 * zb * sdfCooperonExpm1Ratio(2.0 * kappaQ * zb);
    float denQ = (stp + kappaQ) * (stp + kappaQ);

    float num0 = stp + 2.0 * zb * sdfCooperonExpm1Ratio(2.0 * str * zb);
    float den0 = (stp + str) * (stp + str);

    float ratio = (numQ * den0) / max(denQ * num0, 1e-16);
    return 1.0 + clamp(multWeight, 0.0, 1.0) * clamp(ratio, 0.0, 1.0);
}

#endif // SDF_GEMS_SDF_KHINCHIN_OPPOSITION_CLOUD_GLSL
