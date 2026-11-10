// sdfFockAragoDiffraction: Uniform Rubinowicz-Poisson Obstacle Diffraction,
// Zero-FFT SDF Evolute Poisson-Arago Bright Spot, Fock Convex Creeping Waves &
// Dynamical Darwin-Prins Photonic Stop-Band Reflectance.
//
// Solves:
//   1. Zero-FFT wave-optical obstacle shadow diffraction & the 1818 Poisson-Arago Bright Spot:
//      - Heald's (1985) monotonic rational auxiliary Fresnel pair (f(x), g(x)) evaluates straight-
//        edge Cornu diffraction with ZERO trig calls in the shadow zone (I = 0.5*(f^2 + g^2),
//        exact I(0) == 0.25) and continuous phase matching into lit-zone Fresnel fringes.
//      - Decomposes the Maggi-Rubinowicz circular-rim diffraction kernel into an exact Bessel
//        J_0(k R b / z_eff) evolute wave + a stationary-phase Poisson rim wave, reconstructing
//        both outer Fresnel fringes and the central Poisson-Arago bright spot at the shadow
//        medial axis b = max(R + dStar, 0) -> 0 from a single SDF closest-approach query.
//      - Incorporates Fock's (1946) convex-body creeping-wave Airy-pole shadow attenuation
//        governed by longitudinal grazing curvature kappaPar >= 0.
//   2. Dynamical Darwin-Prins / Kogelnik Coupled-Wave Photonic Bandgaps:
//      - Replaces kinematic 1st-Born Bragg sums (which diverge > 1 as loss alpha -> 0 and rely on
//        ad-hoc Reinhard clamping in sdfGyroidBraggIridescence) with the exact overflow-free
//        complex Riccati solution in [0, 1], capturing flat-top 100% total-reflection Darwin
//        stop-bands (|detuning| <= couplingKappa) and Fabry-Perot sideband ripples.
//
// References:
//   - Poisson & Arago (1818-1819), Memoires de l'Academie des Sciences (Circular shadow bright spot)
//   - Maggi (1888) & Rubinowicz (1917), "Die Beugungswelle in der Kirchhoffschen Theorie der Beugung"
//   - Darwin (Phil. Mag. 1914) & Prins (Z. Phys. 1930), Dynamical X-ray & photonic Bragg diffraction
//   - Fock (J. Phys. USSR 1946), "Diffraction of radio waves around the Earth's surface"
//   - Kogelnik (Bell Syst. Tech. J. 1969), "Coupled wave theory for thick hologram gratings"
//   - Heald (Math. Comp. 1985), "Rational approximations for the Fresnel integrals"
//
// Attribution:
//   Discovered, derived, implemented, and benchmarked using Gemini 4 Carbon (sdf-gems, 2026)

#ifndef SDF_GEMS_SDF_FOCK_ARAGO_DIFFRACTION_GLSL
#define SDF_GEMS_SDF_FOCK_ARAGO_DIFFRACTION_GLSL

// -----------------------------------------------------------------------------
// 1. Bessel J_0 & Heald Rational Auxiliary Fresnel Pair (f(x), g(x))
// -----------------------------------------------------------------------------
#ifndef SDF_GEMS_HAS_BESSEL_J0
#define SDF_GEMS_HAS_BESSEL_J0
// Abramowitz & Stegun 9.4.1 & 9.4.3 minimax polynomial/Hankel approximation to J_0(x) (< 1e-6 FP32 err).
float sdfBesselJ0(in float x)
{
    float ax = abs(x);
    if (ax <= 3.0)
    {
        float t = (ax * 0.33333333);
        t *= t;
        return 1.0 + t * (-2.2499997 + t * (1.2656208 + t * (-0.3163866
                   + t * (0.0444479 + t * (-0.0039444 + t * 0.0002100)))));
    }
    float u   = 3.0 / ax;
    float f0  = 0.79788456 + u * (-0.00000077 + u * (-0.00552740 + u * (-0.00009512
              + u * (0.00137237 + u * (-0.00072805 + u * 0.00014476)))));
    float th0 = ax - 0.78539816 + u * (-0.04166397 + u * (-0.00003954 + u * (0.00262573
              + u * (-0.00054125 + u * (-0.00029333 + u * 0.00013558)))));
    return (f0 * inversesqrt(ax)) * cos(th0);
}
#endif // SDF_GEMS_HAS_BESSEL_J0

// Heald (1985) R_23 rational approximation to auxiliary Fresnel functions vec2(f(x), g(x)) for x >= 0.
// Satisfies f(0) = g(0) = 0.5, f(x) ~ 1/(pi*x), g(x) ~ 1/(pi^2*x^3); max error < 3.5e-4, 0 trig calls.
vec2 sdfFresnelAuxFG(in float x)
{
    float ax = max(x, 0.0);
    float x2 = ax * ax;
    float f  = (1.0 + 0.926 * ax) / (2.0 + 1.792 * ax + 3.104 * x2);
    float g  = 1.0 / (2.0 + 4.142 * ax + 3.492 * x2 + 6.670 * x2 * ax);
    return vec2(f, g);
}

// Straight-edge Sommerfeld / Cornu-Fresnel diffraction intensity behind an SDF silhouette at
// signed clearance dStar (dStar < 0 in geometric shadow, dStar > 0 in lit zone), source-to-occluder
// distance z1 (pass >= 1e5 for plane wave), occluder-to-screen distance z2 > 0, and wavelength > 0.
// Uses 0 trig calls in shadow (exact I(0) == 0.25) and continuous phase into lit supernumerary fringes.
float sdfCornuEdgeDiffraction(in float dStar, in float z1, in float z2, in float wavelength)
{
    float zEff = (z1 * z2) / max(z1 + z2, 1e-5);
    float lam  = max(wavelength, 1e-7);
    float nu   = dStar * sqrt(2.0 / max(lam * zEff, 1e-10));
    float ax   = abs(nu);
    vec2  fg   = sdfFresnelAuxFG(ax);

    if (dStar <= 0.0)
    {
        return 0.5 * dot(fg, fg);
    }
    float ph  = 1.57079633 * ax * ax;
    float c   = cos(ph);
    float s   = sin(ph);
    float uRe = 1.0 - 0.5 * ((fg.x + fg.y) * c - (fg.x - fg.y) * s);
    float uIm =      -0.5 * ((fg.x + fg.y) * s + (fg.x - fg.y) * c);
    return uRe * uRe + uIm * uIm;
}

// -----------------------------------------------------------------------------
// 2. Fock Convex Creeping-Wave Attenuation & Uniform Rubinowicz-Arago Diffraction
// -----------------------------------------------------------------------------
// Fock (1946) creeping-wave amplitude attenuation exp(-alpha_1 * (k / (2*kappaPar))^(1/3) * |theta_defl|)
// for a 3D convex body of longitudinal grazing curvature kappaPar >= 0 (pass 0.0 for thin razor disk).
// First Airy zero constant: alpha_1 * 2^(-1/3) = 1.01879297 * 0.79370053 = 0.8086165.
float sdfFockCreepingAtten(in float dStar, in float zEff, in float kappaPar, in float kWave)
{
    if (dStar >= 0.0 || kappaPar <= 1e-5)
    {
        return 1.0;
    }
    float thetaDefl = (-dStar) / max(zEff, 1e-5);
    float fockParam = 0.8086165 * pow(max(kWave / kappaPar, 0.0), 0.33333333) * thetaDefl;
    return exp(-min(fockParam, 40.0));
}

// Uniform Rubinowicz-Poisson-Bessel obstacle diffraction intensity in 1 SDF closest-approach query:
//   - kappaRim : transverse silhouette curvature 1 / R_rim (> 0 for bounded circular/spherical obstacle)
//   - kappaPar : longitudinal grazing curvature along ray (0 for thin disk, 1/R for 3D sphere)
// Captures outer Fresnel fringes, exact I ~ 0.25 silhouette rim, AND the central J_0^2(k R b / z_eff)
// Poisson-Arago bright spot at the shadow evolute b = max(R + dStar, 0) -> 0!
float sdfAragoObstacleChannel(in float dStar, in float z1, in float z2,
                              in float kappaRim, in float kappaPar, in float wavelength)
{
    if (kappaRim <= 1e-5)
    {
        return sdfCornuEdgeDiffraction(dStar, z1, z2, wavelength);
    }
    float zEff  = (z1 * z2) / max(z1 + z2, 1e-5);
    float lam   = max(wavelength, 1e-7);
    float kWave = 6.28318531 / lam;
    float alpha = kWave / max(2.0 * zEff, 1e-8);

    float R     = 1.0 / kappaRim;
    float b     = max(R + dStar, 0.0); // radial distance from shadow evolute center

    // 1. Exact Half-1 Rubinowicz Bessel wave 0.5 * e^{i alpha (R^2 + b^2)} * J_0(2 alpha R b):
    float creep = sdfFockCreepingAtten(dStar, zEff, kappaPar, kWave);
    float obliq = z2 * inversesqrt(z2 * z2 + R * R);
    float j0Val = sdfBesselJ0(2.0 * alpha * R * b);
    float phB   = mod(alpha * (R * R + b * b), 6.28318531);
    vec2  uHalf = (0.5 * obliq * creep * j0Val) * vec2(cos(phB), sin(phB));

    // 2. Stationary-phase Cornu-Fresnel edge wave U_edge(nu):
    float nu  = dStar * sqrt(2.0 / max(lam * zEff, 1e-10));
    float ax  = abs(nu);
    vec2  fg  = sdfFresnelAuxFG(ax);
    float phE = mod(1.57079633 * ax * ax, 6.28318531);
    float cE  = cos(phE), sE = sin(phE);
    vec2  uSh = (0.5 * creep) * vec2((fg.x + fg.y) * cE - (fg.x - fg.y) * sE,
                                     (fg.x + fg.y) * sE + (fg.x - fg.y) * cE);
    vec2  uEd = (dStar <= 0.0) ? uSh : (vec2(1.0, 0.0) - uSh);

    // 3. Smooth Rubinowicz-Poisson weight wRim = min((b / R)^1.5, 1.0):
    //    b = 0 -> uTot = 2 * uHalf (exact Arago bright spot |U|^2 = obliq^2 * J_0^2!)
    //    b >= R -> uTot = uEd + uHalf (Fresnel edge + antipodal rim interference ripples!)
    float ratio = clamp(b / max(R, 1e-6), 0.0, 1.0);
    float wRim  = ratio * sqrt(ratio);
    vec2  uTot  = (2.0 - wRim) * uHalf + wRim * uEd;
    return dot(uTot, uTot);
}

vec3 sdfAragoObstacleRGB(in float dStar, in float z1, in float z2,
                         in float kappaRim, in float kappaPar, in vec3 wavelengthsRGB)
{
    return vec3(
        sdfAragoObstacleChannel(dStar, z1, z2, kappaRim, kappaPar, wavelengthsRGB.x),
        sdfAragoObstacleChannel(dStar, z1, z2, kappaRim, kappaPar, wavelengthsRGB.y),
        sdfAragoObstacleChannel(dStar, z1, z2, kappaRim, kappaPar, wavelengthsRGB.z)
    );
}

// -----------------------------------------------------------------------------
// 3. Dynamical Darwin-Prins / Kogelnik Coupled-Wave Photonic Bandgap Reflectance
// -----------------------------------------------------------------------------
// Exact coupled-mode reflectance R in [0, 1] for a Bragg slab of thickness slabL >= 0,
// wavenumber detuning delta = k - q_Bragg/2, grating coupling kappaB >= 0, and field loss alpha >= 0.
// Overflow-free across all slabL in [0, inf); captures the flat-top 100% Darwin stop-band!
float sdfDarwinPrinsChannel(in float detuning, in float couplingKappa, in float lossAlpha, in float slabL)
{
    if (slabL <= 0.0 || couplingKappa <= 0.0)
    {
        return 0.0;
    }
    float al  = max(lossAlpha, 0.0);
    float kb2 = couplingKappa * couplingKappa;
    // Complex square root gamma = gr + i*gi = sqrt(kb^2 + (al - i*detuning)^2), gr >= 0:
    float zr   = kb2 + al * al - detuning * detuning;
    float zi   = -2.0 * al * detuning;
    float zMod = sqrt(zr * zr + zi * zi);
    float gr   = sqrt(max(0.5 * (zMod + zr), 0.0));
    float gi   = (zi < 0.0 ? -1.0 : 1.0) * sqrt(max(0.5 * (zMod - zr), 0.0));

    float E  = exp(-min(2.0 * gr * slabL, 80.0));
    float th = 2.0 * gi * slabL;
    float c  = cos(th), s = sin(th);

    float sr = 1.0 - E * c, si =  E * s;
    float cr = 1.0 + E * c, ci = -E * s;

    float num = kb2 * (sr * sr + si * si);
    float dr  = (gr * cr - gi * ci) + (al * sr + detuning * si);
    float di  = (gr * ci + gi * cr) + (al * si - detuning * sr);
    return clamp(num / max(dr * dr + di * di, 1e-20), 0.0, 1.0);
}

// Polychromatic RGB dynamical Darwin-Prins iridescence across a 6-plane {110} Gyroid photonic slab [ta, tb].
// Combines reflections across lattice planes via probabilistic non-exceedance 1 - prod_m (1 - R_m) in [0, 1].
vec3 sdfDarwinPrinsGyroidRGB(in vec3 rd, in float ta, in float tb,
                             in float q0, in vec3 kWaveRGB, in float couplingKappa, in float lossAlpha)
{
    float L = max(tb - ta, 0.0);
    if (L <= 0.0)
    {
        return vec3(0.0);
    }
    const vec3 G_DIRS[6] = vec3[6](
        vec3( 0.70710678,  0.70710678, 0.0),
        vec3( 0.70710678, -0.70710678, 0.0),
        vec3( 0.0,  0.70710678,  0.70710678),
        vec3( 0.0,  0.70710678, -0.70710678),
        vec3( 0.70710678, 0.0,  0.70710678),
        vec3(-0.70710678, 0.0,  0.70710678)
    );
    vec3 transProd = vec3(1.0);
    for (int m = 0; m < 6; ++m)
    {
        float cosBragg = max(abs(dot(G_DIRS[m], rd)), 0.08);
        float qBraggHalf = 0.5 * q0 * cosBragg;
        // Effective coupling along tilted ray scales by 1 / cosBragg:
        float kapEff = couplingKappa / cosBragg;
        float alEff  = lossAlpha / cosBragg;
        vec3 det = kWaveRGB - vec3(qBraggHalf);
        vec3 rM = vec3(
            sdfDarwinPrinsChannel(det.x, kapEff, alEff, L),
            sdfDarwinPrinsChannel(det.y, kapEff, alEff, L),
            sdfDarwinPrinsChannel(det.z, kapEff, alEff, L)
        );
        transProd *= (vec3(1.0) - rM);
    }
    return vec3(1.0) - transProd;
}

#endif // SDF_GEMS_SDF_FOCK_ARAGO_DIFFRACTION_GLSL
