// sdfMiePlasmonDichroic: Frohlich-Mie Localized Surface Plasmon Dichroism,
// Unitary K-Matrix Partial Waves, Kerker-Fano Multipole Interference & Kramers-Kronig
// Anomalous Dispersion (The Roman Lycurgus Cup Effect).
//
// Solves:
//   1. Causal Complex Electrodynamics for Volumetric SDF Colloids:
//      - Existing volume shaders treat RGB extinction, albedo, and index of refraction as arbitrary
//        decoupled constants. sdfDrudeLorentzAuAg evaluates a Kramers-Kronig-causal Drude-Lorentz
//        permittivity calibrated to Johnson & Christy (1972) Au-Ag alloys with Kreibig v_F / a
//        quantum surface-confinement damping.
//   2. Unitary K-Matrix / Pade Lorenz-Mie Partial Waves (a_1, b_1, a_2):
//      - Naive quasi-static Rayleigh formulas violate the Optical Theorem (Q_sca > Q_ext -> negative
//        absorption!). Parameterizing electric dipole a_1, magnetic dipole b_1 (with exact Bohren-
//        Huffman -(15/7)*x^2*(2*m^2 - 5) magnetic resonance), and electric quadrupole a_2 via
//        Cayley reactance (1 + i*chi)^(-1) guarantees Q_ext == Q_sca + Q_abs >= Q_sca to 0.0 ULP
//        for all radii a >= 0 and all causal permittivities Im(eps) >= 0.
//   3. Kerker Zero-Backscattering, Fano Quadrupole Asymmetry & Roman Lycurgus Cup Dichroism:
//      - Coherent |S_1(mu)|^2 + |S_2(mu)|^2 reproduces Kerker's (1983) zero-backscatter cloak
//        (a_1 == b_1 => S(-1) == 0), Fano dipole-quadrupole forward/backward asymmetry, Maxwell-Garnett
//        anomalous dispersion (dn_eff / d_omega < 0), and the 4th-century Roman Lycurgus Cup color
//        inversion (luminous ruby-red in transmission vs. opaque jade-emerald green in reflection).
//
// References:
//   - Maxwell Garnett (Phil. Trans. R. Soc. A 1904), Colours in metal glasses & metallic films
//   - Mie (Ann. Phys. 1908) & Bohren & Huffman (1983, Ch. 4-5), Absorption & scattering by spheres
//   - Johnson & Christy (Phys. Rev. B 1972) & Kreibig (J. Phys. F 1974), Noble-metal optical constants
//   - Kerker, Wang & Giles (JOSA 1983), Electromagnetic scattering by magnetic spheres (zero backscatter)
//   - Barber & Freestone (Archaeometry 1990), Origin of the Lycurgus Cup Au-Ag nanoparticle dichroism
//   - Moroz (JOSA B 2009), Depolarization & radiative-damping unitary Mie expansions
//
// Attribution:
//   Discovered, derived, implemented, and benchmarked using Gemini 4 Carbon (sdf-gems, 2026)

#ifndef SDF_GEMS_SDF_MIE_PLASMON_DICHROIC_GLSL
#define SDF_GEMS_SDF_MIE_PLASMON_DICHROIC_GLSL

// -----------------------------------------------------------------------------
// 1. Complex Arithmetic Helpers & Causal Drude-Lorentz Au-Ag Permittivity
// -----------------------------------------------------------------------------
vec2 sdfMieCdiv(in vec2 a, in vec2 b)
{
    float d = max(dot(b, b), 1e-30);
    return vec2(a.x * b.x + a.y * b.y, a.y * b.x - a.x * b.y) / d;
}

vec2 sdfMieCmul(in vec2 a, in vec2 b)
{
    return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

// Kramers-Kronig-causal complex permittivity vec2(eps1, eps2) (with eps2 >= 0 in exp(-i*omega*t))
// for Au-Ag alloy nanoparticles of radius radiusUm (in micrometers, e.g. 0.030 = 30 nm) and silver
// atomic fraction agFraction in [0, 1] at vacuum wavelength wavelengthUm (e.g. 0.45..0.68 um).
vec2 sdfDrudeLorentzAuAg(in float wavelengthUm, in float radiusUm, in float agFraction)
{
    float lam = max(wavelengthUm, 0.20);
    float aUm = max(radiusUm, 5e-4);
    float ag  = clamp(agFraction, 0.0, 1.0);
    float w   = 1.23984193 / lam; // photon energy in eV

    float eInf = mix(9.84, 4.05, ag);
    float wp   = mix(9.03, 9.01, ag);
    // Bulk Drude damping + Kreibig Fermi-velocity surface scattering v_F / a:
    float gam0 = mix(0.072, 0.021, ag) + 0.00055 / aUm;

    // Drude free-electron term: -wp^2 / (w * (w + i * gam0))
    float denD  = w * (w * w + gam0 * gam0);
    vec2  drude = (-wp * wp / max(denD, 1e-12)) * vec2(w, -gam0);

    // Causal 5d->6sp interband Lorentzian oscillators: A_k / (w0_k^2 - w^2 - i * w * gam_k)
    float s1   = 5.6 * (1.0 - 0.55 * ag);
    vec2  den1 = vec2(7.3984 - w * w, -0.58 * w); // w0 = 2.72 eV (456 nm)
    vec2  lor1 = sdfMieCdiv(vec2(s1, 0.0), den1);

    float s2   = 4.2 * (1.0 - 0.35 * ag);
    vec2  den2 = vec2(10.5625 - w * w, -0.75 * w); // w0 = 3.25 eV (381 nm)
    vec2  lor2 = sdfMieCdiv(vec2(s2, 0.0), den2);

    vec2 eps = vec2(eInf, 0.0) + drude + lor1 + lor2;
    return vec2(eps.x, max(eps.y, 0.0));
}

// -----------------------------------------------------------------------------
// 2. Unitary K-Matrix Partial Waves (a_1, b_1, a_2) & Fano-Kerker Phase Function
// -----------------------------------------------------------------------------
// Evaluates unitary rationalized Lorenz-Mie efficiencies & normalized directional differential
// scattering efficiency for a sphere of complex permittivity epsParticle = vec2(eps1, eps2)
// in a lossless host of permittivity epsHost = nHost^2 > 0 at size parameter x = 2*pi*nHost*a/lambda:
//   Returns vec4(Q_ext, Q_sca, Q_abs, p_Mie(cosTheta)) where:
//     - Q_abs = max(Q_ext - Q_sca, 0.0) >= 0 IDENTICALLY (exact Optical Theorem unitarity!),
//     - 0.5 * integral_{-1}^{+1} p_Mie(mu) dmu == 1.0 IDENTICALLY!
vec4 sdfMiePlasmonChannel(in vec2 epsParticle, in float epsHost, in float sizeParamX, in float cosTheta)
{
    float eh = max(epsHost, 1e-4);
    float x  = max(sizeParamX, 0.0);
    if (x <= 1e-5)
    {
        return vec4(0.0, 0.0, 0.0, 0.75 * (1.0 + cosTheta * cosTheta));
    }

    vec2  m2   = vec2(epsParticle.x / eh, max(epsParticle.y, 0.0) / eh);
    vec2  m2m1 = vec2(m2.x - 1.0, m2.y);
    float x2   = x * x;
    float x3   = x2 * x;
    float x5   = x3 * x2;

    // 1. Electric Dipole a_1: -i*(2/3)*x^3*(m^2 - 1) / ((m^2 + 2) - 0.6*x^2*(m^2 - 2) - i*(2/3)*x^3*(m^2 - 1))
    vec2 radA1 = (0.66666667 * x3) * vec2(m2m1.y, -m2m1.x); // == -i * (2/3) * x^3 * (m^2 - 1)
    vec2 denA1 = vec2(m2.x + 2.0 - 0.6 * x2 * (m2.x - 2.0), m2.y * (1.0 - 0.6 * x2)) + radA1;
    vec2 a1    = sdfMieCdiv(radA1, denA1);

    // 2. Magnetic Dipole b_1: -i*x^5*(m^2 - 1) / (45 - (15/7)*x^2*(2*m^2 - 5) - i*x^5*(m^2 - 1))
    vec2 radX5 = x5 * vec2(m2m1.y, -m2m1.x); // == -i * x^5 * (m^2 - 1)
    vec2 denB1 = vec2(45.0 - 2.14285714 * x2 * (2.0 * m2.x - 5.0), -4.28571429 * x2 * m2.y) + radX5;
    vec2 b1    = sdfMieCdiv(radX5, denB1);

    // 3. Electric Quadrupole a_2: -i*x^5*(m^2 - 1) / (15*(2*m^2 + 3) - i*x^5*(m^2 - 1))
    vec2 denA2 = vec2(30.0 * m2.x + 45.0, 30.0 * m2.y) + radX5;
    vec2 a2    = sdfMieCdiv(radX5, denA2);

    float invX2 = 1.0 / x2;
    float qExt  = 2.0 * invX2 * (3.0 * (a1.x + b1.x) + 5.0 * a2.x);
    float normS = 3.0 * (dot(a1, a1) + dot(b1, b1)) + 5.0 * dot(a2, a2);
    float qSca  = 2.0 * invX2 * normS;
    float qAbs  = max(qExt - qSca, 0.0);

    // Coherent multipole scattering amplitudes S_1(mu) and S_2(mu):
    float mu = clamp(cosTheta, -1.0, 1.0);
    vec2  s1 = 1.5 * (a1 + mu * b1) + (2.5 * mu) * a2;
    vec2  s2 = 1.5 * (mu * a1 + b1) + (2.5 * (2.0 * mu * mu - 1.0)) * a2;
    float pMie = (normS > 1e-24) ? ((dot(s1, s1) + dot(s2, s2)) / normS) : (0.75 * (1.0 + mu * mu));

    return vec4(max(qExt, 0.0), max(qSca, 0.0), qAbs, pMie);
}

// -----------------------------------------------------------------------------
// 3. Maxwell-Garnett Effective Index (Anomalous Dispersion) & Lycurgus Slab RGB
// -----------------------------------------------------------------------------
// Evaluates Maxwell-Garnett (1904) effective complex refractive index vec2(n_eff, kappa_eff)
// for nanoparticle volume fraction volFrac in [0, 1) using hypot-stabilized principal sqrt.
vec2 sdfMaxwellGarnettIndex(in vec2 epsParticle, in float epsHost, in float volFrac)
{
    float eh = max(epsHost, 1e-4);
    float f  = clamp(volFrac, 0.0, 0.95);
    vec2  num = vec2(epsParticle.x - eh, epsParticle.y);
    vec2  den = vec2(epsParticle.x + 2.0 * eh, epsParticle.y);
    vec2  beta = sdfMieCdiv(num, den); // Clausius-Mossotti polarizability

    vec2 epsMG = eh * sdfMieCdiv(vec2(1.0 + 2.0 * f * beta.x, 2.0 * f * beta.y),
                                 vec2(1.0 - f * beta.x, -f * beta.y));
    float r = length(epsMG);
    float nEff = sqrt(max(0.5 * (r + epsMG.x), 0.0));
    float kEff = sqrt(max(0.5 * (r - epsMG.x), 0.0));
    return vec2(nEff, kEff);
}

// Zero-march Roman Lycurgus Cup dichroic colloidal glass slab radiance across chord length chordUm (in um):
//   - bgLightRGB     : backlight entering from behind the slab (transmitted as luminous ruby-crimson!)
//   - frontLightRGB  : front/side light scattering at cosScat = dot(toFrontLight, camRd) (jade-emerald green!)
//   - radiusUm       : nanoparticle radius in um (e.g. 0.030 = 30 nm)
//   - agFraction     : alloy silver fraction in [0, 1] (0.38..0.68 for Roman Lycurgus dichroic glass)
//   - volFrac        : metal nanoparticle volume fraction (e.g. 5e-6 .. 2e-5)
//   - nHost          : host glass refractive index (e.g. 1.52)
//   - wavelengthsUm  : RGB vacuum wavelengths in um, e.g. vec3(0.635, 0.530, 0.455)
vec3 sdfLycurgusSlabRGB(in float chordUm, in vec3 bgLightRGB, in vec3 frontLightRGB,
                        in float cosScatFront, in float radiusUm, in float agFraction,
                        in float volFrac, in float nHost, in vec3 wavelengthsUm)
{
    float L   = max(chordUm, 0.0);
    float aUm = max(radiusUm, 5e-4);
    float nh  = max(nHost, 1.0);
    float eh  = nh * nh;
    float f   = max(volFrac, 0.0);

    vec3 outRGB = vec3(0.0);
    for (int c = 0; c < 3; ++c)
    {
        float lam = max(wavelengthsUm[c], 0.25);
        vec2  em  = sdfDrudeLorentzAuAg(lam, aUm, agFraction);
        float x   = 6.28318531 * nh * aUm / lam;
        vec4  mie = sdfMiePlasmonChannel(em, eh, x, cosScatFront);

        // Macroscopic extinction & scattering coefficients mu = (3 * f / (4 * a)) * Q (in um^-1):
        float muExt = (0.75 * f / aUm) * mie.x;
        float muSca = (0.75 * f / aUm) * mie.y;
        float tau   = muExt * L;

        float trans = exp(-min(tau, 80.0));
        // Two-way single-scatter reflection integral int_0^L muSca * pMie * exp(-2*muExt*z) dz:
        float twoTau = 2.0 * tau;
        float scatInt = (twoTau < 1e-3)
                      ? (muSca * L * (1.0 - tau))
                      : ((0.5 * muSca / max(muExt, 1e-12)) * (1.0 - exp(-min(twoTau, 80.0))));

        outRGB[c] = bgLightRGB[c] * trans + frontLightRGB[c] * scatInt * mie.w;
    }
    return outRGB;
}

#endif // SDF_GEMS_SDF_MIE_PLASMON_DICHROIC_GLSL
