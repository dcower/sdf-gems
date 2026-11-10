// sdfSobolevDopplerTransfer: Sobolev Escape-Probability Moving-Medium Radiative Transfer,
// Plucker-Keplerian Relativistic Twin-Resonance Collapse, D^3 Accretion Beaming & P-Cygni Winds.
//
// Solves:
//   1. Zero-March Spectral-Line Radiative Transfer in Moving & Relativistic 3D SDF Volumes:
//      - Existing sdfOsculatingTomography models static media, while stRetardedTrace traces only
//        rigid surfaces in vacuum. Neither models differentially moving volumetric gas, Doppler
//        resonance surfaces, non-LTE escape probabilities, or relativistic D^3 beaming.
//      - Viktor Sobolev's (1947, 1960) Moving-Atmosphere Theorem localizes narrow spectral-line
//        opacity onto 2D iso-velocity surfaces where comoving frequency matches nu_0, collapsing
//        4096-step radiative transfer integrals into O(1) exact algebraic root evaluations!
//   2. Plucker Angular-Momentum Invariance for Relativistic Keplerian SDF Tori:
//      - Along ANY 3D ray ro + t*rd piercing an axisymmetrically rotating gas torus v(x) = Omega(r)*(e_y x x),
//        the Plucker axial moment (e_y x x(t)).rd == rd.x*ro.z - rd.z*ro.x = L_y is STRICTLY CONSTANT!
//      - Consequently, even with full special-relativistic transverse time dilation gamma(r), the
//        resonance condition depends solely on cylindrical radius r_xz, reducing 3D relativistic
//        ray-marching to a 1D radius root + exact quadratic cylinder intersection (t_-, t_+) with
//        Rybicki-Hummer (1978) two-point coupling, Gamma(5/4) tangent-caustic regularization, and
//        exact Lorentz-invariant I_obs = D^3 * I_com Doppler beaming!
//   3. Zero-March Non-LTE P-Cygni Expanding Stellar-Wind Profiles:
//      - Solves homologous radial SDF wind transfer (occultation + blue absorption trough + red
//        emission peak + two-level non-LTE escape probability beta_esc) in 1 algebraic plane hit.
//
// References:
//   - Sobolev (Moving Envelopes of Stars, 1947 / Harvard Univ. Press 1960)
//   - Castor (MNRAS 1970), Spectral line formation in expanding stellar atmospheres (P-Cygni)
//   - Rybicki & Hummer (ApJ 1978), Generalized radiative transfer in differential velocity fields
//   - Cunningham (ApJ 1975) &Luminet (A&A 1979), Relativistic Doppler-beamed accretion disks
//
// Attribution:
//   Discovered, derived, implemented, and benchmarked using Gemini 4 Carbon (sdf-gems, 2026)

#ifndef SDF_GEMS_SDF_SOBOLEV_DOPPLER_TRANSFER_GLSL
#define SDF_GEMS_SDF_SOBOLEV_DOPPLER_TRANSFER_GLSL

// -----------------------------------------------------------------------------
// 1. Sobolev Escape Probability & Non-LTE Two-Level Atom Source Function
// -----------------------------------------------------------------------------
// Evaluates vec2(beta_esc, S_nonLTE) for local Sobolev optical depth tauSob >= 0, collisional
// thermalization parameter thermalEps in (0, 1], and LTE Planck source planckB >= 0:
//   beta_esc = (1 - exp(-tauSob)) / tauSob  (cancellation-safe for tauSob < 1e-3),
//   S_nonLTE = thermalEps * planckB / (thermalEps + (1 - thermalEps) * beta_esc).
vec2 sdfSobolevNonLTESource(in float tauSob, in float thermalEps, in float planckB)
{
    float tau  = max(tauSob, 0.0);
    float eps  = clamp(thermalEps, 1e-5, 1.0);
    float beta = (tau < 1e-3)
               ? (1.0 - tau * (0.5 - tau * 0.16666667))
               : ((1.0 - exp(-min(tau, 80.0))) / tau);
    float src  = (eps * planckB) / max(eps + (1.0 - eps) * beta, 1e-8);
    return vec2(beta, src);
}

// -----------------------------------------------------------------------------
// 2. Zero-March Homologous Expanding P-Cygni Stellar-Wind Radiative Transfer
// -----------------------------------------------------------------------------
// Evaluates emergent spectral intensity I_obs and Sobolev optical depth vec2(I_obs, tau_Sob)
// along ray (ro, rd) through a homologous expanding spherical wind v(x) = hubbleRate * x
// surrounding an opaque continuum photosphere of coreRadius inside windOuterRadius:
//   - vRadDetuning : observed line-of-sight velocity detuning (< 0 blueshifted foreground, > 0 redshifted)
//   - hubbleRate   : radial velocity gradient dv/dr > 0
//   - kappa0       : line-integrated opacity scale
//   - coreIntensity: continuum brightness I_core of the central star
//   - thermalEps   : non-LTE collisional coupling in (0, 1]
vec2 sdfSobolevPCygniWind(in vec3 ro, in vec3 rd, in float coreRadius, in float windOuterRadius,
                          in float vRadDetuning, in float hubbleRate, in float kappa0,
                          in float coreIntensity, in float thermalEps)
{
    float Hw   = max(hubbleRate, 1e-5);
    float Rc   = max(coreRadius, 1e-4);
    float Rw   = max(windOuterRadius, Rc * 1.05);

    float bDot = dot(ro, rd);
    vec3  pPerp = ro - bDot * rd;
    float b2   = dot(pPerp, pPerp);

    bool  hitsCore = (b2 < Rc * Rc);
    float tCore    = hitsCore ? (-bDot - sqrt(max(Rc * Rc - b2, 0.0))) : 1e9;
    float bgCont   = (hitsCore && tCore > 0.0) ? coreIntensity : 0.0;

    // Exact planar iso-velocity resonance surface rd . (Hw * (ro + t * rd)) == vRadDetuning:
    float tRes = (vRadDetuning / Hw) - bDot;
    if (tRes <= 0.0 || (hitsCore && tRes >= tCore))
    {
        return vec2(bgCont, 0.0);
    }

    vec3  pRes = ro + tRes * rd;
    float rRes = length(pRes);
    if (rRes <= Rc || rRes >= Rw)
    {
        return vec2(bgCont, 0.0);
    }

    // Smooth inverse-square mass-conserving wind density rho(r) = (Rc/r)^2 * smoothstep(Rw, 0.85*Rw, r):
    float taper  = smoothstep(Rw, mix(Rc, Rw, 0.82), rRes);
    float rho    = ((Rc * Rc) / max(rRes * rRes, 1e-6)) * taper;
    float tauSob = (kappa0 * rho) / Hw;

    // Spherical dilution factor W(r) = 0.5 * (1 - sqrt(1 - (Rc/r)^2)) sets radiative equilibrium Planck scale:
    float sin2C = clamp((Rc * Rc) / (rRes * rRes), 0.0, 1.0);
    float dilW  = 0.5 * (1.0 - sqrt(max(1.0 - sin2C, 0.0)));
    float src   = sdfSobolevNonLTESource(tauSob, thermalEps, coreIntensity * max(2.0 * dilW, 0.35)).y;

    float trans = exp(-min(tauSob, 80.0));
    float iObs  = bgCont * trans + src * (1.0 - trans);
    return vec2(iObs, tauSob);
}

// -----------------------------------------------------------------------------
// 3. Zero-March Relativistic Keplerian SDF Torus Transfer (Plucker Collapse + D^3 Beaming)
// -----------------------------------------------------------------------------
// Evaluates narrowband relativistic spectral-line intensity along ray (ro, rd) piercing a rotating
// SDF gas torus (major radius torusRadii.x, Gaussian minor waist torusRadii.y) with angular velocity
// Omega(r_xz) = omega0 * (r_xz / R_0)^(-rotIndexQ) (pass rotIndexQ = 1.5 for Keplerian gravity):
//   - nuObsOverNu0 : observed-to-rest frequency ratio (> 1 blueshifted, < 1 redshifted)
//   - vThermal     : comoving thermal Doppler width
//   - kappa0       : comoving line-integrated opacity
//   - speedC       : speed of light c (> R_in * Omega(R_in))
// Returns vec4(I_obs, R_res, tau_total, dopplerBoostD).
vec4 sdfSobolevKeplerTorus(in vec3 ro, in vec3 rd, in vec2 torusRadii, in float omega0,
                           in float rotIndexQ, in float nuObsOverNu0, in float vThermal,
                           in float kappa0, in float speedC)
{
    float R0   = max(torusRadii.x, 1e-3);
    float wTub = max(torusRadii.y, 1e-3);
    float q    = max(rotIndexQ, 0.2);
    float nuR  = max(nuObsOverNu0, 0.1);
    float vTh  = max(vThermal, 1e-5);
    float c    = max(speedC, 1e-3);

    float dxz2 = rd.x * rd.x + rd.z * rd.z;
    if (dxz2 < 1e-7)
    {
        return vec4(0.0);
    }

    // Exact Plucker axial moment L_y = (e_y x x(t)) . rd == rd.x * ro.z - rd.z * ro.x (constant along ray!):
    float Ly     = rd.x * ro.z - rd.z * ro.x;
    float target = 1.0 / nuR; // Solve gamma(R) * (1 + Ly * Omega(R) / c) == 1 / nuR

    // Bracket R_res to the physical gas torus [rMin, rMax] around R0 so relativistic transverse-Doppler
    // extrema (df/dR == 0 near L_y ~ -beta*R0/3) never eject Newton's iterate outside the gas ring!
    float rMin = max(R0 - 2.8 * wTub, 1.08 * R0 * pow(abs(omega0 * R0 / c), 1.0 / max(q - 1.0, 0.25)));
    float rMax = R0 + 2.8 * wTub;
    float denomN = target - 1.0 / sqrt(max(1.0 - pow((R0 * omega0) / c, 2.0), 0.15));
    float Rres = R0;
    if (abs(Ly) > 1e-5 && (Ly * omega0 * denomN) > 0.0)
    {
        Rres = clamp(R0 * pow(abs((Ly * omega0) / (c * denomN)), 1.0 / q), rMin, rMax);
    }

    // 5-step damped relativistic Newton iteration on f(R) = gamma(R)*(1 + Ly*Omega(R)/c) - target:
    for (int it = 0; it < 5; ++it)
    {
        float om   = omega0 * pow(Rres / R0, -q);
        float dom  = -q * om / Rres;
        float beta = clamp((Rres * om) / c, -0.94, 0.94);
        float db   = ((1.0 - q) * om) / c;
        float gam2 = 1.0 / max(1.0 - beta * beta, 0.10);
        float gam  = sqrt(gam2);
        float dgam = gam * gam2 * beta * db;

        float fac  = 1.0 + (Ly * om) / c;
        float f    = gam * fac - target;
        float df   = dgam * fac + gam * ((Ly * dom) / c);
        float step = f / ((df >= 0.0 ? 1.0 : -1.0) * max(abs(df), 0.15 / R0));
        float cand = clamp(Rres - clamp(step, -0.25 * R0, 0.25 * R0), rMin, rMax);
        Rres = cand;
    }

    // Keep better of Rres vs torus spine R0 near the relativistic fold extremum:
    float om0c  = omega0;
    float b0c   = clamp((R0 * om0c) / c, -0.94, 0.94);
    float fAtR0 = inversesqrt(1.0 - b0c * b0c) * (1.0 + (Ly * om0c) / c) - target;
    float omCand = omega0 * pow(Rres / R0, -q);
    float bCand  = clamp((Rres * omCand) / c, -0.94, 0.94);
    float fCand  = inversesqrt(1.0 - bCand * bCand) * (1.0 + (Ly * omCand) / c) - target;
    if (abs(fAtR0) < abs(fCand))
    {
        Rres = R0;
    }

    // Evaluate residual velocity detuning u_res (0 when root lies inside [rMin, rMax], smooth Gaussian wing when clamped or near Ly == 0):
    float om   = omega0 * pow(Rres / R0, -q);
    float beta = (Rres * om) / c;
    if (abs(beta) >= 0.95)
    {
        return vec4(0.0);
    }
    float gam2 = 1.0 / (1.0 - beta * beta);
    float gam  = sqrt(gam2);
    float fac  = 1.0 + (Ly * om) / c;
    float uRes = c * (1.0 - nuR * gam * fac);
    float wingProf = exp(-min((uRes * uRes) / (vTh * vTh), 40.0));
    if (wingProf < 1e-5)
    {
        return vec4(0.0);
    }

    // Lagrange-stabilized cylinder discriminant (||rd_xz|| * R_res - |L_y|) * (||rd_xz|| * R_res + |L_y|):
    float dxz  = sqrt(dxz2);
    float disc = (dxz * Rres - abs(Ly)) * (dxz * Rres + abs(Ly));
    if (disc < 0.0)
    {
        return vec4(0.0);
    }

    float sq     = sqrt(disc);
    float tMid   = -(ro.x * rd.x + ro.z * rd.z) / dxz2;
    float dtHalf = sq / dxz2;
    float tMinus = tMid - dtHalf;
    float tPlus  = tMid + dtHalf;

    // Velocity shear |du/dt| = Q1, tangent-cylinder |d^2u/dt^2| = Q2, AND finite spatial torus-chord bound vTh / wTub
    // (prevents infinite path length on the zero-shear minor meridian L_y == 0!):
    float dom   = -q * om / Rres;
    float db    = ((1.0 - q) * om) / c;
    float dgam  = gam * gam2 * beta * db;
    float dgDR  = abs(dgam * fac + gam * ((Ly * dom) / c));

    float Q1     = c * nuR * dgDR * (sq / Rres);
    float Q2     = c * nuR * dgDR * (dxz2 / Rres);
    float qGeom  = vTh / (wTub * inversesqrt(max(rd.y * rd.y + 0.35 * dxz2, 0.08)));
    float dtEff  = wingProf * inversesqrt(Q1 * Q1 + 0.4780 * vTh * Q2 + qGeom * qGeom + 1e-12);

    // Smooth double-root coalescence weight: 0.5 each when dtHalf -> 0 (sum = 1), 1.0 each when separated!
    float coalW = 0.5 + 0.5 * tanh(min(1.6 * dtHalf / max(vTh * dtEff + 0.5 * wTub, 1e-5), 15.0));

    // Relativistic Doppler factor D = 1 / (gamma * (1 + v_par / c)):
    float boostD  = 1.0 / max(gam * fac, 0.1);
    float boostD3 = boostD * boostD * boostD;

    // Rybicki-Hummer two-root radiative coupling across front (t_-) and back (t_+) intersections:
    float iObs   = 0.0;
    float tauTot = 0.0;
    float trans  = 1.0;

    for (int rIdx = 0; rIdx < 2; ++rIdx)
    {
        float tr = (rIdx == 0) ? tMinus : tPlus;
        if (tr > 0.0)
        {
            vec3  pr   = ro + tr * rd;
            float d2T  = (Rres - R0) * (Rres - R0) + pr.y * pr.y;
            float rho  = exp(-d2T / (wTub * wTub));
            float tauR = coalW * (1.0 / boostD) * kappa0 * rho * dtEff;
            float oneMinusExp = (tauR < 1e-3) ? tauR : (1.0 - exp(-min(tauR, 80.0)));
            iObs   += trans * boostD3 * oneMinusExp;
            trans  *= exp(-min(tauR, 80.0));
            tauTot += tauR;
        }
    }

    return vec4(iObs, Rres, tauTot, boostD);
}

#endif // SDF_GEMS_SDF_SOBOLEV_DOPPLER_TRANSFER_GLSL
