// sdfBeltramiPhotoelastic: Beltrami-Eikonal Zero-FEM 3D Residual Stress Tensors,
// Non-Abelian SU(2) Magnus-Rabi Twisting Polarimetry & Aben Characteristic Optical Rotation.
//
// Solves:
//   1. Zero-FEM equilibrated 3D residual stress inside arbitrary smooth SDF solids (d <= 0):
//      By Beltrami's (1892) curl-curl representation sigma_B = (Lap psi) I - nabla^2 psi with
//      quench potential psi(x) = f(d(x)), f'(d) = -sigmaSkin * d * exp(d / quenchDepth):
//      - Static equilibrium div(sigma_B) == 0 holds identically across R^3 without a 3D mesh.
//      - Boundary traction sigma_B * n |_{d=0} == 0 vanishes identically on every smooth SDF
//        via the Eikonal identity H * grad(d) = 0 and f'(0) = 0, self-equilibrating a biaxial
//        compressive outer skin (|d| < quenchDepth) against a tensile core (Prince Rupert's drop).
//   2. Fixes two hard polarimetry bugs in sdfPauliBirefringence (glsl/volume/sdfPauliSpinorSoliton.glsl):
//      - Replaces double-squared Bloch azimuth (which produced an unphysical 8-arm sin^2(4*theta)
//        cross) with exact SU(2) Stokes-Jones projection (4-arm Brewster sin^2(2*theta) cross).
//      - Restores circular/Faraday magneto-optic rotation gamma_3 (which sdfPauliBirefringence
//        dropped to literal 0.0) and separates Plane vs Circular Polariscope readouts.
//   3. Solves Aben's (1966) Non-Abelian Torsion Paradox via 4th-order 2-node Gauss-Legendre
//      Magnus SU(2) propagation (Omega_4 = kh/2 (g_- + g_+) - sqrt(3)/6 k^2 h^2 (g_- x g_+))
//      and exact Gooch-Tarry / Rabi helical propagators: twisting principal stress axes generate
//      non-Abelian Pancharatnam-Berry optical rotation along the circular pole where Abelian ray
//      integrals cancel to zero.
//
// References:
//   - Beltrami (1892), "Osservazioni sulla nota precedente" (Curl-curl stress tensor potentials)
//   - Wertheim (1854), "Memoire sur la double refraction temporairement produite dans les corps isotropes"
//   - Aben (1966, 1979), "Optical phenomena in photoelastic models by the rotation of principal axes"
//     & "Integrated Photoelasticity" (Characteristic directions & optical rotation)
//   - Gooch & Tarry (J. Phys. D 1975), "The optical properties of twisted nematic liquid crystal structures"
//   - Iserles & Norsett (Phil. Trans. R. Soc. A 1999), "On the solution of linear differential equations in Lie groups"
//   - Sharafutdinov (1994), "Integral Geometry of Tensor Fields" (Transverse ray transform)
//
// Attribution:
//   Discovered, derived, implemented, and benchmarked using Gemini 4 Carbon (sdf-gems, 2026)

#ifndef SDF_GEMS_SDF_BELTRAMI_PHOTOELASTIC_GLSL
#define SDF_GEMS_SDF_BELTRAMI_PHOTOELASTIC_GLSL

// -----------------------------------------------------------------------------
// 1. Beltrami-Eikonal Traction-Free 3D Stress Projected to Transverse Ray Frame
// -----------------------------------------------------------------------------
// Given interior SDF value d <= 0, unit gradient g = nabla d, Hessian H = nabla^2 d,
// unit ray direction rd, unit transverse horizontal polarizer axis e1 (perp to rd),
// thermal quench skin stress sigmaSkin, quench skin thickness quenchDepth > 0,
// optional torsion axis twistAxis & twist shear tauTwist, and Faraday term faradayRate:
// Projects the divergence-free, traction-free 3D stress tensor onto (e1, e2 = cross(rd, e1)),
// returning the SU(2) Bloch birefringence generator vec3(gamma1_0_90, gamma2_pm45, gamma3_faraday).
vec3 sdfBeltramiTransverseStress(in float d, in vec3 g, in mat3 H,
                                 in vec3 rd, in vec3 e1,
                                 in float sigmaSkin, in float quenchDepth,
                                 in vec3 twistAxis, in float tauTwist, in float faradayRate)
{
    if (d >= 0.0)
    {
        return vec3(0.0, 0.0, faradayRate);
    }
    vec3  e2  = cross(rd, e1);
    float lq  = max(quenchDepth, 1e-4);
    float u   = max(d / lq, -80.0); // d <= 0 -> u <= 0 (zero overflow)
    float ex  = exp(u);

    // Quench potential derivatives: f'(d) = -sigmaSkin * d * exp(d/lq), f''(d) = -sigmaSkin * (1 + d/lq) * exp(d/lq)
    float fp  = -sigmaSkin * d * ex;
    float fpp = -sigmaSkin * (1.0 + u) * ex;

    float g1  = dot(g, e1);
    float g2  = dot(g, e2);
    vec3  He1 = H * e1;
    vec3  He2 = H * e2;
    float h11 = dot(e1, He1);
    float h22 = dot(e2, He2);
    float h12 = dot(e1, He2);

    // Isotropic trace (Lap psi) I_3 cancels identically in (sigma_11 - sigma_22) and sigma_12:
    float gamma1 = -0.5 * fpp * (g1 * g1 - g2 * g2) - 0.5 * fp * (h11 - h22);
    float gamma2 = -fpp * g1 * g2 - fp * h12;

    // Boundary-vanishing level-set-tangent torsion shear (a (x) s_tw + s_tw (x) a), s_tw = -tau * d * e^(d/lq) * (g x a):
    if (abs(tauTwist) > 1e-7)
    {
        vec3  sTw = (-tauTwist * d * ex) * cross(g, twistAxis);
        float a1  = dot(twistAxis, e1);
        float a2  = dot(twistAxis, e2);
        float s1  = dot(sTw, e1);
        float s2  = dot(sTw, e2);
        gamma1 += a1 * s1 - a2 * s2;
        gamma2 += a1 * s2 + a2 * s1;
    }

    return vec3(gamma1, gamma2, faradayRate);
}

// -----------------------------------------------------------------------------
// 2. Unitary SU(2) Exponential, Composition, Rabi Helical Slab & 4th-Order Magnus
// -----------------------------------------------------------------------------
// Maps generator vector omega = vec3(theta1, theta2, theta3) to unit SU(2) quaternion
// q = (w, x, y, z) = (cos ||omega||, sinc(||omega||) * omega), isomorphic to
// U = w * I_2 - i * (x * sigma_3 + y * sigma_1 + z * sigma_2).
vec4 sdfExpSU2(in vec3 omega)
{
    float n2 = dot(omega, omega);
    if (n2 < 1e-8)
    {
        float sinc = 1.0 - n2 * 0.16666667;
        return normalize(vec4(1.0 - 0.5 * n2, sinc * omega));
    }
    float n = sqrt(n2);
    return vec4(cos(n), (sin(n) / n) * omega);
}

// Exact Hamilton composition qOut = qLate * qEarly corresponding to Jones product U_late * U_early.
vec4 sdfComposeSU2(in vec4 qLate, in vec4 qEarly)
{
    return vec4(
        qLate.x * qEarly.x - dot(qLate.yzw, qEarly.yzw),
        qLate.x * qEarly.yzw + qEarly.x * qLate.yzw + cross(qLate.yzw, qEarly.yzw)
    );
}

// Exact closed-form Gooch-Tarry / Rabi SU(2) propagator for a uniformly twisting & Faraday-active
// birefringent slab of length chordL (principal physical axis theta(s) = startAzimuth + twistRate * s,
// optical wavenumber kWave, linear retardance rate linRetardRate, Faraday rate faradayRate).
vec4 sdfGoochTarryTwistedSU2(in float linRetardRate, in float startAzimuth,
                             in float twistRate, in float faradayRate,
                             in float chordL, in float kWave)
{
    float endAz = startAzimuth + twistRate * chordL;
    vec4 qStart = vec4(cos(startAzimuth), 0.0, 0.0, -sin(startAzimuth));
    vec4 qCore  = sdfExpSU2(chordL * vec3(kWave * linRetardRate, 0.0, kWave * faradayRate - twistRate));
    vec4 qEnd   = vec4(cos(endAz), 0.0, 0.0, sin(endAz));
    return sdfComposeSU2(qEnd, sdfComposeSU2(qCore, qStart));
}

// 4th-order O(h^5) local / O(h^4) global unitary 2-node Gauss-Legendre Magnus SU(2) step
// across slab thickness h from Bloch samples gammaMinus = gamma(t_mid - h/(2*sqrt(3)))
// and gammaPlus = gamma(t_mid + h/(2*sqrt(3))). Preserves ||q|| == 1 identically!
vec4 sdfMagnus4StepSU2(in vec3 gammaMinus, in vec3 gammaPlus, in float h, in float kWave)
{
    float kh  = kWave * h;
    vec3  avg = (0.5 * kh) * (gammaMinus + gammaPlus);
    // Non-Abelian Lie commutator -(sqrt(3)/6) * (k*h)^2 * (gamma_- x gamma_+):
    vec3  com = (-0.2886751346 * kh * kh) * cross(gammaMinus, gammaPlus);
    return sdfExpSU2(avg + com);
}

// -----------------------------------------------------------------------------
// 3. Aben Characteristic Polariscope Readout & Polychromatic RGB Polariscopes
// -----------------------------------------------------------------------------
// Extracts exact physical observables from unit SU(2) propagator q = (w, x, y, z):
//   - .x: Crossed Plane Polariscope transmittance |<V| U |H>|^2 = y^2 + z^2 in [0, 1]
//         (4-arm Brewster sin^2(2*theta) isoclinics + Faraday & non-Abelian twist transmission!)
//   - .y: Dark-Field Circular Polariscope transmittance |<L| U |R>|^2 = x^2 + y^2 in [0, 1]
//         (eliminates both azimuthal isoclinics and Faraday rotation -> pure stress fringes!)
//   - .z: Aben characteristic optical rotation angle atan(z, w) in [-pi, pi]
//   - .w: Characteristic phase retardance half-angle acos(sqrt(w^2 + z^2)) in [0, pi/2]
vec4 sdfAbenPolariscopeReadout(in vec4 qSU2)
{
    float w = qSU2.x, x = qSU2.y, y = qSU2.z, z = qSU2.w;
    float iPlane = clamp(y * y + z * z, 0.0, 1.0);
    float iCirc  = clamp(x * x + y * y, 0.0, 1.0);
    float rotAng = atan(z, w);
    float retard = acos(clamp(sqrt(w * w + z * z), 0.0, 1.0));
    return vec4(iPlane, iCirc, rotAng, retard);
}

// Evaluates polychromatic RGB Plane Polariscope (col 0) & Circular Polariscope (col 1) transmission
// across a chord of total length chordL partitioned into two 4th-order Magnus sub-slabs using 4
// interior Gauss-Legendre stress samples (g0, g1 on [0, L/2] and g2, g3 on [L/2, L]).
mat2x3 sdfBeltramiPolariscopeRGB(in vec3 g0, in vec3 g1, in vec3 g2, in vec3 g3,
                                 in float chordL, in vec3 kWaveRGB)
{
    if (chordL <= 0.0)
    {
        return mat2x3(vec3(0.0), vec3(0.0));
    }
    float h = 0.5 * chordL;
    vec3 planeRGB = vec3(0.0);
    vec3 circRGB  = vec3(0.0);
    for (int c = 0; c < 3; ++c)
    {
        float kw = kWaveRGB[c];
        vec4 qEarly = sdfMagnus4StepSU2(g0, g1, h, kw);
        vec4 qLate  = sdfMagnus4StepSU2(g2, g3, h, kw);
        vec4 qTot   = sdfComposeSU2(qLate, qEarly);
        vec4 ro     = sdfAbenPolariscopeReadout(qTot);
        planeRGB[c] = ro.x;
        circRGB[c]  = ro.y;
    }
    return mat2x3(planeRGB, circRGB);
}

#endif // SDF_GEMS_SDF_BELTRAMI_PHOTOELASTIC_GLSL
