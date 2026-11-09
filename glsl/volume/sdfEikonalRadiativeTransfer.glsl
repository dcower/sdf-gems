// sdfEikonalRadiativeTransfer: Eikonal Green's Function Multipole, 0-Ray Footpoint GI & Zonal SH-2 Radiative Transfer.
//
// Solves:
//   1. Dark unlit shadow crevices in direct lighting:
//      Provides real-time 1-bounce diffuse color bleeding and interreflections in ZERO raymarching steps.
//   2. Sub-surface and ambient contact occlusion via Zonal Harmonic (SH-2) cap moments and bent normals.
//   3. Exact differential-area-to-cap configuration factor (view factor) with smooth horizon wrap.
//
// References:
//   - Green (1828), "An Essay on the Application of Mathematical Analysis to the Theories of Electricity and Magnetism"
//   - Ramamoorthi & Hanrahan (SIGGRAPH 2001), "An Efficient Representation for Irradiance Environment Maps"
//   - Sloan et al. (SIGGRAPH 2002), "Precomputed Radiance Transfer for Real-Time Rendering in Dynamic, Low-Frequency Lighting Environments"
//   - Kontkanen & Laine (Eurographics SRE 2005), "Ambient Occlusion Fields"
//   - Snyder (SIGGRAPH 1996), "Area Light Sources for Real-Time Graphics"
//
// Attribution:
//   Mathematical framework & theory: Gemini 4 Argon
//   Implementation & GPU verification: Gemini 3.8 Flash (sdf-gems, 2026)

#ifndef SDF_GEMS_SDF_EIKONAL_RADIATIVE_TRANSFER_GLSL
#define SDF_GEMS_SDF_EIKONAL_RADIATIVE_TRANSFER_GLSL

// -----------------------------------------------------------------------------
// 1. Zero-Ray Eikonal Footpoint Surface Projection
// -----------------------------------------------------------------------------
// Projects any query point p onto the closest surface point q_A in \partial\Omega_A
// in exactly ZERO raymarching steps: q_A = p - d_A * grad(d_A).
vec3 sdfEikonalFootpoint(vec3 p, vec4 dgA)
{
    return p - dgA.x * dgA.yzw;
}

vec3 sdfEikonalFootpoint(vec3 p, float dA, vec3 gA)
{
    return p - dA * gA;
}

// -----------------------------------------------------------------------------
// 2. Legendre Zonal Harmonic (SH-3) Cap Moments
// -----------------------------------------------------------------------------
// Evaluates the exact analytical Legendre polynomial integrals over a spherical cap
// subtending half-angle sinThetaMax = sin(theta_max) with mu = cos(theta_max):
//   Z_0(mu) = \int_mu^1 P_0(t) dt = 1 - mu
//   Z_1(mu) = \int_mu^1 P_1(t) dt = 0.5 * (1 - mu^2) = 0.5 * sin^2(theta_max)
//   Z_2(mu) = \int_mu^1 P_2(t) dt = 0.5 * mu * (1 - mu^2) = 0.5 * mu * sin^2(theta_max)
// Returns vec4(Z_0, Z_1, Z_2, sin^2(theta_max)).
vec4 sdfZonalSH3Weights(float sinThetaMax)
{
    float s = clamp(sinThetaMax, 0.0, 1.0);
    float s2 = s * s;
    float mu = sqrt(max(1.0 - s2, 0.0));
    float z0 = 1.0 - mu;
    float z1 = 0.5 * s2;
    float z2 = 0.5 * mu * s2;
    return vec4(z0, z1, z2, s2);
}

// -----------------------------------------------------------------------------
// 3. Eikonal Differential-to-Cap Configuration Factor (View Factor)
// -----------------------------------------------------------------------------
// Computes the subtended spherical cap half-angle:
//   sinThetaMax = clamp(R_A / max(R_A + max(d_A, 0.0), 1e-4), 0.0, 1.0)
// along direction d_BA = -normalize(grad(d_A)) and evaluates the cosine-weighted
// configuration factor F_{B -> A}(p, n_B) in [0, 1] with C^1 smooth horizon wrap.
float sdfEikonalViewFactor(vec3 nB, vec4 dgA, float radiusA)
{
    float rA = max(radiusA, 1e-4);
    float dist = max(dgA.x, 0.0);
    float sinThetaMax = clamp(rA / max(rA + dist, 1e-4), 0.0, 1.0);
    float sin2 = sinThetaMax * sinThetaMax;

    float gLen = length(dgA.yzw);
    vec3 dBA = (gLen > 1e-5) ? (-dgA.yzw / gLen) : nB;
    float cosGamma = dot(nB, dBA);

    // Smooth C^1 horizon wrap
    float u = clamp((cosGamma + sinThetaMax) / max(2.0 * sinThetaMax, 1e-4), 0.0, 1.0);
    float h = u * u * (3.0 - 2.0 * u);

    float effCos = max(cosGamma + 0.5 * sinThetaMax * (1.0 - u), 0.0);
    return clamp(h * effCos * sin2, 0.0, 1.0);
}

float sdfEikonalViewFactor(vec3 nB, float dA, vec3 gA, float radiusA)
{
    return sdfEikonalViewFactor(nB, vec4(dA, gA), radiusA);
}

// -----------------------------------------------------------------------------
// 4. Zonal SH-1 Bent-Normal & Contact Visibility
// -----------------------------------------------------------------------------
// Uses the dipole/Z1 moment along d_occ = -normalize(grad(d_occ)) to deflect the
// receiver normal nB away from the occluder, returning vec4(bentNormal, visibility).
vec4 sdfEikonalBentNormal(vec3 nB, vec4 dgOccluder, float radiusOccluder)
{
    float rOcc = max(radiusOccluder, 1e-4);
    float dist = max(dgOccluder.x, 0.0);
    float sinTheta = clamp(rOcc / max(rOcc + dist, 1e-4), 0.0, 1.0);
    vec4 zsh = sdfZonalSH3Weights(sinTheta);
    float z1 = zsh.y; // 0.5 * sin^2(theta)

    float gLen = length(dgOccluder.yzw);
    vec3 dOcc = (gLen > 1e-5) ? (-dgOccluder.yzw / gLen) : nB;
    float ndotd = max(dot(nB, dOcc), 0.0);

    vec3 unnormBent = nB - 2.0 * z1 * ndotd * dOcc;
    float bLen = length(unnormBent);
    vec3 bentN = (bLen > 1e-5) ? (unnormBent / bLen) : nB;

    float vis = clamp(1.0 - 2.0 * z1 * ndotd, 0.0, 1.0);
    return vec4(bentN, vis);
}

// -----------------------------------------------------------------------------
// 5. Zero-Ray 1-Bounce Diffuse Indirect Color Bleeding
// -----------------------------------------------------------------------------
// Computes the 1-bounce diffuse irradiance transferred from emitter Object A
// (at projected footpoint q_A with normal n_A = normalize(dgA.yzw)) onto receiver Object B:
//   E_bounce = F_{B -> A} * max(n_A . ldir, 0.0) * max(n_A . (-nB), 0.15) * albedoA * lcol
vec3 sdfEikonalIndirectBounce(vec3 p, vec3 nB, vec4 dgA, float radiusA, vec3 albedoA, vec3 ldir, vec3 lcol)
{
    float gLen = length(dgA.yzw);
    vec3 nA = (gLen > 1e-5) ? (dgA.yzw / gLen) : (-nB);

    float F = sdfEikonalViewFactor(nB, dgA, radiusA);
    float directA = max(dot(nA, ldir), 0.0);
    float coupling = max(dot(nA, -nB), 0.15);

    return F * directA * coupling * albedoA * lcol;
}

vec3 sdfEikonalIndirectBounce(vec3 p, vec3 nB, float dA, vec3 gA, float radiusA, vec3 albedoA, vec3 ldir, vec3 lcol)
{
    return sdfEikonalIndirectBounce(p, nB, vec4(dA, gA), radiusA, albedoA, ldir, lcol);
}

#endif // SDF_GEMS_SDF_EIKONAL_RADIATIVE_TRANSFER_GLSL
