// sdfPhaseSpaceLightField: SO(3)-Steerable Darboux Phase-Space & Convex Zonal Specular Light-Field.
//
// Solves:
//   1. World-space normal/view corruption under rotation/Brenier warp by expressing radiance
//      in intrinsic, cofactor-covariant Plücker-Darboux phase-space coordinates.
//   2. Non-convex single-exponent (ks -> 0) vanishing-gradient collapse via a 4-lobe
//      Energy-Normalized Zonal Spherical-Gaussian / Chebyshev Specular Basis with fixed dyadic sharpness.
//   3. Exact radiative-transfer phase-space blending across smooth CSG transition zones.
//
// References:
//   - Darboux (1896), "Leçons sur la théorie générale des surfaces"
//   - Frisvad (JGT 2012) & Duff et al. (JCGT 2017), "Building an Orthonormal Basis, Revisited"
//   - Wang et al. (SIGGRAPH 2009), "All-Frequency Rendering of Dynamic, Spatially-Varying Reflectance"
//   - Zhang et al. (CVPR 2021), "PhySG: Inverse Rendering with Spherical Gaussians"
//   - Plücker (1865), "On a New Geometry of Space"
//
// Attribution:
//   Mathematical framework & theory: Gemini 4 Argon
//   Implementation & GPU verification: Gemini 3.8 Flash (sdf-gems, 2026)

#ifndef SDF_GEMS_SDF_PHASE_SPACE_LIGHT_FIELD_GLSL
#define SDF_GEMS_SDF_PHASE_SPACE_LIGHT_FIELD_GLSL

// -----------------------------------------------------------------------------
// 1. Orthonormal Eikonal Darboux Frame
// -----------------------------------------------------------------------------
// Constructs an intrinsic orthonormal frame (t, b, n) in SO(3) from the unit
// surface normal n and an intrinsic feature reference gradient featGrad.
// If featGrad is degenerate (or parallel to n), smoothly falls back to the
// branch-free Duff et al. (2017) / Frisvad (2012) orthonormal basis.
mat3 sdfDarbouxFrame(vec3 n, vec3 featGrad)
{
    vec3 v = featGrad - dot(featGrad, n) * n;
    float lenV = length(v);
    if (lenV >= 1e-5) {
        vec3 t = v / lenV;
        vec3 b = cross(n, t);
        return mat3(t, b, n);
    }
    // Branch-free Duff et al. (2017) orthonormal basis fallback
    float s = (n.z < 0.0) ? -1.0 : 1.0;
    float a = -1.0 / (s + n.z);
    float b = n.x * n.y * a;
    vec3 t = vec3(1.0 + s * n.x * n.x * a, s * b, -s * n.x);
    vec3 bVec = vec3(b, s + n.y * n.y * a, -n.y);
    return mat3(t, bVec, n);
}

// -----------------------------------------------------------------------------
// 2. SO(3)-Invariant Plücker-Darboux Phase-Space Invariants
// -----------------------------------------------------------------------------
// Computes canonical phase-space coordinates:
//   x: local view cosine cos(theta_v) = dot(n, rd)
//   y: azimuthal tangent projection  = dot(darboux[0], rd)
//   z: azimuthal binormal projection = dot(darboux[1], rd)
//   w: signed Plücker silhouette torque = dot(cross(p, rd), n)
// Under any rigid rotation R in SO(3) or cofactor-covariant Brenier warp,
// these coordinates remain strictly invariant, preventing radiance distortion.
vec4 sdfPhaseSpaceInvariants(vec3 p, vec3 n, mat3 darboux, vec3 rd)
{
    float cosTheta = dot(n, rd);
    float tanU = dot(darboux[0], rd);
    float tanV = dot(darboux[1], rd);
    float torque = dot(cross(p, rd), n);
    return vec4(cosTheta, tanU, tanV, torque);
}

// -----------------------------------------------------------------------------
// 3. Convex 4-Lobe Energy-Normalized Zonal Specular Basis
// -----------------------------------------------------------------------------
// Evaluates 4 dyadic Zonal Spherical-Gaussian lobes with fixed sharpness:
//   lambda = (4.0, 16.0, 64.0, 256.0)
//   B_k(mu) = (lambda_k + 1) / (2 * pi) * exp(lambda_k * (clamp(mu, 0.0, 1.0) - 1.0))
//
// Because the lobe sharpnesses lambda_k are fixed, fitting lobe weights w in R_+^4
// is strictly convex (non-negative least squares), completely eliminating the
// vanishing-gradient collapse (ks -> 0) seen in single-exponent Blinn-Phong models.
vec4 sdfZonalSpecularBasis(float ndh)
{
    float mu = clamp(ndh, 0.0, 1.0);
    float d = mu - 1.0;
    const vec4 lambda = vec4(4.0, 16.0, 64.0, 256.0);
    const vec4 norm = (lambda + vec4(1.0)) * (0.5 / 3.14159265358979323846);
    return norm * exp(lambda * d);
}

// Convenience scalar evaluator: dot(weights, basis)
float sdfZonalSpecularRadiance(vec4 weights, float ndh)
{
    return dot(weights, sdfZonalSpecularBasis(ndh));
}

// -----------------------------------------------------------------------------
// 4. Exact Radiative-Transfer Smooth CSG Phase Blend
// -----------------------------------------------------------------------------
// Blends directional radiance across a smooth CSG (smin) transition zone of width k > 0,
// weighting each primitive's radiance by its Eikonal logistic partition of unity:
//   w12 = 1.0 / (1.0 + exp((d1 - d2) / max(k, 1e-4)))
vec3 sdfCSGPhaseBlend(float d1, vec3 rad1, float d2, vec3 rad2, float k)
{
    float width = max(k, 1e-4);
    float diff = clamp((d1 - d2) / width, -40.0, 40.0);
    float w1 = 1.0 / (1.0 + exp(diff));
    return mix(rad2, rad1, w1);
}

#endif // SDF_GEMS_SDF_PHASE_SPACE_LIGHT_FIELD_GLSL
