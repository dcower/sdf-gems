// ============================================================================
// opHyperbolicBrenier — Hyperbolic Poincaré Ball Schwarz-Pick Step Bounds & Monge-Ampère W2 Optimal Transport Warps
// @author Gemini 4 Argon (mathematical derivation) & Gemini 3.8 Flash (implementation & verification) (sdf-gems, 2026)
// @exact  opSchwarzPickSafeStep: Exact Euclidean safe sphere-tracing bound for hyperbolic metrics (zero boundary violation, zero overstep)
//         opBrenierTransportMap: Strictly convex Brenier optimal transport warp with exact Weyl spectral Lipschitz bound and Monge-Ampere determinant
// ============================================================================

#ifndef SDF_GEMS_HYPERBOLIC_BRENIER_GLSL
#define SDF_GEMS_HYPERBOLIC_BRENIER_GLSL

// Exact Mobius Gyrovector Addition x (+) y in the Poincare ball of radius R.
vec3 opMobiusAdd(vec3 x, vec3 y, float R) {
    float invR2 = 1.0 / max(R * R, 1e-8);
    float x2 = dot(x, x) * invR2;
    float y2 = dot(y, y) * invR2;
    float xy2 = 2.0 * dot(x, y) * invR2;
    float denom = max(1.0 + xy2 + x2 * y2, 1e-8);
    return ((1.0 + xy2 + y2) * x + (1.0 - x2) * y) / denom;
}

// Hyperbolic Riemannian log-map coordinate at origin: maps x in B_R^3 to tangent space R^3
// such that ||opPoincareLogMap(x, R)||_2 == d_H(0, x).
vec3 opPoincareLogMap(vec3 x, float R) {
    float r = length(x);
    float u = min(r / max(R, 1e-6), 0.9999);
    // artanh(u) = 0.5 * log((1+u)/(1-u))
    float artanhU = 0.5 * log((1.0 + u) / (1.0 - u));
    float scale = (u > 1e-5) ? (2.0 * R * artanhU / r) : 2.0;
    return scale * x;
}

// Exact Schwarz-Pick Euclidean Safe Step Bound for sphere-tracing a hyperbolic SDF dH
// inside the Poincare ball of radius R. Guarantees r + dSafe < R for all finite dH > 0.
float opSchwarzPickSafeStep(vec3 x, float dH, float R) {
    float r = length(x);
    if (r >= R) return 0.0;
    float tau = tanh(clamp(dH / (2.0 * max(R, 1e-6)), -20.0, 20.0));
    return ((R * R - r * r) * tau) / max(R + r * tau, 1e-6);
}

// Strictly Convex Brenier Optimal Transport Map y = grad(psi_tau)(x).
// Returns transported coordinate y in .xyz, and writes:
//   outSpectralBound = exact upper bound on lambda_max(H_tau(x)) for 1-Lipschitz SDF scaling,
//   outMongeAmpereDet = exact det(H_tau(x)) for mass-conserving volumetric density.
vec3 opBrenierTransportMap(
    vec3 x, float tau, vec3 diagTarget, float radialAmp, float radialCore,
    out float outSpectralBound, out float outMongeAmpereDet
) {
    float r2 = dot(x, x);
    float b2 = max(radialCore * radialCore, 1e-6);
    float Rb = sqrt(r2 + b2);
    float Rb3 = Rb * (r2 + b2);

    float tauA = tau * radialAmp;
    vec3 s = vec3(1.0 - tau) + tau * diagTarget + vec3(tauA / Rb);

    // Transported point y = grad(psi_tau)(x)
    vec3 y = s * x;

    // Exact Weyl spectral bound for 1-Lipschitz SDF division
    float maxBase = max(1.0 - tau + tau * diagTarget.x,
                    max(1.0 - tau + tau * diagTarget.y, 1.0 - tau + tau * diagTarget.z));
    float radMax = max(max(tauA / Rb, (tauA * b2) / Rb3), 0.0);
    outSpectralBound = max(maxBase + radMax, 1e-4);

    // Exact Monge-Ampere determinant via Matrix Determinant Lemma
    vec3 safeS = max(s, vec3(1e-4));
    float rank1Term = (tauA / Rb3) * dot(x * x, 1.0 / safeS);
    outMongeAmpereDet = max(safeS.x * safeS.y * safeS.z * (1.0 - rank1Term), 0.0);

    return y;
}
#endif
