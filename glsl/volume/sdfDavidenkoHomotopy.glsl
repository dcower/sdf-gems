#ifndef SDF_GEMS_DAVIDENKO_HOMOTOPY_GLSL
#define SDF_GEMS_DAVIDENKO_HOMOTOPY_GLSL

// 2nd-Order Screen-Space Davidenko Surface Continuation Predictor & Arnold Fold Detector.
// Given reference hit (t0, d0) with unit normal g0 = gradF(p0) and Hessian H0 = grad^2 F(p0),
// predicts hit depth tPred along neighbor ray d1.
// Returns vec3(tPred, foldDiscriminant, isValidSheet) where isValidSheet == 1.0 if no fold jump occurred.
vec3 sdfDavidenkoSurfacePredict(
    float t0, vec3 d0, vec3 d1,
    vec3 g0, mat3 H0, float silhouetteEps
) {
    vec3 dd = d1 - d0;
    float mu0 = dot(g0, d0);
    float mu1 = dot(g0, d1);
    float gdd = dot(g0, dd);

    // Normal curvature along reference ray
    float kappaN = dot(d0, H0 * d0);
    // Approximate SDF value at p0 + t0 * dd
    float fShift = t0 * gdd + 0.5 * t0 * t0 * dot(dd, H0 * dd);
    float foldDisc = mu1 * mu1 - 2.0 * kappaN * fShift;

    // Detect Arnold A2 fold bifurcation (silhouette jump)
    bool valid = (mu0 < -silhouetteEps) && (mu1 < -silhouetteEps) && (foldDisc >= 0.0);
    float safeMu0 = min(mu0, -max(silhouetteEps, 1e-4));

    // 1st-order Davidenko tangent predictor
    float dt1 = -(t0 * gdd) / safeMu0;
    vec3 v1 = dt1 * d0 + t0 * dd;

    // 2nd-order Davidenko curvature correction
    float vHv = dot(v1, H0 * v1);
    float dt2 = (2.0 * dt1 * dt1) / max(t0, 1e-4) - vHv / safeMu0;

    // Clamp 2nd-order correction near fold singularities for robustness
    float maxCorr = 0.5 * abs(dt1) + 0.1 * t0;
    float halfDt2 = clamp(0.5 * dt2, -maxCorr, maxCorr);

    float tPred = max(t0 + dt1 + halfDt2, 0.0);
    return vec3(tPred, foldDisc, valid ? 1.0 : 0.0);
}

// 1-Step Eikonal-Newton Surface Corrector given predicted depth tPred and sdg(o + tPred * d1) = (F, gradF).
float sdfDavidenkoSurfaceCorrect(float tPred, vec3 d1, vec4 dgPred) {
    float slope = dot(dgPred.yzw, d1);
    // If slope is negative, take damped Newton step; otherwise take conservative sphere-tracing step
    float stepT = (slope < -1e-3)
        ? (-dgPred.x / slope)
        : dgPred.x;
    return max(tPred + clamp(stepT, -2.0 * abs(dgPred.x), 2.0 * abs(dgPred.x)), 0.0);
}

// Screen-Space Davidenko Volumetric Closest-Approach (Tangential Silhouette) Predictor.
// Tracks the parameter t_crit where dot(gradF(o + t*d), d) == 0 from ray d0 to neighbor ray d1.
// Returns vec2(tCritPred, fMinPred).
vec2 sdfDavidenkoClosestApproachPredict(
    float tCrit0, float fMin0,
    vec3 d0, vec3 d1,
    vec3 gCrit, mat3 HCrit
) {
    vec3 dd = d1 - d0;
    float kappaN = max(dot(d0, HCrit * d0), 1e-4);
    float rhs = dot(gCrit, dd) + tCrit0 * dot(d0, HCrit * dd);
    float dtCrit = -rhs / kappaN;

    float tCritPred = max(tCrit0 + dtCrit, 0.0);
    vec3 dp = dtCrit * d0 + tCrit0 * dd;
    float fMinPred = fMin0 + dot(gCrit, dp) + 0.5 * dot(dp, HCrit * dp);
    return vec2(tCritPred, fMinPred);
}
#endif
