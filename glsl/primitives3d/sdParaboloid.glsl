// Exact 3D Solid Capped Paraboloid & Parabolic Reflector Dish SDF
// Discovered using Gemini 4 Argon (sdf-gems, 2026)
//
// Computes the exact 3D Euclidean signed distance to:
//   1. sdParaboloid(p, w, h): A solid capped paraboloid (bullet / nose-cone / dome)
//      of circular base radius w at y = 0 and apex at y = h, bounded by the revolved
//      parabolic mantle Y(u) = h * (1 - u^2 / w^2) for u = length(p.xz) in [0, w]
//      and the flat circular base disk y = 0, u in [0, w].
//   2. sdParaboloidDish(p, w, h, th): A hollow parabolic reflector dish / satellite
//      dish of aperture radius w, depth h, and shell half-thickness th.
//
// Mathematical Derivation:
//   In cylindrical coordinates (u0, y0) = (length(p.xz), p.y), the squared Euclidean
//   distance from (u0, y0) to a point (u, h - u^2 / ik) on the parabola (where
//   ik = w^2 / h) has derivative:
//     (1/2) d(D^2)/du = (u - u0) - (2 u / ik) * (h - u^2 / ik - y0) = 0
//   which is directly a depressed cubic in u (zero quadratic term!):
//     u^3 - 3 * pc * u - 2 * qc = 0
//   with pc = ik * (h - y0 - 0.5 * ik) / 3 and qc = 0.25 * ik^2 * u0 >= 0.
//   Solving this depressed cubic via a singularity-free Cardano/trigonometric solver
//   (guarded against 0/0 NaNs at the focal center of curvature u0 = 0, y0 = h - 0.5*ik)
//   and clamping u to [0, w] gives the exact closest point on the parabolic mantle.

#ifndef SDF_GEMS_SD_PARABOLOID_GLSL
#define SDF_GEMS_SD_PARABOLOID_GLSL

#include "../common/common.glsl"

// Internal helper: solves u^3 - 3*pc*u - 2*qc = 0 for u >= 0 without 0/0 NaNs
float sdParaboloidRoot(float u0, float dy, float ik) {
    float pc = ik * (dy - 0.5 * ik) / 3.0;
    float qc = 0.25 * ik * ik * u0;
    float h = qc * qc - pc * pc * pc;
    if (h >= 0.0) {
        float r = pow(qc + sqrt(h), 1.0 / 3.0);
        return (r > 0.0) ? max(r + pc / r, 0.0) : 0.0;
    }
    float r = sqrt(pc);
    return 2.0 * r * cos(acos(clamp(qc / (pc * r), -1.0, 1.0)) / 3.0);
}

// Solid capped 3D paraboloid of base radius w (at y = 0) and apex height h (at y = h)
float sdParaboloid(vec3 p, float w, float h) {
    float u0 = length(p.xz);
    float ik = w * w / h;
    float u = min(sdParaboloidRoot(u0, h - p.y, ik), w);
    float dm2 = dot2(vec2(u0 - u, p.y - (h - u * u / ik)));
    float db2 = dot2(vec2(max(u0 - w, 0.0), p.y));
    float s = (p.y > 0.0 && ik * (p.y - h) + u0 * u0 < 0.0) ? -1.0 : 1.0;
    return s * sqrt(min(dm2, db2));
}

// Open parabolic reflector dish of aperture radius w, depth h, and half-thickness th
float sdParaboloidDish(vec3 p, float w, float h, float th) {
    float u0 = length(p.xz);
    float ik = w * w / h;
    float u = min(sdParaboloidRoot(u0, p.y, ik), w);
    return length(vec2(u0 - u, p.y - u * u / ik)) - th;
}

#endif // SDF_GEMS_SD_PARABOLOID_GLSL
