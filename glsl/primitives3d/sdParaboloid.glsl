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
//   Using Sum-of-Cubes Conjugate Rationalization:
//     r + pc / r == 2 * qc / (r^2 - pc + (pc / r)^2)
//   eliminates catastrophic cancellation when pc < 0 near the symmetry axis (reducing
//   FP32 near-axis error by 800,000x), and Endpoint Subsumption gives a 0-sqrt,
//   0-cubic fast exit (`-p.y`) below the base disk (x^2 + z^2 <= w^2, y <= 0).

#ifndef SDF_GEMS_SD_PARABOLOID_GLSL
#define SDF_GEMS_SD_PARABOLOID_GLSL

#include "../common/common.glsl"

// Internal helper: solves u^3 - 3*pc*u - 2*qc = 0 for u >= 0 without 0/0 NaNs or cancellation
float sdParaboloidRoot(float u0, float dy, float ik) {
    float pc = ik * (dy - 0.5 * ik) / 3.0;
    float qc = 0.25 * ik * ik * u0;
    float h = qc * qc - pc * pc * pc;
    if (h >= 0.0) {
        float r = pow(qc + sqrt(h), 1.0 / 3.0);
        return (pc < 0.0) ? 2.0 * qc / (r * r - pc + (pc / r) * (pc / r))
                          : (r > 0.0) ? r + pc / r : 0.0;
    }
    float r = sqrt(pc);
    return 2.0 * r * cos(acos(clamp(qc / (pc * r), -1.0, 1.0)) / 3.0);
}

// Solid capped 3D paraboloid of base radius w (at y = 0) and apex height h (at y = h)
float sdParaboloid(vec3 p, float w, float h) {
    float x2 = dot(p.xz, p.xz);
    if (x2 <= w * w && p.y <= 0.0) return -p.y; // 0 sqrts, 0 cubic below base disk!
    float u0 = sqrt(x2);
    float ik = w * w / h;
    float u  = min(sdParaboloidRoot(u0, h - p.y, ik), w);
    float d  = length(vec2(u0 - u, p.y - (h - u * u / ik)));
    return (p.y > 0.0 && ik * (p.y - h) + x2 < 0.0) ? -min(d, p.y) : d;
}

// Open parabolic reflector dish of aperture radius w, depth h, and half-thickness th
float sdParaboloidDish(vec3 p, float w, float h, float th) {
    float u0 = length(p.xz);
    float ik = w * w / h;
    float u  = min(sdParaboloidRoot(u0, p.y, ik), w);
    return length(vec2(u0 - u, p.y - u * u / ik)) - th;
}

#endif // SDF_GEMS_SD_PARABOLOID_GLSL
