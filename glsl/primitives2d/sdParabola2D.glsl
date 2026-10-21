// Singularity-Free Exact 2D Parabola & Capped Parabola Segment SDF
// Original author: Inigo Quilez (https://iquilezles.org/articles/distfunctions2d/)
// Singularity & y-axis sign-inversion fixes: Discovered using Gemini 4 Argon (sdf-gems, 2026)
//
// Improvements over IQ's sdParabola (ws3GD7) & sdParabola segment (3lSczz):
//   1. Fixes the 0.0 / 0.0 NaN at the center of curvature (0, 1/(2k)) where pc = 0
//      and qc = 0 cause IQ's else-branch `acos(q / (p * r))` to evaluate `acos(0/0) = NaN`.
//   2. Fixes IQ's y-axis sign-inversion and zero-collapse bug in `sdParabola(pos, k)`:
//      IQ uses `sign(pos.x - x)` for the exterior/interior sign. On the symmetry axis
//      `pos.x == 0` below the center of curvature (`pos.y <= 0.5/k`), the exact root is
//      `x = 0`, so `sign(0.0 - 0.0) = 0.0` collapses the entire negative y-axis to 0
//      (or in FP32, `r = pow(sqrt(-p^3), 1/3)` rounds slightly above `-p`, making
//      `x = r + p/r > 0` and flipping `sign(0 - x) = -1` on the entire exterior y < 0!).
//      Using the algebraic implicit sign `sign(k * pos.x * pos.x - pos.y)` is exact
//      everywhere in R^2 and decoupled from root-finding rounding error.

#ifndef SDF_GEMS_SD_PARABOLA2D_GLSL
#define SDF_GEMS_SD_PARABOLA2D_GLSL

#include "../common/common.glsl"

// Exact signed distance to the infinite 2D parabola y = k * x^2 (negative inside y > k*x^2)
float sdParabola2D(vec2 pos, float k) {
    pos.x = abs(pos.x);
    float ik = 1.0 / k;
    float p = ik * (pos.y - 0.5 * ik) / 3.0;
    float q = 0.25 * ik * ik * pos.x;
    float h = q * q - p * p * p;
    float x;
    if (h >= 0.0) {
        float r = pow(q + sqrt(h), 1.0 / 3.0);
        x = (r > 0.0) ? max(r + p / r, 0.0) : 0.0;
    } else {
        float r = sqrt(p);
        x = 2.0 * r * cos(acos(clamp(q / (p * r), -1.0, 1.0)) / 3.0);
    }
    float s = (k * pos.x * pos.x > pos.y) ? 1.0 : -1.0;
    return length(pos - vec2(x, k * x * x)) * s;
}

// Exact signed distance to the 2D parabola segment y = he - (he/wi^2)*x^2 for |x| <= wi
float sdParabolaSeg2D(vec2 pos, float wi, float he) {
    pos.x = abs(pos.x);
    float ik = wi * wi / he;
    float p = ik * (he - pos.y - 0.5 * ik) / 3.0;
    float q = 0.25 * ik * ik * pos.x;
    float h = q * q - p * p * p;
    float x;
    if (h >= 0.0) {
        float r = pow(q + sqrt(h), 1.0 / 3.0);
        x = (r > 0.0) ? max(r + p / r, 0.0) : 0.0;
    } else {
        float r = sqrt(p);
        x = 2.0 * r * cos(acos(clamp(q / (p * r), -1.0, 1.0)) / 3.0);
    }
    x = min(x, wi);
    float s = (ik * (pos.y - he) + pos.x * pos.x > 0.0) ? 1.0 : -1.0;
    return length(pos - vec2(x, he - x * x / ik)) * s;
}

#endif // SDF_GEMS_SD_PARABOLA2D_GLSL
