// ============================================================================
// sdRoundedX2D & sdRoundedXPrism — Exact 0/1-Sqrt 2D Rounded X & 1-Sqrt 3D Prism
// @author Inigo Quilez (1-sqrt 2D) — 0-sqrt arm reduction & 1-sqrt 3D prism discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`0` sqrts in central diamond & arm strips p.x + p.y <= w, `1` in tip caps)
// ============================================================================
//
// Mathematical Derivation:
//   In IQ's sdRoundedX (Shadertoy 7dSGWK):
//     p = abs(p); return length(p - min(p.x + p.y, w) * 0.5) - r;
//   Whenever p.x + p.y <= w (covering the central intersection and all 4 arm normal
//   strips), min(p.x + p.y, w) = p.x + p.y, so:
//     p - (p.x + p.y) * 0.5 = 0.5 * vec2(p.x - p.y, p.y - p.x),
//   whose Euclidean norm is identically:
//     ||0.5 * vec2(p.x - p.y, p.y - p.x)||_2 == abs(p.x - p.y) * sqrt(0.5),
//   requiring ZERO sqrts! A `sqrt` (`length(p - w * 0.5)`) is only evaluated in the
//   4 circular tip caps where p.x + p.y > w.

#ifndef SDF_GEMS_SD_ROUNDED_X2D_GLSL
#define SDF_GEMS_SD_ROUNDED_X2D_GLSL

// 1. Exact 0/1-sqrt 2D Rounded X of arm span w and radius r
float sdRoundedX2D(vec2 p, float w, float r) {
    p = abs(p);
    return ((p.x + p.y <= w) ? abs(p.x - p.y) * 0.70710678 : length(p - w * 0.5)) - r;
}

// 2. Exact 1-sqrt 3D Extruded Rounded-X Prism of arm span w, radius r, and half-height hz
float sdRoundedXPrism(vec3 p, float w, float r, float hz) {
    p = abs(p);
    float d2d = ((p.x + p.y <= w) ? abs(p.x - p.y) * 0.70710678 : length(p.xy - w * 0.5)) - r;
    vec2  d   = vec2(d2d, p.z - hz);
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
}

#endif // SDF_GEMS_SD_ROUNDED_X2D_GLSL
