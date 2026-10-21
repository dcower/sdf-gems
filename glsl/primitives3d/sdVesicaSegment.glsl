// ============================================================================
// sdVesicaSegment — Exact 1-to-2 Sqrt 3D Revolved Vesica / Spindle
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` in tip cones & 2 in body, down from 3 `sqrt`s)
// @requires dot2 (from common.glsl)
// ============================================================================

float sdVesicaSegment( vec3 p, vec3 a, vec3 b, float w )
{
    vec3  ba = b - a;
    vec3  pc = p - (a + b) * 0.5;
    float l2 = dot(ba, ba);
    float u  = dot(pc, ba);
    float d  = (l2 - 4.0 * w * w) / (8.0 * w);
    float x2 = dot2(pc - ba * (u / l2));
    float z  = d * (2.0 * abs(u) - l2);
    return (z > 0.0 && l2 * l2 * x2 < z * z)
         ? length(pc - ba * (u >= 0.0 ? 0.5 : -0.5))
         : sqrt(dot2(pc) + d * d + 2.0 * d * sqrt(x2)) - (d + w);
}
