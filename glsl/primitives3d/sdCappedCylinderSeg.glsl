// ============================================================================
// sdCappedCylinderSeg — Exact 3D Arbitrary-Segment Capped Cylinder (FP16 Overflow Fixed)
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdCappedCylinderSeg( vec3 p, vec3 a, vec3 b, float r )
{
    vec3  ba = b - a, pa = p - a;
    float L  = length(ba);
    vec3  n  = ba / L;
    float y  = dot(pa, n);
    float x  = length(pa - n * y) - r;
    vec2  d  = vec2(x, abs(y - 0.5 * L) - 0.5 * L);
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
}
