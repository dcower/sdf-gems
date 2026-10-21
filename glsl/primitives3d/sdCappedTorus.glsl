// ============================================================================
// sdCappedTorus — Exact 3D Capped Torus (Cancellation-Free Core Circle)
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 0% NaN rate)
// ============================================================================

float sdCappedTorus( vec3 p, vec2 sc, float ra, float rb )
{
    p.x = abs(p.x);
    vec2 q = (sc.y*p.x > sc.x*p.y) ? p.xy - ra*sc : vec2(length(p.xy) - ra, 0.0);
    return length(vec3(q, p.z)) - rb;
}
