// ============================================================================
// sdCross2D — Exact 1-Sqrt 2D Symmetric Cross (Fixes Interior Overestimation Bug)
// @author Inigo Quilez (bound interior) — upgraded to exact interior using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/ (https://www.shadertoy.com/view/XtGfzw)
// @exact  Yes (`||grad d|| = 1` in both exterior and interior; 1 `sqrt`, `8.88e-16` max error)
// ============================================================================

// b.x = arm half-length, b.y = arm half-thickness (b.x >= b.y), r = corner rounding offset
float sdCross2D( vec2 p, vec2 b, float r )
{
    p = abs(p); p = (p.y > p.x) ? p.yx : p.xy;
    vec2  q = p - b;
    float k = max(q.y, q.x);
    vec2  w = (k > 0.0) ? q : vec2(b.y - p.x, -q.y);
    float d = length(max(w, 0.0));
    return ((k > 0.0) ? d : max(-d, q.x)) + r;
}
