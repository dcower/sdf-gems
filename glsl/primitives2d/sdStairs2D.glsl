// ============================================================================
// sdStairs2D — Exact 1-Sqrt Rotation-Free 2D Staircase (Fixes Tall-Step Sign Bug)
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/ (https://www.shadertoy.com/view/7tKSWt)
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` down from 2, 0 `mat2` rotations, fixes `wh.x < wh.y` sign/Voronoi bug)
// @requires dot2 (from common.glsl)
// ============================================================================

float sdStairs2D( vec2 p, vec2 wh, float n )
{
    vec2 ba = wh * n;
    if( wh.x < wh.y ) { p = ba.yx - p.yx; wh = wh.yx; ba = ba.yx; }
    float d = min(dot2(p - vec2(clamp(p.x, 0.0, ba.x), 0.0)),
                  dot2(p - vec2(ba.x, clamp(p.y, 0.0, ba.y))));
    float s = sign(max(-p.y, p.x - ba.x));

    float id = clamp(round(dot(p, wh) / dot2(wh)), 0.0, n - 1.0);
    p -= id * wh;

    float hh = wh.y * 0.5;
    p.y -= hh;
    if( p.y > hh * sign(p.x) ) s = 1.0;
    p = (id < 0.5 || p.x > 0.0) ? p : -p;

    d = min(d, dot2(p - vec2(0.0, clamp(p.y, -hh, hh))));
    d = min(d, dot2(p - vec2(clamp(p.x, 0.0, wh.x), hh)));
    return sqrt(d) * s;
}
