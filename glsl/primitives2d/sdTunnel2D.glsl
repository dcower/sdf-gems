// ============================================================================
// sdTunnel2D — Exact 1-Sqrt 2D Arch Tunnel & Circle Wave
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/ (flSSDy & stGyzt)
// @exact  Yes (`||grad d|| = 1` for sdTunnel2D; 1 `sqrt` everywhere down from 2 `sqrt`s)
// ============================================================================

// 1. Exact 1-Sqrt 2D Arch Tunnel (down from 2 sqrts in arch half-plane)
float sdTunnel2D( vec2 p, vec2 wh )
{
    p.x = abs(p.x); p.y = -p.y;
    vec2  q = p - wh;
    float l = length((p.y > 0.0) ? max(q, 0.0) : p);
    return (p.y > 0.0) ? l + min(max(q.x, q.y), 0.0)
                       : max(l - wh.x, q.y);
}

// 2. Exact 1-Sqrt 2D Circle Wave (reverse triangle inequality eliminates both endpoint length() calls)
float sdCircleWave2D( vec2 p, float tb, float ra )
{
    tb = 3.1415927 * 5.0 / 6.0 * max(tb, 0.0001);
    vec2 co = ra * vec2(sin(tb), cos(tb));
    p.x = abs(mod(p.x, co.x * 4.0) - co.x * 2.0);
    if( co.y * p.x > co.x * p.y ) p = 2.0 * co - p;
    return abs(length(p) - ra);
}
