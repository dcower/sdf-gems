// ============================================================================
// sdTunnel2D — Exact 0/1-Sqrt 2D Arch Tunnel, Analytical Gradient, & 3D Prism
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/ (flSSDy)
// @exact  Yes (`||grad d|| = 1` everywhere; 0 sqrts in rectangular strips, 1 sqrt in arch & corners; down from 2 sqrts)
// ============================================================================

// 1. Exact 0/1-Sqrt 2D Arch Tunnel (0 sqrts in rectangular interior & wall strips, 1 sqrt in arch & bottom corner)
float sdTunnel2D( vec2 p, vec2 wh )
{
    p.x = abs(p.x); p.y = -p.y;
    vec2 q = p - wh;
    if( p.y > 0.0 )
        return (q.x > 0.0 && q.y > 0.0) ? length(q) : max(q.x, q.y);
    return max(length(p) - wh.x, q.y);
}

// 2. Analytical 2D Distance and Unit Gradient vec3(d, grad.xy) (0/1-inversesqrt)
vec3 sdgTunnel2D( vec2 p, vec2 wh )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x); p.y = -p.y;
    vec2 q = p - wh;
    if( p.y > 0.0 )
    {
        if( q.x > 0.0 && q.y > 0.0 )
        {
            float q2   = dot(q, q);
            float invQ = inversesqrt(q2);
            return vec3(q2 * invQ, sx * q.x * invQ, -q.y * invQ);
        }
        return (q.x > q.y) ? vec3(q.x, sx, 0.0) : vec3(q.y, 0.0, -1.0); // 0 SFU ops!
    }
    float p2   = dot(p, p);
    float invP = inversesqrt(max(p2, 1e-30));
    float da   = p2 * invP - wh.x;
    return (da >= q.y) ? vec3(da, sx * p.x * invP, -p.y * invP) : vec3(q.y, 0.0, -1.0);
}
