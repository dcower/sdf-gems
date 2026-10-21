// ============================================================================
// sdgRevolution — Exact Full & Limited-Angle 3D Revolution with Gradients
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

// Full 360-degree revolution gradient finalize for d2g = sdgShape2D(vec2(length(p.xz) - o, p.y))
vec4 sdgRevolutionFinalize( vec3 p, vec3 d2g )
{
    float l = length(p.xz);
    return vec4(d2g.x, d2g.y * p.x / l, d2g.z, d2g.y * p.z / l);
}

// Limited-angle partial revolution gradient finalize for q = opRevolutionLim(p, o, sc) and d2g = sdgShape2D(q.xy)
vec4 sdgRevolutionLimFinalize( vec3 p, vec2 sc, vec3 d2g )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    vec2  q  = vec2(abs(p.x), p.z);
    float u  = q.x * sc.x + q.y * sc.y;
    float w  = q.x * sc.y - q.y * sc.x;
    float l  = length(q);
    vec3  gw = vec3(sx * sc.y, 0.0, -sc.x);
    if( w > 0.0 ) {
        if( d2g.x <= 0.0 ) return vec4(w, gw);
        float d  = length(vec2(d2g.x, w));
        vec3  gu = vec3(sx * d2g.y * sc.x, d2g.z, d2g.y * sc.y);
        return vec4(d, (d2g.x * gu + w * gw) / d);
    }
    if( d2g.x < 0.0 && u > 0.0 && w > d2g.x ) return vec4(w, gw);
    return vec4(d2g.x, sx * d2g.y * q.x / l, d2g.z, d2g.y * q.y / l);
}
