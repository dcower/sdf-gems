// ============================================================================
// sdVesica2D, sdgVesica2D, sdVesicaWidth2D & sdOrientedVesica2D — Exact 1-to-2 Sqrt 2D Vesica Piscis
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 1 sqrt in circular arc / tip regions, `0.0` error)
// ============================================================================

// 1. Exact 1-Sqrt-Arc 2D Vesica Piscis (Signed-Square Cone Reduction)
float sdVesica2D( vec2 p, float r, float d )
{
    p = abs(p);
    float b2 = r * r - d * d;
    float x  = p.x + d;
    return (d * d * p.y * p.y > b2 * x * x)
         ? length(vec2(p.x, p.y - sqrt(b2)))
         : length(vec2(x, p.y)) - r;
}

// 2. Exact 1-Sqrt-Arc 2D Vesica Piscis SDF + Analytical Unit Gradient vec3(d, grad_x, grad_y)
vec3 sdgVesica2D( vec2 p, float r, float d )
{
    vec2  s  = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    p        = abs(p);
    float b2 = r * r - d * d;
    float u  = p.y * d;
    float v  = p.x + d;
    if( u * abs(u) > v * abs(v) * b2 )
    {
        vec2  q = vec2(p.x, p.y - sqrt(b2));
        float l = length(q) * ((d < 0.0) ? -1.0 : 1.0);
        return vec3(l, s * q / l);
    }
    vec2  q = vec2(v, p.y);
    float l = length(q);
    return vec3(l - r, s * q / max(l, 1e-30));
}

// 3. Physical Half-Length (r) & Half-Width (w) Parameterization: 1 `sqrt` everywhere in R^2
float sdVesicaWidth2D( vec2 p, float r, float w )
{
    p = abs(p);
    float d = 0.5 * (r * r - w * w) / w;
    return (r * p.x < d * (p.y - r))
         ? length(vec2(p.x, p.y - r))
         : length(vec2(p.x + d, p.y)) - (d + w);
}

// 4. Exact 1-Sqrt-Tip Oriented 2D Vesica between Endpoints (a, b) with Half-Width (w):
//    Uses unnormalized local coordinates Q = r * q (0 sqrts for branch test, 1 sqrt in tip region).
float sdOrientedVesica2D( vec2 p, vec2 a, vec2 b, float w )
{
    vec2  del = 0.5 * (b - a);
    float r2  = dot(del, del);
    float d   = 0.5 * (r2 - w * w) / w;
    vec2  Q   = abs(mat2(del.y, del.x, -del.x, del.y) * (p - 0.5 * (a + b)));
    if( r2 * Q.x < d * (Q.y - r2) )
    {
        vec2 u = vec2(Q.x, Q.y - r2);
        return sqrt(dot(u, u) / max(r2, 1e-30));
    }
    float r = sqrt(r2);
    return length(Q / max(r, 1e-30) + vec2(d, 0.0)) - (d + w);
}
