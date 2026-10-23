// ============================================================================
// sdRing2D, sdgRing2D, sdArc2D & sdgArc2D — Exact 0/1-Sqrt 2D Ring Sector & Circular Arc + Gradients
// @author Inigo Quilez — upgraded to 0/1-sqrt via Parabolic v^2-Cancellation using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 0 sqrts in cap Voronoi region, 1 SFU elsewhere, 0 NaNs)
// ============================================================================

#ifndef SDF_GEMS_SD_RING2D_GLSL
#define SDF_GEMS_SD_RING2D_GLSL

// 1. Exact 0/1-Sqrt 2D Annular Ring Sector (0 sqrts in exterior & interior cap Voronoi region via
//    r1 * (r1 - 2*v) < u^2 < r2 * (r2 + 2*v), 1 sqrt elsewhere)
float sdRing2D( in vec2 p, in vec2 n, in float r, in float th )
{
    p.x = abs(p.x);
    float v  = n.x * p.x - n.y * p.y, u = n.y * p.x + n.x * p.y;
    float h  = 0.5 * th, w = abs(u - r) - h, r1 = r - h, r2 = r + h, u2 = u * u;
    if( (v >= 0.0) ? (w <= 0.0) : (u >= r1 && h + v > 0.0 && r1 * (r1 - 2.0 * v) < u2 && u2 < r2 * (r2 + 2.0 * v)) )
        return v;
    float l = length((v > 0.0) ? vec2(v, max(w, 0.0)) : p);
    return (v > 0.0) ? l : abs(l - r) - h;
}

// 2. Exact 0/1-SFU 2D Annular Ring Sector Distance & Analytical Unit Gradient
vec3 sdgRing2D( in vec2 p, in vec2 n, in float r, in float th )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    float v  = n.x * p.x - n.y * p.y, u = n.y * p.x + n.x * p.y;
    float h  = 0.5 * th, w = abs(u - r) - h, r1 = r - h, r2 = r + h, u2 = u * u;
    if( (v >= 0.0) ? (w <= 0.0) : (u >= r1 && h + v > 0.0 && r1 * (r1 - 2.0 * v) < u2 && u2 < r2 * (r2 + 2.0 * v)) )
        return vec3(v, sx * n.x, -n.y);
    if( v > 0.0 )
    {
        float du = ((u < r) ? -1.0 : 1.0) * max(w, 0.0);
        float l2 = v * v + du * du, il = inversesqrt(max(l2, 1e-30));
        return vec3(l2 * il, sx * (n.x * v + n.y * du) * il, (-n.y * v + n.x * du) * il);
    }
    float l2 = dot(p, p), il = inversesqrt(max(l2, 1e-30)), l = l2 * il;
    float sr = (l < r) ? -1.0 : 1.0;
    return vec3(abs(l - r) - h, sx * sr * p.x * il, sr * p.y * il);
}

// 3. Exact 1-Sqrt 2D Round-Capped Circular Arc of radius ra, cap radius rb, and sc = vec2(sin(theta), cos(theta))
float sdArc2D( in vec2 p, in vec2 sc, in float ra, in float rb )
{
    p.x = abs(p.x);
    return ((sc.y * p.x > sc.x * p.y) ? length(p - sc * ra) : abs(length(p) - ra)) - rb;
}

// 4. Exact 1-SFU 2D Round-Capped Circular Arc Distance & Analytical Unit Gradient
#ifndef SDF_GEMS_HAS_SDG_ARC2D
#define SDF_GEMS_HAS_SDG_ARC2D
vec3 sdgArc2D( in vec2 p, in vec2 sc, in float ra, in float rb )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    if( sc.y * p.x > sc.x * p.y )
    {
        vec2  q  = p - sc * ra;
        float l2 = dot(q, q);
        if( l2 > 0.0 )
        {
            float il = inversesqrt(l2);
            return vec3(l2 * il - rb, sx * q.x * il, q.y * il);
        }
        return vec3(-rb, sx * sc.y, -sc.x);
    }
    float l2 = dot(p, p);
    if( l2 > 0.0 )
    {
        float il = inversesqrt(l2), l = l2 * il;
        float sr = (l < ra) ? -1.0 : 1.0;
        return vec3(abs(l - ra) - rb, sx * sr * p.x * il, sr * p.y * il);
    }
    return vec3(ra - rb, 0.0, -1.0);
}
#endif // SDF_GEMS_HAS_SDG_ARC2D

#endif // SDF_GEMS_SD_RING2D_GLSL
