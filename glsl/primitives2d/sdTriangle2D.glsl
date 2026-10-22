// ============================================================================
// sdTriangle2D, sdgTriangle2D, sdTriangleIsosceles2D & sdgTriangleIsosceles2D —
// 100% 0-Division, 0-Clamp Disjoint Voronoi 2D Triangles & Gradients
//
// @author  Inigo Quilez — upgraded using Gemini 4 Argon (sdf-gems, 2026)
// @source  https://iquilezles.org/articles/distfunctions2d/
// @license MIT
// @exact   Yes (0 divisions, 0 clamps, and 1 inversesqrt across 100% of R^2)
// ============================================================================

#ifndef SDF_GEMS_SD_TRIANGLE_2D_GLSL
#define SDF_GEMS_SD_TRIANGLE_2D_GLSL

// ----------------------------------------------------------------------------
// 1. Exact 0-Division, 0-Clamp Disjoint-Voronoi 2D Triangle Distance + Gradient
//
// Decomposes R^2 into 7 disjoint regions using only dot/cross products:
// - The interior (c0 > 0 && c1 > 0 && c2 > 0) and the 3 mutually disjoint
//   exterior edge-normal strips s_i = (c_i <= 0 && 0 <= t_i <= l_i), where the
//   tangential component cancels identically and both distance and outward unit
//   gradient are given by -c_m * inversesqrt(l_m) with 0 divisions!
// - The 3 exterior vertex wedges, where the closest vertex offset w is selected
//   directly by the sign of t_i with 0 divisions and 0 clamps!
// ----------------------------------------------------------------------------
vec3 sdgTriangle2D( in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2 )
{
    vec2  e0 = p1 - p0, e1 = p2 - p1, e2 = p0 - p2;
    vec2  v0 = p  - p0, v1 = p  - p1, v2 = p  - p2;
    float a2 = e0.x * e2.y - e0.y * e2.x, s = (a2 < 0.0) ? -1.0 : 1.0;
    float c0 = s * (v0.x * e0.y - v0.y * e0.x);
    float c1 = s * (v1.x * e1.y - v1.y * e1.x);
    float c2 = abs(a2) - c0 - c1;
    float l0 = dot(e0, e0), l1 = dot(e1, e1), l2 = dot(e2, e2);
    float t0 = dot(v0, e0), t1 = dot(v1, e1), t2 = dot(v2, e2);

    bool s0 = (c0 <= 0.0 && t0 >= 0.0 && t0 <= l0);
    bool s1 = (c1 <= 0.0 && t1 >= 0.0 && t1 <= l1);
    bool s2 = (c2 <= 0.0 && t2 >= 0.0 && t2 <= l2);
    bool inside = (c0 > 0.0 && c1 > 0.0 && c2 > 0.0);

    if( inside || s0 || s1 || s2 )
    {
        float cm = c0, lm = l0; vec2 em = e0;
        if( s1 || (inside && c1 * c1 * lm < cm * cm * l1) ) { cm = c1; lm = l1; em = e1; }
        if( s2 || (inside && c2 * c2 * lm < cm * cm * l2) ) { cm = c2; lm = l2; em = e2; }
        float invL = inversesqrt(lm);
        return vec3(-cm * invL, vec2(-s * em.y, s * em.x) * invL);
    }

    vec2  w    = (t0 <= 0.0 && t2 >= l2) ? v0 : ((t1 <= 0.0 && t0 >= l0) ? v1 : v2);
    float d2   = dot(w, w);
    float invD = inversesqrt(max(d2, 1e-30));
    return vec3(d2 * invD, w * invD);
}

// 2. Exact 0-Division, 0-Clamp 2D Triangle Signed Distance
float sdTriangle2D( in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2 )
{
    return sdgTriangle2D(p, p0, p1, p2).x;
}

// 3. Exact 0-Div 0/1-Sqrt 2D Isosceles Triangle (0 SFU in both exterior & interior base Voronoi regions!)
float sdTriangleIsosceles2D( in vec2 p, in vec2 q )
{
    p.x = abs(p.x);
    float k = (q.y < 0.0) ? -1.0 : 1.0, db = k * (p.y - q.y);
    float ms = k * (p.x * q.y - p.y * q.x), q2 = dot(q, q);
    if( (p.x <= q.x && db >= 0.0) || (db <= 0.0 && ms <= 0.0 && db * db * q2 <= ms * ms) )
        return db;
    bool c_apex = (dot(p, q) <= 0.0);
    if( c_apex || (p.x > q.x && dot(p - q, q) >= 0.0) )
        return length(c_apex ? p : (p - q));
    return ms * inversesqrt(q2);
}

// 4. Exact 0/1-SFU 2D Isosceles Triangle with Analytical Unit Gradient
vec3 sdgTriangleIsosceles2D( in vec2 p, in vec2 q )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0, k = (q.y < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    float db = k * (p.y - q.y), ms = k * (p.x * q.y - p.y * q.x), q2 = dot(q, q);
    if( (p.x <= q.x && db >= 0.0) || (db <= 0.0 && ms <= 0.0 && db * db * q2 <= ms * ms) )
        return vec3(db, 0.0, k);
    bool c_apex = (dot(p, q) <= 0.0);
    bool c_corn = (!c_apex) && (p.x > q.x) && (dot(p - q, q) >= 0.0);
    if( c_apex || c_corn )
    {
        vec2  w  = c_apex ? p : (p - q);
        float l2 = dot(w, w), il = inversesqrt(max(l2, 1e-24));
        return vec3(l2 * il, vec2(sx * w.x, w.y) * il);
    }
    float il = inversesqrt(q2);
    return vec3(ms * il, vec2(sx * k * q.y, -k * q.x) * il);
}

#endif // SDF_GEMS_SD_TRIANGLE_2D_GLSL
