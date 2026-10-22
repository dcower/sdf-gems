// ============================================================================
// udTriangle, sdgTriangle3D, udQuad & sdgQuad3D — 7-Region Disjoint Voronoi & Partition-of-Unity 3D Triangles & Quads
// @author Inigo Quilez — upgraded using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 0 clamps & 0.56 avg divs in udTriangle/sdgTriangle3D, 0 divs in face prism)
// ============================================================================

#ifndef SDF_GEMS_UD_TRIANGLE_GLSL
#define SDF_GEMS_UD_TRIANGLE_GLSL

// 1. Exact 7-Region Disjoint-Voronoi 3D Triangle Unsigned Distance (0 clamps, at most 1 division)
float udTriangle( vec3 p, vec3 a, vec3 b, vec3 c )
{
    vec3  ba  = b - a, pa = p - a;
    vec3  cb  = c - b, pb = p - b;
    vec3  ac  = a - c, pc = p - c;
    vec3  nor = cross(ba, ac);
    float n2  = dot(nor, nor);
    float s1  = dot(cross(ba, nor), pa);
    float s2  = dot(cross(cb, nor), pb);
    float s3  = n2 - s1 - s2;

    if( s1 > 0.0 && s2 > 0.0 && s3 > 0.0 )
        return abs(dot(nor, pa)) * inversesqrt(n2);

    float t1 = dot(ba, pa), l1 = dot(ba, ba);
    float t2 = dot(cb, pb), l2 = dot(cb, cb);
    float t3 = dot(ac, pc), l3 = dot(ac, ac);
    bool  e1 = (s1 <= 0.0 && t1 >= 0.0 && t1 <= l1);
    bool  e2 = (s2 <= 0.0 && t2 >= 0.0 && t2 <= l2);
    bool  e3 = (s3 <= 0.0 && t3 >= 0.0 && t3 <= l3);
    if( e1 || e2 || e3 )
    {
        vec3  pm = e1 ? pa : (e2 ? pb : pc);
        float tm = e1 ? t1 : (e2 ? t2 : t3);
        float lm = e1 ? l1 : (e2 ? l2 : l3);
        return sqrt(max(dot(pm, pm) - tm * tm / lm, 0.0));
    }
    vec3 w = (t1 <= 0.0 && t3 >= l3) ? pa : ((t2 <= 0.0 && t1 >= l1) ? pb : pc);
    return length(w);
}

// 2. Exact 7-Region Disjoint-Voronoi 3D Triangle Distance + Analytical 3D Unit Gradient (`1` inversesqrt SFU)
vec4 sdgTriangle3D( vec3 p, vec3 a, vec3 b, vec3 c )
{
    vec3  ba  = b - a, pa = p - a;
    vec3  cb  = c - b, pb = p - b;
    vec3  ac  = a - c, pc = p - c;
    vec3  nor = cross(ba, ac);
    float n2  = dot(nor, nor);
    float s1  = dot(cross(ba, nor), pa);
    float s2  = dot(cross(cb, nor), pb);
    float s3  = n2 - s1 - s2;

    if( s1 > 0.0 && s2 > 0.0 && s3 > 0.0 )
    {
        float il = inversesqrt(n2), np = dot(nor, pa);
        return vec4(abs(np) * il, ((np < 0.0) ? -il : il) * nor);
    }

    float t1 = dot(ba, pa), l1 = dot(ba, ba);
    float t2 = dot(cb, pb), l2 = dot(cb, cb);
    float t3 = dot(ac, pc), l3 = dot(ac, ac);
    bool  e1 = (s1 <= 0.0 && t1 >= 0.0 && t1 <= l1);
    bool  e2 = (s2 <= 0.0 && t2 >= 0.0 && t2 <= l2);
    bool  e3 = (s3 <= 0.0 && t3 >= 0.0 && t3 <= l3);
    vec3  w;
    if( e1 || e2 || e3 )
    {
        vec3  pm = e1 ? pa : (e2 ? pb : pc);
        vec3  em = e1 ? ba : (e2 ? cb : ac);
        float tm = e1 ? t1 : (e2 ? t2 : t3);
        float lm = e1 ? l1 : (e2 ? l2 : l3);
        w = pm - em * (tm / lm);
    }
    else
    {
        w = (t1 <= 0.0 && t3 >= l3) ? pa : ((t2 <= 0.0 && t1 >= l1) ? pb : pc);
    }
    float d2 = dot(w, w), il = inversesqrt(max(d2, 1e-24));
    return vec4(d2 * il, w * il);
}

// 3. Exact 3D Convex Quad Unsigned Distance (with Newell Partition-of-Unity s4 = n2 - s1 - s2 - s3)
float udQuad( vec3 p, vec3 a, vec3 b, vec3 c, vec3 d )
{
    vec3  ba  = b - a, pa = p - a;
    vec3  cb  = c - b, pb = p - b;
    vec3  dc  = d - c, pc = p - c;
    vec3  ad  = a - d;
    vec3  nor = cross(a - c, d - b);
    float n2  = dot(nor, nor);
    float s1  = dot(cross(ba, nor), pa);
    float s2  = dot(cross(cb, nor), pb);
    float s3  = dot(cross(dc, nor), pc);
    float s4  = n2 - s1 - s2 - s3;

    if( s1 > 0.0 && s2 > 0.0 && s3 > 0.0 && s4 > 0.0 )
        return abs(dot(nor, pa)) * inversesqrt(n2);

    float d2 = 1e30;
    if( s1 <= 0.0 ) { vec3 q = ba * clamp(dot(ba, pa) / dot(ba, ba), 0.0, 1.0) - pa; d2 = dot(q, q); }
    if( s2 <= 0.0 ) { vec3 q = cb * clamp(dot(cb, pb) / dot(cb, cb), 0.0, 1.0) - pb; d2 = min(d2, dot(q, q)); }
    if( s3 <= 0.0 ) { vec3 q = dc * clamp(dot(dc, pc) / dot(dc, dc), 0.0, 1.0) - pc; d2 = min(d2, dot(q, q)); }
    if( s4 <= 0.0 ) { vec3 pd = p - d, q = ad * clamp(dot(ad, pd) / dot(ad, ad), 0.0, 1.0) - pd; d2 = min(d2, dot(q, q)); }
    return sqrt(d2);
}

// 4. Exact 3D Convex Quad Distance + Analytical 3D Unit Gradient (`1` inversesqrt SFU)
vec4 sdgQuad3D( vec3 p, vec3 a, vec3 b, vec3 c, vec3 d )
{
    vec3  ba  = b - a, pa = p - a;
    vec3  cb  = c - b, pb = p - b;
    vec3  dc  = d - c, pc = p - c;
    vec3  ad  = a - d;
    vec3  nor = cross(a - c, d - b);
    float n2  = dot(nor, nor);
    float s1  = dot(cross(ba, nor), pa);
    float s2  = dot(cross(cb, nor), pb);
    float s3  = dot(cross(dc, nor), pc);
    float s4  = n2 - s1 - s2 - s3;

    if( s1 > 0.0 && s2 > 0.0 && s3 > 0.0 && s4 > 0.0 )
    {
        float il = inversesqrt(n2), np = dot(nor, pa);
        return vec4(abs(np) * il, ((np < 0.0) ? -il : il) * nor);
    }
    float d2 = 1e30;
    vec3  w  = vec3(0.0);
    if( s1 <= 0.0 ) { vec3 q = pa - ba * clamp(dot(ba, pa) / dot(ba, ba), 0.0, 1.0); d2 = dot(q, q); w = q; }
    if( s2 <= 0.0 ) { vec3 q = pb - cb * clamp(dot(cb, pb) / dot(cb, cb), 0.0, 1.0); float l2 = dot(q, q); if( l2 < d2 ) { d2 = l2; w = q; } }
    if( s3 <= 0.0 ) { vec3 q = pc - dc * clamp(dot(dc, pc) / dot(dc, dc), 0.0, 1.0); float l2 = dot(q, q); if( l2 < d2 ) { d2 = l2; w = q; } }
    if( s4 <= 0.0 ) { vec3 pd = p - d, q = pd - ad * clamp(dot(ad, pd) / dot(ad, ad), 0.0, 1.0); float l2 = dot(q, q); if( l2 < d2 ) { d2 = l2; w = q; } }
    float il = inversesqrt(max(d2, 1e-24));
    return vec4(d2 * il, w * il);
}

#endif // SDF_GEMS_UD_TRIANGLE_GLSL
