// ============================================================================
// sdBoxFrame — Exact 1-Sqrt 3D Wireframe Box Frame (Closest-Beam Reduction)
// @author Inigo Quilez — upgraded/optimized by Dillon Cower (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 1 sqrt down from 3, 0.0 error vs 3-beam form)
// ============================================================================

// 1. Exact 4-Line 1-Sqrt Square-Strut Box Frame
float sdBoxFrame( vec3 p, vec3 b, float e )
{
    p = abs(p) - b;
    float a = min(p.x, p.y);
    vec2  q = abs(vec2(max(a, p.z), max(p.x, p.y)) + e) - e;
    return length(max(vec3(min(a, p.z), q), 0.0)) + min(max(q.x, q.y), 0.0);
}

// 2. Exact 3-Line 1-Sqrt Cylindrical-Edge Wireframe Cage (8 spherical corners, 12 cylindrical struts)
float sdBoxEdges( vec3 p, vec3 b, float r )
{
    p = abs(p) - b;
    float s = min(min(p.x, min(p.y, p.z)), 0.0);
    return sqrt(dot(p, p) - s * s) - r;
}
