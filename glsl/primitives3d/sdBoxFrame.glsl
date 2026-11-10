// ============================================================================
// sdBoxFrame — Exact 1-Sqrt 3D Wireframe Box Frame (Closest-Beam Reduction)
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 1 sqrt down from 3; exterior identical to the 3-beam form, interior exact
//         including the 8 corner cubes, where the 3-beam form is only a bound)
// ============================================================================

// 1. Exact 1-Sqrt Square-Strut Box Frame (struts occupy |p_i| in [b_i - 2e, b_i] on two axes)
//    Outside the solid: closest-beam reduction, length(max(vec3(min w, q), 0)).
//    Inside the solid: the nearest non-solid point is either across an outer face (-max w) or in the hollow,
//    which needs TWO coordinates below b - 2e: sqrt of the two smallest t_i^2, t = max(w + 2e, 0).
//    (The former interior term min(max(q.x, q.y), 0) treated the inner planes w_j = -2e as surface inside the
//    8 corner cubes, where they are interior to the solid: it under-estimated the depth there by up to
//    (sqrt(2) - 1) * 2e. Fixed with Claude Fable 5.2, sdf-gems 2026.)
float sdBoxFrame( vec3 p, vec3 b, float e )
{
    p = abs(p) - b;
    float a = min(p.x, p.y);
    vec2  q = abs(vec2(max(a, p.z), max(p.x, p.y)) + e) - e;
    float l = length(max(vec3(min(a, p.z), q), 0.0));
    if( l > 0.0 ) return l;
    // An axis with b_i <= 2e has no hollow side: its struts span the whole box.
    vec3  t  = vec3((b.x > 2.0 * e) ? max(p.x + 2.0 * e, 0.0) : 1e15,
                    (b.y > 2.0 * e) ? max(p.y + 2.0 * e, 0.0) : 1e15,
                    (b.z > 2.0 * e) ? max(p.z + 2.0 * e, 0.0) : 1e15);
    float lo = min(t.x, t.y), hi = max(t.x, t.y);
    float h  = length(vec2(min(lo, t.z), max(lo, min(hi, t.z))));   // two smallest of t
    return -min(-max(p.x, max(p.y, p.z)), h);
}

// 2. Exact 3-Line 1-Sqrt Cylindrical-Edge Wireframe Cage (8 spherical corners, 12 cylindrical struts)
float sdBoxEdges( vec3 p, vec3 b, float r )
{
    p = abs(p) - b;
    float s = min(min(p.x, min(p.y, p.z)), 0.0);
    return sqrt(dot(p, p) - s * s) - r;
}
