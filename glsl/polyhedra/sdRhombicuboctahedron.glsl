// ============================================================================
// sdRhombicuboctahedron — Exact 0/1-Sqrt 3D Archimedean Small Rhombicuboctahedron (26 Faces)
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & triangle face prism, 1 sqrt in octagon chain)
// ============================================================================
//
// Mathematical Foundation:
//   The Archimedean Small Rhombicuboctahedron has 26 faces (6 axial squares
//   |x| <= s, 12 diagonal squares |x| + |y| <= sqrt(2)*s, and 8 corner
//   equilateral triangles |x| + |y| + |z| <= (2*sqrt(2) - 1)*s), 48 edges, and
//   24 vertices at permutations of (+/-s, +/-(sqrt(2) - 1)*s, +/-(sqrt(2) - 1)*s).
//   Standard 13-plane intersection bounds underestimate exterior Euclidean
//   distance by up to 11.31% near vertices and fail under `opRound` (sd(p) - r).
//
//   1. Weyl-Chamber Reduction & Interior / Triangle-Face Early Exit (0 Sqrts):
//      After folding p into the fundamental Weyl chamber B_3 of O_h via
//      p = abs(p) and sorting p.x >= p.y >= p.z >= 0 (order 48), all 26 faces,
//      48 edges, and 24 vertices reduce to:
//        - Axial square face F_1:    d_1 = x - s
//        - Diagonal square face F_2: d_2 = (x + y - sqrt(2)*s) / sqrt(2)
//        - Triangular face F_3:      d_3 = (x + y + z - (2*sqrt(2) - 1)*s) / sqrt(3)
//        - 2 Edges E_1, E_2 and 1 Vertex v_0 = (s, (sqrt(2) - 1)*s, (sqrt(2) - 1)*s).
//      Whenever p lies inside the polyhedron (max(d_1, d_2, d_3) <= 0) OR in
//      the exterior triangular face prism (d_3 > 0 and x + y - 2*z <= (2 - sqrt(2))*s),
//      the exact signed distance is max(d_1, max(d_2, d_3)) with ZERO sqrts!
//
//   2. 2D Octagon Chain Unification Theorem:
//      Outside the triangular face Voronoi region, all 5 remaining exterior
//      regions (the axial square F_1, diagonal square F_2, both edges E_1 and
//      E_2, and vertex v_0) decouple into z-clamping at (sqrt(2) - 1)*s and a
//      2D clipped octagon chain in (x, y) with a single clamped parameter:
//        t = clamp(0.5*(x + y - sqrt(2)*s), 0.0, max(y - (sqrt(2) - 1)*s, 0.0)),
//      giving the exact 3D Euclidean difference vector:
//        d = vec3(max(t, x - s), t, max(z - (sqrt(2) - 1)*s, 0.0)).

// 1. Exact 0/1-Sqrt Archimedean Rhombicuboctahedron (0 sqrts in interior & triangle face prism, 1 sqrt in octagon chain)
float sdRhombicuboctahedron( vec3 p, float s )
{
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;

    float d1 = p.x - s;
    float d2 = (p.x + p.y - 1.41421356 * s) * 0.70710678;
    float d3 = (p.x + p.y + p.z - 1.82842712 * s) * 0.57735027;
    float dm = max(d1, max(d2, d3));

    // 0-sqrt early exit for 100% of interior and the exterior triangular face prism
    if( dm <= 0.0 || (d3 > 0.0 && p.x + p.y - 2.0 * p.z <= 0.58578644 * s) )
        return dm;

    // Unified 2D octagon chain in (x, y) + z-clamping for all remaining faces, edges, and vertex v_0
    float t = clamp(0.5 * (p.x + p.y - 1.41421356 * s), 0.0, max(p.y - 0.41421356 * s, 0.0));
    return length(vec3(max(t, d1), t, max(p.z - 0.41421356 * s, 0.0)));
}

// 2. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) (0/1-inversesqrt)
vec4 sdgRhombicuboctahedron( vec3 p, float s )
{
    vec3 sg = vec3((p.x < 0.0) ? -1.0 : 1.0,
                   (p.y < 0.0) ? -1.0 : 1.0,
                   (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) p.xy = p.yx;
    bool s2 = (p.y < p.z); if( s2 ) p.yz = p.zy;
    bool s3 = (p.x < p.y); if( s3 ) p.xy = p.yx;

    float d1 = p.x - s;
    float d2 = (p.x + p.y - 1.41421356 * s) * 0.70710678;
    float d3 = (p.x + p.y + p.z - 1.82842712 * s) * 0.57735027;
    float dm = max(d1, max(d2, d3));
    float d;
    vec3  g;
    if( dm <= 0.0 || (d3 > 0.0 && p.x + p.y - 2.0 * p.z <= 0.58578644 * s) )
    {
        d = dm;
        g = (d1 >= d2 && d1 >= d3) ? vec3(1.0, 0.0, 0.0)
          : (d2 >= d3)             ? vec3(0.70710678, 0.70710678, 0.0)
                                   : vec3(0.57735027);
    }
    else
    {
        float t    = clamp(0.5 * (p.x + p.y - 1.41421356 * s), 0.0, max(p.y - 0.41421356 * s, 0.0));
        vec3  v    = vec3(max(t, d1), t, max(p.z - 0.41421356 * s, 0.0));
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = v * invV;
    }

    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4(d, sg * g);
}
