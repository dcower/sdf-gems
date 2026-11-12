// ============================================================================
// sdTruncatedCube — Exact 0/1-Sqrt 3D Archimedean Truncated Cube (14 Faces)
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & triangle face prism, 1 sqrt in octagon chain)
// ============================================================================
//
// Mathematical Foundation:
//   The Archimedean Truncated Cube has 14 faces (6 regular octagons |x| <= s,
//   |y| <= s, |z| <= s and 8 equilateral triangles |x| + |y| + |z| <= (1 + sqrt(2))*s),
//   36 edges, and 24 vertices at permutations of (+/-s, +/-s, +/-(sqrt(2) - 1)*s).
//   Standard Shadertoy/hg_sdf implementations use the 7-plane intersection bound
//   max(max(|x|,|y|,|z|) - s, (|x|+|y|+|z| - (1+sqrt(2))*s)/sqrt(3)), which
//   underestimates exterior Euclidean distance by up to 29.29% near vertices and
//   produces sharp planes instead of rounded edges under `opRound` (sd(p) - r).
//
//   1. Weyl-Chamber Reduction & Interior / Triangle-Face Early Exit (0 Sqrts):
//      After folding p into the fundamental Weyl chamber B_3 of the octahedral
//      group O_h via p = abs(p) and sorting p.x >= p.y >= p.z >= 0 (order 48),
//      all 14 faces, 36 edges, and 24 vertices reduce to:
//        - 1 Octagonal face F_1:  x = s
//        - 1 Triangular face F_2: x + y + z = (1 + sqrt(2))*s
//        - 2 Edges E_1, E_2 and 1 Vertex v_0 = (s, s, (sqrt(2) - 1)*s).
//      Whenever p lies inside the polyhedron (max(d_1, d_2) <= 0) OR in the
//      exterior triangular face prism (p.x <= s or 2*p.x - p.y - p.z <= (2 - sqrt(2))*s),
//      the exact signed distance is simply max(d_1, d_2) with ZERO sqrts!
//
//   2. 2D Octagon Chain Unification Theorem:
//      Outside the triangular face Voronoi region, all remaining features in
//      the Weyl chamber (the octagonal face F_1, both edges E_1 and E_2, and
//      vertex v_0) lie in the plane x = s and form a 2D clipped octagon chain
//      in (y, z) with a single clamped parameter:
//        t = clamp(0.5*(y + z - sqrt(2)*s), 0.0, max(z - (sqrt(2) - 1)*s, 0.0)),
//      giving the exact 3D Euclidean difference vector:
//        d = vec3(x - s, max(t, y - s), t).

// 1. Exact 0/1-Sqrt Archimedean Truncated Cube (0 sqrts in interior & triangle face prism, 1 sqrt in octagon chain)
float sdTruncatedCube( vec3 p, float s )
{
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;

    // 0-sqrt early exit for 100% of interior and the exterior triangular face prism
    if( min(p.x - s, 2.0 * p.x - p.y - p.z - 0.58578644 * s) <= 0.0 )
        return max(p.x - s, (p.x + p.y + p.z - 2.41421356 * s) * 0.57735027);

    // Unified 2D octagon chain for octagonal face, both edges, and vertex v_0
    float t = clamp(0.5 * (p.y + p.z - 1.41421356 * s), 0.0, max(p.z - 0.41421356 * s, 0.0));
    return length(vec3(p.x - s, max(t, p.y - s), t));
}

// 2. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) (0/1-inversesqrt)
vec4 sdgTruncatedCube( vec3 p, float s )
{
    vec3 sg = vec3((p.x < 0.0) ? -1.0 : 1.0,
                   (p.y < 0.0) ? -1.0 : 1.0,
                   (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) p.xy = p.yx;
    bool s2 = (p.y < p.z); if( s2 ) p.yz = p.zy;
    bool s3 = (p.x < p.y); if( s3 ) p.xy = p.yx;

    float d1 = p.x - s;
    float d2 = (p.x + p.y + p.z - 2.41421356 * s) * 0.57735027;
    float d;
    vec3  g;
    if( min(d1, 2.0 * p.x - p.y - p.z - 0.58578644 * s) <= 0.0 )
    {
        d = max(d1, d2);
        g = (d1 >= d2) ? vec3(1.0, 0.0, 0.0) : vec3(0.57735027);
    }
    else
    {
        float t    = clamp(0.5 * (p.y + p.z - 1.41421356 * s), 0.0, max(p.z - 0.41421356 * s, 0.0));
        vec3  v    = vec3(d1, max(t, p.y - s), t);
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
