// ============================================================================
// sdTetrakisHexahedron — Exact 0/1-Sqrt 3D Tetrakis Hexahedron (Catalan Dual of Truncated Octahedron)
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prism, 1 sqrt in edge/vertex cone)
// ============================================================================
//
// Mathematical Foundation:
//   The Tetrakis Hexahedron is the 24-faced Catalan dual of the Archimedean
//   Truncated Octahedron (obtained by erecting a square pyramid onto each of
//   the 6 faces of a cube, and a classic crystal habit of fluorite). It has
//   24 congruent isosceles triangular faces, 36 edges (12 cube-like edges and
//   24 pyramid hip edges), and 14 vertices (6 four-valent pyramid apices at
//   (+/-s, 0, 0) and 8 six-valent cube corners at (+/-2s/3, +/-2s/3, +/-2s/3)).
//
//   1. Catalan Single-Face Weyl-Chamber Theorem:
//      After folding p into the fundamental Weyl chamber B_3 of the octahedral
//      group O_h via p = abs(p) and sorting p.x >= p.y >= p.z >= 0 (order 48),
//      all 24 faces collapse to a SINGLE face:
//        2*x + y = 2*s,   unit normal n = (2, 1, 0) / sqrt(5),
//        inradius h = (2 / sqrt(5)) * s.
//      Because n = (2, 1, 0) is orthogonal to the chamber wall W_3 (z = 0),
//      any point in the sorted chamber projecting onto the face automatically
//      satisfies z' = z >= 0. Furthermore, the unconstrained face projection
//      lies strictly inside the chamber iff:
//        m := 2*x + y - 2*s <= 5.0 * min(x - y, y - z).
//      This holds across 100% of the polyhedron interior (m <= 0) and the
//      entire exterior face prism (0 sqrts!).
//
//   2. Catalan Wedge-Bisector Unification & Pythagorean Weyl-Basis Form:
//      In the exterior edge/vertex region (m > 5.0*min(x - y, y - z)), the
//      homogeneous bisector x + z = 2*y (w_1 = w_2) 100% separates the Voronoi
//      region of the cube-like edge E_1 (x + z <= 2*y) from the pyramid hip
//      edge E_2 (x + z > 2*y), so only ONE edge is ever evaluated!
//      Furthermore, because ||(2,1,0)||^2 = 5, ||(2,-4,5)||^2 = 45, and
//      ||(-1,2,2)||^2 = 9 form an orthogonal rational frame with integer weights
//      (9, 1, 5) / 45, the distance also admits a branchless 6-line Pythagorean
//      Weyl-basis form!

// 1. Exact 0/1-Sqrt Tetrakis Hexahedron (0 sqrts in interior & face prism, 1 sqrt in edge/vertex region)
float sdTetrakisHexahedron( vec3 p, float s )
{
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;

    float m = 2.0 * p.x + p.y - 2.0 * s;
    if( m <= 5.0 * min(p.x - p.y, p.y - p.z) )
        return m * 0.44721360; // 1.0 / sqrt(5.0): 0 sqrts in interior & face prism!

    float c = (2.0 / 3.0) * s;
    float u = clamp((s - p.x + 2.0 * p.y + 2.0 * p.z) * (1.0 / 9.0), 0.0, (1.0 / 3.0) * s);
    vec3  d = (p.x + p.z <= 2.0 * p.y)
            ? vec3(p.x - c, p.y - c, max(p.z - c, 0.0))
            : vec3(p.x - s + u, p.y - 2.0 * u, p.z - 2.0 * u);
    return length(d);
}

// 2. Ultra-Compact Branchless Pythagorean Weyl-Basis Form (6 lines, 0 sorting branches, 1 sqrt)
float sdTetrakisHexahedronBranchless( vec3 p, float s )
{
    p = abs(p);
    float x = max(p.x, max(p.y, p.z)), z = min(p.x, min(p.y, p.z));
    float y = p.x + p.y + p.z - x - z, m = 2.0 * x + y - 2.0 * s;
    float a = max(m - 5.0 * min(x - y, y - z), 0.0);
    float b = (x + z <= 2.0 * y) ? max(3.0 * z - 2.0 * s, 0.0)
                                 : max(abs(2.0 * y + 2.0 * z - x - 0.5 * s) - 1.5 * s, 0.0);
    return ((m < 0.0) ? -1.0 : 1.0) * sqrt((9.0 * m * m + a * a + 5.0 * b * b) / 45.0);
}

// 3. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) (0/1-inversesqrt)
vec4 sdgTetrakisHexahedron( vec3 p, float s )
{
    vec3 sg = vec3((p.x < 0.0) ? -1.0 : 1.0,
                   (p.y < 0.0) ? -1.0 : 1.0,
                   (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) p.xy = p.yx;
    bool s2 = (p.y < p.z); if( s2 ) p.yz = p.zy;
    bool s3 = (p.x < p.y); if( s3 ) p.xy = p.yx;

    float m = 2.0 * p.x + p.y - 2.0 * s;
    float d;
    vec3  g;
    if( m <= 5.0 * min(p.x - p.y, p.y - p.z) )
    {
        d = m * 0.44721360;
        g = vec3(0.89442719, 0.44721360, 0.0);
    }
    else
    {
        float c    = (2.0 / 3.0) * s;
        float u    = clamp((s - p.x + 2.0 * p.y + 2.0 * p.z) * (1.0 / 9.0), 0.0, (1.0 / 3.0) * s);
        vec3  v    = (p.x + p.z <= 2.0 * p.y)
                   ? vec3(p.x - c, p.y - c, max(p.z - c, 0.0))
                   : vec3(p.x - s + u, p.y - 2.0 * u, p.z - 2.0 * u);
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
