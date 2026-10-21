// ============================================================================
// sdTriakisOctahedron — Exact 0/1-Sqrt 3D Triakis Octahedron (Catalan Dual of Truncated Cube)
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prism, 1 sqrt in edge/vertex cone)
// ============================================================================
//
// Mathematical Foundation:
//   The Triakis Octahedron is the 24-faced Catalan dual of the Archimedean
//   Truncated Cube (obtained by erecting a low triangular pyramid onto each of
//   the 8 faces of a regular octahedron). It has 24 congruent isosceles
//   triangular faces, 36 edges, and 14 vertices (6 eight-valent octahedron
//   corners at (+/-s, 0, 0) and 8 three-valent pyramid apices at (+/-a, +/-a, +/-a)
//   where a = s / (1 + sqrt(2)) = (sqrt(2) - 1)*s).
//
//   1. Catalan Single-Face Weyl-Chamber Theorem:
//      In the sorted octahedral Weyl chamber p = abs(p), p.x >= p.y >= p.z >= 0,
//      all 24 faces collapse to a SINGLE face:
//        x + y + (sqrt(2) - 1)*z = s,   unit normal n = (1, 1, k) / sqrt(5 - 2*sqrt(2)),
//      where k = sqrt(2) - 1 approx 0.41421356 and 1/||N|| = 0.67859834.
//      Because N = (1, 1, k) is orthogonal to the chamber wall W_1 (x = y),
//      any point in the sorted chamber projecting onto the face automatically
//      satisfies x' >= y'. The unconstrained face projection lies strictly
//      inside the chamber iff:
//        m := x + y + k*z - s <= min((3 + sqrt(2)/2)*(y - z), (3*sqrt(2) + 1)*z)
//             = min(3.70710678*(y - z), 5.24264069*z),
//      which holds across 100% of the interior (m <= 0) and the exterior face prism (0 sqrts!).
//
//   2. Catalan Wedge-Bisector Unification:
//      In the exterior edge/vertex region, the homogeneous bisector z <= k*y
//      (w_1 <= w_2) 100% separates the Voronoi region of the octahedral edge E_1
//      (z <= k*y) from the pyramid hip edge E_2 (z > k*y, whose tangent (-sqrt(2), 1, 1)
//      has exact integer squared norm 2 + 1 + 1 = 4!), so only ONE edge is ever evaluated!

// 1. Exact 0/1-Sqrt Triakis Octahedron (0 sqrts in interior & face prism, 1 sqrt in edge/vertex region)
float sdTriakisOctahedron( vec3 p, float s )
{
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;

    float m = p.x + p.y + 0.41421356 * p.z - s;
    if( m <= min(5.24264069 * p.z, 3.70710678 * (p.y - p.z)) )
        return m * 0.67859834; // 0 sqrts in interior & face prism!

    float t1 = max(0.5 * (s - p.x + p.y), 0.0);
    float t2 = clamp(0.25 * (1.41421356 * (s - p.x) + p.y + p.z), 0.0, 0.41421356 * s);
    vec3  d  = (p.z <= 0.41421356 * p.y)
             ? vec3(p.x - s + t1, p.y - t1, p.z)
             : vec3(p.x - s + 1.41421356 * t2, p.y - t2, p.z - t2);
    return length(d);
}

// 2. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) (0/1-inversesqrt)
vec4 sdgTriakisOctahedron( vec3 p, float s )
{
    vec3 sg = vec3((p.x < 0.0) ? -1.0 : 1.0,
                   (p.y < 0.0) ? -1.0 : 1.0,
                   (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) p.xy = p.yx;
    bool s2 = (p.y < p.z); if( s2 ) p.yz = p.zy;
    bool s3 = (p.x < p.y); if( s3 ) p.xy = p.yx;

    float m = p.x + p.y + 0.41421356 * p.z - s;
    float d;
    vec3  g;
    if( m <= min(5.24264069 * p.z, 3.70710678 * (p.y - p.z)) )
    {
        d = m * 0.67859834;
        g = vec3(0.67859834, 0.67859834, 0.28108464);
    }
    else
    {
        float t1   = max(0.5 * (s - p.x + p.y), 0.0);
        float t2   = clamp(0.25 * (1.41421356 * (s - p.x) + p.y + p.z), 0.0, 0.41421356 * s);
        vec3  v    = (p.z <= 0.41421356 * p.y)
                   ? vec3(p.x - s + t1, p.y - t1, p.z)
                   : vec3(p.x - s + 1.41421356 * t2, p.y - t2, p.z - t2);
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
