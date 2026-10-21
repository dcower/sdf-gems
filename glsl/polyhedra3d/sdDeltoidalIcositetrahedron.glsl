// ============================================================================
// sdDeltoidalIcositetrahedron — Exact 0/1-Sqrt 3D Deltoidal Icositetrahedron (Catalan Dual of Rhombicuboctahedron)
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prism, 1 sqrt in edge/vertex cone)
// ============================================================================
//
// Mathematical Foundation:
//   The Deltoidal Icositetrahedron (trapezikositetrahedron) is the 24-faced
//   Catalan dual of the Archimedean Rhombicuboctahedron, and is a classic
//   crystal habit of garnet and leucite. It has 24 congruent kite (deltoid)
//   faces, 48 edges, and 26 vertices:
//     - 6 four-valent axial vertices at (+/-s, 0, 0),
//     - 12 four-valent face-diagonal vertices at (+/-s/sqrt(2), +/-s/sqrt(2), 0),
//     - 8 three-valent body-diagonal vertices at (+/-a, +/-a, +/-a) where
//       a = (2*sqrt(2) - 1)/3 * s.
//
//   1. Catalan Single-Face Weyl-Chamber Theorem:
//      In the sorted octahedral Weyl chamber p = abs(p), p.x >= p.y >= p.z >= 0,
//      all 24 kite faces collapse to a SINGLE face:
//        x + (sqrt(2) - 1)*(y + z) = s,   unit normal n = (1, k, k) / sqrt(7 - 4*sqrt(2)),
//      where k = sqrt(2) - 1 approx 0.41421356 and 1/||N|| = 0.86285621.
//      Because N = (1, k, k) is orthogonal to the chamber wall W_2 (y = z),
//      any point in the sorted chamber projecting onto the face automatically
//      satisfies y' >= z'. The unconstrained face projection lies strictly
//      inside the chamber iff:
//        m := x + k*(y + z) - s <= min(3.24264069*z, 2.29289322*(x - y)),
//      which holds across 100% of the interior (m <= 0) and the exterior face prism (0 sqrts!).
//
//   2. Catalan Wedge-Bisector Unification:
//      Since (1 - k)/k = sqrt(2), the homogeneous wedge bisector w_1 <= w_2
//      simplifies identically to sqrt(2)*z <= x - y, which 100% separates the
//      lower kite edge E_1 (on z = 0) from the upper kite edge E_2 (on x = y),
//      and both edges share the same tangent squared norm 1 + k^2 = 4 - 2*sqrt(2)
//      = 1.0 / 0.85355339!

// 1. Exact 0/1-Sqrt Deltoidal Icositetrahedron (0 sqrts in interior & face prism, 1 sqrt in edge/vertex region)
float sdDeltoidalIcositetrahedron( vec3 p, float s )
{
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;

    float m = p.x + 0.41421356 * (p.y + p.z) - s;
    if( m <= min(3.24264069 * p.z, 2.29289322 * (p.x - p.y)) )
        return m * 0.86285621; // 0 sqrts in interior & face prism!

    float t1 = clamp((0.41421356 * (s - p.x) + p.y) * 0.85355339, 0.0, 0.70710678 * s);
    float t2 = clamp((0.29289322 * (1.41421356 * s - p.x - p.y) + p.z) * 0.85355339, 0.0, 0.54691816 * s);
    float xy = 0.70710678 * s - 0.29289322 * t2;
    vec3  d  = (1.41421356 * p.z <= p.x - p.y)
             ? vec3(p.x - s + 0.41421356 * t1, p.y - t1, p.z)
             : vec3(p.x - xy, p.y - xy, p.z - t2);
    return length(d);
}

// 2. Analytical 3D Distance and Unit Gradient vec4(d, grad.xyz) (0/1-inversesqrt)
vec4 sdgDeltoidalIcositetrahedron( vec3 p, float s )
{
    vec3 sg = vec3((p.x < 0.0) ? -1.0 : 1.0,
                   (p.y < 0.0) ? -1.0 : 1.0,
                   (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) p.xy = p.yx;
    bool s2 = (p.y < p.z); if( s2 ) p.yz = p.zy;
    bool s3 = (p.x < p.y); if( s3 ) p.xy = p.yx;

    float m = p.x + 0.41421356 * (p.y + p.z) - s;
    float d;
    vec3  g;
    if( m <= min(3.24264069 * p.z, 2.29289322 * (p.x - p.y)) )
    {
        d = m * 0.86285621;
        g = vec3(0.86285621, 0.35740674, 0.35740674);
    }
    else
    {
        float t1   = clamp((0.41421356 * (s - p.x) + p.y) * 0.85355339, 0.0, 0.70710678 * s);
        float t2   = clamp((0.29289322 * (1.41421356 * s - p.x - p.y) + p.z) * 0.85355339, 0.0, 0.54691816 * s);
        float xy   = 0.70710678 * s - 0.29289322 * t2;
        vec3  v    = (1.41421356 * p.z <= p.x - p.y)
                   ? vec3(p.x - s + 0.41421356 * t1, p.y - t1, p.z)
                   : vec3(p.x - xy, p.y - xy, p.z - t2);
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
