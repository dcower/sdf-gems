// ============================================================================
// sdTruncatedIcosahedron & sdgTruncatedIcosahedron — Exact 0/1-Sqrt 3D Archimedean Truncated Icosahedron (C60 Soccer Ball)
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prisms, 1 sqrt in edge/vertex cone)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for the Archimedean
//          truncated icosahedron (12 regular pentagons + 20 regular hexagons, 60 vertices,
//          90 edges) via 3-step unrolled I_h Weyl chamber folding and 1-vertex Voronoi partition.
//
// Parameters:
//   p : 3D query point
//   r : inradius to the 12 regular pentagonal faces
//       (hexagon inradius = sqrt((47 + 21*sqrt(5))/100) * r = 0.97415421 * r,
//        edge length a = (sqrt(5) - 1)/3 * r = 0.42965690 * r,
//        circumradius R = 1.06453182 * r)
// ============================================================================

float sdTruncatedIcosahedron( vec3 p, float r )
{
    // 1. Fold R^3 into the fundamental Weyl chamber of the icosahedral group I_h (order 120)
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    // 2. Signed distances to the pentagon face (np) and hexagon face (nt) in the chamber
    const vec2 np = vec2(0.52573111, 0.85065081);
    const vec2 nt = vec2(0.35682209, 0.93417236);
    float dp = dot(p.xz, np) - r;
    float dt = dot(p.yz, nt) - 0.97415421 * r; // sqrt((47 + 21*sqrt(5))/100) * r

    // 0-sqrt interior and exterior face-prism exit (cos(dihedral) = sqrt(5)/3 = 0.74535599 -> 0.79465447)
    if( max(dp, dt) <= 0.0 || dt <= 0.79465447 * dp || (dp <= 0.79465447 * dt && p.y >= dt * nt.x) )
        return max(dp, dt);

    // 3. Single-sqrt projection onto the pentagon-hexagon edge (e_pt) or hexagon-hexagon edge (e_tt)
    float x0  = 0.21482845 * r, z0 = 1.04279922 * r;
    float dn3 = dot(p, n3);
    vec3  d   = (p.x <= dn3)
              ? vec3(max(p.x - x0, 0.0), p.y, p.z - z0)
              : p - vec3(x0, 0.0, z0) + max(x0 - dn3, 0.0) * n3;
    return length(d);
}

// ============================================================================
// sdgTruncatedIcosahedron — Exact 0/1-SFU Analytical 3D Distance & Gradient
// Returns vec4(d, grad.x, grad.y, grad.z) with 0 SFUs in interior/face prisms,
// 1 inversesqrt in the exterior edge/vertex cone.
// ============================================================================

vec4 sdgTruncatedIcosahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    vec3 s0 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r0 = (dot(p, n3) < 0.0); if( r0 ) p -= 2.0 * dot(p, n3) * n3;
    vec3 s1 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r1 = (dot(p, n3) < 0.0); if( r1 ) p -= 2.0 * dot(p, n3) * n3;
    vec3 s2 = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p); bool r2 = (dot(p, n3) < 0.0); if( r2 ) p -= 2.0 * dot(p, n3) * n3;
    vec2 s3 = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    p.xy = abs(p.xy);

    const vec2 np = vec2(0.52573111, 0.85065081);
    const vec2 nt = vec2(0.35682209, 0.93417236);
    float dp = dot(p.xz, np) - r;
    float dt = dot(p.yz, nt) - 0.97415421 * r;
    float d;
    vec3  g;
    if( max(dp, dt) <= 0.0 || dt <= 0.79465447 * dp || (dp <= 0.79465447 * dt && p.y >= dt * nt.x) )
    {
        d = max(dp, dt);
        g = (dp >= dt) ? vec3(np.x, 0.0, np.y) : vec3(0.0, nt.x, nt.y);
    }
    else
    {
        float x0   = 0.21482845 * r, z0 = 1.04279922 * r;
        float dn3  = dot(p, n3);
        vec3  v    = (p.x <= dn3)
                   ? vec3(max(p.x - x0, 0.0), p.y, p.z - z0)
                   : p - vec3(x0, 0.0, z0) + max(x0 - dn3, 0.0) * n3;
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = v * invV;
    }

    g.xy *= s3;
    if( r2 ) g -= 2.0 * dot(g, n3) * n3; g *= s2;
    if( r1 ) g -= 2.0 * dot(g, n3) * n3; g *= s1;
    if( r0 ) g -= 2.0 * dot(g, n3) * n3; g *= s0;
    return vec4(d, g);
}
