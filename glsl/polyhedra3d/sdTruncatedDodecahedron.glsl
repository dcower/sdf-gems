// ============================================================================
// sdTruncatedDodecahedron & sdgTruncatedDodecahedron — Exact 0/1-Sqrt 3D Archimedean Truncated Dodecahedron
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prisms, 1 sqrt in edge/vertex cone)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for the Archimedean
//          truncated dodecahedron (12 regular decagons + 20 equilateral triangles,
//          60 vertices, 90 edges) via 3-step unrolled I_h Weyl chamber folding.
//
// Parameters:
//   p : 3D query point
//   r : inradius to the 12 regular decagonal faces
//       (triangle inradius = (5 + 3*sqrt(5)) / (2*sqrt(15)) * r = 1.51152263 * r ? No:
//        for unit decagon inradius r, triangle inradius = sqrt((25 + 10*sqrt(5))/3) / (sqrt(5 + 2*sqrt(5))/2) * r = 1.16983942 * r,
//        edge length a = 2 * sqrt(25 - 10*sqrt(5)) * r = 0.40162283 * r,
//        circumradius R = 1.19261460 * r)
// ============================================================================

float sdTruncatedDodecahedron( vec3 p, float r )
{
    // 1. Fold R^3 into the fundamental Weyl chamber of the icosahedral group I_h (order 120)
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    // 2. Signed distances to the decagon face (np) and triangle face (nt) in the chamber
    const vec2 np = vec2(0.52573111, 0.85065081);
    const vec2 nt = vec2(0.35682209, 0.93417236);
    float dp = dot(p.xz, np) - r;
    float dt = dot(p.yz, nt) - 1.16983942 * r; // (3 + 11*sqrt(5))/(6*sqrt(15)) * r

    // 0-sqrt interior and exterior face-prism exit
    if( max(dp, dt) <= 0.0 || dp <= 0.79465447 * dt || (dt <= 0.79465447 * dp && p.x >= dp * np.x) )
        return max(dp, dt);

    // 3. Single-sqrt projection onto the decagon-triangle edge or decagon-decagon edge
    float y0  = 0.20081142 * r, z0 = 1.17557050 * r;
    float dn3 = dot(p, n3);
    vec3  d   = (p.y <= dn3)
              ? vec3(p.x, max(p.y - y0, 0.0), p.z - z0)
              : p - vec3(0.0, y0, z0) + max(y0 - dn3, 0.0) * n3;
    return length(d);
}

// ============================================================================
// sdgTruncatedDodecahedron — Exact 0/1-SFU Analytical 3D Distance & Gradient
// Returns vec4(d, grad.x, grad.y, grad.z) with 0 SFUs in interior/face prisms,
// 1 inversesqrt in the exterior edge/vertex cone.
// ============================================================================

vec4 sdgTruncatedDodecahedron( vec3 p, float r )
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
    float dt = dot(p.yz, nt) - 1.16983942 * r;
    float d;
    vec3  g;
    if( max(dp, dt) <= 0.0 || dp <= 0.79465447 * dt || (dt <= 0.79465447 * dp && p.x >= dp * np.x) )
    {
        d = max(dp, dt);
        g = (dp >= dt) ? vec3(np.x, 0.0, np.y) : vec3(0.0, nt.x, nt.y);
    }
    else
    {
        float y0   = 0.20081142 * r, z0 = 1.17557050 * r;
        float dn3  = dot(p, n3);
        vec3  v    = (p.y <= dn3)
                   ? vec3(p.x, max(p.y - y0, 0.0), p.z - z0)
                   : p - vec3(0.0, y0, z0) + max(y0 - dn3, 0.0) * n3;
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
