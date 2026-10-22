// ============================================================================
// sdDeltoidalHexecontahedron & sdgDeltoidalHexecontahedron — Exact 0/1-Sqrt 3D Catalan Deltoidal Hexecontahedron (60 Kites)
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prism, 1 sqrt in edge/vertex cone)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for the 60-kite Catalan
//          deltoidal hexecontahedron (dual of the rhombicosidodecahedron, 60 kites, 120 edges,
//          62 vertices). Because N = (1, 1, phi^3)/sqrt(8*phi + 7) lies on W3 (N.n3 = 0),
//          both chamber edges share vertex V2 = (0, 0, z2) and slope (phi^3, -1) with
//          phi^6 + 1 = 2*sqrt(5)*phi^3, reducing both edges to 1 unified clamp!
//
// Parameters:
//   p : 3D query point
//   r : insphere radius tangent to all 60 kite faces
// ============================================================================

#ifndef SDF_GEMS_SD_DELTOIDAL_HEXECONTAHEDRON_GLSL
#define SDF_GEMS_SD_DELTOIDAL_HEXECONTAHEDRON_GLSL

float sdDeltoidalHexecontahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    float df = (p.x + p.y) * 0.22391898 + p.z * 0.94853602 - r;
    float u  = max(p.x, p.y), v = min(p.x, p.y);
    if( df <= 4.46590102 * v )
        return df;

    float dz = p.z - 1.05425622 * r;
    float tu = clamp(0.94721360 * (u - 0.23606798 * dz), 0.0, (p.x <= p.y) ? 0.36938282 * r : 0.56860747 * r);
    return length(vec3(u - tu, v, dz + 0.23606798 * tu));
}

vec4 sdgDeltoidalHexecontahedron( vec3 p, float r )
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

    float df  = (p.x + p.y) * 0.22391898 + p.z * 0.94853602 - r;
    float u   = max(p.x, p.y), v_m = min(p.x, p.y);
    float d;
    vec3  g;
    if( df <= 4.46590102 * v_m )
    {
        d = df;
        g = vec3(0.22391898, 0.22391898, 0.94853602);
    }
    else
    {
        float dz   = p.z - 1.05425622 * r;
        float tu   = clamp(0.94721360 * (u - 0.23606798 * dz), 0.0, (p.x <= p.y) ? 0.36938282 * r : 0.56860747 * r);
        vec3  v    = (p.x <= p.y) ? vec3(p.x, p.y - tu, dz + 0.23606798 * tu)
                                  : vec3(p.x - tu, p.y, dz + 0.23606798 * tu);
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

#endif // SDF_GEMS_SD_DELTOIDAL_HEXECONTAHEDRON_GLSL
