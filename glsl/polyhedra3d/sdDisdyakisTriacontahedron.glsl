// ============================================================================
// sdDisdyakisTriacontahedron & sdgDisdyakisTriacontahedron — Exact 0/1-Sqrt 3D Catalan Disdyakis Triacontahedron (120 Faces)
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prism, 1 sqrt in edge/vertex cone)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for the 120-face Catalan
//          disdyakis triacontahedron (hexakis icosahedron: 120 scalene triangles, 180 edges,
//          62 vertices — the largest Catalan solid) via the Equi-Dihedral Incircle Theorem:
//          all 120 faces collapse to 1 triangle in C_{I_h} with N.w1 = N.w2 = N.w3 = k,
//          so argmin(p.x, p.y, p.n3) 100% partitions the 3 edges!
//
// Parameters:
//   p : 3D query point
//   r : insphere radius tangent to all 120 scalene triangular faces
// ============================================================================

#ifndef SDF_GEMS_SD_DISDYAKIS_TRIACONTAHEDRON_GLSL
#define SDF_GEMS_SD_DISDYAKIS_TRIACONTAHEDRON_GLSL

float sdDisdyakisTriacontahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    // Single face normal N = (1, 1, 4*phi + 1) / sqrt(31 + 12*sqrt(5)) in C_{I_h}
    float df  = (p.x + p.y) * 0.13149609 + p.z * 0.98255664 - r;
    float dn3 = dot(p, n3);
    float v   = min(p.x, p.y);
    if( df <= 7.60478900 * min(v, dn3) )
        return df;

    // Unified E1/E2 projection (when min(p.x, p.y) <= dn3) vs E3 projection (on p.n3 = 0)
    float u   = max(p.x, p.y);
    float dz2 = p.z - 1.01775303 * r;
    float tu  = clamp(0.98240560 * (u - 0.13383054 * dz2), 0.0, (p.x <= p.y) ? 0.36984125 * r : 0.58095418 * r);
    vec3  q3  = p - vec3(0.0, 0.36984125, 0.96825697) * r;
    const vec3 e3 = vec3(0.84285835, -0.53657207, -0.04099046);
    vec3  d   = (v <= dn3)
              ? vec3(u - tu, v, dz2 + 0.13383054 * tu)
              : q3 - clamp(dot(q3, e3), 0.0, 0.68926669 * r) * e3;
    return length(d);
}

vec4 sdgDisdyakisTriacontahedron( vec3 p, float r )
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

    float df  = (p.x + p.y) * 0.13149609 + p.z * 0.98255664 - r;
    float dn3 = dot(p, n3);
    float v_m = min(p.x, p.y);
    float d;
    vec3  g;
    if( df <= 7.60478900 * min(v_m, dn3) )
    {
        d = df;
        g = vec3(0.13149609, 0.13149609, 0.98255664);
    }
    else
    {
        float u   = max(p.x, p.y);
        float dz2 = p.z - 1.01775303 * r;
        float tu  = clamp(0.98240560 * (u - 0.13383054 * dz2), 0.0, (p.x <= p.y) ? 0.36984125 * r : 0.58095418 * r);
        vec3  v12 = (p.x <= p.y) ? vec3(p.x, p.y - tu, dz2 + 0.13383054 * tu)
                                 : vec3(p.x - tu, p.y, dz2 + 0.13383054 * tu);
        vec3  q3  = p - vec3(0.0, 0.36984125, 0.96825697) * r;
        const vec3 e3 = vec3(0.84285835, -0.53657207, -0.04099046);
        vec3  v   = (v_m <= dn3) ? v12 : q3 - clamp(dot(q3, e3), 0.0, 0.68926669 * r) * e3;
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

#endif // SDF_GEMS_SD_DISDYAKIS_TRIACONTAHEDRON_GLSL
