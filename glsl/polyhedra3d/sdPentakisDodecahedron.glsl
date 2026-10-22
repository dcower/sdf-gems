// ============================================================================
// sdPentakisDodecahedron & sdgPentakisDodecahedron — Exact 0/1-Sqrt 3D Catalan Pentakis Dodecahedron (60 Faces)
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prism, 1 sqrt in edge/vertex cone)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for the 60-face Catalan
//          pentakis dodecahedron (geodesic dome / dual of the C60 truncated icosahedron,
//          60 isosceles triangles, 90 edges, 32 vertices) via I_h Weyl chamber folding,
//          y-axis-aligned base edge reduction, and the Golden-Pythagorean identity
//          (3*phi)^2 + (5 - phi)^2 + 1^2 = 6^2 = 36.
//
// Parameters:
//   p : 3D query point
//   r : insphere radius tangent to all 60 isosceles triangular faces
// ============================================================================

#ifndef SDF_GEMS_SD_PENTAKIS_DODECAHEDRON_GLSL
#define SDF_GEMS_SD_PENTAKIS_DODECAHEDRON_GLSL

float sdPentakisDodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    // Single face normal N = (1, 0, 3*phi) / sqrt(10 + 9*phi) in C_{I_h}
    float df  = p.x * 0.20177411 + p.z * 0.97943209 - r;
    float dn3 = dot(p, n3);
    if( df <= 4.95603732 * min(p.x, dn3) )
        return df;

    // Exterior edges: E1 (along y-axis on x = 0) and E3 (on p.n3 = 0 with rational golden tangent)
    float y3 = 0.38998724 * r, z3 = 1.02099984 * r;
    vec3  q3 = p - vec3(0.0, y3, z3);
    const vec3 e3 = vec3(0.80901699, -0.56366100, -0.16666667); // (3*phi, -(5 - phi), -1) / 6
    vec3  d  = (p.x <= dn3)
             ? vec3(p.x, max(q3.y, 0.0), q3.z)
             : q3 - clamp(dot(q3, e3), 0.0, 0.69188259 * r) * e3;
    return length(d);
}

vec4 sdgPentakisDodecahedron( vec3 p, float r )
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

    float df  = p.x * 0.20177411 + p.z * 0.97943209 - r;
    float dn3 = dot(p, n3);
    float d;
    vec3  g;
    if( df <= 4.95603732 * min(p.x, dn3) )
    {
        d = df;
        g = vec3(0.20177411, 0.0, 0.97943209);
    }
    else
    {
        float y3 = 0.38998724 * r, z3 = 1.02099984 * r;
        vec3  q3 = p - vec3(0.0, y3, z3);
        const vec3 e3 = vec3(0.80901699, -0.56366100, -0.16666667);
        vec3  v  = (p.x <= dn3)
                 ? vec3(p.x, max(q3.y, 0.0), q3.z)
                 : q3 - clamp(dot(q3, e3), 0.0, 0.69188259 * r) * e3;
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

#endif // SDF_GEMS_SD_PENTAKIS_DODECAHEDRON_GLSL
