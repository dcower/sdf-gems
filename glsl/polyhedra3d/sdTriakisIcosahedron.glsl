// ============================================================================
// sdTriakisIcosahedron & sdgTriakisIcosahedron — Exact 0/1-Sqrt 3D Catalan Triakis Icosahedron (60 Faces)
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prism, 1 sqrt in edge/vertex cone)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for the 60-face Catalan
//          triakis icosahedron (dual of the truncated dodecahedron, 60 isosceles triangles,
//          90 edges, 32 vertices) via I_h Weyl chamber folding, x-axis-aligned base edge
//          reduction, and the Golden-Pythagorean identity (phi + 7)^2 + 5^2 + (3*phi - 4)^2 = 10^2.
//
// Parameters:
//   p : 3D query point
//   r : insphere radius tangent to all 60 isosceles triangular faces
// ============================================================================

#ifndef SDF_GEMS_SD_TRIAKIS_ICOSAHEDRON_GLSL
#define SDF_GEMS_SD_TRIAKIS_ICOSAHEDRON_GLSL

float sdTriakisIcosahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    // Single face normal N = (0, 1, 3*phi + 1) / sqrt(11 + 15*phi) in C_{I_h}
    float df  = p.y * 0.16838141 + p.z * 0.98572192 - r;
    float dn3 = dot(p, n3);
    if( df <= 5.93889803 * min(p.y, dn3) )
        return df;

    // Exterior edges: E2 (along x-axis on y = 0) and E3 (on p.n3 = 0 with rational golden tangent)
    float x5 = 0.62698615 * r, z5 = 1.01448490 * r;
    vec3  q3 = p - vec3(x5, 0.0, z5);
    const vec3 e3 = vec3(-0.86180340, 0.5, -0.08541020); // (-(phi + 7), 5, -(3*phi - 4)) / 10
    vec3  d  = (p.y <= dn3)
             ? vec3(max(q3.x, 0.0), p.y, q3.z)
             : q3 - clamp(dot(q3, e3), 0.0, 0.72752805 * r) * e3;
    return length(d);
}

vec4 sdgTriakisIcosahedron( vec3 p, float r )
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

    float df  = p.y * 0.16838141 + p.z * 0.98572192 - r;
    float dn3 = dot(p, n3);
    float d;
    vec3  g;
    if( df <= 5.93889803 * min(p.y, dn3) )
    {
        d = df;
        g = vec3(0.0, 0.16838141, 0.98572192);
    }
    else
    {
        float x5 = 0.62698615 * r, z5 = 1.01448490 * r;
        vec3  q3 = p - vec3(x5, 0.0, z5);
        const vec3 e3 = vec3(-0.86180340, 0.5, -0.08541020);
        vec3  v  = (p.y <= dn3)
                 ? vec3(max(q3.x, 0.0), p.y, q3.z)
                 : q3 - clamp(dot(q3, e3), 0.0, 0.72752805 * r) * e3;
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

#endif // SDF_GEMS_SD_TRIAKIS_ICOSAHEDRON_GLSL
