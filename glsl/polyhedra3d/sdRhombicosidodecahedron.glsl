// ============================================================================
// sdRhombicosidodecahedron & sdgRhombicosidodecahedron — Exact 0/1-Sqrt 3D Archimedean Rhombicosidodecahedron (62 Faces)
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & all 62 face prisms, 1 sqrt in edge/vertex cone)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for the 62-face Archimedean
//          rhombicosidodecahedron (30 squares, 12 regular pentagons, 20 equilateral triangles,
//          120 edges, 60 vertices). By the Vertex-Centered Golden Wedge Theorem, shifting to
//          the chamber vertex v0 = (phi^-3, phi^-3, 1)*r makes both edges axis-aligned along
//          x and y, so the exact distance to all 120 edges and 60 vertices is literally
//          length(vec3(max(q.x, 0.0), max(q.y, 0.0), q.z)) with 0 dot products!
//
// Parameters:
//   p : 3D query point
//   r : inradius to the 30 axial square faces
//       (pentagon inradius = 0.97475909 * r, triangle inradius = 1.01840663 * r,
//        edge length a = 2*(sqrt(5) - 2)*r = 0.47213595 * r, circumradius R = 1.05425622 * r)
// ============================================================================

#ifndef SDF_GEMS_SD_RHOMBICOSIDODECAHEDRON_GLSL
#define SDF_GEMS_SD_RHOMBICOSIDODECAHEDRON_GLSL

float sdRhombicosidodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    vec3  q  = p - vec3(0.23606798, 0.23606798, 1.0) * r;
    float d5 = dot(q.xz, vec2(0.52573111, 0.85065081));
    float d3 = dot(q.yz, vec2(0.35682209, 0.93417236));
    if( q.z <= max(1.61803399 * q.x, 2.61803399 * q.y) || max(q.x, q.y) <= 0.0 )
        return max(q.z, max(d5, d3));

    return length(vec3(max(q.x, 0.0), max(q.y, 0.0), q.z));
}

vec4 sdgRhombicosidodecahedron( vec3 p, float r )
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

    vec3  q  = p - vec3(0.23606798, 0.23606798, 1.0) * r;
    float d5 = dot(q.xz, vec2(0.52573111, 0.85065081));
    float d3 = dot(q.yz, vec2(0.35682209, 0.93417236));
    float d;
    vec3  g;
    if( q.z <= max(1.61803399 * q.x, 2.61803399 * q.y) || max(q.x, q.y) <= 0.0 )
    {
        d = max(q.z, max(d5, d3));
        g = (q.z >= d5 && q.z >= d3) ? vec3(0.0, 0.0, 1.0)
          : (d5 >= d3)               ? vec3(0.52573111, 0.0, 0.85065081)
                                     : vec3(0.0, 0.35682209, 0.93417236);
    }
    else
    {
        vec3  v    = vec3(max(q.x, 0.0), max(q.y, 0.0), q.z);
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

#endif // SDF_GEMS_SD_RHOMBICOSIDODECAHEDRON_GLSL
