// ============================================================================
// sdTruncatedIcosidodecahedron & sdgTruncatedIcosidodecahedron — Exact 0/1-Sqrt 3D Archimedean Truncated Icosidodecahedron (62 Faces, 120 Vertices)
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & all 62 face prisms, 1 sqrt in edge/vertex cone)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for the 62-face Archimedean
//          truncated icosidodecahedron (great rhombicosidodecahedron: 30 squares, 20 regular
//          hexagons, 12 regular decagons, 180 edges, 120 vertices) via the Omni-Truncated
//          Weyl-Chamber Ray-Projection Theorem: sqrt(dot(q, q) - u * u).
//
// Parameters:
//   p : 3D query point
//   r : inradius to the 30 axial square faces
// ============================================================================

#ifndef SDF_GEMS_SD_TRUNCATED_ICOSIDODECAHEDRON_GLSL
#define SDF_GEMS_SD_TRUNCATED_ICOSIDODECAHEDRON_GLSL

float sdTruncatedIcosidodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for( int i = 0; i < 3; i++ )
    {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);

    vec3  q  = p - vec3(0.13383054, 0.13383054, 1.0) * r;
    float d2 = q.z;
    float d5 = dot(q.xz, vec2(0.52573111, 0.85065081));
    float d3 = dot(q.yz, vec2(0.35682209, 0.93417236));
    float dm = max(d2, max(d5, d3));
    if( dm <= 0.0 ||
        max(q.x, q.y) <= 0.0 ||
        (d5 > 0.0 && q.z <= 1.61803399 * q.x && d3 <= 0.79465447 * d5) ||
        (d3 > 0.0 && q.z <= 2.61803399 * q.y && d5 <= 0.79465447 * d3) )
        return dm;

    float u = min(min(q.x, q.y), min(dot(q, n3), 0.0));
    return sqrt(dot(q, q) - u * u);
}

vec4 sdgTruncatedIcosidodecahedron( vec3 p, float r )
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

    vec3  q  = p - vec3(0.13383054, 0.13383054, 1.0) * r;
    float d2 = q.z;
    float d5 = dot(q.xz, vec2(0.52573111, 0.85065081));
    float d3 = dot(q.yz, vec2(0.35682209, 0.93417236));
    float dm = max(d2, max(d5, d3));
    float d;
    vec3  g;
    if( dm <= 0.0 ||
        max(q.x, q.y) <= 0.0 ||
        (d5 > 0.0 && q.z <= 1.61803399 * q.x && d3 <= 0.79465447 * d5) ||
        (d3 > 0.0 && q.z <= 2.61803399 * q.y && d5 <= 0.79465447 * d3) )
    {
        d = dm;
        g = (d2 >= d5 && d2 >= d3) ? vec3(0.0, 0.0, 1.0)
          : (d5 >= d3)             ? vec3(0.52573111, 0.0, 0.85065081)
                                   : vec3(0.0, 0.35682209, 0.93417236);
    }
    else
    {
        float qn3  = dot(q, n3);
        float u    = min(min(q.x, q.y), min(qn3, 0.0));
        vec3  w    = (q.x <= q.y && q.x <= qn3) ? vec3(1.0, 0.0, 0.0)
                   : (q.y <= qn3)               ? vec3(0.0, 1.0, 0.0)
                                                : n3;
        vec3  v    = q - u * w;
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

#endif // SDF_GEMS_SD_TRUNCATED_ICOSIDODECAHEDRON_GLSL
