// ============================================================================
// sdTruncatedCuboctahedron & sdgTruncatedCuboctahedron — Exact 0/1-Sqrt 3D Archimedean Truncated Cuboctahedron (26 Faces, 48 Vertices)
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & all 26 face prisms, 1 sqrt in edge/vertex cone)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for the 26-face Archimedean
//          truncated cuboctahedron (great rhombicuboctahedron: 6 regular octagons, 8 regular
//          hexagons, 12 squares, 72 edges, 48 vertices) via the Omni-Truncated Weyl-Chamber
//          Ray-Projection Theorem in O_h: sqrt(dot(q, q) - u * u).
//
// Parameters:
//   p : 3D query point
//   s : inradius / half-extent to the 6 axial octagon faces
// ============================================================================

#ifndef SDF_GEMS_SD_TRUNCATED_CUBOCTAHEDRON_GLSL
#define SDF_GEMS_SD_TRUNCATED_CUBOCTAHEDRON_GLSL

float sdTruncatedCuboctahedron( vec3 p, float s )
{
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;

    vec3  q  = p - vec3(1.0, 0.63060194, 0.26120387) * s;
    float d1 = q.x;
    float d2 = (q.x + q.y) * 0.70710678;
    float d3 = (q.x + q.y + q.z) * 0.57735027;
    float dm = max(d1, max(d2, d3));
    if( dm <= 0.0 ||
        (d1 > 0.0 && q.y <= 0.0 && q.y + q.z <= 0.0) ||
        (d2 > 0.0 && q.x <= q.y && q.z <= 0.0) ||
        (d3 > 0.0 && q.x + q.y <= 2.0 * q.z && 2.0 * q.x <= q.y + q.z) )
        return dm;

    float u = min(0.70710678 * min(q.x - q.y, q.y - q.z), min(q.z, 0.0));
    return sqrt(dot(q, q) - u * u);
}

vec4 sdgTruncatedCuboctahedron( vec3 p, float s )
{
    vec3 sg = vec3((p.x < 0.0) ? -1.0 : 1.0,
                   (p.y < 0.0) ? -1.0 : 1.0,
                   (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) p.xy = p.yx;
    bool s2 = (p.y < p.z); if( s2 ) p.yz = p.zy;
    bool s3 = (p.x < p.y); if( s3 ) p.xy = p.yx;

    vec3  q  = p - vec3(1.0, 0.63060194, 0.26120387) * s;
    float d1 = q.x;
    float d2 = (q.x + q.y) * 0.70710678;
    float d3 = (q.x + q.y + q.z) * 0.57735027;
    float dm = max(d1, max(d2, d3));
    float d;
    vec3  g;
    if( dm <= 0.0 ||
        (d1 > 0.0 && q.y <= 0.0 && q.y + q.z <= 0.0) ||
        (d2 > 0.0 && q.x <= q.y && q.z <= 0.0) ||
        (d3 > 0.0 && q.x + q.y <= 2.0 * q.z && 2.0 * q.x <= q.y + q.z) )
    {
        d = dm;
        g = (d1 >= d2 && d1 >= d3) ? vec3(1.0, 0.0, 0.0)
          : (d2 >= d3)             ? vec3(0.70710678, 0.70710678, 0.0)
                                   : vec3(0.57735027);
    }
    else
    {
        float w1   = (q.x - q.y) * 0.70710678;
        float w2   = (q.y - q.z) * 0.70710678;
        float w3   = q.z;
        float u    = min(min(w1, w2), min(w3, 0.0));
        vec3  w    = (w1 <= w2 && w1 <= w3) ? vec3(0.70710678, -0.70710678, 0.0)
                   : (w2 <= w3)             ? vec3(0.0, 0.70710678, -0.70710678)
                                            : vec3(0.0, 0.0, 1.0);
        vec3  v    = q - u * w;
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

#endif // SDF_GEMS_SD_TRUNCATED_CUBOCTAHEDRON_GLSL
