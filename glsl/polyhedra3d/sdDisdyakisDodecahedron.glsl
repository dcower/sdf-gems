// ============================================================================
// sdDisdyakisDodecahedron & sdgDisdyakisDodecahedron — Exact 0/1-Sqrt 3D Catalan Disdyakis Dodecahedron (48 Faces)
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 sqrts in interior & face prism, 1 sqrt in edge/vertex cone)
// @novelty Exact 0/1-sqrt Euclidean SDF and analytical gradient for the 48-face Catalan
//          disdyakis dodecahedron (hexakis octahedron: 48 scalene triangles, 72 edges,
//          26 vertices, dual of the truncated cuboctahedron) via the Equi-Dihedral Incircle
//          Theorem in O_h: all 48 faces collapse to 1 triangle with N.w1 = N.w2 = N.w3 = k,
//          so argmin(x - y, y - z, sqrt(2)*z) 100% partitions the 3 edges!
//
// Parameters:
//   p : 3D query point
//   s : axial half-extent to the 6 eight-valent vertices at (+/-s, 0, 0)
//       (insphere radius r = 0.82594259 * s)
// ============================================================================

#ifndef SDF_GEMS_SD_DISDYAKIS_DODECAHEDRON_GLSL
#define SDF_GEMS_SD_DISDYAKIS_DODECAHEDRON_GLSL

float sdDisdyakisDodecahedron( vec3 p, float s )
{
    p = abs(p);
    if( p.x < p.y ) p.xy = p.yx;
    if( p.y < p.z ) p.yz = p.zy;
    if( p.x < p.y ) p.xy = p.yx;

    float m = p.x + 0.63060194 * p.y + 0.26120387 * p.z - s;
    float w = min(p.x - p.y, p.y - p.z);
    if( m <= min(3.96848604 * w, 5.61228668 * p.z) )
        return m * 0.82594259;

    float t3  = clamp((0.63060194 * (s - p.x) + p.y) * 0.71548220, 0.0, 0.61327046 * s);
    float t1  = clamp((0.16018862 * (1.22654092 * s - p.x - p.y) + p.z) * 0.95118446, 0.0, 0.52859548 * s);
    float t2  = clamp((0.89180581 * (s - p.x) + p.y + p.z) * 0.35774110, 0.0, 0.52859548 * s);
    float xy1 = 0.61327046 * s - 0.16018862 * t1;
    vec3  d   = (1.41421356 * p.z <= w)
              ? vec3(p.x - s + 0.63060194 * t3, p.y - t3, p.z)
              : (p.x - p.y <= p.y - p.z)
              ? vec3(p.x - xy1, p.y - xy1, p.z - t1)
              : vec3(p.x - s + 0.89180581 * t2, p.y - t2, p.z - t2);
    return length(d);
}

vec4 sdgDisdyakisDodecahedron( vec3 p, float s )
{
    vec3 sg = vec3((p.x < 0.0) ? -1.0 : 1.0,
                   (p.y < 0.0) ? -1.0 : 1.0,
                   (p.z < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    bool s1 = (p.x < p.y); if( s1 ) p.xy = p.yx;
    bool s2 = (p.y < p.z); if( s2 ) p.yz = p.zy;
    bool s3 = (p.x < p.y); if( s3 ) p.xy = p.yx;

    float m = p.x + 0.63060194 * p.y + 0.26120387 * p.z - s;
    float w = min(p.x - p.y, p.y - p.z);
    float d;
    vec3  g;
    if( m <= min(3.96848604 * w, 5.61228668 * p.z) )
    {
        d = m * 0.82594259;
        g = vec3(0.82594259, 0.52084106, 0.21573952);
    }
    else
    {
        float t3   = clamp((0.63060194 * (s - p.x) + p.y) * 0.71548220, 0.0, 0.61327046 * s);
        float t1   = clamp((0.16018862 * (1.22654092 * s - p.x - p.y) + p.z) * 0.95118446, 0.0, 0.52859548 * s);
        float t2   = clamp((0.89180581 * (s - p.x) + p.y + p.z) * 0.35774110, 0.0, 0.52859548 * s);
        float xy1  = 0.61327046 * s - 0.16018862 * t1;
        vec3  v    = (1.41421356 * p.z <= w)
                   ? vec3(p.x - s + 0.63060194 * t3, p.y - t3, p.z)
                   : (p.x - p.y <= p.y - p.z)
                   ? vec3(p.x - xy1, p.y - xy1, p.z - t1)
                   : vec3(p.x - s + 0.89180581 * t2, p.y - t2, p.z - t2);
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

#endif // SDF_GEMS_SD_DISDYAKIS_DODECAHEDRON_GLSL
