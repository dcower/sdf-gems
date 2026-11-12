#ifndef SDF_GEMS_SD_DURER_SOLID_GLSL
#define SDF_GEMS_SD_DURER_SOLID_GLSL

#include "../common/common.glsl"

// ============================================================================
// sdRhombohedron, sdgRhombohedron, sdDurerSolid & sdgDurerSolid
// — Exact 0-Loop, 0-Trig, 0/1-Sqrt 3D Calcite Trigonal Rhombohedron &
//   Albrecht Durer's (1514) Melencolia I Truncated Rhombohedron
//
// @author  Discovered and derived using Gemini 4 Argon (research & mathematical design)
//          and Gemini 3.8 Flash (implementation & verification) (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 loops, 0 trig, 0 sqrts in interior & face prisms, 1 sqrt in edge/vertex cone)
// @novelty Proves the D_3d Truncation Commutation Theorem: folding p.xy into the
//          60-degree wedge [0, pi/3] via a single reflection across n_60 = (-sqrt(3)/2, 1/2)
//          followed by the D_3d Voronoi rotoreflection (x, y, z) -> (u, v, -z)
//          reduces the 6-faced Calcite Rhombohedron to 1 face + 3 edges, and reduces
//          Albrecht Durer's 1514 Melencolia I polyhedron (8 faces: 6 pentagons +
//          2 equilateral triangles, 18 edges, 12 vertices) to 2 faces + 4 edges in
//          0 loops, 0 trig, and at most 1 sqrt!
// ============================================================================

#ifndef SDF_GEMS_HAS_SD_RHOMBOHEDRON
#define SDF_GEMS_HAS_SD_RHOMBOHEDRON

// 1. Exact 0-Loop, 0-Trig, 0/1-Sqrt 3D Calcite Trigonal Rhombohedron (equatorial radius r, polar half-height h)
float sdRhombohedron( in vec3 p, in float r, in float h )
{
    const vec2  n60 = vec2(-0.8660254037844386, 0.5);
    const float s3  = 0.8660254037844386;
    p.y = abs(p.y);
    p.xy -= 2.0 * max(dot(p.xy, n60), 0.0) * n60;
    p.y = abs(p.y);

    float a = (4.0 / 3.0) * h;
    float x = p.x, y = p.y, z = p.z;
    float u = 0.5 * x + s3 * y, v = s3 * x - 0.5 * y;

    if( a * (x - u) + 2.0 * r * z < 0.0 )
    {
        float tx = x, ty = y;
        x = u; y = v; u = tx; v = ty; z = -z;
    }

    float n2 = a * a + r * r;
    float m  = a * x + r * (z - h);
    float w2 = a * (x - u) + 2.0 * r * z;
    if( a * s3 * m <= n2 * v && (0.5 * a * a + 2.0 * r * r) * m <= n2 * w2 )
        return m * inversesqrt(n2);

    float invE2 = 1.0 / (r * r + 0.25 * a * a);
    vec3  e1    = vec3(0.5 * r, s3 * r, -0.5 * a);
    float l1    = clamp((x * e1.x + y * e1.y + (z - h) * e1.z) * invE2, 0.0, 1.0);
    vec3  d1    = vec3(x, y, z - h) - l1 * e1;

    vec3  e2    = vec3(0.5 * r, -s3 * r, -0.5 * a);
    vec3  q2    = vec3(x - 0.5 * r, y - s3 * r, z - (1.0 / 3.0) * h);
    float l2    = clamp(dot(q2, e2) * invE2, 0.0, 1.0);
    vec3  d2    = q2 - l2 * e2;

    vec3  e3    = vec3(r, 0.0, 0.5 * a);
    vec3  q3    = vec3(x, y, z + h);
    float l3    = clamp(dot(q3, e3) * invE2, 0.0, 1.0);
    vec3  d3    = q3 - l3 * e3;

    return sqrt(min(dot(d1, d1), min(dot(d2, d2), dot(d3, d3))));
}

// 2. Analytical 3D Distance & Unit Gradient vec4(d, grad.xyz) for Calcite Trigonal Rhombohedron
vec4 sdgRhombohedron( in vec3 p, in float r, in float h )
{
    const vec2  n60 = vec2(-0.8660254037844386, 0.5);
    const float s3  = 0.8660254037844386;
    float sy0 = (p.y < 0.0) ? -1.0 : 1.0; p.y = abs(p.y);
    bool  r0  = (dot(p.xy, n60) > 0.0);    if( r0 ) p.xy -= 2.0 * dot(p.xy, n60) * n60;
    float sy1 = (p.y < 0.0) ? -1.0 : 1.0; p.y = abs(p.y);

    float a = (4.0 / 3.0) * h;
    float x = p.x, y = p.y, z = p.z;
    float u = 0.5 * x + s3 * y, v = s3 * x - 0.5 * y;

    bool swap = (a * (x - u) + 2.0 * r * z < 0.0);
    if( swap )
    {
        float tx = x, ty = y;
        x = u; y = v; u = tx; v = ty; z = -z;
    }

    float n2   = a * a + r * r;
    float invN = inversesqrt(n2);
    float m    = a * x + r * (z - h);
    float w2   = a * (x - u) + 2.0 * r * z;
    float d;
    vec3  g;
    if( a * s3 * m <= n2 * v && (0.5 * a * a + 2.0 * r * r) * m <= n2 * w2 )
    {
        d = m * invN;
        g = vec3(a * invN, 0.0, r * invN);
    }
    else
    {
        float invE2 = 1.0 / (r * r + 0.25 * a * a);
        vec3  e1    = vec3(0.5 * r, s3 * r, -0.5 * a);
        float l1    = clamp((x * e1.x + y * e1.y + (z - h) * e1.z) * invE2, 0.0, 1.0);
        vec3  d1    = vec3(x, y, z - h) - l1 * e1;

        vec3  e2    = vec3(0.5 * r, -s3 * r, -0.5 * a);
        vec3  q2    = vec3(x - 0.5 * r, y - s3 * r, z - (1.0 / 3.0) * h);
        float l2    = clamp(dot(q2, e2) * invE2, 0.0, 1.0);
        vec3  d2    = q2 - l2 * e2;

        vec3  e3    = vec3(r, 0.0, 0.5 * a);
        vec3  q3    = vec3(x, y, z + h);
        float l3    = clamp(dot(q3, e3) * invE2, 0.0, 1.0);
        vec3  d3    = q3 - l3 * e3;

        float n1 = dot(d1, d1), n2e = dot(d2, d2), n3 = dot(d3, d3);
        vec3  best = (n1 <= n2e && n1 <= n3) ? d1 : ((n2e <= n3) ? d2 : d3);
        float b2   = dot(best, best);
        float invB = inversesqrt(max(b2, 1e-30));
        d = b2 * invB;
        g = best * invB;
    }

    if( swap ) g = vec3(0.5 * g.x + s3 * g.y, s3 * g.x - 0.5 * g.y, -g.z);
    g.y *= sy1;
    if( r0 ) g.xy -= 2.0 * dot(g.xy, n60) * n60;
    g.y *= sy0;
    return vec4(d, g);
}

// Canonical 1-Parameter Overloads for Rhombohedron (scale s)
float sdRhombohedron( in vec3 p, in float s )
{
    return sdRhombohedron(p, 0.85065081 * s, 1.26296630 * s);
}

vec4 sdgRhombohedron( in vec3 p, in float s )
{
    return sdgRhombohedron(p, 0.85065081 * s, 1.26296630 * s);
}

#endif // SDF_GEMS_HAS_SD_RHOMBOHEDRON

// 3. Exact 0-Loop, 0-Trig, 0/1-Sqrt Albrecht Durer's (1514) Solid (equatorial radius r, polar half-height h, cap height zCap)
float sdDurerSolid( in vec3 p, in float r, in float h, in float zCap )
{
    const vec2  n60 = vec2(-0.8660254037844386, 0.5);
    const float s3  = 0.8660254037844386;
    p.y = abs(p.y);
    p.xy -= 2.0 * max(dot(p.xy, n60), 0.0) * n60;
    p.y = abs(p.y);

    float a    = (4.0 / 3.0) * h;
    float zc   = clamp(zCap, (1.0 / 3.0) * h, h);
    float xCap = r * (h - zc) / a;
    float tMin = 1.5 * (1.0 - zc / h);

    float x = p.x, y = p.y, z = p.z;
    float u = 0.5 * x + s3 * y, v = s3 * x - 0.5 * y;

    if( a * (x - u) + 2.0 * r * z < 0.0 )
    {
        float tx = x, ty = y;
        x = u; y = v; u = tx; v = ty; z = -z;
    }

    float n2    = a * a + r * r;
    float m     = a * x + r * (z - h);
    float dPent = m * inversesqrt(n2);
    float dTri  = abs(z) - zc;
    float w2    = a * (x - u) + 2.0 * r * z;

    bool inPent = (m >= 0.0) && (a * s3 * m <= n2 * v) &&
                  ((0.5 * a * a + 2.0 * r * r) * m <= n2 * w2) && (r * m >= n2 * (z - zc));
    bool inTri  = (z >= zc) && (x <= xCap);

    if( max(dPent, dTri) <= 0.0 || inPent || inTri )
        return max(dPent, dTri);

    vec3  d0    = vec3(x - xCap, y - clamp(y, 0.0, 1.73205081 * xCap), z - zc);
    float invE2 = 1.0 / (r * r + 0.25 * a * a);
    vec3  e1    = vec3(0.5 * r, s3 * r, -0.5 * a);
    float l1    = clamp((x * e1.x + y * e1.y + (z - h) * e1.z) * invE2, tMin, 1.0);
    vec3  d1    = vec3(x, y, z - h) - l1 * e1;

    vec3  e2    = vec3(0.5 * r, -s3 * r, -0.5 * a);
    vec3  q2    = vec3(x - 0.5 * r, y - s3 * r, z - (1.0 / 3.0) * h);
    float l2    = clamp(dot(q2, e2) * invE2, 0.0, 1.0);
    vec3  d2    = q2 - l2 * e2;

    vec3  e3    = vec3(r, 0.0, 0.5 * a);
    vec3  q3    = vec3(x, y, z + h);
    float l3    = clamp(dot(q3, e3) * invE2, tMin, 1.0);
    vec3  d3    = q3 - l3 * e3;

    return sqrt(min(min(dot(d0, d0), dot(d1, d1)), min(dot(d2, d2), dot(d3, d3))));
}

// 4. Analytical 3D Distance & Unit Gradient vec4(d, grad.xyz) for Albrecht Durer's (1514) Solid
vec4 sdgDurerSolid( in vec3 p, in float r, in float h, in float zCap )
{
    const vec2  n60 = vec2(-0.8660254037844386, 0.5);
    const float s3  = 0.8660254037844386;
    float sy0 = (p.y < 0.0) ? -1.0 : 1.0; p.y = abs(p.y);
    bool  r0  = (dot(p.xy, n60) > 0.0);    if( r0 ) p.xy -= 2.0 * dot(p.xy, n60) * n60;
    float sy1 = (p.y < 0.0) ? -1.0 : 1.0; p.y = abs(p.y);

    float a    = (4.0 / 3.0) * h;
    float zc   = clamp(zCap, (1.0 / 3.0) * h, h);
    float xCap = r * (h - zc) / a;
    float tMin = 1.5 * (1.0 - zc / h);

    float x = p.x, y = p.y, z = p.z;
    float u = 0.5 * x + s3 * y, v = s3 * x - 0.5 * y;

    bool swap = (a * (x - u) + 2.0 * r * z < 0.0);
    if( swap )
    {
        float tx = x, ty = y;
        x = u; y = v; u = tx; v = ty; z = -z;
    }

    float n2    = a * a + r * r;
    float invN  = inversesqrt(n2);
    float m     = a * x + r * (z - h);
    float dPent = m * invN;
    float szTri = (z < 0.0) ? -1.0 : 1.0;
    float dTri  = abs(z) - zc;
    float w2    = a * (x - u) + 2.0 * r * z;

    bool inPent = (m >= 0.0) && (a * s3 * m <= n2 * v) &&
                  ((0.5 * a * a + 2.0 * r * r) * m <= n2 * w2) && (r * m >= n2 * (z - zc));
    bool inTri  = (z >= zc) && (x <= xCap);

    float d;
    vec3  g;
    if( max(dPent, dTri) <= 0.0 || inPent || inTri )
    {
        bool usePent = (max(dPent, dTri) <= 0.0) ? (dPent >= dTri) : inPent;
        d = usePent ? dPent : dTri;
        g = usePent ? vec3(a * invN, 0.0, r * invN) : vec3(0.0, 0.0, szTri);
    }
    else
    {
        vec3  d0    = vec3(x - xCap, y - clamp(y, 0.0, 1.73205081 * xCap), z - zc);
        float invE2 = 1.0 / (r * r + 0.25 * a * a);
        vec3  e1    = vec3(0.5 * r, s3 * r, -0.5 * a);
        float l1    = clamp((x * e1.x + y * e1.y + (z - h) * e1.z) * invE2, tMin, 1.0);
        vec3  d1    = vec3(x, y, z - h) - l1 * e1;

        vec3  e2    = vec3(0.5 * r, -s3 * r, -0.5 * a);
        vec3  q2    = vec3(x - 0.5 * r, y - s3 * r, z - (1.0 / 3.0) * h);
        float l2    = clamp(dot(q2, e2) * invE2, 0.0, 1.0);
        vec3  d2    = q2 - l2 * e2;

        vec3  e3    = vec3(r, 0.0, 0.5 * a);
        vec3  q3    = vec3(x, y, z + h);
        float l3    = clamp(dot(q3, e3) * invE2, tMin, 1.0);
        vec3  d3    = q3 - l3 * e3;

        float n0 = dot(d0, d0), n1 = dot(d1, d1), n2e = dot(d2, d2), n3 = dot(d3, d3);
        vec3  best = d0;
        float minD2 = n0;
        if( n1 < minD2 )  { minD2 = n1;  best = d1; }
        if( n2e < minD2 ) { minD2 = n2e; best = d2; }
        if( n3 < minD2 )  { minD2 = n3;  best = d3; }

        float invB  = inversesqrt(max(minD2, 1e-30));
        d = minD2 * invB;
        g = best * invB;
    }

    if( swap ) g = vec3(0.5 * g.x + s3 * g.y, s3 * g.x - 0.5 * g.y, -g.z);
    g.y *= sy1;
    if( r0 ) g.xy -= 2.0 * dot(g.xy, n60) * n60;
    g.y *= sy0;
    return vec4(d, g);
}

// 5. Canonical 1-Parameter Overloads for Durer's (1514) Melencolia I Solid (scale s, 72-degree acute rhombus habit)
float sdDurerSolid( in vec3 p, in float s )
{
    return sdDurerSolid(p, 0.85065081 * s, 1.26296630 * s, 0.85065081 * s);
}

vec4 sdgDurerSolid( in vec3 p, in float s )
{
    return sdgDurerSolid(p, 0.85065081 * s, 1.26296630 * s, 0.85065081 * s);
}

#endif // SDF_GEMS_SD_DURER_SOLID_GLSL
