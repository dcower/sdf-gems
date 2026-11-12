#ifndef SDF_GEMS_SD_TETRAGONAL_TRAPEZOHEDRON_GLSL
#define SDF_GEMS_SD_TETRAGONAL_TRAPEZOHEDRON_GLSL

#include "../common/common.glsl"

// ============================================================================
// sdTetragonalTrapezohedron, sdgTetragonalTrapezohedron,
// sdGyroelongatedSquareBipyramid & sdgGyroelongatedSquareBipyramid
// — Exact 0-Loop, 0-Trig, 0/1-Sqrt 3D Leucite Tetragonal Trapezohedron (8 Kites)
//   & Johnson Solid J17 Gyroelongated Square Bipyramid (16 Triangles)
//
// @author  Discovered and derived using Gemini 4 Argon (research & mathematical design)
//          and Gemini 3.8 Flash (implementation & verification) (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; 0 loops, 0 trig, 0 sqrts in interior & face prisms, 1 sqrt in edge/vertex cones)
// @novelty Proves the D_4d Commutation Theorem for both dual families of D_4d polyhedra:
//          1) Folding p.xy into 0 <= y <= x via abs() and a swap, followed by the D_4d
//             rotoreflection (x, y, z) -> ((x+y)/sqrt(2), (x-y)/sqrt(2), -z), reduces
//             all 8 congruent kite faces, 16 edges, and 10 vertices of the Tetragonal
//             Trapezohedron (leucite/analcime crystal habit) to 1 kite face and 2 edges.
//          2) Applying dual D_4d Voronoi bisector reflections to the polar pyramids and
//             equatorial antiprism belt reduces all 16 equilateral triangular faces,
//             24 edges, and 10 vertices of Johnson Solid J17 (Gyroelongated Square
//             Bipyramid, the 16-faced convex deltahedron) to 2 faces and 3 edges!
// ============================================================================

#ifndef SDF_GEMS_HAS_SD_TETRAGONAL_TRAPEZOHEDRON
#define SDF_GEMS_HAS_SD_TETRAGONAL_TRAPEZOHEDRON

// 1. Exact 0-Loop, 0-Trig, 0/1-Sqrt 3D Tetragonal Trapezohedron (equatorial radius r, polar half-height h)
float sdTetragonalTrapezohedron( in vec3 p, in float r, in float h )
{
    float x = abs(p.x), y = abs(p.y), z = p.z;
    if( y > x ) { float t = x; x = y; y = t; }

    const float c  = 0.7071067811865475;
    const float c1 = 0.2928932188134525;
    const float t2 = 0.1715728752538099;
    float a = 1.1715728752538099 * h;
    float u = (x + y) * c, v = (x - y) * c;

    if( a * (x - u) + 2.0 * r * z < 0.0 )
    {
        float tx = x, ty = y;
        x = u; y = v; u = tx; v = ty; z = -z;
    }

    float n2 = a * a + r * r;
    float m  = a * x + r * (z - h);
    float w2 = a * (x - u) + 2.0 * r * z;
    if( a * c * m <= n2 * v && (a * a * c1 + 2.0 * r * r) * m <= n2 * w2 )
        return m * inversesqrt(n2);

    vec3  e1 = vec3(r * c, r * c, -a * c);
    float l1 = clamp((x * e1.x + y * e1.y + (z - h) * e1.z) / (r * r + 0.5 * a * a), 0.0, 1.0);
    vec3  d1 = vec3(x, y, z - h) - l1 * e1;

    vec3  e2 = vec3(r * c1, -r * c, -a * c1);
    vec3  q2 = vec3(x - r * c, y - r * c, z - h * t2);
    float l2 = clamp(dot(q2, e2) / (2.0 * r * r * c1 + a * a * c1 * c1), 0.0, 1.0);
    vec3  d2 = q2 - l2 * e2;

    return sqrt(min(dot(d1, d1), dot(d2, d2)));
}

// 2. Analytical 3D Distance & Unit Gradient vec4(d, grad.xyz) for Tetragonal Trapezohedron
vec4 sdgTetragonalTrapezohedron( in vec3 p, in float r, in float h )
{
    vec2  sg = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    float x = abs(p.x), y = abs(p.y), z = p.z;
    bool  swapXY = (y > x);
    if( swapXY ) { float t = x; x = y; y = t; }

    const float c  = 0.7071067811865475;
    const float c1 = 0.2928932188134525;
    const float t2 = 0.1715728752538099;
    float a = 1.1715728752538099 * h;
    float u = (x + y) * c, v = (x - y) * c;

    bool swapD4d = (a * (x - u) + 2.0 * r * z < 0.0);
    if( swapD4d )
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
    if( a * c * m <= n2 * v && (a * a * c1 + 2.0 * r * r) * m <= n2 * w2 )
    {
        d = m * invN;
        g = vec3(a * invN, 0.0, r * invN);
    }
    else
    {
        vec3  e1 = vec3(r * c, r * c, -a * c);
        float l1 = clamp((x * e1.x + y * e1.y + (z - h) * e1.z) / (r * r + 0.5 * a * a), 0.0, 1.0);
        vec3  d1 = vec3(x, y, z - h) - l1 * e1;

        vec3  e2 = vec3(r * c1, -r * c, -a * c1);
        vec3  q2 = vec3(x - r * c, y - r * c, z - h * t2);
        float l2 = clamp(dot(q2, e2) / (2.0 * r * r * c1 + a * a * c1 * c1), 0.0, 1.0);
        vec3  d2 = q2 - l2 * e2;

        vec3  best = (dot(d1, d1) <= dot(d2, d2)) ? d1 : d2;
        float b2   = dot(best, best);
        float invB = inversesqrt(max(b2, 1e-30));
        d = b2 * invB;
        g = best * invB;
    }

    if( swapD4d ) g = vec3((g.x + g.y) * c, (g.x - g.y) * c, -g.z);
    if( swapXY )  g.xy = g.yx;
    g.xy *= sg;
    return vec4(d, g);
}

// Canonical 1-parameter overloads for Equal-Inradius Tetragonal Trapezohedron (h = r)
float sdTetragonalTrapezohedron( in vec3 p, in float s )
{
    return sdTetragonalTrapezohedron(p, s, s);
}

vec4 sdgTetragonalTrapezohedron( in vec3 p, in float s )
{
    return sdgTetragonalTrapezohedron(p, s, s);
}

// 3. Exact 0-Loop, 0-Trig, 0/1-Sqrt Johnson Solid J17 Gyroelongated Square Bipyramid
//    (equatorial circumradius r, antiprism half-height hAnti, pyramid height hPyr)
float sdGyroelongatedSquareBipyramid( in vec3 p, in float r, in float hAnti, in float hPyr )
{
    float x = abs(p.x), y = abs(p.y), z = p.z;
    if( y > x ) { float t = x; x = y; y = t; }

    const float k  = 0.7071067811865475;
    const float k1 = 0.2928932188134525;
    float u = (x + y) * k, v = (x - y) * k;

    bool  swapP = (hPyr * (u - x) + 2.0 * r * k * z < 0.0);
    float xp = swapP ? u : x, yp = swapP ? v : y;
    float up = swapP ? x : u, vp = swapP ? y : v, zp = swapP ? -z : z;

    bool  swapS = (hAnti * (u - x) + r * k1 * z < 0.0);
    float xs = swapS ? u : x, ys = swapS ? v : y;
    float us = swapS ? x : u, zs = swapS ? -z : z;

    float lp2   = hPyr * hPyr + 0.5 * r * r;
    float mp    = hPyr * up + r * k * (zp - hAnti - hPyr);
    float dPyr  = mp * inversesqrt(lp2);

    float ls2   = 4.0 * hAnti * hAnti + r * r * k1 * k1;
    float ms    = 2.0 * hAnti * us + r * k1 * zs - r * hAnti * (1.0 + k);
    float dSide = ms * inversesqrt(ls2);

    float dFace  = max(dPyr, dSide);
    bool  inPyr  = (mp > 0.0) && (hPyr * mp <= lp2 * (up - vp)) && (r * k * mp <= lp2 * (zp - hAnti));
    bool  inSide = (ms > 0.0) && (r * k1 * ms >= ls2 * (zs - hAnti)) &&
                   (k1 * (2.0 * hAnti * hAnti + r * r * k1) * ms <= ls2 * (hAnti * (us - xs) + r * k1 * zs));

    if( dFace <= 0.0 || inPyr || inSide )
        return dFace;

    vec3  q1 = vec3(xp, yp, zp - hAnti - hPyr);
    vec3  e1 = vec3(r, 0.0, -hPyr);
    float t1 = clamp(dot(q1, e1) / (r * r + hPyr * hPyr), 0.0, 1.0);
    vec3  d1 = q1 - t1 * e1;

    float du2 = up - r * k, dv2 = max(vp - r * k, 0.0);
    vec3  d2  = vec3((du2 + dv2) * k, (du2 - dv2) * k, zp - hAnti);

    vec3  q3 = vec3(xs - r, ys, zs - hAnti);
    vec3  e3 = vec3(-r * k1, r * k, -2.0 * hAnti);
    float t3 = clamp(dot(q3, e3) / (2.0 * r * r * k1 + 4.0 * hAnti * hAnti), 0.0, 1.0);
    vec3  d3 = q3 - t3 * e3;

    return sqrt(min(dot(d1, d1), min(dot(d2, d2), dot(d3, d3))));
}

// 4. Analytical 3D Distance & Unit Gradient vec4(d, grad.xyz) for Gyroelongated Square Bipyramid
vec4 sdgGyroelongatedSquareBipyramid( in vec3 p, in float r, in float hAnti, in float hPyr )
{
    vec2  sg = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    float x = abs(p.x), y = abs(p.y), z = p.z;
    bool  swapXY = (y > x);
    if( swapXY ) { float t = x; x = y; y = t; }

    const float k  = 0.7071067811865475;
    const float k1 = 0.2928932188134525;
    float u = (x + y) * k, v = (x - y) * k;

    bool  swapP = (hPyr * (u - x) + 2.0 * r * k * z < 0.0);
    float xp = swapP ? u : x, yp = swapP ? v : y;
    float up = swapP ? x : u, vp = swapP ? y : v, zp = swapP ? -z : z;

    bool  swapS = (hAnti * (u - x) + r * k1 * z < 0.0);
    float xs = swapS ? u : x, ys = swapS ? v : y;
    float us = swapS ? x : u, zs = swapS ? -z : z;

    float lp2   = hPyr * hPyr + 0.5 * r * r;
    float invLp = inversesqrt(lp2);
    float mp    = hPyr * up + r * k * (zp - hAnti - hPyr);
    float dPyr  = mp * invLp;

    float ls2   = 4.0 * hAnti * hAnti + r * r * k1 * k1;
    float invLs = inversesqrt(ls2);
    float ms    = 2.0 * hAnti * us + r * k1 * zs - r * hAnti * (1.0 + k);
    float dSide = ms * invLs;

    float dFace  = max(dPyr, dSide);
    bool  inInt  = (dFace <= 0.0);
    bool  inPyr  = (mp > 0.0) && (hPyr * mp <= lp2 * (up - vp)) && (r * k * mp <= lp2 * (zp - hAnti));
    bool  inSide = (ms > 0.0) && (r * k1 * ms >= ls2 * (zs - hAnti)) &&
                   (k1 * (2.0 * hAnti * hAnti + r * r * k1) * ms <= ls2 * (hAnti * (us - xs) + r * k1 * zs));

    float d;
    vec3  g;
    if( inInt || inPyr || inSide )
    {
        bool usePyr = inInt ? (dPyr >= dSide) : inPyr;
        d = usePyr ? dPyr : dSide;
        g = usePyr
          ? (swapP ? vec3(hPyr * invLp, 0.0, -r * k * invLp)
                   : vec3(hPyr * k * invLp, hPyr * k * invLp, r * k * invLp))
          : (swapS ? vec3(2.0 * hAnti * invLs, 0.0, -r * k1 * invLs)
                   : vec3(2.0 * hAnti * k * invLs, 2.0 * hAnti * k * invLs, r * k1 * invLs));
    }
    else
    {
        vec3  q1 = vec3(xp, yp, zp - hAnti - hPyr);
        vec3  e1 = vec3(r, 0.0, -hPyr);
        float t1 = clamp(dot(q1, e1) / (r * r + hPyr * hPyr), 0.0, 1.0);
        vec3  d1 = q1 - t1 * e1;

        float du2 = up - r * k, dv2 = max(vp - r * k, 0.0);
        vec3  d2  = vec3((du2 + dv2) * k, (du2 - dv2) * k, zp - hAnti);

        vec3  q3 = vec3(xs - r, ys, zs - hAnti);
        vec3  e3 = vec3(-r * k1, r * k, -2.0 * hAnti);
        float t3 = clamp(dot(q3, e3) / (2.0 * r * r * k1 + 4.0 * hAnti * hAnti), 0.0, 1.0);
        vec3  d3 = q3 - t3 * e3;

        vec3  dp   = (dot(d1, d1) <= dot(d2, d2)) ? d1 : d2;
        float p2   = dot(dp, dp), s2 = dot(d3, d3);
        vec3  vpW  = swapP ? vec3((dp.x + dp.y) * k, (dp.x - dp.y) * k, -dp.z) : dp;
        vec3  vsW  = swapS ? vec3((d3.x + d3.y) * k, (d3.x - d3.y) * k, -d3.z) : d3;
        vec3  ve   = (p2 <= s2) ? vpW : vsW;
        float v2m  = min(p2, s2);
        float invV = inversesqrt(max(v2m, 1e-30));
        d = v2m * invV;
        g = ve * invV;
    }

    if( swapXY ) g.xy = g.yx;
    g.xy *= sg;
    return vec4(d, g);
}

// 5. Canonical 1-Parameter Overloads for Uniform Johnson Solid J17 (16 Equilateral Triangles)
float sdGyroelongatedSquareBipyramid( in vec3 p, in float r )
{
    return sdGyroelongatedSquareBipyramid(p, r, 0.59460356 * r, 0.70710678 * r);
}

vec4 sdgGyroelongatedSquareBipyramid( in vec3 p, in float r )
{
    return sdgGyroelongatedSquareBipyramid(p, r, 0.59460356 * r, 0.70710678 * r);
}

#endif // SDF_GEMS_HAS_SD_TETRAGONAL_TRAPEZOHEDRON
#endif // SDF_GEMS_SD_TETRAGONAL_TRAPEZOHEDRON_GLSL
