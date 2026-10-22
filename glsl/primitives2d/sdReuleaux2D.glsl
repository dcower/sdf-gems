// ============================================================================
// sdReuleauxTriangle2D, sdgReuleauxTriangle2D, sdReuleauxPentagon2D,
// sdgReuleauxPentagon2D & sdReuleauxPolygon2D — Exact 1-Sqrt Constant-Width 2D Reuleaux Polygons
//
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (0 trig, 0 divs, 1 sqrt across 100% of R^2 for N = 3 and N = 5)
// ============================================================================

#ifndef SDF_GEMS_SD_REULEAUX_2D_GLSL
#define SDF_GEMS_SD_REULEAUX_2D_GLSL

// ----------------------------------------------------------------------------
// sdReuleauxTriangle2D(p, w) — Exact 0-Trig, 1-Sqrt 2D Reuleaux Triangle
//
// A Reuleaux triangle of constant width w is the intersection of 3 disks of
// radius w centered at the vertices of an equilateral triangle of edge length w
// and circumradius w / sqrt(3).
//
// Naive CSG max(d0, max(d1, d2)) evaluates 3 sqrts and underestimates exterior
// corner distances by 1 - cos(pi/6) = 13.40%, producing pointed 120-degree
// corners when offset!
//
// By folding p into the 6-fold D_3 Weyl chamber via 1 abs + 1 reflection:
//   p.x = abs(p.x);
//   if( p.x + k * p.y > 0.0 ) p = vec2(abs(p.x - k * p.y), -k * p.x - p.y) * 0.5;
// all 3 circular arcs collapse to 1 arc centered at V0 = (0, w / sqrt(3)) and
// all 3 vertices collapse to 1 corner V1 = (w / 2, -w / (2 * sqrt(3))),
// separated by the exact normal ray sqrt(3)*p.x + p.y > w / sqrt(3).
// Evaluates in 0 trig, 0 divisions, and 1 sqrt across 100% of R^2!
// ----------------------------------------------------------------------------
float sdReuleauxTriangle2D( vec2 p, float w )
{
    const float k = 1.73205081, invK = 0.57735027;
    p.x = abs(p.x);
    if( p.x + k * p.y > 0.0 ) p = vec2(abs(p.x - k * p.y), -k * p.x - p.y) * 0.5;
    bool c = (k * p.x + p.y > w * invK);
    return length(p - (c ? vec2(0.5, -0.5 * invK) : vec2(0.0, invK)) * w) - (c ? 0.0 : w);
}

// ----------------------------------------------------------------------------
// sdgReuleauxTriangle2D(p, w) — Exact 0-Trig, 1-SFU Distance + 2D Unit Gradient
// ----------------------------------------------------------------------------
vec3 sdgReuleauxTriangle2D( vec2 p, float w )
{
    const float k = 1.73205081, invK = 0.57735027;
    float sx0 = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    bool rf = (p.x + k * p.y > 0.0);
    if( rf ) p = vec2(p.x - k * p.y, -k * p.x - p.y) * 0.5;
    float sx1 = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);

    bool  c    = (k * p.x + p.y > w * invK);
    vec2  v    = p - (c ? vec2(0.5, -0.5 * invK) : vec2(0.0, invK)) * w;
    float v2   = dot(v, v);
    float invV = inversesqrt(max(v2, 1e-30));
    float d    = v2 * invV - (c ? 0.0 : w);
    vec2  g    = vec2(v.x * sx1, v.y) * invV;
    if( rf ) g = vec2(g.x - k * g.y, -k * g.x - g.y) * 0.5;
    g.x *= sx0;
    return vec3(d, g);
}

// ----------------------------------------------------------------------------
// sdReuleauxPentagon2D(p, w) — Exact 0-Trig, 1-Sqrt 2D Reuleaux Pentagon
//
// Uses the 2-reflection 0-trig D_5 dihedral fold into the 36-degree fundamental
// sector where the active arc center is (0, -R_c) and the corner vertex is
// R_c * (sin 36, cos 36) with R_c = w / (2 * cos(pi/10)) = 0.52573111 * w.
// ----------------------------------------------------------------------------
float sdReuleauxPentagon2D( vec2 p, float w )
{
    const vec2  n1 = vec2(-0.80901699, 0.58778525);
    const vec2  n2 = vec2( 0.80901699, 0.58778525);
    const float rc = 0.52573111; // 1.0 / (2.0 * cos(pi / 10.0))
    p.x = abs(p.x);
    p -= 2.0 * min(dot(n1, p), 0.0) * n1;
    p -= 2.0 * min(dot(n2, p), 0.0) * n2;
    p.x = abs(p.x);
    float r = rc * w;
    bool  c = (0.95105652 * p.x - 0.30901699 * (p.y + r) > 0.0);
    return length(p - (c ? vec2(0.58778525, 0.80901699) * r : vec2(0.0, -r))) - (c ? 0.0 : w);
}

// ----------------------------------------------------------------------------
// sdgReuleauxPentagon2D(p, w) — Exact 0-Trig, 1-SFU Distance + 2D Unit Gradient
// ----------------------------------------------------------------------------
vec3 sdgReuleauxPentagon2D( vec2 p, float w )
{
    const vec2  n1 = vec2(-0.80901699, 0.58778525);
    const vec2  n2 = vec2( 0.80901699, 0.58778525);
    const float rc = 0.52573111;
    float sx0 = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    float d1 = dot(n1, p); bool r1 = (d1 < 0.0); if( r1 ) p -= 2.0 * d1 * n1;
    float d2 = dot(n2, p); bool r2 = (d2 < 0.0); if( r2 ) p -= 2.0 * d2 * n2;
    float sx1 = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);

    float r    = rc * w;
    bool  c    = (0.95105652 * p.x - 0.30901699 * (p.y + r) > 0.0);
    vec2  v    = p - (c ? vec2(0.58778525, 0.80901699) * r : vec2(0.0, -r));
    float v2   = dot(v, v);
    float invV = inversesqrt(max(v2, 1e-30));
    float d    = v2 * invV - (c ? 0.0 : w);
    vec2  g    = vec2(v.x * sx1, v.y) * invV;
    if( r2 ) g -= 2.0 * dot(n2, g) * n2;
    if( r1 ) g -= 2.0 * dot(n1, g) * n1;
    g.x *= sx0;
    return vec3(d, g);
}

// ----------------------------------------------------------------------------
// sdReuleauxPolygon2D(p, n, w) — General Odd-N Constant-Width Reuleaux Polygon
// ----------------------------------------------------------------------------
float sdReuleauxPolygon2D( vec2 p, float n, float w )
{
    float an  = 3.14159265 / n;
    vec2  csh = vec2(cos(0.5 * an), sin(0.5 * an));
    float rc  = 0.5 * w / csh.x;
    float rad = length(p);
    float ang = mod(atan(p.y, p.x), 2.0 * an);
    ang = an - abs(ang - an);
    vec2  q   = rad * vec2(cos(ang), sin(ang));
    bool  c   = ((q.x - rc) * csh.y - q.y * csh.x > 0.0);
    vec2  ctr = c ? vec2(rc, 0.0) : -rc * vec2(cos(an), sin(an));
    return length(q - ctr) - (c ? 0.0 : w);
}

#endif // SDF_GEMS_SD_REULEAUX_2D_GLSL
