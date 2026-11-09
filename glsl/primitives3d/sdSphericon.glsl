#ifndef SDF_GEMS_HAS_SD_SPHERICON
#define SDF_GEMS_HAS_SD_SPHERICON

#include "../common/common.glsl"

// ============================================================================
// sdSphericon, sdOctasphericon & opSphericonCoords — Exact Developable Rollers
// @author Gemini 4 Argon (research & mathematical design) and Gemini 3.8 Flash (experiments & implementation)
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt`, 0 loops)
// ============================================================================

// S4 Rotoreflection domain operator: folds (x, y, z) into the upper half-space
// and returns exact 2D meridian coordinates (u, rho) = (|q.x|, length(q.yz)).
vec2 opSphericonCoords( vec3 p )
{
    vec3 q = (p.z < 0.0) ? vec3(p.y, p.x, -p.z) : p;
    return vec2(abs(q.x), length(q.yz));
}

// ----------------------------------------------------------------------------
// sdSphericon: Exact Euclidean signed distance to a Sphericon of radius r
// ----------------------------------------------------------------------------
float sdSphericon( vec3 p, float r )
{
    vec2 m = opSphericonCoords(p);
    float u = m.x;
    float rho = m.y;

    if (abs(u - rho) <= r)
    {
        return (u + rho - r) * 0.7071067811865475;
    }

    return (u > rho) ? length(vec2(u - r, rho)) : length(vec2(u, rho - r));
}

// ----------------------------------------------------------------------------
// sdgSphericon: Exact Euclidean signed distance and analytical unit gradient
// ----------------------------------------------------------------------------
vec4 sdgSphericon( vec3 p, float r )
{
    bool flip = (p.z < 0.0);
    vec3 q = flip ? vec3(p.y, p.x, -p.z) : p;

    float sx = (q.x < 0.0) ? -1.0 : 1.0;
    float u = abs(q.x);
    float rho2 = q.y * q.y + q.z * q.z;
    float rho = sqrt(rho2);

    vec2 dir_yz = (rho > 1e-12) ? (q.yz / rho) : vec2(0.0, 1.0);

    const float INV_SQRT2 = 0.7071067811865475;

    float d;
    vec2 gu_gr;

    if (abs(u - rho) <= r)
    {
        d = (u + rho - r) * INV_SQRT2;
        gu_gr = vec2(INV_SQRT2, INV_SQRT2);
    }
    else if (u - rho > r)
    {
        vec2 v = vec2(u - r, rho);
        d = length(v);
        gu_gr = (d > 0.0) ? (v / d) : vec2(1.0, 0.0);
    }
    else
    {
        vec2 v = vec2(u, rho - r);
        d = length(v);
        gu_gr = (d > 0.0) ? (v / d) : vec2(0.0, 1.0);
    }

    vec3 gq = vec3(gu_gr.x * sx, gu_gr.y * dir_yz);
    vec3 gp = flip ? vec3(gq.y, gq.x, -gq.z) : gq;

    return vec4(d, gp);
}

// ----------------------------------------------------------------------------
// sdOctasphericon: Exact Euclidean signed distance to an Octasphericon of radius r
// ----------------------------------------------------------------------------
float sdOctasphericon( vec3 p, float r )
{
    vec2 m = opSphericonCoords(p);
    float u = max(m.x, m.y);
    float v = min(m.x, m.y);

    const float cos8 = 0.9238795325112867; // cos(pi/8)
    const float sin8 = 0.3826834323650898; // sin(pi/8)

    float h = r * cos8;
    float L = r * sin8;

    float dn = u * cos8 + v * sin8 - h;
    float dt = -u * sin8 + v * cos8;
    float dt_excess = max(abs(dt) - L, 0.0);

    return (dn <= 0.0) ? dn : length(vec2(dn, dt_excess));
}

// ----------------------------------------------------------------------------
// sdgOctasphericon: Exact Euclidean signed distance and analytical unit gradient
// ----------------------------------------------------------------------------
vec4 sdgOctasphericon( vec3 p, float r )
{
    bool flip = (p.z < 0.0);
    vec3 q = flip ? vec3(p.y, p.x, -p.z) : p;

    float sx = (q.x < 0.0) ? -1.0 : 1.0;
    float u = abs(q.x);
    float rho2 = q.y * q.y + q.z * q.z;
    float rho = sqrt(rho2);

    vec2 dir_yz = (rho > 1e-12) ? (q.yz / rho) : vec2(0.0, 1.0);

    const float cos8 = 0.9238795325112867; // cos(pi/8)
    const float sin8 = 0.3826834323650898; // sin(pi/8)
    const vec2 n = vec2(cos8, sin8);
    const vec2 tau = vec2(-sin8, cos8);

    bool sw = (u < rho);
    float a = sw ? rho : u;
    float b = sw ? u : rho;

    float h = r * cos8;
    float L = r * sin8;

    float dn = a * cos8 + b * sin8 - h;
    float dt = -a * sin8 + b * cos8;
    float dt_excess = max(abs(dt) - L, 0.0);

    float d;
    vec2 gab;

    if (dn <= 0.0 || dt_excess == 0.0)
    {
        d = dn;
        gab = n;
    }
    else
    {
        d = length(vec2(dn, dt_excess));
        float st = (dt < 0.0) ? -1.0 : 1.0;
        gab = (dn * n + (st * dt_excess) * tau) / d;
    }

    vec2 gu_gr = sw ? gab.yx : gab;
    vec3 gq = vec3(gu_gr.x * sx, gu_gr.y * dir_yz);
    vec3 gp = flip ? vec3(gq.y, gq.x, -gq.z) : gq;

    return vec4(d, gp);
}

#endif // SDF_GEMS_HAS_SD_SPHERICON
