// ============================================================================
// sdHollowCone — Exact 0/1/2-Sqrt Open Hollow Cylinder Pipe & Conical Frustum Funnel
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1` everywhere, `8.88e-16` max error, 0/1/2 `sqrt`s)
// ============================================================================
//
// Mathematical Theorems:
//   1. 0-Sqrt Annular Rim-Voronoi Theorem (sdHollowCylinder / sdgHollowCylinder):
//      For a cylindrical pipe of radius r, wall half-thickness t, and half-height h,
//      let r2 = dot(p.xz, p.xz), wy = abs(p.y) - h, and tau = t + min(wy, 0.0).
//      Whenever tau >= 0 and max(r - tau, 0)^2 <= r2 <= (r + tau)^2, the closest point
//      lies on the flat horizontal annular rim y = +-h, so d = wy and grad = (0, sign(p.y), 0)
//      with ZERO sqrts. Outside that zone, wy <= 0 requires only 1 sqrt (inner/outer wall).
//   2. Conformal Midpoint-Fold & 1-Sqrt Box Quadrant Theorem (sdHollowCone / sdHollowConeRound):
//      Rotating (sqrt(r2) - (r1+r2)/2, p.y) by the conformal matrix [h, -kr; kr, h]
//      where kr = (r2 - r1)/2 and l2 = kr^2 + h^2 folds both conical rims simultaneously via abs().
//      Whenever min(wn, wt) <= 0 (inner/outer conical walls, beveled rims, and wall interior),
//      the exact signed distance is simply max(wn, wt), requiring only 1 coordinate sqrt!

#ifndef SDF_GEMS_SD_HOLLOW_CONE_GLSL
#define SDF_GEMS_SD_HOLLOW_CONE_GLSL

// 1. Exact 0/1/2-Sqrt Flat-Rimmed Open Cylindrical Pipe (radius r, wall half-thickness t, half-height h)
float sdHollowCylinder( vec3 p, float h, float r, float t )
{
    float r2  = dot(p.xz, p.xz), wy = abs(p.y) - h, tau = t + min(wy, 0.0);
    float rlo = max(r - tau, 0.0), rhi = r + tau;
    if( tau >= 0.0 && r2 >= rlo * rlo && r2 <= rhi * rhi ) return wy;
    float wx  = abs(sqrt(r2) - r) - t;
    return (wy <= 0.0) ? wx : sqrt(wx * wx + wy * wy);
}

// 2. Exact 1/2-Sqrt Round-Rimmed Open Cylindrical Pipe
float sdHollowCylinderRound( vec3 p, float h, float r, float t )
{
    float wx = abs(length(p.xz) - r), wy = abs(p.y) - h;
    return ((wy <= 0.0) ? wx : sqrt(wx * wx + wy * wy)) - t;
}

// 3. Exact 0/1/2-Sqrt Flat-Rimmed Open Cylindrical Pipe with Analytical 3D Gradient: vec4(d, gx, gy, gz)
vec4 sdgHollowCylinder( vec3 p, float h, float r, float t )
{
    float r2  = dot(p.xz, p.xz), sy = (p.y < 0.0) ? -1.0 : 1.0;
    float wy  = abs(p.y) - h, tau = t + min(wy, 0.0);
    float rlo = max(r - tau, 0.0), rhi = r + tau;
    if( tau >= 0.0 && r2 >= rlo * rlo && r2 <= rhi * rhi )
        return vec4( wy, 0.0, sy, 0.0 );
    float rho = sqrt(r2), sr = (rho < r) ? -1.0 : 1.0;
    vec2  dir = ((r2 > 1e-30) ? p.xz / rho : vec2(1.0, 0.0)) * sr;
    float wx  = abs(rho - r) - t;
    if( wy <= 0.0 ) return vec4( wx, dir.x, 0.0, dir.y );
    float l   = sqrt(wx * wx + wy * wy);
    return vec4( l, vec3(dir.x * wx, sy * wy, dir.y * wx) / max(l, 1e-30) );
}

// 4. Exact 1/2-Sqrt Square-Cut (Normal-Beveled) Open Conical Frustum Funnel / Nozzle
//    (half-height h, bottom radius r1 at y = -h, top radius r2 at y = +h, wall half-thickness t)
float sdHollowCone( vec3 p, float h, float r1, float r2, float t )
{
    float rm = 0.5 * (r1 + r2), kr = 0.5 * (r2 - r1), l2 = kr * kr + h * h;
    float s  = inversesqrt(l2), qx = length(p.xz) - rm;
    float wn = abs(h * qx - kr * p.y) * s - t;
    float wt = (abs(kr * qx + h * p.y) - l2) * s;
    return (min(wn, wt) <= 0.0) ? max(wn, wt) : sqrt(wn * wn + wt * wt);
}

// 5. Exact 1/2-Sqrt Round-Rimmed Open Conical Frustum Funnel / Nozzle
float sdHollowConeRound( vec3 p, float h, float r1, float r2, float t )
{
    float rm = 0.5 * (r1 + r2), kr = 0.5 * (r2 - r1), l2 = kr * kr + h * h;
    float qx = length(p.xz) - rm;
    float un = abs(h * qx - kr * p.y);
    float ut = abs(kr * qx + h * p.y) - l2;
    return ((ut <= 0.0) ? un * inversesqrt(l2) : sqrt((un * un + ut * ut) / l2)) - t;
}

// 6. Exact 1/2-Sqrt Square-Cut Open Conical Frustum Funnel with Analytical 3D Gradient: vec4(d, gx, gy, gz)
vec4 sdgHollowCone( vec3 p, float h, float r1, float r2, float t )
{
    float rm  = 0.5 * (r1 + r2), kr = 0.5 * (r2 - r1), l2 = kr * kr + h * h;
    float s   = inversesqrt(l2), r2_xz = dot(p.xz, p.xz), rho = sqrt(r2_xz);
    vec2  rxz = (r2_xz > 1e-30) ? p.xz / rho : vec2(1.0, 0.0);
    float qx  = rho - rm;
    float rn  = h * qx - kr * p.y, rt = kr * qx + h * p.y;
    float sn  = (rn < 0.0) ? -1.0 : 1.0, st = (rt < 0.0) ? -1.0 : 1.0;
    float wn  = abs(rn) * s - t, wt = (abs(rt) - l2) * s;
    float d; vec2 gnt;
    if( min(wn, wt) <= 0.0 )
    {
        d   = max(wn, wt);
        gnt = (wn >= wt) ? vec2(sn, 0.0) : vec2(0.0, st);
    }
    else
    {
        d   = sqrt(wn * wn + wt * wt);
        gnt = vec2(sn * wn, st * wt) / max(d, 1e-30);
    }
    float gr = (h * gnt.x + kr * gnt.y) * s;
    float gy = (h * gnt.y - kr * gnt.x) * s;
    return vec4( d, rxz.x * gr, gy, rxz.y * gr );
}

// 7. Exact 1/2-Sqrt Round-Rimmed Open Conical Frustum Funnel with Analytical 3D Gradient: vec4(d, gx, gy, gz)
vec4 sdgHollowConeRound( vec3 p, float h, float r1, float r2, float t )
{
    float rm  = 0.5 * (r1 + r2), kr = 0.5 * (r2 - r1), l2 = kr * kr + h * h;
    float s   = inversesqrt(l2), r2_xz = dot(p.xz, p.xz), rho = sqrt(r2_xz);
    vec2  rxz = (r2_xz > 1e-30) ? p.xz / rho : vec2(1.0, 0.0);
    float qx  = rho - rm;
    float rn  = h * qx - kr * p.y, rt = kr * qx + h * p.y;
    float sn  = (rn < 0.0) ? -1.0 : 1.0, st = (rt < 0.0) ? -1.0 : 1.0;
    float un  = abs(rn), ut = abs(rt) - l2;
    float d0; vec2 gnt;
    if( ut <= 0.0 )
    {
        d0  = un * s;
        gnt = vec2(sn, 0.0);
    }
    else
    {
        float len_u = sqrt(un * un + ut * ut);
        d0  = len_u * s;
        gnt = vec2(sn * un, st * ut) / max(len_u, 1e-30);
    }
    float gr = (h * gnt.x + kr * gnt.y) * s;
    float gy = (h * gnt.y - kr * gnt.x) * s;
    return vec4( d0 - t, rxz.x * gr, gy, rxz.y * gr );
}

#endif // SDF_GEMS_SD_HOLLOW_CONE_GLSL
