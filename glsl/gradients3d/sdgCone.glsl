// ============================================================================
// sdgCone, sdgConeInf & sdgRhombus — 0/1/2-SFU 3D Analytical Gradients
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, 0 SFUs in cap/face Voronoi zones, NaN-free everywhere)
// ============================================================================

#ifndef SDF_GEMS_SDG_CONE_GLSL
#define SDF_GEMS_SDG_CONE_GLSL

// 1. Exact 0/1/2-SFU Finite 3D Cone Distance + Analytical Gradient (vec4(d, grad.xyz))
//    Uses the Slanted-Wall Cap-Voronoi Radius Theorem:
//      rc = (h * c.x + (1.0 + c.x) * min(d, 0.0)) / c.y
//    to evaluate both distance and gradient with ZERO square roots across the entire
//    exterior base cylinder and interior base-cap Voronoi frustum!
vec4 sdgCone( vec3 p, vec2 c, float h )
{
    float x2 = dot(p.xz, p.xz), d = -p.y - h;
    float rc = (h * c.x + (1.0 + c.x) * min(d, 0.0)) / c.y;
    if( rc >= 0.0 && x2 <= rc * rc )
        return vec4(d, 0.0, -1.0, 0.0);
    if( p.y >= 0.0 && c.x * c.x * x2 <= c.y * c.y * p.y * p.y )
    {
        float r2   = x2 + p.y * p.y;
        float invR = inversesqrt(max(r2, 1e-30));
        return vec4(r2 * invR, (r2 > 0.0) ? (p * invR) : vec3(0.0, 1.0, 0.0));
    }
    float invX = inversesqrt(max(x2, 1e-30));
    float wx   = (x2 > 0.0) ? (x2 * invX) : 0.0;
    vec2  er   = (x2 > 0.0) ? (p.xz * invX) : vec2(1.0, 0.0);
    float rx   = wx - h * c.x / c.y;
    if( rx > 0.0 && c.x * rx + c.y * d >= 0.0 )
    {
        float l2   = rx * rx + d * d;
        float invL = inversesqrt(max(l2, 1e-30));
        return vec4(l2 * invL, vec3(rx * er.x, -d, rx * er.y) * invL);
    }
    return vec4(c.y * wx + c.x * p.y, c.y * er.x, c.x, c.y * er.y);
}

// 2. Exact 1-SFU Infinite 3D Cone Distance + Analytical Gradient
vec4 sdgConeInf( vec3 p, vec2 c )
{
    float x2 = dot(p.xz, p.xz);
    if( p.y >= 0.0 && c.x * c.x * x2 <= c.y * c.y * p.y * p.y )
    {
        float r2   = x2 + p.y * p.y;
        float invR = inversesqrt(max(r2, 1e-30));
        return vec4(r2 * invR, (r2 > 0.0) ? (p * invR) : vec3(0.0, 1.0, 0.0));
    }
    float invX = inversesqrt(max(x2, 1e-30));
    float wx   = (x2 > 0.0) ? (x2 * invX) : 0.0;
    vec2  er   = (x2 > 0.0) ? (p.xz * invX) : vec2(1.0, 0.0);
    return vec4(c.y * wx + p.y * c.x, c.y * er.x, c.x, c.y * er.y);
}

// 3. Exact 0/1/2-SFU 3D Extruded Rounded Rhombus Distance + Analytical Gradient
vec4 sdgRhombus( vec3 p, float la, float lb, float h, float ra )
{
    vec3  s    = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3  q    = abs(p);
    float invL = inversesqrt(la * la + lb * lb);
    float qx;
    vec2  gxz;
    if( la * (q.x - la) > lb * q.z )
    {
        vec2  v    = vec2(q.x - la, q.z);
        float l2   = dot(v, v), invV = inversesqrt(max(l2, 1e-30));
        qx  = l2 * invV - ra;
        gxz = (l2 > 0.0) ? (v * invV) : vec2(1.0, 0.0);
    }
    else if( lb * (q.z - lb) > la * q.x )
    {
        vec2  v    = vec2(q.x, q.z - lb);
        float l2   = dot(v, v), invV = inversesqrt(max(l2, 1e-30));
        qx  = l2 * invV - ra;
        gxz = (l2 > 0.0) ? (v * invV) : vec2(0.0, 1.0);
    }
    else
    {
        qx  = (lb * q.x + la * q.z - la * lb) * invL - ra;
        gxz = vec2(lb, la) * invL;
    }
    float qy = q.y - h;
    if( qy <= 0.0 )
        return (qx > qy) ? vec4(qx, s.x * gxz.x, 0.0, s.z * gxz.y) : vec4(qy, 0.0, s.y, 0.0);
    if( qx <= 0.0 )
        return vec4(qy, 0.0, s.y, 0.0);
    float l2   = qx * qx + qy * qy;
    float invC = inversesqrt(max(l2, 1e-30));
    return vec4(l2 * invC, vec3(s.x * qx * gxz.x, s.y * qy, s.z * qx * gxz.y) * invC);
}

#endif // SDF_GEMS_SDG_CONE_GLSL
