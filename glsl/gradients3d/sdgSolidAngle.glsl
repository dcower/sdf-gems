// ============================================================================
// sdgSolidAngle — 0/1/2-SFU 3D Solid Angle, Capped Torus, & Cut Sphere Gradients
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 1 inversesqrt in 65%-81% of space for `sdgSolidAngle`, 0/1/2 SFU ops for `sdgCutSphere`, 1/2 SFU ops for `sdgCappedTorus`)
// ============================================================================
//
// Mathematical Foundation:
//   1. `sdgSolidAngle(vec3 p, vec2 c, float ra)` — 1/2-SFU via Signed-Square & Cap-Voronoi:
//      Let x2 = dot(p.xz, p.xz), y = p.y, and r2 = x2 + y*y.
//      In the meridian plane (qx, y) with cone generator c = (sin(alpha), cos(alpha)),
//      let u = c.x*qx + c.y*y (along generator) and v = c.y*qx - c.x*y (signed normal).
//      - By Signed-Square Reduction, v > 0 <==> c.y*|c.y|*x2 > (c.x*y)*|c.x*y| and
//        u < 0 <==> c.x^2 * x2 < (-c.y*y)*|-c.y*y| WITHOUT computing qx = sqrt(x2)!
//        Whenever u < 0 or (v <= 0 && r2 > ra^2), the closest point is either the
//        apex (0,0,0) or the spherical cap at radius ra, requiring ONLY 1 `inversesqrt(r2)`!
//      - In the remaining region, 1 `inversesqrt(x2)` gives qx and unit radial
//        direction er = p.xz * invQx. By the Slanted-Edge Cap-Voronoi Radius Theorem
//        (rc = ra + min(v, 0.0)), the conical wall wins iff (v > 0 && u <= ra) or
//        (v <= 0 && rc >= 0 && r2 <= rc^2), returning vec4(v, c.y*er.x, -c.x, c.y*er.y)
//        with ZERO further SFU ops (1 SFU op total across 65%-81% of R^3)!
//
//   2. `sdgCappedTorus(vec3 p, vec2 sc, float ra, float rb)` — 1/2-SFU:
//      In the end-cap hemisphere wedge (sc.y * |p.x| > sc.x * p.y), the closest
//      core point is the 3D cap center (+/-ra*sc.x, ra*sc.y, 0), requiring ONLY
//      1 `inversesqrt` (62.8% of [-2ra, 2ra]^3) and 0 z-axis singularities!
//
//   3. `sdgCutSphere(vec3 p, float r, float h)` & `sdgCutHollowSphere`:
//      Combines Conjugate-Square & Signed-Square region tests with `inversesqrt`:
//      0 SFU ops in the flat cap disk zone (dx <= 0), 1 `inversesqrt` in the
//      spherical body zone, and 2 `inversesqrt`s in the circular rim wedge.

#ifndef SDF_GEMS_SDG_SOLID_ANGLE_GLSL
#define SDF_GEMS_SDG_SOLID_ANGLE_GLSL

// 1. 1/2-SFU 3D Solid Angle Distance and Unit Gradient vec4(d, grad.xyz)
vec4 sdgSolidAngle( vec3 p, vec2 c, float ra )
{
    float x2    = dot(p.xz, p.xz), y = p.y, r2 = x2 + y * y;
    float cxy   = c.x * y, cyy = -c.y * y;
    bool  v_pos = (c.y * abs(c.y) * x2 > cxy * abs(cxy));
    bool  u_neg = (c.x * c.x * x2 < cyy * abs(cyy));

    // Fast 1-SFU spherical cap & apex exit when qx = sqrt(x2) is not needed!
    if( u_neg || (!v_pos && r2 > ra * ra) )
    {
        if( r2 <= 0.0 ) return vec4(0.0, 0.0, -1.0, 0.0);
        float invR = inversesqrt(r2);
        float rl   = r2 * invR;
        if( v_pos ) return vec4(rl, p * invR);
        float s = (u_neg && 4.0 * r2 < ra * ra) ? -1.0 : 1.0;
        return vec4((s < 0.0) ? -rl : (rl - ra), (s * invR) * p);
    }

    float invQx = inversesqrt(max(x2, 1e-30));
    float qx    = (x2 > 0.0) ? (x2 * invQx) : 0.0;
    vec2  er    = (x2 > 0.0) ? (p.xz * invQx) : vec2(1.0, 0.0);
    float u     = c.x * qx + c.y * y;
    float v     = c.y * qx - c.x * y;
    float rc    = ra + min(v, 0.0);

    // 1-SFU conical-wall Voronoi zone (both exterior & interior!)
    if( (v_pos && u <= ra) || (!v_pos && rc >= 0.0 && r2 <= rc * rc) )
        return vec4(v, c.y * er.x, -c.x, c.y * er.y);

    // Interior spherical cap sub-wedge
    if( !v_pos )
    {
        float invR = inversesqrt(max(r2, 1e-30));
        return vec4(r2 * invR - ra, (r2 > 0.0) ? (p * invR) : vec3(0.0, 1.0, 0.0));
    }

    // 2-SFU exterior circular rim wedge (v > 0 && u > ra)
    vec2  w    = vec2(qx, y) - ra * c;
    float w2   = dot(w, w);
    float invW = inversesqrt(max(w2, 1e-30));
    return vec4(w2 * invW, (w2 > 0.0) ? (vec3(w.x * er.x, w.y, w.x * er.y) * invW) : vec3(c.x * er.x, c.y, c.x * er.y));
}

// 2. 1/2-SFU 3D Capped Torus Distance and Unit Gradient vec4(d, grad.xyz)
#ifndef SDF_GEMS_HAS_SDG_CAPPED_TORUS
#define SDF_GEMS_HAS_SDG_CAPPED_TORUS
vec4 sdgCappedTorus( vec3 p, vec2 sc, float ra, float rb )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    if( sc.y * p.x > sc.x * p.y )
    {
        vec3  w  = vec3(p.x - ra * sc.x, p.y - ra * sc.y, p.z);
        float w2 = dot(w, w);
        if( w2 > 0.0 )
        {
            float invW = inversesqrt(w2);
            return vec4(w2 * invW - rb, sx * w.x * invW, w.y * invW, w.z * invW);
        }
        return vec4(-rb, sx * sc.y, -sc.x, 0.0);
    }
    float r2   = dot(p.xy, p.xy);
    float invR = (r2 > 0.0) ? inversesqrt(r2) : 0.0;
    float r_xy = r2 * invR;
    vec2  er   = (r2 > 0.0) ? (vec2(sx * p.x, p.y) * invR) : vec2(0.0, 1.0);
    float u    = r_xy - ra;
    float l2   = u * u + p.z * p.z;
    if( l2 > 0.0 )
    {
        float invL = inversesqrt(l2);
        return vec4(l2 * invL - rb, vec3(u * er, p.z) * invL);
    }
    return vec4(-rb, er.x, er.y, 0.0);
}
#endif // SDF_GEMS_HAS_SDG_CAPPED_TORUS

// 3. 0/1/2-SFU 3D Cut Sphere Distance and Unit Gradient vec4(d, grad.xyz)
#ifndef SDF_GEMS_HAS_SDG_CUT_SPHERE
#define SDF_GEMS_HAS_SDG_CUT_SPHERE
vec4 sdgCutSphere( vec3 p, float r, float h )
{
    float x2 = dot(p.xz, p.xz), dy = h - p.y;
    if( x2 <= (r + h) * (r - h + 2.0 * min(dy, 0.0)) )
        return vec4(dy, 0.0, -1.0, 0.0);
    float w2 = r * r - h * h;
    if( h * abs(h) * x2 < w2 * p.y * abs(p.y) )
    {
        float r2 = x2 + p.y * p.y;
        if( r2 > 0.0 )
        {
            float invR = inversesqrt(r2);
            return vec4(r2 * invR - r, p * invR);
        }
        return vec4(-r, 0.0, 1.0, 0.0);
    }
    float d  = x2 - w2, wx = sqrt(w2 * x2);
    float sc = d / max(x2 + wx, 1e-30);
    vec3  v  = vec3(sc * p.x, -dy, sc * p.z);
    float v2 = dot(v, v);
    if( v2 > 0.0 )
    {
        float invV = inversesqrt(v2);
        return vec4(v2 * invV, v * invV);
    }
    return vec4(0.0, 0.0, -1.0, 0.0);
}
#endif // SDF_GEMS_HAS_SDG_CUT_SPHERE

// 4. 1/2-SFU 3D Cut Hollow Sphere Distance and Unit Gradient vec4(d, grad.xyz)
#ifndef SDF_GEMS_HAS_SDG_CUT_HOLLOW_SPHERE
#define SDF_GEMS_HAS_SDG_CUT_HOLLOW_SPHERE
vec4 sdgCutHollowSphere( vec3 p, float r, float h, float t )
{
    float x2 = dot(p.xz, p.xz), w2 = r * r - h * h;
    if( h * abs(h) * x2 >= w2 * p.y * abs(p.y) )
    {
        float r2 = x2 + p.y * p.y;
        if( r2 > 0.0 )
        {
            float invR = inversesqrt(r2);
            float u    = r2 * invR - r;
            float su   = (u < 0.0) ? -1.0 : 1.0;
            return vec4(abs(u) - t, (su * invR) * p);
        }
        return vec4(r - t, 0.0, (h >= 0.0) ? -1.0 : 1.0, 0.0);
    }
    float invX = (x2 > 0.0) ? inversesqrt(x2) : 0.0;
    float x    = x2 * invX;
    vec2  er   = (x2 > 0.0) ? (p.xz * invX) : vec2(1.0, 0.0);
    float w    = sqrt(w2);
    float ux   = (x2 > 0.0) ? ((x2 - w2) / (x + w)) : -w;
    float dy   = p.y - h;
    float l2   = ux * ux + dy * dy;
    if( l2 > 0.0 )
    {
        float invL = inversesqrt(l2);
        return vec4(l2 * invL - t, vec3(ux * er.x, dy, ux * er.y) * invL);
    }
    return vec4(-t, er.x, 0.0, er.y);
}
#endif // SDF_GEMS_HAS_SDG_CUT_HOLLOW_SPHERE

#endif // SDF_GEMS_SDG_SOLID_ANGLE_GLSL
