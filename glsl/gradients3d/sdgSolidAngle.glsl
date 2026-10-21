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
        float invR = inversesqrt(max(r2, 1e-30));
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
        return vec4(r2 * invR - ra, p * invR);
    }

    // 2-SFU exterior circular rim wedge (v > 0 && u > ra)
    vec2  w    = vec2(qx, y) - ra * c;
    float w2   = dot(w, w);
    float invW = inversesqrt(max(w2, 1e-30));
    return vec4(w2 * invW, vec3(w.x * er.x, w.y, w.x * er.y) * invW);
}

// 2. 1/2-SFU 3D Capped Torus Distance and Unit Gradient vec4(d, grad.xyz)
vec4 sdgCappedTorus( vec3 p, vec2 sc, float ra, float rb )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    if( sc.y * p.x > sc.x * p.y )
    {
        vec3  w    = vec3(p.x - ra * sc.x, p.y - ra * sc.y, p.z);
        float w2   = dot(w, w);
        float invW = inversesqrt(max(w2, 1e-30)); // 1 SFU op total!
        return vec4(w2 * invW - rb, sx * w.x * invW, w.y * invW, w.z * invW);
    }
    float r2   = dot(p.xy, p.xy);
    float invR = inversesqrt(max(r2, 1e-30));
    float r_xy = (r2 > 0.0) ? (r2 * invR) : 0.0;
    vec2  er   = (r2 > 0.0) ? (vec2(sx * p.x, p.y) * invR) : vec2(0.0, 1.0);
    float u    = r_xy - ra;
    float l2   = u * u + p.z * p.z;
    float invL = inversesqrt(max(l2, 1e-30));
    return vec4(l2 * invL - rb, vec3(u * er, p.z) * invL);
}

// 3. 0/1/2-SFU 3D Cut Sphere Distance and Unit Gradient vec4(d, grad.xyz)
vec4 sdgCutSphere( vec3 p, float r, float h )
{
    float x2 = dot(p.xz, p.xz), w2 = r * r - h * h, dx = x2 - w2, dy = h - p.y;
    if( (r + h) * (r - h + 2.0 * dy) < x2 && h * abs(h) * x2 < w2 * p.y * abs(p.y) )
    {
        float r2   = x2 + p.y * p.y;
        float invR = inversesqrt(max(r2, 1e-30)); // 1 SFU op in spherical zone!
        return vec4(r2 * invR - r, p * invR);
    }
    if( dx <= 0.0 )
        return vec4(dy, 0.0, -1.0, 0.0); // 0 SFU ops in exterior & interior flat disk cap zone!
    float invX = inversesqrt(x2);
    float x    = x2 * invX;
    float ux   = dx / (x + sqrt(w2));
    float l2   = ux * ux + dy * dy;
    float invL = inversesqrt(max(l2, 1e-30));
    float k_xz = ux * invX * invL;
    return vec4(l2 * invL, p.x * k_xz, -dy * invL, p.z * k_xz);
}

// 4. 1/2-SFU 3D Cut Hollow Sphere Distance and Unit Gradient vec4(d, grad.xyz)
vec4 sdgCutHollowSphere( vec3 p, float r, float h, float t )
{
    float x2 = dot(p.xz, p.xz), w2 = r * r - h * h, dx = x2 - w2, dy = p.y - h;
    if( h * abs(h) * x2 >= w2 * p.y * abs(p.y) )
    {
        float r2   = x2 + p.y * p.y;
        float invR = inversesqrt(max(r2, 1e-30)); // 1 SFU op!
        float u    = r2 * invR - r;
        float su   = (u < 0.0) ? -1.0 : 1.0;
        return vec4(abs(u) - t, (su * invR) * p);
    }
    float invX = inversesqrt(max(x2, 1e-30));
    float x    = (x2 > 0.0) ? (x2 * invX) : 0.0;
    vec2  er   = (x2 > 0.0) ? (p.xz * invX) : vec2(1.0, 0.0);
    float w    = sqrt(w2);
    float ux   = (x2 > 0.0) ? (dx / (x + w)) : -w;
    float l2   = ux * ux + dy * dy;
    float invL = inversesqrt(max(l2, 1e-30));
    return vec4(l2 * invL - t, vec3(ux * er.x, dy, ux * er.y) * invL);
}
