// ============================================================================
// sdBarrel, sdgBarrel, sdBarrel2D & sdgBarrel2D — Exact 0/1/2-Sqrt 3D & 2D Barrels and Spools
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1` everywhere; 0 sqrts in flat cap cylinder, 1/2 sqrts in walls/rim)
// ============================================================================

#ifndef SDF_GEMS_SD_BARREL_GLSL
#define SDF_GEMS_SD_BARREL_GLSL

// 1. Exact 0/1/2-Sqrt 3D Convex Barrel (r2 > r1) & Concave Waisted Spool (r2 < r1)
float sdBarrel( vec3 p, float r1, float r2, float h )
{
    float x2 = dot(p.xz, p.xz), qy = abs(p.y);
    float b  = r2 - r1;
    if( abs(b) < 1e-5 * h )
    {
        vec2 d = vec2(sqrt(x2) - r1, qy - h);
        return length(max(d, 0.0)) + min(max(d.x, d.y), 0.0);
    }
    float sb = (b < 0.0) ? -1.0 : 1.0;
    float d  = 0.5 * (h * h - b * b) / b;
    float R  = d + b, xc = r1 - d;
    float w  = h * xc + d * qy;
    if( x2 <= r1 * r1 && (qy >= h || sb * (h * h * x2 - w * abs(w)) <= 0.0) )
        return qy - h; // 0 sqrts in exterior & interior flat-cap cylinder!

    float qx = sqrt(x2);
    float u  = qx - xc;
    float s  = max(sb * (d * qy - h * u), sb * ((R - h) * (R - h + 2.0 * qy) - u * u));
    if( s >= 0.0 && qx <= r1 )
        return qy - h; // 1 sqrt in remaining parabolic cap sub-wedge

    vec2  v  = (s < 0.0) ? vec2(u, qy) : vec2(qx - r1, qy - h);
    float vl = length(v);
    return (s < 0.0) ? (sb * vl - R) : vl;
}

// 2. Analytical 3D Distance & Unit Gradient vec4(d, grad.xyz) for 3D Barrel & Spool
vec4 sdgBarrel( vec3 p, float r1, float r2, float h )
{
    float x2 = dot(p.xz, p.xz), qy = abs(p.y);
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float b  = r2 - r1;
    if( abs(b) < 1e-5 * h )
    {
        float invQx = inversesqrt(max(x2, 1e-30));
        float qx    = (x2 > 0.0) ? (x2 * invQx) : 0.0;
        vec2  er    = (x2 > 0.0) ? (p.xz * invQx) : vec2(1.0, 0.0);
        vec2  d     = vec2(qx - r1, qy - h);
        if( max(d.x, d.y) <= 0.0 )
            return (d.x > d.y) ? vec4(d.x, er.x, 0.0, er.y) : vec4(d.y, 0.0, sy, 0.0);
        vec2  vp   = max(d, 0.0);
        float vl   = length(vp);
        vec2  g2   = (vl > 0.0) ? (vp / vl) : vec2(1.0, 0.0);
        return vec4(vl, g2.x * er.x, sy * g2.y, g2.x * er.y);
    }
    float sb = (b < 0.0) ? -1.0 : 1.0;
    float d0 = 0.5 * (h * h - b * b) / b;
    float R  = d0 + b, xc = r1 - d0;
    float w  = h * xc + d0 * qy;
    if( x2 <= r1 * r1 && (qy >= h || sb * (h * h * x2 - w * abs(w)) <= 0.0) )
        return vec4(qy - h, 0.0, sy, 0.0); // 0 SFU ops!

    float invQx = inversesqrt(max(x2, 1e-30));
    float qx    = (x2 > 0.0) ? (x2 * invQx) : 0.0;
    vec2  er    = (x2 > 0.0) ? (p.xz * invQx) : vec2(1.0, 0.0);
    float u     = qx - xc;
    float s     = max(sb * (d0 * qy - h * u), sb * ((R - h) * (R - h + 2.0 * qy) - u * u));
    if( s >= 0.0 && qx <= r1 )
        return vec4(qy - h, 0.0, sy, 0.0);

    bool  in_arc = (s < 0.0);
    vec2  v      = in_arc ? vec2(u, qy) : vec2(qx - r1, qy - h);
    float v2     = dot(v, v);
    float invV   = inversesqrt(max(v2, 1e-30));
    float vl     = v2 * invV;
    float d      = in_arc ? (sb * vl - R) : vl;
    vec2  g2     = ((in_arc ? sb : 1.0) * invV) * v;
    return vec4(d, g2.x * er.x, sy * g2.y, g2.x * er.y);
}

// 3. Exact 0/1-Sqrt 2D Barrel & Waisted Spool
float sdBarrel2D( vec2 p, float r1, float r2, float h )
{
    vec2  q  = abs(p);
    float b  = r2 - r1;
    if( abs(b) < 1e-5 * h )
    {
        vec2 d = q - vec2(r1, h);
        return length(max(d, 0.0)) + min(max(d.x, d.y), 0.0);
    }
    float sb = (b < 0.0) ? -1.0 : 1.0;
    float d  = 0.5 * (h * h - b * b) / b;
    float R  = d + b;
    float u  = q.x - (r1 - d);
    float s  = max(sb * (d * q.y - h * u), sb * ((R - h) * (R - h + 2.0 * q.y) - u * u));
    if( s >= 0.0 && q.x <= r1 )
        return q.y - h; // 0 sqrts in flat cap zone!

    vec2  v  = (s < 0.0) ? vec2(u, q.y) : (q - vec2(r1, h));
    float vl = length(v);
    return (s < 0.0) ? (sb * vl - R) : vl;
}

// 4. Analytical 2D Distance & Unit Gradient vec3(d, grad.xy) for 2D Barrel & Spool
vec3 sdgBarrel2D( vec2 p, float r1, float r2, float h )
{
    vec2  sg = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    vec2  q  = abs(p);
    float b  = r2 - r1;
    if( abs(b) < 1e-5 * h )
    {
        vec2 d = q - vec2(r1, h);
        if( max(d.x, d.y) <= 0.0 )
            return (d.x > d.y) ? vec3(d.x, sg.x, 0.0) : vec3(d.y, 0.0, sg.y);
        vec2  vp = max(d, 0.0);
        float vl = length(vp);
        vec2  g  = (vl > 0.0) ? (vp / vl) : vec2(1.0, 0.0);
        return vec3(vl, sg * g);
    }
    float sb = (b < 0.0) ? -1.0 : 1.0;
    float d0 = 0.5 * (h * h - b * b) / b;
    float R  = d0 + b;
    float u  = q.x - (r1 - d0);
    float s  = max(sb * (d0 * q.y - h * u), sb * ((R - h) * (R - h + 2.0 * q.y) - u * u));
    if( s >= 0.0 && q.x <= r1 )
        return vec3(q.y - h, 0.0, sg.y);

    bool  in_arc = (s < 0.0);
    vec2  v      = in_arc ? vec2(u, q.y) : (q - vec2(r1, h));
    float v2     = dot(v, v);
    float invV   = inversesqrt(max(v2, 1e-30));
    float vl     = v2 * invV;
    float d      = in_arc ? (sb * vl - R) : vl;
    vec2  g      = ((in_arc ? sb : 1.0) * invV) * v;
    return vec3(d, sg * g);
}

#endif // SDF_GEMS_SD_BARREL_GLSL
