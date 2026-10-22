// ============================================================================
// sdgCross3D — Exact 0/1-SFU 3D Gradients for sdCross3D, sdInfCross3D & sdJack3D
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1` in both exterior and interior; 0 SFUs on faces/arms, 1 in corners/jacks)
// ============================================================================

#ifndef SDF_GEMS_SDG_CROSS3D_GLSL
#define SDF_GEMS_SDG_CROSS3D_GLSL

// ----------------------------------------------------------------------------
// 1. sdgCross3D(p, b) — 0/1-SFU 3D Orthogonal Cross Distance + Unit Gradient
//    b.x = arm half-length, b.y = arm square half-thickness (b.x >= b.y)
// ----------------------------------------------------------------------------
vec4 sdgCross3D( vec3 p, vec2 b )
{
    vec3 s  = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3 q  = abs(p);
    bool s1 = (q.x < q.y); if( s1 ) q.xy = q.yx;
    bool s2 = (q.y < q.z); if( s2 ) q.yz = q.zy;
    bool s3 = (q.x < q.y); if( s3 ) q.xy = q.yx;

    vec3  w = q - b.xyy;
    float d;
    vec3  g;
    if( max(w.x, w.y) > 0.0 )
    {
        if( min(w.x, w.y) <= 0.0 && w.z <= 0.0 )
        {
            d = max(w.x, w.y);
            g = (w.x > w.y) ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
        }
        else
        {
            vec3  v    = max(w, 0.0);
            float l2   = dot(v, v), invL = inversesqrt(l2);
            d = l2 * invL;
            g = v * invL;
        }
    }
    else
    {
        float ux = b.y - q.x, uy = -w.y;
        if( ux <= 0.0 || w.x * w.x <= ux * ux + uy * uy )
        {
            d = max(w.x, w.y);
            g = (w.x > w.y) ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
        }
        else
        {
            float l2   = ux * ux + uy * uy, invL = inversesqrt(l2);
            d = -l2 * invL;
            g = vec3(ux * invL, uy * invL, 0.0);
        }
    }

    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4( d, s * g );
}

// ----------------------------------------------------------------------------
// 2. sdgInfCross3D(p, b) — 0/1-SFU Infinite Menger Cross Distance + Gradient
// ----------------------------------------------------------------------------
vec4 sdgInfCross3D( vec3 p, float b )
{
    vec3 s  = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3 q  = abs(p);
    bool s1 = (q.x < q.y); if( s1 ) q.xy = q.yx;
    bool s2 = (q.y < q.z); if( s2 ) q.yz = q.zy;
    bool s3 = (q.x < q.y); if( s3 ) q.xy = q.yx;

    float wy = q.y - b, wz = q.z - b;
    float d;
    vec3  g;
    if( wy > 0.0 )
    {
        if( wz <= 0.0 )
        {
            d = wy;
            g = vec3(0.0, 1.0, 0.0);
        }
        else
        {
            float l2 = wy * wy + wz * wz, invL = inversesqrt(l2);
            d = l2 * invL;
            g = vec3(0.0, wy * invL, wz * invL);
        }
    }
    else
    {
        float ux = b - q.x, uy = -wy;
        if( ux <= 0.0 )
        {
            d = wy;
            g = vec3(0.0, 1.0, 0.0);
        }
        else
        {
            float l2 = ux * ux + uy * uy, invL = inversesqrt(l2);
            d = -l2 * invL;
            g = vec3(ux * invL, uy * invL, 0.0);
        }
    }

    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4( d, s * g );
}

// ----------------------------------------------------------------------------
// 3. sdgJack3D(p, w, r) — 1-SFU 0-Swap 8-Pointed Body-Diagonal Jack + Gradient
// ----------------------------------------------------------------------------
vec4 sdgJack3D( vec3 p, float w, float r )
{
    vec3  s  = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3  q  = abs(p);
    vec3  v  = q - min(q.x + q.y + q.z, w) * (1.0 / 3.0);
    float l2 = dot(v, v);
    if( l2 > 0.0 )
    {
        float invL = inversesqrt(l2);
        return vec4( l2 * invL - r, s * v * invL );
    }
    return vec4( -r, s * 0.57735027 );
}

// ----------------------------------------------------------------------------
// 4. sdgAxisJack3D(p, w, r) — 1-SFU 0-Swap 6-Pointed Coordinate-Axis Jack + Gradient
// ----------------------------------------------------------------------------
vec4 sdgAxisJack3D( vec3 p, float w, float r )
{
    vec3  s  = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3  q  = abs(p);
    bool  ix = (q.x >= q.y) && (q.x >= q.z);
    bool  iy = (!ix) && (q.y >= q.z);
    vec3  v  = vec3( ix ? (q.x - min(q.x, w)) : q.x,
                     iy ? (q.y - min(q.y, w)) : q.y,
                     (!ix && !iy) ? (q.z - min(q.z, w)) : q.z );
    float l2 = dot(v, v);
    if( l2 > 0.0 )
    {
        float invL = inversesqrt(l2);
        return vec4( l2 * invL - r, s * v * invL );
    }
    return vec4( -r, s * (ix ? vec3(1.0, 0.0, 0.0) : (iy ? vec3(0.0, 1.0, 0.0) : vec3(0.0, 0.0, 1.0))) );
}

// ----------------------------------------------------------------------------
// 5. sdgFaceJack3D(p, w, r) — 1-SFU 12-Pointed Cuboctahedral Jack + Gradient
// ----------------------------------------------------------------------------
vec4 sdgFaceJack3D( vec3 p, float w, float r )
{
    vec3 s  = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3 q  = abs(p);
    bool s1 = (q.x < q.y); if( s1 ) q.xy = q.yx;
    bool s2 = (q.y < q.z); if( s2 ) q.yz = q.zy;
    bool s3 = (q.x < q.y); if( s3 ) q.xy = q.yx;

    float t  = min(q.x + q.y, w) * 0.5;
    vec3  v  = vec3(q.x - t, q.y - t, q.z);
    float l2 = dot(v, v);
    float d  = -r;
    vec3  g  = vec3(0.0, 0.0, 1.0);
    if( l2 > 0.0 )
    {
        float invL = inversesqrt(l2);
        d = l2 * invL - r;
        g = v * invL;
    }
    if( s3 ) g.xy = g.yx;
    if( s2 ) g.yz = g.zy;
    if( s1 ) g.xy = g.yx;
    return vec4( d, s * g );
}

#endif // SDF_GEMS_SDG_CROSS3D_GLSL
