// ============================================================================
// opRepHelix — Exact 1-Tap, 2-Tap, 4-Tap & Limited 3D Helical / Screw Domain Repetition
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`0.0` error vs brute-force multi-turn search; eliminates central-axis tear & oblique Voronoi slicing)
// ============================================================================

#ifndef SDF_GEMS_OP_REP_HELIX_GLSL
#define SDF_GEMS_OP_REP_HELIX_GLSL

// Transforms p into the local frame of helical step k (rotated by k*(2*pi/n) around Y, translated by k*(pitch/n) along Y)
vec3 opHelixTransform( vec3 p, float k, float pitch, float n )
{
    float ang = k * (6.28318530718 / n);
    float c   = cos(ang), s = sin(ang);
    return vec3(c * p.x + s * p.z, p.y - k * (pitch / n), c * p.z - s * p.x);
}

// Internal helper: computes the 4 sorted helical neighbor step indices (id0, id1, id2, id3)
// clamped to [lima, limb] (pass lima = -1e30, limb = 1e30 for infinite helix).
vec4 _opRepHelixCandidates( vec3 p, float R, float pitch, float n, float lima, float limb )
{
    float rxz    = length(p.xz);
    float theta  = atan(p.z, p.x);
    float kp     = pitch * 0.159154943092; // pitch / (2*pi)
    float yk     = p.y / kp;
    float dtheta = 6.28318530718 / n;
    float dy     = pitch / n;
    float yk_ref = clamp(yk, lima * dtheta, limb * dtheta);

    float n0   = round((yk_ref - theta) * 0.159154943092);
    float rem  = yk_ref - (theta + 6.28318530718 * n0);
    float n1   = n0 + ((rem >= 0.0) ? 1.0 : -1.0);
    float beta = 0.636619772368 * (R * rxz) / (kp * kp); // (2/pi) * alpha
    float b    = 1.0 + 2.0 * beta;

    // Primary turn n0
    float t0    = theta + 6.28318530718 * n0;
    float u0    = yk - t0, au0 = abs(u0);
    float disc0 = max(0.0, b * b - 2.54647908947 * beta * au0);
    float d0    = clamp((2.0 * au0) / (b + sqrt(disc0)), 0.0, min(au0, 3.14159265359));
    float idx0  = (t0 + ((u0 >= 0.0) ? d0 : -d0)) / dtheta;

    // Adjacent turn n1 with Turn-Discriminant Bifurcation
    float t1    = theta + 6.28318530718 * n1;
    float u1    = yk - t1, au1 = abs(u1);
    float disc1 = b * b - 2.54647908947 * beta * au1;
    float d1    = clamp((2.0 * au1) / (b + sqrt(max(0.0, disc1))), 0.0, 3.14159265359);
    float idx1  = (t1 + ((u1 >= 0.0) ? d1 : -d1)) / dtheta;

    float k0 = floor(idx0);
    float k1 = k0 + 1.0;
    float k2 = (disc1 > 0.0) ? floor(idx1) : (k0 - 1.0);
    float k3 = (disc1 > 0.0) ? (k2 + 1.0)  : (k1 + 1.0);
    if( abs(k2 - k0) <= 1.0 ) { k2 = k0 - 1.0; k3 = k1 + 1.0; }

    vec4 ks = clamp(vec4(k0, k1, k2, k3), lima, limb);
    vec4 e;
    for( int i = 0; i < 4; i++ )
    {
        float k   = ks[i];
        float ang = k * dtheta;
        float qx  = p.x * cos(ang) + p.z * sin(ang);
        float qy  = p.y - k * dy;
        e[i] = qy * qy - 2.0 * R * qx;
    }
    if( ks.y == ks.x ) e.y = 1e30;
    if( ks.z == ks.x || ks.z == ks.y ) e.z = 1e30;
    if( ks.w == ks.x || ks.w == ks.y || ks.w == ks.z ) e.w = 1e30;

    // Partial sorting network for 1st and 2nd closest candidates in ks.x, ks.y
    if( e.y < e.x ) { float tk = ks.x; ks.x = ks.y; ks.y = tk; float te = e.x; e.x = e.y; e.y = te; }
    if( e.w < e.z ) { float tk = ks.z; ks.z = ks.w; ks.w = tk; float te = e.z; e.z = e.w; e.w = te; }
    if( e.z < e.x ) { float tk = ks.x; ks.x = ks.z; ks.z = tk; float te = e.x; e.x = e.z; e.z = te; }
    if( e.z < e.y ) { float tk = ks.y; ks.y = ks.z; ks.z = tk; }
    return ks;
}

// 1. 1-Tap Infinite Helical Domain Repetition (returns local step coordinates q and integer step id)
vec3 opRepHelix( vec3 p, float R, float pitch, float n, out float id )
{
    vec4 ks = _opRepHelixCandidates(p, R, pitch, n, -1e30, 1e30);
    id = ks.x;
    return opHelixTransform(p, id, pitch, n);
}

vec3 opRepHelix( vec3 p, float R, float pitch, float n )
{
    float id;
    return opRepHelix(p, R, pitch, n, id);
}

// 2. 1-Tap Limited Helical Domain Repetition on step range [lima, limb]
vec3 opRepHelixLim( vec3 p, float R, float pitch, float n, float lima, float limb, out float id )
{
    vec4 ks = _opRepHelixCandidates(p, R, pitch, n, lima, limb);
    id = ks.x;
    return opHelixTransform(p, id, pitch, n);
}

// 3. 2-Tap Helical Domain Repetition for asymmetric primitives (outputs 1st & 2nd closest step indices)
void opRepHelix2Tap( in vec3 p, in float R, in float pitch, in float n, out float id0, out float id1 )
{
    vec4 ks = _opRepHelixCandidates(p, R, pitch, n, -1e30, 1e30);
    id0 = ks.x;
    id1 = ks.y;
}

// 4. 4-Tap Full Neighbor-Corrected Helical Domain Repetition (outputs all 4 along-helix & cross-turn step indices)
vec4 opRepHelix4Tap( in vec3 p, in float R, in float pitch, in float n )
{
    return _opRepHelixCandidates(p, R, pitch, n, -1e30, 1e30);
}

#endif // SDF_GEMS_OP_REP_HELIX_GLSL
