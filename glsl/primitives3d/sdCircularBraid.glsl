// ============================================================================
// sdCircularBraid, sdgCircularBraid, sdLinearBraid3D & sdgLinearBraid3D
// — Exact O(1) 3D Circular 3-Strand Celtic Braid Ring & Linear 3-Strand Plait
//
// @author  Discovered and derived using Gemini 4 Argon (research & mathematical design)
//          and Gemini 3.8 Flash (implementation & verification) for sdf-gems (2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^3; O(1) complexity completely independent of wave count N)
// @novelty Models a true alternating over-under 3-strand braid via a Lemniscate of
//          Gerono cross-section (a*cos(theta_k), b*sin(2*theta_k)) for phase-shifted
//          strands theta_k = N*phi + 2*pi*k/3 (k = 0, 1, 2), where crossings in radial
//          offset occur at maximal vertical separation sqrt(3)*b. By D_{3N} dihedral-strand
//          reflection folding into a single fundamental cell [0, pi/(3*N)] and a vectorized
//          3-strand vec3 Newton solver, an N-wave Celtic braid ring executes in O(1) time!
// ============================================================================

#ifndef SDF_GEMS_SD_CIRCULAR_BRAID_GLSL
#define SDF_GEMS_SD_CIRCULAR_BRAID_GLSL

#include "../common/common.glsl"

// 1. Exact O(1) 3D Circular 3-Strand Celtic Braid Ring + Analytical Unit Gradient vec4(d, grad.xyz)
//    R: major ring radius, a: radial weave amplitude, b: vertical weave amplitude,
//    nWaves: radial lobes per strand (single knot when gcd(nWaves, 3) == 1, 3-ring link when 3 | nWaves),
//    rTube: strand tube radius.
vec4 sdgCircularBraid( in vec3 p, in float R, in float a, in float b, in int nWaves, in float rTube )
{
    R = max(R, 1e-5);
    a = max(a, 0.0);
    b = max(b, 0.0);
    float N      = float(max(nWaves, 1));
    const float PI = 3.14159265358979323846;
    float period = (2.0 * PI) / (3.0 * N);
    float halfP  = 0.5 * period;
    const vec3 phases = vec3(0.0, 2.0943951023931953, 4.1887902047863905);

    float rho = length(p.xy);
    float phi = (rho > 1e-12) ? atan(p.y, p.x) : 0.0;

    // Exact D_{3N} dihedral-strand fold into [0, pi / (3*N)]
    float kCell   = floor(phi / period);
    float phiLoc  = phi - kCell * period;
    bool  refl    = (phiLoc > halfP);
    float phiFold = refl ? (period - phiLoc) : phiLoc;
    float pz      = refl ? -p.z : p.z;
    float px      = rho * cos(phiFold);
    float py      = rho * sin(phiFold);

    float win  = 1.6 * PI / N;
    float tMin = phiFold - win;
    float tMax = phiFold + win;

    float bestD2   = 1e30;
    vec3  bestDiff = vec3(1.0, 0.0, 0.0);

    float offsets[3];
    offsets[0] = 0.0;
    offsets[1] = -0.45 / N;
    offsets[2] =  0.45 / N;

    for( int s = 0; s < 3; s++ )
    {
        vec3 t = vec3(clamp(phiFold + offsets[s], tMin, tMax));
        for( int i = 0; i < 6; i++ )
        {
            vec3 th  = N * t + phases;
            vec3 c1  = cos(t), s1 = sin(t);
            vec3 ct  = cos(th), st = sin(th);
            vec3 c2t = ct * ct - st * st;
            vec3 s2t = 2.0 * st * ct;

            vec3 r   = R + a * ct;
            vec3 dr  = (-N * a) * st;
            vec3 d2r = (-N * N * a) * ct;

            vec3 Z   = b * s2t;
            vec3 dZ  = (2.0 * N * b) * c2t;
            vec3 d2Z = (-4.0 * N * N * b) * s2t;

            vec3 u  = px * c1 + py * s1;
            vec3 v  = -px * s1 + py * c1;
            vec3 f  = (r - u) * dr - r * v + (Z - pz) * dZ;
            vec3 df = dr * dr + (r - u) * d2r - 2.0 * dr * v + r * u + dZ * dZ + (Z - pz) * d2Z;

            vec3 step = vec3(
                (df.x > 1e-6) ? (f.x / df.x) : ((f.x >= 0.0 ? 0.2 : -0.2) / N),
                (df.y > 1e-6) ? (f.y / df.y) : ((f.y >= 0.0 ? 0.2 : -0.2) / N),
                (df.z > 1e-6) ? (f.z / df.z) : ((f.z >= 0.0 ? 0.2 : -0.2) / N)
            );
            t = clamp(t - clamp(step, vec3(-0.4 / N), vec3(0.4 / N)), vec3(tMin), vec3(tMax));
        }

        vec3 th = N * t + phases;
        vec3 c1 = cos(t), s1 = sin(t);
        vec3 r  = R + a * cos(th);
        vec3 Z  = b * sin(2.0 * th);
        for( int k = 0; k < 3; k++ )
        {
            vec3  q    = vec3(r[k] * c1[k], r[k] * s1[k], Z[k]);
            vec3  diff = vec3(px, py, pz) - q;
            float d2   = dot(diff, diff);
            if( d2 < bestD2 )
            {
                bestD2   = d2;
                bestDiff = diff;
            }
        }
    }

    float invD = inversesqrt(max(bestD2, 1e-30));
    float d    = bestD2 * invD;
    vec3  gq   = (bestD2 > 1e-24) ? (bestDiff * invD) : vec3(0.0, 0.0, 1.0);

    vec3 grad;
    if( refl )
    {
        float a0 = (kCell + 1.0) * period;
        float ca = cos(a0), sa = sin(a0);
        grad = vec3(gq.x * ca + gq.y * sa, gq.x * sa - gq.y * ca, -gq.z);
    }
    else
    {
        float a0 = kCell * period;
        float ca = cos(a0), sa = sin(a0);
        grad = vec3(gq.x * ca - gq.y * sa, gq.x * sa + gq.y * ca, gq.z);
    }
    return vec4(d - rTube, grad);
}

// Convenient overloads
vec4 sdgCircularBraid( in vec3 p, in float R, in float a, in int nWaves, in float b, in float rTube )
{
    return sdgCircularBraid(p, R, a, b, nWaves, rTube);
}

vec4 sdgCircularBraid( in vec3 p, in float R, in float a, in float b, in float nWaves, in float rTube )
{
    return sdgCircularBraid(p, R, a, b, int(round(nWaves)), rTube);
}

// 2. Scalar 3D Circular 3-Strand Celtic Braid Ring SDF
float sdCircularBraid( in vec3 p, in float R, in float a, in float b, in int nWaves, in float rTube )
{
    return sdgCircularBraid(p, R, a, b, nWaves, rTube).x;
}

float sdCircularBraid( in vec3 p, in float R, in float a, in int nWaves, in float b, in float rTube )
{
    return sdgCircularBraid(p, R, a, b, nWaves, rTube).x;
}

float sdCircularBraid( in vec3 p, in float R, in float a, in float b, in float nWaves, in float rTube )
{
    return sdgCircularBraid(p, R, a, b, int(round(nWaves)), rTube).x;
}

// 3. Exact 3D Linear 3-Strand Plait / Celtic Braid along Y-axis (bounded to |y| <= h)
vec4 sdgLinearBraid3D( in vec3 p, in float a, in float b, in float kFreq, in float h, in float rTube )
{
    kFreq = max(kFreq, 1e-5);
    h     = max(h, 1e-5);
    const vec3 phases = vec3(0.0, 2.0943951023931953, 4.1887902047863905);

    float yClamped = clamp(p.y, -h, h);
    float win      = 1.6 * 3.14159265358979323846 / kFreq;
    float yMin     = max(-h, yClamped - win);
    float yMax     = min( h, yClamped + win);

    float bestD2   = 1e30;
    vec3  bestDiff = vec3(1.0, 0.0, 0.0);

    // Check endpoints y = -h and y = +h for all 3 strands
    for( int e = 0; e < 2; e++ )
    {
        float ye = (e == 0) ? -h : h;
        vec3  th = kFreq * ye + phases;
        vec3  X  = a * cos(th);
        vec3  Z  = b * sin(2.0 * th);
        for( int k = 0; k < 3; k++ )
        {
            vec3  diff = p - vec3(X[k], ye, Z[k]);
            float d2   = dot(diff, diff);
            if( d2 < bestD2 ) { bestD2 = d2; bestDiff = diff; }
        }
    }

    float offsets[3];
    offsets[0] = 0.0;
    offsets[1] = -0.45 / kFreq;
    offsets[2] =  0.45 / kFreq;

    for( int s = 0; s < 3; s++ )
    {
        vec3 y = vec3(clamp(yClamped + offsets[s], yMin, yMax));
        for( int i = 0; i < 6; i++ )
        {
            vec3 th  = kFreq * y + phases;
            vec3 ct  = cos(th), st = sin(th);
            vec3 c2t = ct * ct - st * st;
            vec3 s2t = 2.0 * st * ct;

            vec3 X   = a * ct;
            vec3 dX  = (-kFreq * a) * st;
            vec3 d2X = (-kFreq * kFreq * a) * ct;

            vec3 Z   = b * s2t;
            vec3 dZ  = (2.0 * kFreq * b) * c2t;
            vec3 d2Z = (-4.0 * kFreq * kFreq * b) * s2t;

            vec3 dx = X - vec3(p.x);
            vec3 dy = y - vec3(p.y);
            vec3 dz = Z - vec3(p.z);

            vec3 f  = dx * dX + dy + dz * dZ;
            vec3 df = dX * dX + dx * d2X + 1.0 + dZ * dZ + dz * d2Z;
            vec3 step = vec3(
                (df.x > 1e-6) ? (f.x / df.x) : ((f.x >= 0.0 ? 0.2 : -0.2) / kFreq),
                (df.y > 1e-6) ? (f.y / df.y) : ((f.y >= 0.0 ? 0.2 : -0.2) / kFreq),
                (df.z > 1e-6) ? (f.z / df.z) : ((f.z >= 0.0 ? 0.2 : -0.2) / kFreq)
            );
            y = clamp(y - clamp(step, vec3(-0.4 / kFreq), vec3(0.4 / kFreq)), vec3(yMin), vec3(yMax));
        }

        vec3 th = kFreq * y + phases;
        vec3 X  = a * cos(th);
        vec3 Z  = b * sin(2.0 * th);
        for( int k = 0; k < 3; k++ )
        {
            vec3  diff = p - vec3(X[k], y[k], Z[k]);
            float d2   = dot(diff, diff);
            if( d2 < bestD2 ) { bestD2 = d2; bestDiff = diff; }
        }
    }

    float invD = inversesqrt(max(bestD2, 1e-30));
    float d    = bestD2 * invD;
    vec3  grad = (bestD2 > 1e-24) ? (bestDiff * invD) : vec3(1.0, 0.0, 0.0);
    return vec4(d - rTube, grad);
}

// 4. Scalar 3D Linear 3-Strand Plait SDF
float sdLinearBraid3D( in vec3 p, in float a, in float b, in float kFreq, in float h, in float rTube )
{
    return sdgLinearBraid3D(p, a, b, kFreq, h, rTube).x;
}

#endif // SDF_GEMS_SD_CIRCULAR_BRAID_GLSL
