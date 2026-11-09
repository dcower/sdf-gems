// ============================================================================
// sdConchoidNicomedes2D, sdgConchoidNicomedes2D, sdConchoidSpindle3D & sdgConchoidSpindle3D
// — Exact 0-Trig 2D Conchoid of Nicomedes Nodal Loop & 3D Conchoid Teardrop Spindle
//
// @author  Discovered and derived using Gemini 4 Argon (research & mathematical design)
//          and Gemini 3.8 Flash (implementation & verification) (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^2 / R^3; 0 trig in solver & 0-trig 0-sqrt polynomial interior test)
// @novelty Eliminates all trigonometric evaluations from Nicomedes' (c. 200 BC)
//          conchoid nodal loop r(theta) = a - b*sec(theta) (a > b > 0) via the
//          scaled half-angle rational parameterization s = tan(theta/2)/tan(theta_max/2) in [0, 1]:
//            q_x(s) = L * (1 - s^2) / (1 + k * s^2),
//            q_y(s) = 2 * sqrt(k) * L * s * (1 - s^2) / (1 - k^2 * s^4),
//          where L = a - b and k = (a - b) / (a + b), paired with the exact
//          0-trig 0-sqrt factored quartic interior test (x + b)^2 * y^2 < x^2 * (L - x) * (x + a + b).
// ============================================================================

#ifndef SDF_GEMS_SD_CONCHOID_NICOMEDES_2D_GLSL
#define SDF_GEMS_SD_CONCHOID_NICOMEDES_2D_GLSL

#include "../common/common.glsl"

// 1. Exact 0-Trig 2D Conchoid of Nicomedes Nodal Loop + Analytical Unit Gradient vec3(d, grad.xy)
//    Spans x in [0, a - b] with nodal crunode at (0, 0) and smooth rounded apex at (a - b, 0).
vec3 sdgConchoidNicomedes2D( in vec2 p, in float a, in float b )
{
    a = max(a, 1e-5);
    b = clamp(b, 1e-5, a - 1e-5);
    float L  = a - b;
    float k  = L / (a + b);
    float sk = sqrt(k);

    float px = p.x;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ay = abs(p.y);

    // Endpoints s = 0 -> apex (L, 0) and s = 1 -> crunode (0, 0)
    vec2  bestQ  = vec2(L, 0.0);
    float bestD2 = (px - L) * (px - L) + ay * ay;

    float d2_1 = px * px + ay * ay;
    if( d2_1 < bestD2 )
    {
        bestD2 = d2_1;
        bestQ  = vec2(0.0, 0.0);
    }

    // 4 SIMD analytical seeds in [0, 1]:
    // s0 (x-inversion seed), s1 (polar-secant seed), s2 (mid-body seed), s3 (crunode-neck seed)
    float xc = clamp(px, 0.0, L);
    float s0 = sqrt(clamp((L - xc) / max(L + k * xc, 1e-12), 0.0, 1.0));
    float rp = sqrt(px * px + ay * ay);
    float s1 = clamp(ay / max(sk * (max(px, 0.0) + rp), 1e-12), 0.0, 1.0);
    vec4  s  = vec4(s0, s1, 0.50, 0.85);

    for( int i = 0; i < 6; i++ )
    {
        vec4 s2   = s * s;
        vec4 ks2  = k * s2;
        vec4 invP = 1.0 / (1.0 + ks2);
        vec4 invM = 1.0 / max(vec4(1.0) - ks2, vec4(1e-7));

        vec4 u   = (1.0 - s2) * invP;
        vec4 du  = (-2.0 * (1.0 + k)) * s * (invP * invP);
        vec4 d2u = (2.0 * (1.0 + k)) * (3.0 * ks2 - 1.0) * (invP * invP * invP);

        vec4 v   = 2.0 * s * invM;
        vec4 dv  = 2.0 * (1.0 + ks2) * (invM * invM);
        vec4 d2v = (4.0 * k) * s * (3.0 + ks2) * (invM * invM * invM);

        vec4 qx   = L * u;
        vec4 dqx  = L * du;
        vec4 d2qx = L * d2u;

        vec4 qy   = (sk * L) * (u * v);
        vec4 dqy  = (sk * L) * (du * v + u * dv);
        vec4 d2qy = (sk * L) * (d2u * v + 2.0 * du * dv + u * d2v);

        vec4 dx = qx - vec4(px);
        vec4 dy = qy - vec4(ay);
        vec4 f  = dx * dqx + dy * dqy;
        vec4 df = dqx * dqx + dqy * dqy + dx * d2qx + dy * d2qy;

        vec4 step = vec4(
            (df.x > 1e-12) ? (f.x / df.x) : ((f.x >= 0.0) ? 0.15 : -0.15),
            (df.y > 1e-12) ? (f.y / df.y) : ((f.y >= 0.0) ? 0.15 : -0.15),
            (df.z > 1e-12) ? (f.z / df.z) : ((f.z >= 0.0) ? 0.15 : -0.15),
            (df.w > 1e-12) ? (f.w / df.w) : ((f.w >= 0.0) ? 0.15 : -0.15)
        );
        s = clamp(s - clamp(step, vec4(-0.25), vec4(0.25)), vec4(0.0), vec4(1.0));
    }

    for( int j = 0; j < 4; j++ )
    {
        float sj  = s[j];
        float s2j = sj * sj;
        float ksj = k * s2j;
        float uj  = (1.0 - s2j) / (1.0 + ksj);
        float vj  = 2.0 * sj / max(1.0 - ksj, 1e-7);
        vec2  qj  = vec2(L * uj, sk * L * uj * vj);
        vec2  dv  = vec2(px, ay) - qj;
        float d2  = dot(dv, dv);
        if( d2 < bestD2 )
        {
            bestD2 = d2;
            bestQ  = qj;
        }
    }

    // Exact 0-trig 0-sqrt factored quartic interior test:
    bool inside = (px > 0.0) && (px < L) &&
                  ((px + b) * (px + b) * (ay * ay) < (px * px) * (L - px) * (px + a + b));
    float sgn  = inside ? -1.0 : 1.0;
    float invD = inversesqrt(max(bestD2, 1e-30));
    float d    = bestD2 * invD;
    vec2  diff = vec2(px - bestQ.x, sy * (ay - bestQ.y));
    vec2  grad = (bestD2 > 1e-24) ? (sgn * diff * invD) : vec2(1.0, 0.0);
    return vec3(sgn * d, grad);
}

// 2. Scalar 2D Conchoid of Nicomedes Nodal Loop SDF
float sdConchoidNicomedes2D( in vec2 p, in float a, in float b )
{
    return sdgConchoidNicomedes2D(p, a, b).x;
}

// 3. Exact 3D Conchoid Teardrop Spindle of Revolution + Analytical Unit Gradient vec4(d, grad.xyz)
//    Revolves the Nicomedes conchoid loop around its symmetry axis (the x-axis).
vec4 sdgConchoidSpindle3D( in vec3 p, in float a, in float b )
{
    float rho = length(p.yz);
    vec3  r2d = sdgConchoidNicomedes2D(vec2(p.x, rho), a, b);
    vec2  er  = (rho > 1e-12) ? (p.yz / rho) : vec2(1.0, 0.0);
    return vec4(r2d.x, r2d.y, r2d.z * er.x, r2d.z * er.y);
}

// 4. Scalar 3D Conchoid Teardrop Spindle SDF
float sdConchoidSpindle3D( in vec3 p, in float a, in float b )
{
    return sdgConchoidNicomedes2D(vec2(p.x, length(p.yz)), a, b).x;
}

#endif // SDF_GEMS_SD_CONCHOID_NICOMEDES_2D_GLSL
