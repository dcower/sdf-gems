#ifndef SDF_GEMS_HAS_OP_HOPF_PULLBACK
#define SDF_GEMS_HAS_OP_HOPF_PULLBACK

#include "../common/common.glsl"

/**
 * @name opHopfMap
 * @category deform
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 * @summary Stereographic projection from R3 to S3, followed by the Hopf fibration to S2.
 */
vec3 opHopfMap(vec3 p, float scale) {
    p /= scale;
    float r2 = dot(p, p);
    float d = 1.0 - r2;
    vec2 z1 = 2.0 * p.xy;
    vec2 z2 = vec2(2.0 * p.z, d);
    float d2 = (1.0 + r2) * (1.0 + r2);
    // n = (2 Re(z1 z2*), 2 Im(z1 z2*), |z1|^2 - |z2|^2)
    float nx = 2.0 * (z1.x * z2.x + z1.y * z2.y);
    float ny = 2.0 * (z1.y * z2.x - z1.x * z2.y);
    float nz = dot(z1, z1) - dot(z2, z2);
    return vec3(nx, ny, nz) / d2;
}

/**
 * @name opHopfPullback
 * @category deform
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 * @summary Pulls back a 1-Lipschitz spherical SDF on S2 into a conservative Euclidean SDF in R3 via the Hopf fibration.
 */
float opHopfPullback(vec3 p, float dS2, float rTube, float scale) {
    p /= scale;
    float r2 = dot(p, p);
    float D_naive = dS2 * (1.0 + r2) * 0.25;
    float d = D_naive / (1.0 + D_naive * inversesqrt(1.0 + r2));
    return (d - rTube) * scale;
}

/**
 * @name sdHopfFiber
 * @category deform
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 * @summary Exact Euclidean SDF to a Hopf fiber (a Villarceau circle in R3) over a base point q on S2.
 */
float sdHopfFiber(vec3 p, vec3 q, float rTube, float scale) {
    p /= scale;
    float k = 1.0 / max(1.0 + q.z, 1e-7);
    vec3 c = vec3(k * q.y, -k * q.x, 0.0);
    float r = sqrt(2.0 * k);
    vec3 n = vec3(-q.x, -q.y, 1.0 + q.z) * sqrt(k * 0.5);
    
    vec3 pc = p - c;
    float dot_n = dot(pc, n);
    vec3 proj = pc - dot_n * n;
    
    vec2 q2 = vec2(length(proj) - r, dot_n);
    return (length(q2) - rTube) * scale;
}

/**
 * @name sdgHopfFiber
 * @category deform
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 * @summary Analytical gradient for sdHopfFiber.
 */
vec4 sdgHopfFiber(vec3 p, vec3 q, float rTube, float scale) {
    p /= scale;
    float k = 1.0 / max(1.0 + q.z, 1e-7);
    vec3 c = vec3(k * q.y, -k * q.x, 0.0);
    float r = sqrt(2.0 * k);
    vec3 n = vec3(-q.x, -q.y, 1.0 + q.z) * sqrt(k * 0.5);
    
    vec3 pc = p - c;
    float dot_n = dot(pc, n);
    vec3 proj = pc - dot_n * n;
    
    float l = length(proj);
    vec2 q2 = vec2(l - r, dot_n);
    float d = length(q2);
    
    vec3 d_proj;
    if (l > 1e-15) {
        d_proj = proj / l;
    } else {
        vec3 t = (abs(n.z) < 0.999) ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
        d_proj = normalize(cross(n, t));
    }
    
    vec3 grad;
    if (d > 0.0) {
        grad = (q2.x * d_proj + q2.y * n) / d;
    } else {
        grad = vec3(0.0, 1.0, 0.0); // inside the 1D circle? Not possible unless rTube > 0. But d is distance to circle (1D), so d >= 0 always.
    }
    
    return vec4((d - rTube) * scale, grad);
}

/**
 * @name sdHopfLink
 * @category deform
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 * @summary Exact Euclidean SDF of linked Hopf Villarceau circles (a Hopf link).
 */
float sdHopfLink(vec3 p, float sep, float rTube, float scale) {
    float sa = sin(sep);
    float ca = cos(sep);
    float d1 = sdHopfFiber(p, vec3(sa, 0.0, ca), rTube, scale);
    float d2 = sdHopfFiber(p, vec3(-sa, 0.0, ca), rTube, scale);
    return min(d1, d2);
}

/**
 * @name sdHopfRingFibration
 * @category deform
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 * @summary Exact Euclidean SDF to a bundle of N linked Villarceau tori arranged symmetrically on the Hopf equator.
 */
float sdHopfRingFibration(vec3 p, float nRings, float rTube, float scale) {
    vec3 n = opHopfMap(p, scale);
    float angle = atan(n.y, n.x);
    float sector = 6.28318530718 / max(nRings, 1.0);
    float id = floor(angle / sector + 0.5);
    
    float d = 1e20;
    for (float i = -1.0; i <= 1.0; i += 1.0) {
        float a = (id + i) * sector;
        vec3 q = vec3(cos(a), sin(a), 0.0);
        d = min(d, sdHopfFiber(p, q, rTube, scale));
        d = min(d, sdHopfFiber(p, -q, rTube, scale));
    }
    return d;
}

#endif
