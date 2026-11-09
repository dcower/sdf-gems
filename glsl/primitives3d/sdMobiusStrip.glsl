/**
 * Exact 3D Euclidean Ruled Möbius Strip & n-Half-Twist Torus Ribbon
 * 
 * @author Gemini 4 Argon (research & mathematical design)
 * @author Gemini 3.1 Pro (experiments & implementation)
 */

#ifndef SDF_GEMS_HAS_MOBIUS_STRIP
#define SDF_GEMS_HAS_MOBIUS_STRIP

#include "../common/common.glsl"

vec4 sdgTwistedRibbon(vec3 p, float R, float w, float th, float nHalfTwists) {
    float k = 0.5 * nHalfTwists;
    float phi = atan(p.y, p.x);
    vec4 v = phi + vec4(0.0, -0.8, 0.8, 3.141592653589793);
    
    // Vectorized 4-seed Newton refinement using the Cylindrical Envelope-Theorem Ruling Reduction
    for (int i = 0; i < 4; i++) {
        vec4 cv = cos(v), sv = sin(v);
        vec4 alpha = k * v;
        vec4 ca = cos(alpha), sa = sin(alpha);
        
        vec4 P = p.x * cv + p.y * sv - R;
        vec4 Q = -p.x * sv + p.y * cv;
        vec4 U = P * ca + p.z * sa;
        vec4 N = -P * sa + p.z * ca;
        
        vec4 ustar = clamp(U, -w, w);
        vec4 Up = Q * ca + k * N;
        vec4 Upp = -(P + R) * ca - 2.0 * k * Q * sa - (k * k) * U;
        
        vec4 g = -R * Q - ustar * Up;
        vec4 U_less_w = step(abs(U), vec4(w));
        vec4 gp = R * (P + R) - ustar * Upp - U_less_w * Up * Up;
        
        vec4 den = max(abs(gp), vec4(1e-12)) * sign(gp + step(abs(gp), vec4(0.0)));
        vec4 step_v = clamp(g / den, -0.5, 0.5);
        v -= step_v;
    }
    
    // Evaluate exact distances for the refined seeds
    vec4 cv = cos(v), sv = sin(v);
    vec4 alpha = k * v;
    vec4 ca = cos(alpha), sa = sin(alpha);
    vec4 P = p.x * cv + p.y * sv - R;
    vec4 U = P * ca + p.z * sa;
    vec4 ustar = clamp(U, -w, w);
    
    vec4 qx = (R + ustar * ca) * cv;
    vec4 qy = (R + ustar * ca) * sv;
    vec4 qz = ustar * sa;
    
    vec4 dx = p.x - qx;
    vec4 dy = p.y - qy;
    vec4 dz = p.z - qz;
    
    vec4 d2 = dx*dx + dy*dy + dz*dz;
    
    // Find the seed that converged to the global minimum distance
    float min_d2 = d2.x;
    int min_idx = 0;
    if (d2.y < min_d2) { min_d2 = d2.y; min_idx = 1; }
    if (d2.z < min_d2) { min_d2 = d2.z; min_idx = 2; }
    if (d2.w < min_d2) { min_d2 = d2.w; min_idx = 3; }
    
    float best_v = (min_idx == 0) ? v.x : ((min_idx == 1) ? v.y : ((min_idx == 2) ? v.z : v.w));
    float best_ustar = (min_idx == 0) ? ustar.x : ((min_idx == 1) ? ustar.y : ((min_idx == 2) ? ustar.z : ustar.w));
    vec3 best_q = vec3(
        (min_idx == 0) ? qx.x : ((min_idx == 1) ? qx.y : ((min_idx == 2) ? qx.z : qx.w)),
        (min_idx == 0) ? qy.x : ((min_idx == 1) ? qy.y : ((min_idx == 2) ? qy.z : qy.w)),
        (min_idx == 0) ? qz.x : ((min_idx == 1) ? qz.y : ((min_idx == 2) ? qz.z : qz.w))
    );
    
    float d = sqrt(min_d2);
    vec3 grad;
    if (d > 1e-7) {
        grad = (p - best_q) / d;
    } else {
        // Analytical surface normal when exactly on the surface
        float b_ca = cos(k * best_v);
        float b_sa = sin(k * best_v);
        float b_cv = cos(best_v);
        float b_sv = sin(best_v);
        vec3 Su = vec3(b_ca * b_cv, b_ca * b_sv, b_sa);
        vec3 Sv = vec3(-(R + best_ustar * b_ca) * b_sv - k * best_ustar * b_sa * b_cv,
                        (R + best_ustar * b_ca) * b_cv - k * best_ustar * b_sa * b_sv,
                        k * best_ustar * b_ca);
        grad = normalize(cross(Su, Sv));
    }
    
    return vec4(d - th, grad);
}

float sdTwistedRibbon(vec3 p, float R, float w, float th, float nHalfTwists) {
    return sdgTwistedRibbon(p, R, w, th, nHalfTwists).x;
}

vec4 sdgMobiusStrip(vec3 p, float R, float w, float th) {
    return sdgTwistedRibbon(p, R, w, th, 1.0);
}

float sdMobiusStrip(vec3 p, float R, float w, float th) {
    return sdTwistedRibbon(p, R, w, th, 1.0);
}

#endif // SDF_GEMS_HAS_MOBIUS_STRIP
