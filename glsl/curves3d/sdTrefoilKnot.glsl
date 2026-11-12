#ifndef SDF_GEMS_HAS_TREFOIL_KNOT
#define SDF_GEMS_HAS_TREFOIL_KNOT

#include "../common/common.glsl"

/**
 * Exact 3D Euclidean SDF and Analytical Unit Gradient of the general (p, q) Torus Knot.
 * @author Gemini 4 Argon & Gemini 3.1 Pro
 */
vec4 sdgTorusKnotExact(vec3 p, float R, float r0, float rTube, int pWinds, int qWinds) {
    float rho = length(p.xy);
    float phi = (rho > 1e-12) ? atan(p.y, p.x) : 0.0;

    float p_w = float(pWinds);
    float q_w = float(qWinds);
    
    // Safety clamp for bounds in shader
    int n_strands = clamp(pWinds, 1, 4);

    float best_d2 = 1e30;
    float best_t = 0.0;

    float shift = 3.141592653589793 / max(q_w, 1.0);
    float inv_p = 1.0 / p_w;
    
    // We check n_strands base seeds, plus shift and -shift for each.
    for (int k = 0; k < 4; ++k) {
        if (k >= n_strands) break;
        
        float base_seed = (phi + 6.283185307179586 * float(k)) * inv_p;
        
        float seeds[3];
        seeds[0] = base_seed;
        seeds[1] = base_seed + shift;
        seeds[2] = base_seed - shift;
        
        for (int s = 0; s < 3; ++s) {
            float t = seeds[s];
            
            for (int i = 0; i < 5; ++i) {
                float qt = q_w * t;
                float al = p_w * t - phi;
                float cqt = cos(qt), sqt = sin(qt);
                float cal = cos(al), sal = sin(al);

                float W = R + r0 * cqt;
                float dW = -q_w * r0 * sqt;
                float d2W = -q_w * q_w * r0 * cqt;

                float g = -q_w * r0 * (R * sqt + p.z * cqt)
                          - rho * (dW * cal - p_w * W * sal);
                float gp = -q_w * q_w * r0 * (R * cqt - p.z * sqt)
                           - rho * ((d2W - p_w * p_w * W) * cal - 2.0 * p_w * dW * sal);
                           
                float denom = max(abs(gp), 1e-12);
                float step = clamp(g / denom, -0.35 / max(q_w, p_w), 0.35 / max(q_w, p_w));
                t -= step;
            }
            
            float qt_final = q_w * t;
            float pt_final = p_w * t;
            float W_final = R + r0 * cos(qt_final);
            vec3 q_pt = vec3(W_final * cos(pt_final), W_final * sin(pt_final), r0 * sin(qt_final));
            
            vec3 diff = p - q_pt;
            float d2 = dot(diff, diff);
            if (d2 < best_d2) {
                best_d2 = d2;
                best_t = t;
            }
        }
    }

    float qt_best = q_w * best_t;
    float pt_best = p_w * best_t;
    float W_best = R + r0 * cos(qt_best);
    vec3 q_best = vec3(W_best * cos(pt_best), W_best * sin(pt_best), r0 * sin(qt_best));
    vec3 diff = p - q_best;
    
    float d = max(sqrt(best_d2), 1e-15);
    // sdg format is vec4(d, grad.x, grad.y, grad.z)
    return vec4(d - rTube, diff / d);
}

float sdTorusKnotExact(vec3 p, float R, float r0, float rTube, int pWinds, int qWinds) {
    return sdgTorusKnotExact(p, R, r0, rTube, pWinds, qWinds).x;
}

vec4 sdgTrefoilKnot(vec3 p, float R, float r0, float rTube) {
    return sdgTorusKnotExact(p, R, r0, rTube, 2, 3);
}

float sdTrefoilKnot(vec3 p, float R, float r0, float rTube) {
    return sdgTrefoilKnot(p, R, r0, rTube).x;
}

#endif
