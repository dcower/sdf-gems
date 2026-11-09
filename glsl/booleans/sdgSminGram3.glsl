#ifndef SDF_GEMS_HAS_SDGSMINGRAM3
#define SDF_GEMS_HAS_SDGSMINGRAM3

#include "../common/common.glsl"
#include "sminVaradhan.glsl"
#include "../gradients3d/sdhJet3D.glsl"

// @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)

#ifndef SDF_GEMS_HAS_SDFJET3D
#define SDF_GEMS_HAS_SDFJET3D
// 2-Jet encoding signed distance, gradient, and 2 ambient principal curvatures.
struct SdfJet3D {
    vec4 dg; // dg.x = distance, dg.yzw = gradient
    vec2 k;  // k.x = k1, k.y = k2 (principal curvatures at p, k1 <= k2)
};
#endif

vec4 sdgSminGram( in vec4 a, in vec4 b, in float r );
vec4 sdgSmaxGram( in vec4 a, in vec4 b, in float r );

// Exact 3D Euclidean rolling-ball smooth union of three 1-jets at arbitrary dihedral angles.
vec4 sdgSminGram3(vec4 a, vec4 b, vec4 c, float r) {
    float ua = r - a.x;
    float ub = r - b.x;
    float uc = r - c.x;

    if (ua <= 0.0 && ub <= 0.0 && uc <= 0.0) {
        vec4 m = (a.x < b.x) ? a : b;
        return (m.x < c.x) ? m : c;
    }

    float cab = dot(a.yzw, b.yzw);
    float cbc = dot(b.yzw, c.yzw);
    float cca = dot(c.yzw, a.yzw);
    float detG = 1.0 + 2.0 * cab * cbc * cca - (cab * cab + cbc * cbc + cca * cca);

    if (detG > 1e-6) {
        float wa_num = (1.0 - cbc * cbc) * ua + (cca * cbc - cab) * ub + (cab * cbc - cca) * uc;
        float wb_num = (cca * cbc - cab) * ua + (1.0 - cca * cca) * ub + (cab * cca - cbc) * uc;
        float wc_num = (cab * cbc - cca) * ua + (cab * cca - cbc) * ub + (1.0 - cab * cab) * uc;

        if (wa_num > 0.0 && wb_num > 0.0 && wc_num > 0.0) {
            float inv_det = 1.0 / detG;
            float alpha = wa_num * inv_det;
            float beta  = wb_num * inv_det;
            float gamma = wc_num * inv_det;
            float rho2 = alpha * ua + beta * ub + gamma * uc;
            float rho = sqrt(max(rho2, 1e-14));
            vec3 g = (alpha * a.yzw + beta * b.yzw + gamma * c.yzw) / rho;
            return vec4(r - rho, g);
        }
    }

    vec4 ab = sdgSminGram(a, b, r);
    vec4 bc = sdgSminGram(b, c, r);
    vec4 ca = sdgSminGram(c, a, r);
    
    vec4 m = (ab.x < bc.x) ? ab : bc;
    return (m.x < ca.x) ? m : ca;
}

// Exact 3D Euclidean rolling-ball smooth intersection of three 1-jets at arbitrary dihedral angles.
vec4 sdgSmaxGram3(vec4 a, vec4 b, vec4 c, float r) {
    float ua = r + a.x;
    float ub = r + b.x;
    float uc = r + c.x;

    if (ua <= 0.0 && ub <= 0.0 && uc <= 0.0) {
        vec4 m = (a.x > b.x) ? a : b;
        return (m.x > c.x) ? m : c;
    }

    float cab = dot(a.yzw, b.yzw);
    float cbc = dot(b.yzw, c.yzw);
    float cca = dot(c.yzw, a.yzw);
    float detG = 1.0 + 2.0 * cab * cbc * cca - (cab * cab + cbc * cbc + cca * cca);

    if (detG > 1e-6) {
        float wa_num = (1.0 - cbc * cbc) * ua + (cca * cbc - cab) * ub + (cab * cbc - cca) * uc;
        float wb_num = (cca * cbc - cab) * ua + (1.0 - cca * cca) * ub + (cab * cca - cbc) * uc;
        float wc_num = (cab * cbc - cca) * ua + (cab * cca - cbc) * ub + (1.0 - cab * cab) * uc;

        if (wa_num > 0.0 && wb_num > 0.0 && wc_num > 0.0) {
            float inv_det = 1.0 / detG;
            float alpha = wa_num * inv_det;
            float beta  = wb_num * inv_det;
            float gamma = wc_num * inv_det;
            float rho2 = alpha * ua + beta * ub + gamma * uc;
            float rho = sqrt(max(rho2, 1e-14));
            vec3 g = (alpha * a.yzw + beta * b.yzw + gamma * c.yzw) / rho;
            return vec4(rho - r, g);
        }
    }

    vec4 ab = sdgSmaxGram(a, b, r);
    vec4 bc = sdgSmaxGram(b, c, r);
    vec4 ca = sdgSmaxGram(c, a, r);
    
    vec4 m = (ab.x > bc.x) ? ab : bc;
    return (m.x > ca.x) ? m : ca;
}

// Exact 2-jet Weingarten curvature-propagating 2-way rolling-ball fillet.
SdfJet3D sdhSminGram(SdfJet3D a, SdfJet3D b, float r) {
    float c = dot(a.dg.yzw, b.dg.yzw);
    float sin2 = 1.0 - c * c;
    if (sin2 > 1e-6) {
        float ua = r - a.dg.x;
        float ub = r - b.dg.x;
        if ((ua - ub * c > 0.0) && (ub - ua * c > 0.0)) {
            float du = 2.0 + 2.0 * c;
            float dv = 2.0 - 2.0 * c;
            float su = a.dg.x + b.dg.x - 2.0 * r;
            float sv = a.dg.x - b.dg.x;
            float rho2 = (su * su) / du + (sv * sv) / dv;
            float rho = sqrt(max(rho2, 1e-14));
            float inv_rho = 1.0 / rho;
            float cu = (su / du) * inv_rho;
            float cv = (sv / dv) * inv_rho;
            float wa = -(cu + cv);
            float wb = -(cu - cv);
            
            vec3 grad = wa * a.dg.yzw + wb * b.dg.yzw;
            
            float k1 = -inv_rho + wa * a.k.x + wb * b.k.x;
            float k2 = wa * a.k.y + wb * b.k.y;
            vec2 k = vec2(min(k1, k2), max(k1, k2));
            
            return SdfJet3D(vec4(r - rho, grad), k);
        }
    }
    return (a.dg.x < b.dg.x) ? a : b;
}

// Exact 2-jet Weingarten curvature-propagating 3-way trihedral rolling-ball fillet.
SdfJet3D sdhSminGram3(SdfJet3D a, SdfJet3D b, SdfJet3D c, float r) {
    float ua = r - a.dg.x;
    float ub = r - b.dg.x;
    float uc = r - c.dg.x;

    if (ua <= 0.0 && ub <= 0.0 && uc <= 0.0) {
        SdfJet3D m = (a.dg.x < b.dg.x) ? a : b;
        return (m.dg.x < c.dg.x) ? m : c;
    }

    float cab = dot(a.dg.yzw, b.dg.yzw);
    float cbc = dot(b.dg.yzw, c.dg.yzw);
    float cca = dot(c.dg.yzw, a.dg.yzw);
    float detG = 1.0 + 2.0 * cab * cbc * cca - (cab * cab + cbc * cbc + cca * cca);

    if (detG > 1e-6) {
        float wa_num = (1.0 - cbc * cbc) * ua + (cca * cbc - cab) * ub + (cab * cbc - cca) * uc;
        float wb_num = (cca * cbc - cab) * ua + (1.0 - cca * cca) * ub + (cab * cca - cbc) * uc;
        float wc_num = (cab * cbc - cca) * ua + (cab * cca - cbc) * ub + (1.0 - cab * cab) * uc;

        if (wa_num > 0.0 && wb_num > 0.0 && wc_num > 0.0) {
            float inv_det = 1.0 / detG;
            float alpha = wa_num * inv_det;
            float beta  = wb_num * inv_det;
            float gamma = wc_num * inv_det;
            float rho2 = alpha * ua + beta * ub + gamma * uc;
            float rho = sqrt(max(rho2, 1e-14));
            float inv_rho = 1.0 / rho;
            
            float wa = alpha * inv_rho;
            float wb = beta * inv_rho;
            float wc = gamma * inv_rho;
            
            vec3 grad = wa * a.dg.yzw + wb * b.dg.yzw + wc * c.dg.yzw;
            
            float k_base = -inv_rho;
            float k_prop_1 = wa * a.k.x + wb * b.k.x + wc * c.k.x;
            float k_prop_2 = wa * a.k.y + wb * b.k.y + wc * c.k.y;
            float k1 = k_base + k_prop_1;
            float k2 = k_base + k_prop_2;
            vec2 k = vec2(min(k1, k2), max(k1, k2));
            
            return SdfJet3D(vec4(r - rho, grad), k);
        }
    }

    SdfJet3D ab = sdhSminGram(a, b, r);
    SdfJet3D bc = sdhSminGram(b, c, r);
    SdfJet3D ca = sdhSminGram(c, a, r);
    
    SdfJet3D m = (ab.dg.x < bc.dg.x) ? ab : bc;
    return (m.dg.x < ca.dg.x) ? m : ca;
}

#endif
