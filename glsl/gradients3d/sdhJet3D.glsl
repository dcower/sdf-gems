#ifndef SDF_GEMS_HAS_SDHJET3D
#define SDF_GEMS_HAS_SDHJET3D

#include "../common/common.glsl"

// @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)

#ifndef SDF_GEMS_HAS_SDFJET3D
#define SDF_GEMS_HAS_SDFJET3D
// 2-Jet encoding signed distance, gradient, and 2 ambient principal curvatures.
struct SdfJet3D {
    vec4 dg; // dg.x = distance, dg.yzw = gradient
    vec2 k;  // k.x = k1, k.y = k2 (principal curvatures at p, k1 <= k2)
};
#endif

float sdhMeanCurvature(SdfJet3D j) {
    return 0.5 * (j.k.x + j.k.y);
}

float sdhGaussianCurvature(SdfJet3D j) {
    return j.k.x * j.k.y;
}

SdfJet3D sdhSphere(vec3 p, float r) {
    float l = length(p);
    vec3 grad = l > 0.0 ? p / l : vec3(1.0, 0.0, 0.0);
    float k = l > 1e-10 ? 1.0 / l : 1e10;
    return SdfJet3D(vec4(l - r, grad), vec2(k, k));
}

SdfJet3D sdhPlane(vec3 p, vec3 n, float h) {
    return SdfJet3D(vec4(dot(p, n) - h, n), vec2(0.0));
}

SdfJet3D sdhBox(vec3 p, vec3 b) {
    vec3 q = abs(p) - b;
    vec3 w = max(q, 0.0);
    float l = length(w);
    float mind = min(max(q.x, max(q.y, q.z)), 0.0);
    float dist = l + mind;
    
    vec3 grad;
    vec2 k = vec2(0.0);
    vec3 s = sign(p);
    s = mix(s, vec3(1.0), equal(s, vec3(0.0)));
    
    if (mind < 0.0) {
        if (q.x >= q.y && q.x >= q.z) {
            grad = vec3(s.x, 0.0, 0.0);
        } else if (q.y >= q.x && q.y >= q.z) {
            grad = vec3(0.0, s.y, 0.0);
        } else {
            grad = vec3(0.0, 0.0, s.z);
        }
    } else {
        grad = s * (w / (l > 0.0 ? l : 1.0));
        int m = int(q.x > 0.0) + int(q.y > 0.0) + int(q.z > 0.0);
        if (m == 2) {
            k = vec2(0.0, 1.0 / max(l, 1e-10));
        } else if (m == 3) {
            k = vec2(1.0 / max(l, 1e-10));
        }
    }
    return SdfJet3D(vec4(dist, grad), k);
}

SdfJet3D sdhRoundBox(vec3 p, vec3 b, float r) {
    SdfJet3D j = sdhBox(p, b);
    j.dg.x -= r;
    return j;
}

SdfJet3D sdhCylinder(vec3 p, float r, float h) {
    float rho = length(p.xz);
    vec2 d = vec2(rho - r, abs(p.y) - h);
    float dist = min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
    
    vec3 grad;
    vec2 k = vec2(0.0);
    vec3 s = vec3(sign(p.x), sign(p.y), sign(p.z));
    s = mix(s, vec3(1.0), equal(s, vec3(0.0)));
    
    if (d.x <= 0.0 && d.y <= 0.0) {
        if (d.x > d.y) {
            grad = vec3((rho > 0.0 ? p.x / rho : 1.0), 0.0, (rho > 0.0 ? p.z / rho : 0.0));
            k = vec2(0.0, rho > 0.0 ? 1.0 / rho : 1e10);
        } else {
            grad = vec3(0.0, s.y, 0.0);
            k = vec2(0.0);
        }
    } else if (d.x > 0.0 && d.y <= 0.0) { // Lateral barrel
        float inv_rho = (rho > 1e-15) ? (1.0 / rho) : 0.0;
        grad = (rho > 1e-15) ? vec3(p.x * inv_rho, 0.0, p.z * inv_rho) : vec3(1.0, 0.0, 0.0);
        k = vec2(0.0, (rho > 1e-10) ? inv_rho : 1e10);
    } else if (d.x <= 0.0 && d.y > 0.0) { // Top/bottom cap
        grad = vec3(0.0, s.y, 0.0);
        k = vec2(0.0);
    } else { // Circular rim
        float drim = max(length(d), 1e-15);
        float inv_rho = (rho > 1e-15) ? (1.0 / rho) : 0.0;
        grad = vec3(p.x * inv_rho * (d.x / drim), s.y * (d.y / drim), p.z * inv_rho * (d.x / drim));
        k = vec2(d.x * inv_rho / drim, 1.0 / max(drim, 1e-10));
    }
    
    // Ensure k1 <= k2
    if (k.x > k.y) k = k.yx;
    
    return SdfJet3D(vec4(dist, grad), k);
}

SdfJet3D sdhTorus(vec3 p, vec2 t) {
    float rho = length(p.xz);
    vec2 q = vec2(rho - t.x, p.y);
    float u = length(q);
    float dist = u - t.y;
    
    vec2 dir_xz = (rho > 1e-15) ? (p.xz / rho) : vec2(1.0, 0.0);
    vec2 dir_q  = (u > 1e-15)   ? (q / u)       : vec2(1.0, 0.0);
    vec3 grad   = vec3(dir_xz.x * dir_q.x, dir_q.y, dir_xz.y * dir_q.x);
    
    float k_merid = 1.0 / max(u, 1e-10);
    float k_paral = (rho > 1e-10) ? ((rho - t.x) / (rho * max(u, 1e-10))) : 1e10;
    
    return SdfJet3D(vec4(dist, grad), vec2(min(k_merid, k_paral), max(k_merid, k_paral)));
}

SdfJet3D sdhOffset(SdfJet3D j, float r) {
    j.dg.x -= r;
    return j;
}

SdfJet3D sdhScale(SdfJet3D j, float s) {
    j.dg.x *= s;
    j.k /= s;
    return j;
}

SdfJet3D sdhUnion(SdfJet3D a, SdfJet3D b) {
    return a.dg.x < b.dg.x ? a : b;
}

SdfJet3D sdhSmoothUnion(SdfJet3D a, SdfJet3D b, float k) {
    float h = clamp(0.5 + 0.5 * (b.dg.x - a.dg.x) / max(k, 1e-10), 0.0, 1.0);
    float d = mix(b.dg.x, a.dg.x, h) - k * h * (1.0 - h);
    
    vec3 gm = mix(b.dg.yzw, a.dg.yzw, h);
    float l_gm = length(gm);
    vec3 grad = (l_gm > 1e-15) ? (gm / l_gm) : a.dg.yzw;
    
    // The blending fillet introduces negative normal curvature
    vec3 diff = b.dg.yzw - a.dg.yzw;
    float blend_k = -dot(diff, diff) / (2.0 * max(k, 1e-10));
    
    // For a smooth union, we approximate the principal curvatures.
    // In the blending region, we take the mix of the base curvatures,
    // but along the gradient difference direction, there is a strong negative curvature.
    // As a rough approximation for now:
    vec2 base_k = mix(b.k, a.k, h);
    
    // To be conservative and approximate the saddle shape:
    // One principal curvature becomes very negative (the blend_k),
    // the other is roughly the max of the base curvatures.
    // Outside the blend region (h=0 or 1), blend_k is 0 but we want the exact base.
    float weight = 4.0 * h * (1.0 - h); // 1 at center of blend, 0 at edges
    float k1 = min(base_k.x, blend_k * weight + mix(b.k.x, a.k.x, h));
    float k2 = mix(b.k.y, a.k.y, h);
    
    return SdfJet3D(vec4(d, grad), vec2(k1, k2));
}

float sdhOsculatingStep(SdfJet3D j, vec3 rd) {
    float d = j.dg.x;
    if (d <= 0.0) return d;
    vec3 n = j.dg.yzw;
    float vn = dot(rd, n);
    if (vn >= -1e-6) return d; // ray moving away or parallel to level set

    // Supporting-plane step (provably conservative for all convex regions and k == 0):
    float t_plane = d / (-vn);

    // If curvature is non-positive (flat plane / face region), step directly to the supporting plane:
    float k_max = j.k.y;
    if (k_max < 0.0) return d; // saddle/concave blend region: safe fallback to d
    if (k_max <= 1e-6) return min(t_plane, 4.0 * d);

    // Exact ray-osculating-sphere intersection via Citardauq's formula (0 division by zero as k -> 0):
    float kd = k_max * d;
    float c_term = d * max(1.0 - 0.5 * kd, 0.5);
    float disc = vn * vn - 2.0 * k_max * c_term;
    if (disc < 0.0) return min(t_plane, 3.0 * d);

    float t_osc = (2.0 * c_term) / (-vn + sqrt(disc));
    return clamp(min(t_osc, t_plane), d, 8.0 * d);
}
#endif
