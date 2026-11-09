/**
 * Exact Catenoid, Ruled Helicoid, and Hyperbolic-Secant Scherk Saddle Tower
 * @author Gemini 3.1 Pro (experiments & implementation) & Gemini 4 Argon (research & mathematical design)
 */

#ifndef SDF_GEMS_HAS_MINIMAL_SURFACE
#define SDF_GEMS_HAS_MINIMAL_SURFACE

#include "../common/common.glsl"

// --- 1. Exact Bounded Catenoid ---

float sdCatenoid(vec3 p, float c, float h, float th) {
    float rho = length(p.xy);
    float z = abs(p.z);
    
    float c_safe = max(c, 1e-6);
    float u = clamp(asinh(max(rho, z) / c_safe), 0.0, h / c_safe);
    
    for (int i = 0; i < 4; ++i) {
        float sinh_u = sinh(u);
        float cosh_u = cosh(u);
        float g = c * (u - z / c_safe) + (c * cosh_u - rho) * sinh_u;
        float g_prime = c * (1.0 + cosh(2.0 * u)) - rho * cosh_u;
        u = clamp(u - g / g_prime, 0.0, h / c_safe);
    }
    
    float u0 = 0.0;
    float uh = h / c_safe;
    
    float dx_root = c * cosh(u) - rho;
    float dz_root = c * u - z;
    float d2_root = dx_root * dx_root + dz_root * dz_root;
    
    float dx_0 = c * cosh(u0) - rho;
    float dz_0 = c * u0 - z;
    float d2_0 = dx_0 * dx_0 + dz_0 * dz_0;
    
    float dx_h = c * cosh(uh) - rho;
    float dz_h = c * uh - z;
    float d2_h = dx_h * dx_h + dz_h * dz_h;
    
    float d2 = min(min(d2_root, d2_0), d2_h);
    return sqrt(d2) - th;
}

vec4 sdgCatenoid(vec3 p, float c, float h, float th) {
    float rho = length(p.xy);
    float z = abs(p.z);
    
    float c_safe = max(c, 1e-6);
    float u = clamp(asinh(max(rho, z) / c_safe), 0.0, h / c_safe);
    
    for (int i = 0; i < 4; ++i) {
        float sinh_u = sinh(u);
        float cosh_u = cosh(u);
        float g = c * (u - z / c_safe) + (c * cosh_u - rho) * sinh_u;
        float g_prime = c * (1.0 + cosh(2.0 * u)) - rho * cosh_u;
        u = clamp(u - g / g_prime, 0.0, h / c_safe);
    }
    
    vec3 u_c = vec3(u, 0.0, h / c_safe);
    vec3 dx_c = c * cosh(u_c) - vec3(rho);
    vec3 dz_c = c * u_c - vec3(z);
    vec3 d2_c = dx_c * dx_c + dz_c * dz_c;
    
    float d2 = d2_c.x;
    float u_best = u_c.x;
    if (d2_c.y < d2) { d2 = d2_c.y; u_best = u_c.y; }
    if (d2_c.z < d2) { d2 = d2_c.z; u_best = u_c.z; }
    
    float drho = rho - c * cosh(u_best);
    float dz = z - c * u_best;
    float d = sqrt(d2);
    
    vec2 n_rz;
    if (d > 1e-7) {
        n_rz = vec2(drho, dz) / d;
    } else {
        n_rz = vec2(1.0, -sinh(u_best)) / cosh(u_best);
    }
    
    vec3 grad = vec3(
        (rho > 1e-7 ? p.x / rho : 1.0) * n_rz.x,
        (rho > 1e-7 ? p.y / rho : 0.0) * n_rz.x,
        (p.z >= 0.0 ? 1.0 : -1.0) * n_rz.y
    );
    
    return vec4(d - th, grad);
}

// --- 2. Exact Ruled Helicoid ---

float sdHelicoid(vec3 p, float c, float R, float th) {
    float c_safe = max(c, 1e-6);
    float v0 = p.z / c_safe;
    float v_atan = atan(p.y, p.x);
    float pi = 3.14159265359;
    
    float v_seed1 = v_atan + pi * round((v0 - v_atan) / pi);
    float v_seed2 = v0;
    float v_seed3 = v_seed1 + pi * (v0 > v_seed1 ? 1.0 : -1.0);
    
    vec3 v = vec3(v_seed1, v_seed2, v_seed3);
    
    for (int i = 0; i < 4; ++i) {
        vec3 sin_v = sin(v);
        vec3 cos_v = cos(v);
        vec3 w = p.x * cos_v + p.y * sin_v;
        vec3 q = -p.x * sin_v + p.y * cos_v;
        vec3 u_star = clamp(w, -R, R);
        
        vec3 g = -u_star * q + c * c * (v - v0);
        vec3 du_dv = vec3(
            abs(w.x) < R ? -q.x : 0.0,
            abs(w.y) < R ? -q.y : 0.0,
            abs(w.z) < R ? -q.z : 0.0
        );
        
        vec3 g_prime = -du_dv * q + u_star * w + c * c;
        v -= g / g_prime;
        v = clamp(v, v0 - pi, v0 + pi);
    }
    
    vec3 sin_v = sin(v);
    vec3 cos_v = cos(v);
    vec3 w = p.x * cos_v + p.y * sin_v;
    vec3 u_star = clamp(w, -R, R);
    
    vec3 dx = vec3(p.x) - u_star * cos_v;
    vec3 dy = vec3(p.y) - u_star * sin_v;
    vec3 dz = c * (vec3(v0) - v);
    
    vec3 d2 = dx*dx + dy*dy + dz*dz;
    return sqrt(min(min(d2.x, d2.y), d2.z)) - th;
}

vec4 sdgHelicoid(vec3 p, float c, float R, float th) {
    float c_safe = max(c, 1e-6);
    float v0 = p.z / c_safe;
    float v_atan = atan(p.y, p.x);
    float pi = 3.14159265359;
    
    float v_seed1 = v_atan + pi * round((v0 - v_atan) / pi);
    float v_seed2 = v0;
    float v_seed3 = v_seed1 + pi * (v0 > v_seed1 ? 1.0 : -1.0);
    
    vec3 v = vec3(v_seed1, v_seed2, v_seed3);
    
    for (int i = 0; i < 4; ++i) {
        vec3 sin_v = sin(v);
        vec3 cos_v = cos(v);
        vec3 w = p.x * cos_v + p.y * sin_v;
        vec3 q = -p.x * sin_v + p.y * cos_v;
        vec3 u_star = clamp(w, -R, R);
        
        vec3 g = -u_star * q + c * c * (v - v0);
        vec3 du_dv = vec3(
            abs(w.x) < R ? -q.x : 0.0,
            abs(w.y) < R ? -q.y : 0.0,
            abs(w.z) < R ? -q.z : 0.0
        );
        
        vec3 g_prime = -du_dv * q + u_star * w + c * c;
        v -= g / g_prime;
        v = clamp(v, v0 - pi, v0 + pi);
    }
    
    vec3 sin_v = sin(v);
    vec3 cos_v = cos(v);
    vec3 w = p.x * cos_v + p.y * sin_v;
    vec3 u_star = clamp(w, -R, R);
    
    vec3 dx = vec3(p.x) - u_star * cos_v;
    vec3 dy = vec3(p.y) - u_star * sin_v;
    vec3 dz = c * (vec3(v0) - v);
    
    vec3 d2_c = dx*dx + dy*dy + dz*dz;
    
    float d2 = d2_c.x;
    float v_best = v.x;
    float u_best = u_star.x;
    float dx_best = dx.x;
    float dy_best = dy.x;
    float dz_best = dz.x;
    
    if (d2_c.y < d2) {
        d2 = d2_c.y; v_best = v.y; u_best = u_star.y;
        dx_best = dx.y; dy_best = dy.y; dz_best = dz.y;
    }
    if (d2_c.z < d2) {
        d2 = d2_c.z; v_best = v.z; u_best = u_star.z;
        dx_best = dx.z; dy_best = dy.z; dz_best = dz.z;
    }
    
    float d = sqrt(d2);
    vec3 grad;
    if (d > 1e-7) {
        grad = vec3(dx_best, dy_best, dz_best) / d;
    } else {
        grad = vec3(c * sin(v_best), -c * cos(v_best), u_best);
        grad = normalize(grad);
    }
    
    return vec4(d - th, grad);
}

// --- 3. Hyperbolic-Secant Scherk Saddle Tower ---

float sdScherkTower(vec3 p) {
    vec2 c = clamp(p.xy, -20.0, 20.0);
    float sech_x = 1.0 / cosh(c.x);
    float sech_y = 1.0 / cosh(c.y);
    float g = sin(p.z) * sech_x * sech_y - tanh(c.x) * tanh(c.y);
    return g * 0.919;
}

vec4 sdgScherkTower(vec3 p) {
    vec2 c = clamp(p.xy, -20.0, 20.0);
    float sech_x = 1.0 / cosh(c.x);
    float sech_y = 1.0 / cosh(c.y);
    float tanh_x = tanh(c.x);
    float tanh_y = tanh(c.y);
    float sin_z = sin(p.z);
    
    float g = sin_z * sech_x * sech_y - tanh_x * tanh_y;
    float d = g * 0.919;
    
    vec3 grad = vec3(
        -sin_z * sech_y * sech_x * tanh_x - sech_x * sech_x * tanh_y,
        -sin_z * sech_x * sech_y * tanh_y - tanh_x * sech_y * sech_y,
        cos(p.z) * sech_x * sech_y
    );
    
    return vec4(d, grad * 0.919);
}

#endif // SDF_GEMS_HAS_MINIMAL_SURFACE
