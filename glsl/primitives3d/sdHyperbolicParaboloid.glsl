// ============================================================================
// sdHyperbolicParaboloid, sdHyperboloid1Sheet, sdHyperboloid2Sheet — Exact Quadrics
// @author Gemini 4 Argon (research & mathematical design), Gemini 3.1 Pro (experiments & implementation)
// @exact  Yes (pole-cleared secular cubic for saddle, pole-cleared bisection for hyperboloids)
// ============================================================================

#ifndef SDF_GEMS_HAS_SD_HYPERBOLIC_PARABOLOID
#define SDF_GEMS_HAS_SD_HYPERBOLIC_PARABOLOID

#include "../common/common.glsl"

// Helper: solve roots of depressed cubic t^3 + p t + q = 0 for the saddle
float sdf_internal_solveDepressedCubic(float p, float q) {
    float Delta = (q * q * 0.25) + (p * p * p / 27.0);
    if (Delta >= 0.0) {
        float sqrt_Delta = sqrt(Delta);
        float u = -q * 0.5 + sqrt_Delta;
        float v = -q * 0.5 - sqrt_Delta;
        u = sign(u) * pow(abs(u), 1.0 / 3.0);
        v = sign(v) * pow(abs(v), 1.0 / 3.0);
        return u + v;
    } else {
        float r = sqrt(-p / 3.0);
        float theta = acos(-q / (2.0 * r * r * r));
        return 2.0 * r * cos(theta / 3.0);
    }
}

// 1. Exact 3D Hyperbolic Paraboloid (Pringle Saddle)
// $z = \frac{x^2 - y^2}{2a}$
float sdHyperbolicParaboloid(vec3 p, float a) {
    float px = abs(p.x), py = abs(p.y), pz = p.z;
    float f0 = 2.0 * a * pz - px * px + py * py;
    float sign_d = (f0 >= 0.0) ? 1.0 : -1.0;
    
    float x0 = px, y0 = py, z0 = pz;
    if (f0 < 0.0) {
        x0 = py; y0 = px; z0 = -pz;
    }
    
    float xf, yf, zf;
    if (x0 < 1e-7) { // On or extremely close to the focal plane
        float cond = 2.0 * a * (z0 - a) + 0.25 * y0 * y0;
        if (cond > 0.0) {
            float s = sqrt(cond);
            xf = s; yf = 0.5 * y0;
            zf = (xf * xf - yf * yf) / (2.0 * a);
        } else {
            xf = 0.0;
            float p_c = 2.0 * a * (a + z0);
            float q_c = -2.0 * a * a * y0;
            yf = sdf_internal_solveDepressedCubic(p_c, q_c);
            zf = -yf * yf / (2.0 * a);
        }
    } else {
        float cond = max(2.0 * a * (z0 - a) + 0.25 * y0 * y0, 0.0);
        float s = max(x0, sqrt(cond) + pow(2.0 * a * a * x0, 1.0 / 3.0));
        
        for (int i = 0; i < 4; i++) {
            float denom = max(2.0 * s - x0, 1e-8);
            float y_tmp = (s * y0) / denom;
            float H = s * s * s - (2.0 * a * (z0 - a) + y_tmp * y_tmp) * s - 2.0 * a * a * x0;
            float H_prime = 3.0 * s * s - 2.0 * a * (z0 - a) - y_tmp * y_tmp * (2.0 * s - 3.0 * x0) / denom;
            s = max(x0, s - H / H_prime);
        }
        xf = s;
        yf = (s * y0) / max(2.0 * s - x0, 1e-8);
        zf = (xf * xf - yf * yf) / (2.0 * a);
    }
    
    float d = length(vec3(x0 - xf, y0 - yf, z0 - zf));
    return sign_d * d;
}

vec4 sdgHyperbolicParaboloid(vec3 p, float a) {
    float px = abs(p.x), py = abs(p.y), pz = p.z;
    float f0 = 2.0 * a * pz - px * px + py * py;
    float sign_d = (f0 >= 0.0) ? 1.0 : -1.0;
    
    float x0 = px, y0 = py, z0 = pz;
    if (f0 < 0.0) {
        x0 = py; y0 = px; z0 = -pz;
    }
    
    float xf, yf, zf;
    if (x0 < 1e-7) {
        float cond = 2.0 * a * (z0 - a) + 0.25 * y0 * y0;
        if (cond > 0.0) {
            float s = sqrt(cond);
            xf = s; yf = 0.5 * y0;
            zf = (xf * xf - yf * yf) / (2.0 * a);
        } else {
            xf = 0.0;
            float p_c = 2.0 * a * (a + z0);
            float q_c = -2.0 * a * a * y0;
            yf = sdf_internal_solveDepressedCubic(p_c, q_c);
            zf = -yf * yf / (2.0 * a);
        }
    } else {
        float cond = max(2.0 * a * (z0 - a) + 0.25 * y0 * y0, 0.0);
        float s = max(x0, sqrt(cond) + pow(2.0 * a * a * x0, 1.0 / 3.0));
        
        for (int i = 0; i < 4; i++) {
            float denom = max(2.0 * s - x0, 1e-8);
            float y_tmp = (s * y0) / denom;
            float H = s * s * s - (2.0 * a * (z0 - a) + y_tmp * y_tmp) * s - 2.0 * a * a * x0;
            float H_prime = 3.0 * s * s - 2.0 * a * (z0 - a) - y_tmp * y_tmp * (2.0 * s - 3.0 * x0) / denom;
            s = max(x0, s - H / H_prime);
        }
        xf = s;
        yf = (s * y0) / max(2.0 * s - x0, 1e-8);
        zf = (xf * xf - yf * yf) / (2.0 * a);
    }
    
    vec3 dp = vec3(x0 - xf, y0 - yf, z0 - zf);
    float d = length(dp);
    
    vec3 grad;
    if (d > 1e-6) {
        grad = sign_d * (dp / d);
        if (f0 < 0.0) grad = vec3(grad.y, grad.x, -grad.z);
    } else {
        grad = normalize(vec3(-px / a, py / a, 1.0));
    }
    grad.x *= (p.x < 0.0 ? -1.0 : 1.0);
    grad.y *= (p.y < 0.0 ? -1.0 : 1.0);
    
    return vec4(sign_d * d, grad);
}

// 2. Bounded Circular/Elliptical Pringle Saddle Shell
float sdPringleChip(vec3 p, float a, float r, float th) {
    float d = sdHyperbolicParaboloid(p, a);
    return max(abs(d) - th, length(p.xy) - r);
}

// Helper: Exact 2D distance to generic hyperbola (x^2/a^2 - y^2/c^2 = 1)
float sdf_internal_distHyperbola2D(float x0, float y0, float a, float c) {
    x0 = abs(x0); y0 = abs(y0);
    
    if (x0 < 1e-6) {
        float S = c * y0 / (a * a + c * c);
        float xf = a * sqrt(1.0 + S * S);
        float yf = c * S;
        return length(vec2(x0 - xf, y0 - yf));
    }
    
    float a2 = a * a;
    float c2 = c * c;
    float a2_c2 = a2 + c2;
    float a2_x0 = a2 * x0;
    
    float s_min = max(a, a2 * x0 / a2_c2);
    float s_max = x0 + sqrt((x0 - a) * (x0 - a) + y0 * y0);
    
    float s = 0.5 * (s_min + s_max);
    for (int i = 0; i < 24; i++) {
        float denom = max(s * a2_c2 - a2_x0, 1e-8);
        float yf = s * c2 * y0 / denom;
        float val = s * s - a2 - (a2 / c2) * yf * yf;
        if (val > 0.0) {
            s_max = s;
        } else {
            s_min = s;
        }
        s = 0.5 * (s_min + s_max);
    }
    
    float xf = s;
    float denom = max(s * a2_c2 - a2_x0, 1e-8);
    float yf = s * c2 * y0 / denom;
    return length(vec2(x0 - xf, y0 - yf));
}

vec3 sdf_internal_gradHyperbola2D(float x0, float y0, float a, float c, float sign_d) {
    x0 = abs(x0); y0 = abs(y0);
    float xf, yf;
    
    if (x0 < 1e-6) {
        float S = c * y0 / (a * a + c * c);
        xf = a * sqrt(1.0 + S * S);
        yf = c * S;
    } else {
        float a2 = a * a;
        float c2 = c * c;
        float a2_c2 = a2 + c2;
        float a2_x0 = a2 * x0;
        
        float s_min = max(a, a2 * x0 / a2_c2);
        float s_max = x0 + sqrt((x0 - a) * (x0 - a) + y0 * y0);
        
        float s = 0.5 * (s_min + s_max);
        for (int i = 0; i < 24; i++) {
            float denom = max(s * a2_c2 - a2_x0, 1e-8);
            float yf_tmp = s * c2 * y0 / denom;
            float val = s * s - a2 - (a2 / c2) * yf_tmp * yf_tmp;
            if (val > 0.0) s_max = s;
            else s_min = s;
            s = 0.5 * (s_min + s_max);
        }
        xf = s;
        float denom = max(s * a2_c2 - a2_x0, 1e-8);
        yf = s * c2 * y0 / denom;
    }
    
    vec2 dp = vec2(x0 - xf, y0 - yf);
    float d = length(dp);
    vec2 grad2d;
    if (d > 1e-6) {
        grad2d = sign_d * (dp / d);
    } else {
        grad2d = normalize(vec2(x0 / (a * a), -y0 / (c * c)));
    }
    return vec3(sign_d * d, grad2d);
}

// 3. Exact 3D Hyperboloid of One Sheet
// $\frac{x^2 + y^2}{a^2} - \frac{z^2}{c^2} = 1$
float sdHyperboloid1Sheet(vec3 p, float a, float c) {
    float rho = length(p.xy);
    float z = abs(p.z);
    
    float f0 = (rho * rho) / (a * a) - (z * z) / (c * c) - 1.0;
    float sign_d = (f0 <= 0.0) ? -1.0 : 1.0;
    
    float d = sdf_internal_distHyperbola2D(rho, z, a, c);
    return sign_d * d;
}

vec4 sdgHyperboloid1Sheet(vec3 p, float a, float c) {
    float rho = length(p.xy);
    float z = abs(p.z);
    
    float f0 = (rho * rho) / (a * a) - (z * z) / (c * c) - 1.0;
    float sign_d = (f0 <= 0.0) ? -1.0 : 1.0;
    
    vec3 res = sdf_internal_gradHyperbola2D(rho, z, a, c, sign_d);
    float d = res.x;
    vec2 grad2d = res.yz;
    
    vec3 grad;
    if (rho > 1e-6) {
        grad = vec3((p.x / rho) * grad2d.x, (p.y / rho) * grad2d.x, (p.z < 0.0 ? -1.0 : 1.0) * grad2d.y);
    } else {
        grad = vec3(0.0, 0.0, (p.z < 0.0 ? -1.0 : 1.0) * grad2d.y);
    }
    
    return vec4(d, grad);
}

// 4. Exact 3D Hyperboloid of Two Sheets
// $\frac{z^2}{c^2} - \frac{x^2 + y^2}{a^2} = 1$
float sdHyperboloid2Sheet(vec3 p, float a, float c) {
    float rho = length(p.xy);
    float z = abs(p.z);
    
    float f0 = (z * z) / (c * c) - (rho * rho) / (a * a) - 1.0;
    float sign_d = (f0 <= 0.0) ? -1.0 : 1.0;
    
    float d = sdf_internal_distHyperbola2D(z, rho, c, a);
    return sign_d * d;
}

#endif // SDF_GEMS_HAS_SD_HYPERBOLIC_PARABOLOID
