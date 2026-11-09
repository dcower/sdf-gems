#ifndef SDF_GEMS_HAS_LOG_SPIRAL_2D
#define SDF_GEMS_HAS_LOG_SPIRAL_2D

// ============================================================================
// sdf-gems: sdLogSpiral2D & sdNautilusShell3D
// AI-Discovered by: Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
//
// Exact 2D Euclidean SDF and analytical unit gradient of Jacob Bernoulli's (1692)
// Equiangular / Logarithmic Spiral r(theta) = a * exp(b * theta).
//
// Exact 3D Euclidean SDF and analytical unit gradient of the Self-Similar Nautilus Shell
// (the canal envelope of spheres centered on the logarithmic spiral with self-similar radius
// R_tube(theta) = gamma * a * exp(b * theta)).
// ============================================================================

#define LOG_SPIRAL_PI 3.14159265358979323846

vec3 sdgLogSpiral2D(vec2 p, float a, float b, float th) {
    float rho = length(p);
    if (rho < 1e-15) return vec3(-th, 1.0, 0.0);
    
    float phi = atan(p.y, p.x);
    float k0 = round((log(rho / a) - b * phi) / (2.0 * LOG_SPIRAL_PI * b));
    
    float best_d2 = 1e30;
    vec2 best_diff = p;
    
    for (int i = -2; i <= 2; ++i) {
        float th_k = phi + 2.0 * LOG_SPIRAL_PI * (k0 + float(i));
        float rk = a * exp(clamp(b * th_k, -40.0, 40.0));
        
        // Closed-form initial seed at delta = 0
        float delta = clamp(b * (rho - rk) / max(b * b * rk + rho, 1e-15), -1.5, 1.5);
        
        for (int j = 0; j < 4; ++j) {
            float eb = rk * exp(b * delta);
            float cd = cos(delta);
            float sd = sin(delta);
            float g = b * eb - rho * (b * cd - sd);
            float gp = b * b * eb + rho * (b * sd + cd);
            delta = clamp(delta - g / max(abs(gp), 1e-15), -LOG_SPIRAL_PI, LOG_SPIRAL_PI);
        }
        
        float eb = rk * exp(b * delta);
        float cd = cos(delta);
        float sd = sin(delta);
        
        // Reconstruct S(th_k + delta) by rotating p / rho by delta
        float sc = eb / rho;
        vec2 s_pt = sc * vec2(cd * p.x - sd * p.y, cd * p.y + sd * p.x);
        vec2 diff = p - s_pt;
        float d2 = dot(diff, diff);
        
        if (d2 < best_d2) {
            best_d2 = d2;
            best_diff = diff;
        }
    }
    
    float d = sqrt(best_d2);
    vec2 grad = best_diff / max(d, 1e-15);
    return vec3(d - th, grad);
}

float sdLogSpiral2D(vec2 p, float a, float b, float th) {
    return sdgLogSpiral2D(p, a, b, th).x;
}

vec4 sdgNautilusShell3D(vec3 p, float a, float b, float gamma) {
    float px = p.x, py = p.y, pz = p.z;
    float rho = length(p.xy);
    float phi = rho > 1e-15 ? atan(py, px) : 0.0;
    float r3d = max(length(p), 1e-12);
    
    float k0 = round((log(r3d / a) - b * phi) / (2.0 * LOG_SPIRAL_PI * b));
    
    float best_f = 1e30;
    vec3 best_grad = vec3(0.0, 0.0, 1.0);
    
    for (int i = -2; i <= 2; ++i) {
        float th_k = phi + 2.0 * LOG_SPIRAL_PI * (k0 + float(i));
        float rk = a * exp(clamp(b * th_k, -40.0, 40.0));
        float delta = clamp(b * (rho - rk) / max(b * b * rk + rho, 1e-15), -1.5, 1.5);
        
        for (int j = 0; j < 8; ++j) {
            float r_th = rk * exp(b * delta);
            float cd = cos(delta);
            float sd = sin(delta);
            float diff_r = rho - r_th;
            float D2 = diff_r * diff_r + 2.0 * rho * r_th * (1.0 - cd) + pz * pz;
            float D = max(sqrt(D2), 1e-15);
            
            float num = b * r_th - rho * (b * cd - sd);
            float dD = r_th * num / D;
            float G = num - gamma * b * D;
            float dG = b * b * r_th + rho * (b * sd + cd) - gamma * b * dD;
            
            delta = clamp(delta - clamp(G / max(abs(dG), 1e-15), -0.6, 0.6), -LOG_SPIRAL_PI, LOG_SPIRAL_PI);
        }
        
        float r_th = rk * exp(b * delta);
        float cd = cos(delta);
        float sd = sin(delta);
        float sx, sy;
        
        if (rho > 1e-15) {
            float sc = r_th / rho;
            sx = sc * (cd * px - sd * py);
            sy = sc * (cd * py + sd * px);
        } else {
            float th = th_k + delta;
            sx = r_th * cos(th);
            sy = r_th * sin(th);
        }
        
        vec3 diff = vec3(px - sx, py - sy, pz);
        float D = max(length(diff), 1e-15);
        float f_val = D - gamma * r_th;
        
        if (f_val < best_f) {
            best_f = f_val;
            best_grad = diff / D;
        }
    }
    
    return vec4(best_f, best_grad);
}

float sdNautilusShell3D(vec3 p, float a, float b, float gamma) {
    return sdgNautilusShell3D(p, a, b, gamma).x;
}

#endif
