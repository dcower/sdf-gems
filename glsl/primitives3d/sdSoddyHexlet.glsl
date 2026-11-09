#ifndef SDF_GEMS_HAS_SD_SODDY_HEXLET
#define SDF_GEMS_HAS_SD_SODDY_HEXLET

#include "../common/common.glsl"

/**
 * @brief $O(1)$ Exact 3D Euclidean SDF for an $N$-sphere Soddy necklace / Dupin cyclide.
 * @param p Point in 3D space
 * @param R Outer characteristic scale
 * @param beta Eccentricity (0 <= beta < 1)
 * @param nSpheres Number of kissing spheres (N >= 3)
 * @param phase Porism rotation phase
 * @return Exact Euclidean distance
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 */
float sdSoddyNecklace(vec3 p, float R, float beta, float nSpheres, float phase) {
    float s = sin(SDF_PI / nSpheres);
    float c_pi_N = cos(SDF_PI / nSpheres);
    float W = sqrt(max(0.0, 1.0 - (beta * c_pi_N) * (beta * c_pi_N)));
    float scale_y = (1.0 - beta * beta) / W;
    float sector = SDF_TAU / nSpheres;

    float rho = length(p.xy);
    float seed1 = atan(sqrt(max(0.0, 1.0 - beta * beta)) * p.y, p.x - beta * rho);
    float seed2 = SDF_PI;

    vec2 seeds = vec2(seed1, seed2);
    float min_d = 1e38;

    for (int s_idx = 0; s_idx < 2; ++s_idx) {
        float th = (s_idx == 0) ? seeds.x : seeds.y;
        for (int i = 0; i < 3; ++i) {
            float ct = cos(th), st = sin(th);
            float D = 1.0 + beta * ct;
            float invD = 1.0 / D;
            float invD2 = invD * invD;
            
            float Cx = R * (beta + ct) * invD;
            float Cy = R * scale_y * st * invD;
            float A  = R * scale_y * s * invD;
            
            float dCx = -R * (1.0 - beta * beta) * st * invD2;
            float dCy =  R * scale_y * (ct + beta) * invD2;
            float dA  =  R * scale_y * s * beta * st * invD2;
            
            float d2Cx = -R * (1.0 - beta * beta) * (ct * D + 2.0 * beta * st * st) * invD2 * invD;
            float d2Cy =  R * scale_y * (-st * D + 2.0 * beta * st * (ct + beta)) * invD2 * invD;
            float d2A  =  R * scale_y * s * beta * (ct * D + 2.0 * beta * st * st) * invD2 * invD;
            
            float dx = Cx - p.x;
            float dy = Cy - p.y;
            float dz = -p.z;
            float dist = max(sqrt(dx*dx + dy*dy + dz*dz), 1e-12);
            
            float fp = (dx * dCx + dy * dCy) / dist - dA;
            float dr_dth = (dx * dCx + dy * dCy) / dist;
            float fpp = (dCx*dCx + dCy*dCy + dx*d2Cx + dy*d2Cy - dr_dth*dr_dth) / dist - d2A;
            
            if (fpp > 1e-6) {
                th -= clamp(fp / fpp, -0.6, 0.6);
            } else {
                th -= sign(fp) * 0.4;
            }
        }
        
        float k0 = round((th - phase) / sector);
        for (int dk = -1; dk <= 1; ++dk) {
            float ang = (k0 + float(dk)) * sector + phase;
            float ct = cos(ang), st = sin(ang);
            float inv_denom = R / (1.0 + beta * ct);
            vec3 C = vec3(inv_denom * (beta + ct), inv_denom * scale_y * st, 0.0);
            float A = inv_denom * scale_y * s;
            min_d = min(min_d, length(p - C) - A);
        }
    }

    return min_d;
}

/**
 * @brief $O(1)$ Exact 3D Euclidean SDF and analytical gradient for an $N$-sphere Soddy necklace.
 * @param p Point in 3D space
 * @param R Outer characteristic scale
 * @param beta Eccentricity (0 <= beta < 1)
 * @param nSpheres Number of kissing spheres (N >= 3)
 * @param phase Porism rotation phase
 * @return vec4(d, grad.x, grad.y, grad.z)
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 */
vec4 sdgSoddyNecklace(vec3 p, float R, float beta, float nSpheres, float phase) {
    float s = sin(SDF_PI / nSpheres);
    float c_pi_N = cos(SDF_PI / nSpheres);
    float W = sqrt(max(0.0, 1.0 - (beta * c_pi_N) * (beta * c_pi_N)));
    float scale_y = (1.0 - beta * beta) / W;
    float sector = SDF_TAU / nSpheres;

    float rho = length(p.xy);
    float seed1 = atan(sqrt(max(0.0, 1.0 - beta * beta)) * p.y, p.x - beta * rho);
    float seed2 = SDF_PI;

    vec2 seeds = vec2(seed1, seed2);
    float min_d = 1e38;
    vec3 min_g = vec3(0.0);

    for (int s_idx = 0; s_idx < 2; ++s_idx) {
        float th = (s_idx == 0) ? seeds.x : seeds.y;
        for (int i = 0; i < 3; ++i) {
            float ct = cos(th), st = sin(th);
            float D = 1.0 + beta * ct;
            float invD = 1.0 / D;
            float invD2 = invD * invD;
            
            float Cx = R * (beta + ct) * invD;
            float Cy = R * scale_y * st * invD;
            float A  = R * scale_y * s * invD;
            
            float dCx = -R * (1.0 - beta * beta) * st * invD2;
            float dCy =  R * scale_y * (ct + beta) * invD2;
            float dA  =  R * scale_y * s * beta * st * invD2;
            
            float d2Cx = -R * (1.0 - beta * beta) * (ct * D + 2.0 * beta * st * st) * invD2 * invD;
            float d2Cy =  R * scale_y * (-st * D + 2.0 * beta * st * (ct + beta)) * invD2 * invD;
            float d2A  =  R * scale_y * s * beta * (ct * D + 2.0 * beta * st * st) * invD2 * invD;
            
            float dx = Cx - p.x;
            float dy = Cy - p.y;
            float dz = -p.z;
            float dist = max(sqrt(dx*dx + dy*dy + dz*dz), 1e-12);
            
            float fp = (dx * dCx + dy * dCy) / dist - dA;
            float dr_dth = (dx * dCx + dy * dCy) / dist;
            float fpp = (dCx*dCx + dCy*dCy + dx*d2Cx + dy*d2Cy - dr_dth*dr_dth) / dist - d2A;
            
            if (fpp > 1e-6) {
                th -= clamp(fp / fpp, -0.6, 0.6);
            } else {
                th -= sign(fp) * 0.4;
            }
        }
        
        float k0 = round((th - phase) / sector);
        for (int dk = -1; dk <= 1; ++dk) {
            float ang = (k0 + float(dk)) * sector + phase;
            float ct = cos(ang), st = sin(ang);
            float inv_denom = R / (1.0 + beta * ct);
            vec3 C = vec3(inv_denom * (beta + ct), inv_denom * scale_y * st, 0.0);
            float A = inv_denom * scale_y * s;
            vec3 diff = p - C;
            float dist = length(diff);
            float d = dist - A;
            if (d < min_d) {
                min_d = d;
                min_g = (dist > 1e-12) ? diff / dist : vec3(1.0, 0.0, 0.0);
            }
        }
    }

    return vec4(min_d, min_g);
}

/**
 * @brief Convenience wrapper for Frederick Soddy's canonical 6-sphere Hexlet.
 * @param p Point in 3D space
 * @param R Outer characteristic scale
 * @param beta Eccentricity (0 <= beta < 1)
 * @param phase Porism rotation phase
 * @return Exact Euclidean distance
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 */
float sdSoddyHexlet(vec3 p, float R, float beta, float phase) {
    return sdSoddyNecklace(p, R, beta, 6.0, phase);
}

/**
 * @brief Convenience wrapper for Frederick Soddy's canonical 6-sphere Hexlet with analytical gradient.
 * @param p Point in 3D space
 * @param R Outer characteristic scale
 * @param beta Eccentricity (0 <= beta < 1)
 * @param phase Porism rotation phase
 * @return vec4(d, grad.x, grad.y, grad.z)
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 */
vec4 sdgSoddyHexlet(vec3 p, float R, float beta, float phase) {
    return sdgSoddyNecklace(p, R, beta, 6.0, phase);
}

/**
 * @brief $O(1)$ Exact 2D Euclidean SDF for an $N$-Circle Steiner Porism Chain.
 * @param p Point in 2D space
 * @param R Outer characteristic scale
 * @param beta Eccentricity (0 <= beta < 1)
 * @param nCircles Number of circles
 * @param phase Porism rotation phase
 * @return Exact Euclidean distance
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 */
float sdSteinerChain2D(vec2 p, float R, float beta, float nCircles, float phase) {
    return sdSoddyNecklace(vec3(p, 0.0), R, beta, nCircles, phase);
}

/**
 * @brief $O(1)$ Exact 2D Euclidean SDF and analytical gradient for an $N$-Circle Steiner Porism Chain.
 * @param p Point in 2D space
 * @param R Outer characteristic scale
 * @param beta Eccentricity (0 <= beta < 1)
 * @param nCircles Number of circles
 * @param phase Porism rotation phase
 * @return vec3(d, grad.x, grad.y)
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 */
vec3 sdgSteinerChain2D(vec2 p, float R, float beta, float nCircles, float phase) {
    vec4 g = sdgSoddyNecklace(vec3(p, 0.0), R, beta, nCircles, phase);
    return vec3(g.x, g.yz);
}

/**
 * @brief Exact inner and outer tangent director circles for a Steiner Porism Chain.
 * @param p Point in 2D space
 * @param R Outer characteristic scale
 * @param beta Eccentricity (0 <= beta < 1)
 * @param nCircles Number of circles
 * @return vec2(dInnerCircle, dOuterCircle)
 * @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
 */
vec2 sdSteinerEnvelopes2D(vec2 p, float R, float beta, float nCircles) {
    float s = sin(SDF_PI / nCircles);
    float c_pi_N = cos(SDF_PI / nCircles);
    float W = sqrt(max(0.0, 1.0 - (beta * c_pi_N) * (beta * c_pi_N)));
    
    vec2 c_in = vec2(R * beta * s / W, 0.0);
    float r_in = R * (1.0 - s / W);
    
    vec2 c_out = vec2(-R * beta * s / W, 0.0);
    float r_out = R * (1.0 + s / W);
    
    return vec2(length(p - c_in) - r_in, length(p - c_out) - r_out);
}

#endif // SDF_GEMS_HAS_SD_SODDY_HEXLET
