#ifndef SDF_GEMS_SD_CONICAL_SPIRAL_GLSL
#define SDF_GEMS_SD_CONICAL_SPIRAL_GLSL

#include "../common/common.glsl"

// Exact O(1) 3D Euclidean SDF + Unit Gradient for Pappus's Archimedean Conical Helix
// q(theta) = (a*theta*cos(theta), a*theta*sin(theta), c*theta), theta in [0, h/c]
vec4 sdgConicalHelix(vec3 p, float a, float c, float h, float r) {
    float rho = length(p.xy);
    float phi = atan(p.y, p.x);
    float thMax = h / c;
    float L2 = a * a + c * c;
    float thCone = clamp((a * rho + c * p.z) / L2, 0.0, thMax);
    float k0 = round((thCone - phi) * 0.15915494309);

    vec3 qBest = vec3(0.0);
    float d2Best = dot(p, p);
    vec3 qEnd = vec3(a * thMax * cos(thMax), a * thMax * sin(thMax), h);
    vec3 dEnd = p - qEnd;
    float d2End = dot(dEnd, dEnd);
    if (d2End < d2Best) {
        d2Best = d2End;
        qBest = qEnd;
    }

    for (int dk = -2; dk <= 2; dk++) {
        float k = k0 + float(dk);
        for (int s = 0; s < 2; s++) {
            float th = clamp(phi + 6.28318530718 * k + float(s) * 3.14159265359, 0.0, thMax);
            for (int it = 0; it < 5; it++) {
                float ang = th - phi;
                float ca = cos(ang), sa = sin(ang);
                float g = L2 * th - c * p.z - a * rho * (ca - th * sa);
                float gp = L2 + a * rho * (2.0 * sa + th * ca);
                th = clamp(th - clamp(g / max(gp, 1e-12), -1.0, 1.0), 0.0, thMax);
            }
            vec3 q = vec3(a * th * cos(th), a * th * sin(th), c * th);
            vec3 dv = p - q;
            float d2 = dot(dv, dv);
            if (d2 < d2Best) {
                d2Best = d2;
                qBest = q;
            }
        }
    }
    float d = sqrt(max(d2Best, 1e-30));
    return vec4(d - r, (p - qBest) / d);
}

float sdConicalHelix(vec3 p, float a, float c, float h, float r) {
    return sdgConicalHelix(p, a, c, h, r).x;
}

vec4 sdgConicalSpiral(vec3 p, float a, float c, float h, float r) {
    return sdgConicalHelix(p, a, c, h, r);
}

float sdConicalSpiral(vec3 p, float a, float c, float h, float r) {
    return sdConicalHelix(p, a, c, h, r);
}

// Exact O(1) 3D Euclidean Canal SDF + Unit Gradient for Raup's (1966) Helicocone Gastropod Shell
// Spine C(theta) = a*exp(b*theta)*(cos(theta), sin(theta), h), Aperture Radius R(theta) = gamma*a*exp(b*theta)
vec4 sdgGastropodShell(vec3 p, float a, float b, float h, float gamma) {
    float rho = length(p.xy);
    float phi = atan(p.y, p.x);
    float h2p1 = 1.0 + h * h;

    // Cone-Generator Canal Footpoint Theorem:
    float gamFactor = gamma * inversesqrt(max(h2p1 - gamma * gamma, 1e-12));
    float rStar0 = max((rho + h * p.z + gamFactor * abs(h * rho - p.z)) / h2p1, 1e-6 * a);
    float rStarPi = max((-rho + h * p.z + gamFactor * abs(-h * rho - p.z)) / h2p1, 1e-6 * a);
    float invTwoPiB = 0.15915494309 / b;
    float k0 = round((log(rStar0 / a) - b * phi) * invTwoPiB);
    float kPi = round((log(rStarPi / a) - b * (phi + 3.14159265359)) * invTwoPiB);

    float dBest = length(p);
    vec3 gBest = (dBest > 1e-12) ? (p / dBest) : vec3(0.0, 0.0, 1.0);

    for (int idx = 0; idx < 8; idx++) {
        float k = (idx < 5) ? (k0 + float(idx - 2)) : (kPi + float(idx - 6));
        float thK = phi + 6.28318530718 * k;
        float rK = a * exp(clamp(b * thK, -40.0, 40.0));
        for (int s = 0; s < 2; s++) {
            float delta = float(s) * 3.14159265359;
            for (int it = 0; it < 5; it++) {
                float r = rK * exp(clamp(b * delta, -20.0, 20.0));
                float cd = cos(delta), sd = sin(delta);
                float D2 = max(h2p1 * r * r - 2.0 * (rho * cd + h * p.z) * r + rho * rho + p.z * p.z, 1e-30);
                float D = sqrt(D2);
                float K = b * h2p1 * r - rho * (b * cd - sd) - b * h * p.z;
                float G = K - gamma * b * D;
                float Gp = b * b * h2p1 * r + rho * (b * sd + cd) - gamma * b * (r / D) * K;
                delta -= clamp(G / max(Gp, 1e-12), -0.8, 0.8);
            }
            float r = rK * exp(clamp(b * delta, -20.0, 20.0));
            float cd = cos(delta), sd = sin(delta);
            vec2 cxy = (rho > 1e-12)
                ? (r / rho) * vec2(p.x * cd - p.y * sd, p.y * cd + p.x * sd)
                : r * vec2(cos(thK + delta), sin(thK + delta));
            vec3 dv = p - vec3(cxy, h * r);
            float D = sqrt(max(dot(dv, dv), 1e-30));
            float dVal = D - gamma * r;
            if (dVal < dBest) {
                dBest = dVal;
                gBest = dv / D;
            }
        }
    }
    return vec4(dBest, gBest);
}

float sdGastropodShell(vec3 p, float a, float b, float h, float gamma) {
    return sdgGastropodShell(p, a, b, h, gamma).x;
}

#endif // SDF_GEMS_SD_CONICAL_SPIRAL_GLSL