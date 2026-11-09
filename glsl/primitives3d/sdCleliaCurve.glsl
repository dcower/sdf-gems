#ifndef SDF_GEMS_SD_CLELIA_CURVE_GLSL
#define SDF_GEMS_SD_CLELIA_CURVE_GLSL

#include "../common/common.glsl"

// Exact 3D Euclidean SDF + Unit Gradient for the Spherical Loxodrome (Mercator Rhumb Line)
// q(theta) = R * (sech(a*theta)*cos(theta), sech(a*theta)*sin(theta), tanh(a*theta))
vec4 sdgSphericalLoxodrome(vec3 p, float R, float a, float r) {
    float rho = length(p.xy);
    float phi = atan(p.y, p.x);

    // Asymptotic polar accumulation points theta -> +-inf: (0, 0, +-R)
    float sz = (p.z >= 0.0) ? 1.0 : -1.0;
    vec3 qBest = vec3(0.0, 0.0, sz * R);
    vec3 dVec = p - qBest;
    float d2Best = dot(dVec, dVec);

    float thMerc = asinh(p.z / max(rho, 1e-12)) / a;
    float kMerc = round((thMerc - phi) * 0.15915494309);
    float kEq = round(-phi * 0.15915494309);

    for (int i = 0; i < 6; i++) {
        float k = (i < 3) ? (kMerc + float(i - 1)) : (kEq + float(i - 4));
        float thK = phi + 6.28318530718 * k;
        for (int m = 0; m < 2; m++) {
            float delta;
            if (m == 0) {
                float ath0 = clamp(a * thK, -40.0, 40.0);
                float sech0 = 1.0 / cosh(ath0);
                float tanh0 = tanh(ath0);
                float g0 = a * (rho * tanh0 - p.z * sech0);
                float gp0 = rho * (1.0 + a * a * sech0 * sech0) + a * a * p.z * sech0 * tanh0;
                delta = (gp0 > 0.0) ? clamp(-g0 / max(gp0, 1e-12), -1.5707963, 1.5707963) : 0.0;
            } else {
                delta = 3.14159265359;
            }
            for (int it = 0; it < 5; it++) {
                float th = thK + delta;
                float ath = clamp(a * th, -40.0, 40.0);
                float sech = 1.0 / cosh(ath);
                float tnh = tanh(ath);
                float cd = cos(delta), sd = sin(delta);
                float g = rho * (sd + a * tnh * cd) - a * p.z * sech;
                float gp = rho * ((1.0 + a * a * sech * sech) * cd - a * tnh * sd) + a * a * p.z * sech * tnh;
                delta -= clamp(g / ((gp >= 0.0 ? 1.0 : -1.0) * max(abs(gp), 1e-12)), -0.8, 0.8);
            }
            float th = thK + delta;
            float ath = clamp(a * th, -40.0, 40.0);
            float sech = 1.0 / cosh(ath);
            float tnh = tanh(ath);
            float cd = cos(delta), sd = sin(delta);
            vec2 qxy = (rho > 1e-12)
                ? (R * sech / rho) * vec2(p.x * cd - p.y * sd, p.y * cd + p.x * sd)
                : (R * sech) * vec2(cos(th), sin(th));
            vec3 q = vec3(qxy, R * tnh);
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

float sdSphericalLoxodrome(vec3 p, float R, float a, float r) {
    return sdgSphericalLoxodrome(p, R, a, r).x;
}

// Exact 3D Euclidean SDF + Unit Gradient for Grandi's (1728) Spherical Clelia Curve
// q(t) = R * (cos(t)*cos(m*t), cos(t)*sin(m*t), sin(t))
vec4 sdgCleliaCurve(vec3 p, float R, float m, float r) {
    float rho = length(p.xy);
    float phi = atan(p.y, p.x);
    float d2Best = 1e30;
    vec3 qBest = vec3(R, 0.0, 0.0);
    int nK = int(ceil(0.5 * m)) + 1;

    for (int halfPeriod = 0; halfPeriod < 2; halfPeriod++) {
        float phiEff = (halfPeriod == 0) ? phi : ((m + 1.0) * 3.14159265359 - phi);
        float kCenter = round(-phiEff * 0.15915494309);
        for (int dk = -4; dk <= 4; dk++) {
            if (abs(dk) > nK) continue;
            float k = kCenter + float(dk);
            for (int sub = 0; sub < 2; sub++) {
                float u = clamp((phiEff + 6.28318530718 * k + float(sub) * 3.14159265359) / m, -1.57079632679, 1.57079632679);
                for (int it = 0; it < 5; it++) {
                    float cu = cos(u), su = sin(u);
                    float ang = m * u - phiEff;
                    float ca = cos(ang), sa = sin(ang);
                    float hp = -rho * (su * ca + m * cu * sa) + p.z * cu;
                    float hpp = -rho * ((1.0 + m * m) * cu * ca - 2.0 * m * su * sa) - p.z * su;
                    u = clamp(u - clamp(hp / min(hpp, -1e-12), -0.6, 0.6), -1.57079632679, 1.57079632679);
                }
                float t = (halfPeriod == 0) ? u : (3.14159265359 - u);
                float ct = cos(t), st = sin(t);
                vec3 q = R * vec3(ct * cos(m * t), ct * sin(m * t), st);
                vec3 dv = p - q;
                float d2 = dot(dv, dv);
                if (d2 < d2Best) {
                    d2Best = d2;
                    qBest = q;
                }
            }
        }
    }
    float d = sqrt(max(d2Best, 1e-30));
    return vec4(d - r, (p - qBest) / d);
}

float sdCleliaCurve(vec3 p, float R, float m, float r) {
    return sdgCleliaCurve(p, R, m, r).x;
}

#endif // SDF_GEMS_SD_CLELIA_CURVE_GLSL
