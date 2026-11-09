#ifndef SDF_GEMS_SD_TOROIDAL_COIL_GLSL
#define SDF_GEMS_SD_TOROIDAL_COIL_GLSL

#include "../common/common.glsl"

vec4 sdgToroidalCoil(vec3 p, float R, float r0, float nTurns, float rTube) {
    const float TWO_PI = 6.28318530718;
    const float PI     = 3.14159265359;
    float N    = max(round(nTurns), 1.0);
    float rho  = length(p.xy);
    float phi0 = (rho > 1e-12) ? atan(p.y, p.x) : 0.0;
    float psi0 = atan(p.z, rho - R);

    float k0 = round((N * phi0 - psi0) / TWO_PI);
    float psi1 = psi0 + PI;
    float k1 = round((N * phi0 - psi1) / TWO_PI);
    float k2 = round((N * (phi0 + PI) - PI) / TWO_PI);

    float seeds[8];
    seeds[0] = (psi0 + TWO_PI * (k0 - 2.0)) / N;
    seeds[1] = (psi0 + TWO_PI * (k0 - 1.0)) / N;
    seeds[2] = (psi0 + TWO_PI * k0) / N;
    seeds[3] = (psi0 + TWO_PI * (k0 + 1.0)) / N;
    seeds[4] = (psi0 + TWO_PI * (k0 + 2.0)) / N;
    seeds[5] = (psi1 + TWO_PI * k1) / N;
    seeds[6] = phi0;
    seeds[7] = (PI + TWO_PI * k2) / N;

    float maxStep = PI / N;
    float bestD2 = 1e30;
    vec3  bestDiff = vec3(p.x - (R + r0), p.y, p.z);

    for (int s = 0; s < 8; ++s) {
        float phi = seeds[s];
        for (int i = 0; i < 6; ++i) {
            float c1 = cos(phi), s1 = sin(phi);
            float cn = cos(N * phi), sn = sin(N * phi);
            float u  = p.x * c1 + p.y * s1;
            float v  = -p.x * s1 + p.y * c1;
            float r  = R + r0 * cn;
            float f  = -r * v + N * r0 * (sn * (u - R) - cn * p.z);
            float df = r * u + 2.0 * N * r0 * sn * v + (N * N * r0) * (cn * (u - R) + sn * p.z);
            float denom = (df > 1e-7) ? df : (abs(df) + 1.0);
            phi -= clamp(f / denom, -maxStep, maxStep);
        }
        float c1 = cos(phi), s1 = sin(phi);
        float cn = cos(N * phi), sn = sin(N * phi);
        float r  = R + r0 * cn;
        vec3  diff = p - vec3(r * c1, r * s1, r0 * sn);
        float d2 = dot(diff, diff);
        if (d2 < bestD2) {
            bestD2 = d2;
            bestDiff = diff;
        }
    }

    float invD = inversesqrt(max(bestD2, 1e-30));
    float d    = bestD2 * invD;
    return vec4(d - rTube, bestDiff * invD);
}

float sdToroidalCoil(vec3 p, float R, float r0, float nTurns, float rTube) {
    return sdgToroidalCoil(p, R, r0, nTurns, rTube).x;
}

#endif // SDF_GEMS_SD_TOROIDAL_COIL_GLSL
