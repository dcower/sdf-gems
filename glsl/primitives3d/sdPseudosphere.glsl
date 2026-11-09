#ifndef SDF_GEMS_SD_PSEUDOSPHERE_GLSL
#define SDF_GEMS_SD_PSEUDOSPHERE_GLSL

#include "../common/common.glsl"

// sdPseudosphere / sdTractrix2D
//
// Exact 3D Euclidean SDF and analytical unit gradient of Eugenio Beltrami's (1868)
// Pseudosphere of Constant Negative Gaussian Curvature (K = -1 / R^2) and
// Christiaan Huygens's (1692) Tractrix Curve.
//
// Designed via the Pythagorean Hyperbolic Tangent-Normal Collapse Theorem:
//   G(t) = rho * sech(t) + (R * t - |z|) * tanh(t) - R = 0
//
// Discovered & derived by Gemini 4 Argon (research & mathematical design)
// and Gemini 3.8 Flash (experiments & implementation).

// Exact 2D Euclidean SDF + Unit Gradient for Huygens's (1692) Tractrix Curve
// q(t) = R * (sech(t), t - tanh(t)), t in [0, tMax]
// Evaluates exact distance to the wire/ribbon profile of half-thickness th.
vec3 sdgTractrix2D(vec2 p, float R, float tMax, float th) {
    vec2 sp = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    vec2 ap = abs(p);

    float tRho = clamp(acosh(max(R / max(ap.x, 1e-12), 1.0)), 0.0, tMax);
    float uz = ap.y / R;
    float tZ = clamp(uz + tanh(uz), 0.0, tMax);

    vec2 qBest = vec2(R, 0.0);
    float d2Best = (ap.x - R) * (ap.x - R) + ap.y * ap.y;
    vec2 qEnd = R * vec2(1.0 / cosh(tMax), tMax - tanh(tMax));
    float d2End = dot(ap - qEnd, ap - qEnd);
    if (d2End < d2Best) {
        d2Best = d2End;
        qBest = qEnd;
    }

    float seeds[3] = float[3](tRho, tZ, 0.5 * tMax);
    for (int s = 0; s < 3; s++) {
        float t = max(seeds[s], 1e-5);
        for (int it = 0; it < 6; it++) {
            float ch = cosh(t), sh = sinh(t);
            float sech = 1.0 / ch, tnh = sh / ch;
            float G = ap.x * sech + (R * t - ap.y) * tnh - R;
            float Gp = R * tnh + sech * ((R * t - ap.y) * sech - ap.x * tnh);
            t = clamp(t - clamp(G / max(Gp, 1e-12), -0.6, 0.6), 0.0, tMax);
        }
        vec2 q = R * vec2(1.0 / cosh(t), t - tanh(t));
        vec2 dv = ap - q;
        float d2 = dot(dv, dv);
        if (d2 < d2Best) {
            d2Best = d2;
            qBest = q;
        }
    }
    float d = sqrt(max(d2Best, 1e-30));
    vec2 g = (d > 1e-12) ? ((ap - qBest) / d) : vec2(0.0, 1.0);
    return vec3(d - th, sp * g);
}

float sdTractrix2D(vec2 p, float R, float tMax, float th) {
    return sdgTractrix2D(p, R, tMax, th).x;
}

// Exact 3D Euclidean SDF + Unit Gradient for Beltrami's (1868) Pseudosphere (K = -1 / R^2)
// Capped solid double horn of equatorial radius R and hyperbolic truncation parameter tMax.
vec4 sdgPseudosphere(vec3 p, float R, float tMax) {
    float rho = length(p.xy);
    float sz = (p.z < 0.0) ? -1.0 : 1.0;
    float az = abs(p.z);

    float chMax = cosh(tMax);
    float rhoMax = R / chMax;
    float zMax = R * (tMax - tanh(tMax));

    float tRho = clamp(acosh(max(R / max(rho, 1e-12), 1.0)), 0.0, tMax);
    float uz = az / R;
    float tZ = clamp(uz + tanh(uz), 0.0, tMax);

    vec2 qBest = vec2(R, 0.0);
    float d2Best = (rho - R) * (rho - R) + az * az;

    float rhoCap = min(rho, rhoMax);
    float d2Cap = (rho - rhoCap) * (rho - rhoCap) + (az - zMax) * (az - zMax);
    if (d2Cap < d2Best) {
        d2Best = d2Cap;
        qBest = vec2(rhoCap, zMax);
    }

    float seeds[3] = float[3](tRho, tZ, 0.5 * tMax);
    for (int s = 0; s < 3; s++) {
        float t = max(seeds[s], 1e-5);
        for (int it = 0; it < 6; it++) {
            float ch = cosh(t), sh = sinh(t);
            float sech = 1.0 / ch, tnh = sh / ch;
            float G = rho * sech + (R * t - az) * tnh - R;
            float Gp = R * tnh + sech * ((R * t - az) * sech - rho * tnh);
            t = clamp(t - clamp(G / max(Gp, 1e-12), -0.6, 0.6), 0.0, tMax);
        }
        vec2 q = R * vec2(1.0 / cosh(t), t - tanh(t));
        vec2 dv = vec2(rho, az) - q;
        float d2 = dot(dv, dv);
        if (d2 < d2Best) {
            d2Best = d2;
            qBest = q;
        }
    }

    bool inside = (az < zMax) && (rho < R) && (az < R * acosh(R / max(rho, 1e-12)) - sqrt(max(R * R - rho * rho, 0.0)));
    float sgn = inside ? -1.0 : 1.0;
    float d = sqrt(max(d2Best, 1e-30));
    vec2 g2 = (d > 1e-12) ? (sgn * (vec2(rho, az) - qBest) / d) : vec2(0.0, 1.0);
    vec2 er = (rho > 1e-12) ? (p.xy / rho) : vec2(1.0, 0.0);
    return vec4(sgn * d, g2.x * er.x, g2.x * er.y, sz * g2.y);
}

float sdPseudosphere(vec3 p, float R, float tMax) {
    return sdgPseudosphere(p, R, tMax).x;
}

#endif // SDF_GEMS_SD_PSEUDOSPHERE_GLSL
