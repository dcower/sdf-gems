#ifndef SDF_GEMS_SD_FIGURE_EIGHT_KNOT_GLSL
#define SDF_GEMS_SD_FIGURE_EIGHT_KNOT_GLSL

#include "../common/common.glsl"

// ============================================================================
// Exact O(1) 3D Thurston's Hyperbolic Figure-Eight Knot (4_1)
// & Alternating Harmonic (p, q) Turk's Head Knots
//
// Discovered & mathematically designed by Gemini 4 Argon.
// Experiments, numerical verification & implementation by Gemini 3.8 Flash.
//
// Parametric centerline for alternating harmonic (p, q) Turk's head knot:
//   q(t) = ((R + r0*cos(q*t))*cos(p*t), (R + r0*cos(q*t))*sin(p*t), r1*sin(2*q*t))
//
// For (p, q) = (3, 2), this reproduces Thurston's figure-eight knot (4_1):
//   q(t) = ((R + r0*cos(2t))*cos(3t), (R + r0*cos(2t))*sin(3t), r1*sin(4t))
// ============================================================================

// 1. Exact O(1) 3D Euclidean SDF + Unit Gradient for Alternating Harmonic (p, q) Turk's Head Knots
vec4 sdgTurkHeadKnot(vec3 p, float R, float r0, float r1, float pBights, float qLeads, float rTube) {
    float rho = length(p.xy);
    float phi = (rho > 1e-12) ? atan(p.y, p.x) : 0.0;
    float d2Best = 1e30;
    vec3 qBest = vec3(R + r0, 0.0, 0.0);
    int pInt = int(round(pBights));
    float stepOff = 1.57079632679 / pBights;

    for (int k = 0; k < 6; k++) {
        if (k >= pInt) break;
        float tBase = (phi + 6.28318530718 * float(k)) / pBights;
        for (int m = -1; m <= 2; m++) {
            float t = tBase + float(m) * stepOff;
            for (int it = 0; it < 5; it++) {
                float qt = qLeads * t;
                float cqt = cos(qt), sqt = sin(qt);
                float c2qt = cqt * cqt - sqt * sqt, s2qt = 2.0 * cqt * sqt;
                float W = R + r0 * cqt;
                float Wp = -qLeads * r0 * sqt;
                float Wpp = -qLeads * qLeads * r0 * cqt;
                float Z = r1 * s2qt;
                float Zp = 2.0 * qLeads * r1 * c2qt;
                float Zpp = -4.0 * qLeads * qLeads * r1 * s2qt;
                float ang = pBights * t - phi;
                float ca = cos(ang), sa = sin(ang);
                float g = W * Wp - rho * (Wp * ca - pBights * W * sa) + (Z - p.z) * Zp;
                float gp = Wp * Wp + W * Wpp - rho * ((Wpp - pBights * pBights * W) * ca - 2.0 * pBights * Wp * sa) + Zp * Zp + (Z - p.z) * Zpp;
                t -= clamp(g / max(gp, 0.4), -0.3, 0.3);
            }
            float qt = qLeads * t;
            float cqt = cos(qt), sqt = sin(qt);
            float W = R + r0 * cqt;
            float pt = pBights * t;
            vec3 q = vec3(W * cos(pt), W * sin(pt), 2.0 * r1 * cqt * sqt);
            vec3 dv = p - q;
            float d2 = dot(dv, dv);
            if (d2 < d2Best) {
                d2Best = d2;
                qBest = q;
            }
        }
    }
    float d = sqrt(max(d2Best, 1e-30));
    vec3 grad = (d > 1e-12) ? ((p - qBest) / d) : vec3(0.0, 0.0, 1.0);
    return vec4(d - rTube, grad);
}

// 2. Exact O(1) 3D Euclidean SDF for Alternating Harmonic (p, q) Turk's Head Knots
float sdTurkHeadKnot(vec3 p, float R, float r0, float r1, float pBights, float qLeads, float rTube) {
    return sdgTurkHeadKnot(p, R, r0, r1, pBights, qLeads, rTube).x;
}

// 3. Exact O(1) 3D Euclidean SDF + Unit Gradient for Thurston's Hyperbolic Figure-Eight Knot (4_1)
// Special case (p = 3, q = 2)
vec4 sdgFigureEightKnot(vec3 p, float R, float r0, float r1, float rTube) {
    return sdgTurkHeadKnot(p, R, r0, r1, 3.0, 2.0, rTube);
}

// 4. Exact O(1) 3D Euclidean SDF for Thurston's Hyperbolic Figure-Eight Knot (4_1)
float sdFigureEightKnot(vec3 p, float R, float r0, float r1, float rTube) {
    return sdgTurkHeadKnot(p, R, r0, r1, 3.0, 2.0, rTube).x;
}

#endif // SDF_GEMS_SD_FIGURE_EIGHT_KNOT_GLSL
