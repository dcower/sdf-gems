// sdfBernsteinMomentRay: 1D Bernstein-Bezier ray basis, mat4 Gram quadrature & Christoffel-Darboux Moment-SOS kernels.
// See docs/volume/sdfBernsteinMomentRay.md.
#ifndef SDF_GEMS_SDF_BERNSTEIN_MOMENT_RAY_GLSL
#define SDF_GEMS_SDF_BERNSTEIN_MOMENT_RAY_GLSL

// Exact 4x4 Symmetric Bernstein Gram Matrix G_{3,3}
// (G_{3,3})_{ij} = (1/7) * C(3,i)*C(3,j) / C(6,i+j) = (1/140) * [20,10,4,1; 10,12,9,4; 4,9,12,10; 1,4,10,20]
const mat4 SDF_BERNSTEIN_GRAM_33 = mat4(
    0.142857143, 0.071428571, 0.028571429, 0.007142857,
    0.071428571, 0.085714286, 0.064285714, 0.028571429,
    0.028571429, 0.064285714, 0.085714286, 0.071428571,
    0.007142857, 0.028571429, 0.071428571, 0.142857143
);

// Convert endpoint values & ray derivatives (fa, dfa_dt, fb, dfb_dt) on [ta, tb]
// into cubic Bernstein control points c = (c0, c1, c2, c3).
vec4 sdfRayHermiteToBernstein(float fa, float dfa, float fb, float dfb, float dt) {
    float thirdDt = (1.0 / 3.0) * dt;
    return vec4(fa, fa + thirdDt * dfa, fb - thirdDt * dfb, fb);
}

// Non-negative clamped cubic Bernstein control points (guarantees p(u) >= 0 on [ta, tb]).
vec4 sdfRayHermiteToBernsteinPos(float fa, float dfa, float fb, float dfb, float dt) {
    return max(sdfRayHermiteToBernstein(fa, dfa, fb, dfb, dt), vec4(0.0));
}

// Exact integral of a single cubic Bernstein field over [ta, tb].
float sdfBernsteinIntegral1(vec4 c, float dt) {
    return 0.25 * dt * dot(c, vec4(1.0));
}

// Exact integral of the degree-6 product of two cubic Bernstein fields over [ta, tb].
float sdfBernsteinGramIntegral(vec4 rho, vec4 w, float dt) {
    return dt * dot(rho, SDF_BERNSTEIN_GRAM_33 * w);
}

// Evaluate Christoffel-Darboux SOS sublevel SDF d_CD(p) = sqrt(P_2(p)) - r0
// where P_2(p) = ||L^T * vec4(p, 1.0)||^2.
float sdChristoffelSOS(vec3 p, mat4 L, float r0) {
    vec4 v = transpose(L) * vec4(p, 1.0);
    return length(v) - r0;
}

// Exact 1-atan slab ray integral of Christoffel-Darboux Moment-SOS kernels
// P_2(x)^{-1} (returned in .x) and P_2(x)^{-2} (returned in .y) over t in [ta, tb],
// where P_2(x) = ||L^T * vec4(x, 1.0)||^2.
vec2 sdfChristoffelRayIntegral(vec3 ro, vec3 rd, float ta, float tb, mat4 L) {
    float dt = max(tb - ta, 0.0);
    if (dt <= 0.0) return vec2(0.0);

    vec4 u0 = transpose(L) * vec4(ro, 1.0);
    vec4 u1 = transpose(L) * vec4(rd, 0.0);

    float c2 = max(dot(u1, u1), 1e-12);
    float c1 = dot(u0, u1);
    float tStar = -c1 / c2;

    // Cancellation-free perpendicular residual
    vec4 uPerp = u0 + tStar * u1;
    float h2 = max(dot(uPerp, uPerp), 1e-12);
    float h = sqrt(h2);
    float s = sqrt(c2);

    float da = ta - tStar;
    float db = tb - tStar;
    float xiA = (s * da) / h;
    float xiB = (s * db) / h;

    // Single atan2 call in [0, pi]
    float dTheta = atan(s * h * dt, h2 + c2 * da * db);

    float invSH = 1.0 / (s * h);
    float I1 = invSH * dTheta;

    float rationalTerm = (xiB / (1.0 + xiB * xiB)) - (xiA / (1.0 + xiA * xiA));
    float I2 = (0.5 * invSH / h2) * (dTheta + rationalTerm);

    return vec2(I1, I2);
}

#endif // SDF_GEMS_SDF_BERNSTEIN_MOMENT_RAY_GLSL
