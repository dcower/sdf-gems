// sdfPadeAutoIntSlab: Exact Stieltjes-Bernstein & Padé-Laplace closed-form slab radiative transfer.
// Derived by Gemini 4 Argon; implementation & verification by Gemini 3.8 Flash.
// Closed-form O(1) contraction of cubic volume radiance with cubic Beer-Lambert transmittance.
// See docs/volume/sdfPadeAutoIntSlab.md.
#ifndef SDF_GEMS_SDF_PADE_AUTO_INT_SLAB_GLSL
#define SDF_GEMS_SDF_PADE_AUTO_INT_SLAB_GLSL

// Exact 4x4 Rational Stieltjes-Bernstein Matrix M in column-major layout:
//   M_ij = \int_0^1 B_i^3(u) * (- (B_j^3)'(u)) du
// Column-sum identity: (1, 1, 1, 1) * M = (1, 0, 0, -1) identically,
// guaranteeing that sum(M * T_bern) == T0 - T3 == slab opacity.
const mat4 SDF_STIELTJES_BERNSTEIN_MAT4 = mat4(
     0.50,  0.30,  0.15,  0.05,  // Column 0
    -0.30,  0.00,  0.15,  0.15,  // Column 1
    -0.15, -0.15,  0.00,  0.30,  // Column 2
    -0.05, -0.15, -0.30, -0.50   // Column 3
);

// Numerically robust logistic CDF Phi_s(x) = 1 / (1 + exp(-s x)).
float sdfPadeLogistic(float x, float s) {
    return 1.0 / (1.0 + exp(-clamp(s * x, -35.0, 35.0)));
}

// Constructs the monotone C^1 cubic Bernstein control vector T_bern = (1.0, T1, T2, T3)
// of slab transmittance T(u) = clamp(Phi_s(d(u)) / max(Phi_s(d0), 1e-6), 0.0, 1.0) across [t0, t1].
vec4 sdfTransmittanceBernstein(float d0, float d1, float dotD0, float dotD1, float dt, float s) {
    float c0 = sdfPadeLogistic(d0, s);
    float c1 = sdfPadeLogistic(d1, s);
    float T3 = clamp(c1 / max(c0, 1e-6), 0.0, 1.0);

    float u0 = clamp(s * d0, -35.0, 35.0);
    float dT_dt0 = (s * dotD0) / (1.0 + exp(u0));
    float T1 = clamp(1.0 + (dt * (1.0 / 3.0)) * dT_dt0, T3, 1.0);

    float u1 = clamp(s * d1, -35.0, 35.0);
    float dT_dt1 = T3 * ((s * dotD1) / (1.0 + exp(u1)));
    float T2 = clamp(T3 - (dt * (1.0 / 3.0)) * dT_dt1, T3, T1);

    return vec4(1.0, T1, T2, T3);
}

// Contracts cubic Bernstein transmittance T_bern with cubic Bernstein colour control points
// (c0, c1, c2, c3) via the exact 4x4 Stieltjes-Bernstein matrix M.
// Returns vec4(integrated_rgb, slab_opacity).
vec4 sdfStieltjesSlabIntegral(vec4 T_bern, vec3 c0, vec3 c1, vec3 c2, vec3 c3) {
    vec4 w = SDF_STIELTJES_BERNSTEIN_MAT4 * T_bern;
    vec3 col = w.x * c0 + w.y * c1 + w.z * c2 + w.w * c3;
    float opacity = T_bern.x - T_bern.w;
    return vec4(col, opacity);
}

// Evaluates a monotonic (non-turning) slab with Hermite endpoints.
vec4 _sdfPadeMonotoneSlab(float d0, float d1, float dotD0, float dotD1,
                          vec3 c0, vec3 dc0, vec3 c1, vec3 dc1, float dt, float s) {
    float thirdDt = (1.0 / 3.0) * dt;
    vec3 b0 = c0;
    vec3 b1 = c0 + thirdDt * dc0;
    vec3 b2 = c1 - thirdDt * dc1;
    vec3 b3 = c1;
    vec4 T_bern = sdfTransmittanceBernstein(d0, d1, dotD0, dotD1, dt, s);
    return sdfStieltjesSlabIntegral(T_bern, b0, b1, b2, b3);
}

// Exact O(1) closed-form radiative transfer integral of spatially varying volume radiance
// c(t) against Beer-Lambert absorption across a ray slab of length dt.
// Converts Hermite colour endpoints (c0, dc0, c1, dc1) to cubic Bernstein basis,
// handles interior closest-approach turning points when dotD0 < 0 and dotD1 > 0,
// and returns vec4(integrated_rgb, slab_opacity) in O(1) arithmetic operations with 0 sub-steps.
vec4 sdfPadeLaplaceVolumeSlab(float d0, float d1, float dotD0, float dotD1,
                             vec3 c0, vec3 dc0, vec3 c1, vec3 dc1, float dt, float s) {
    if (dotD0 >= 0.0 || dotD1 <= 0.0) {
        return _sdfPadeMonotoneSlab(d0, d1, dotD0, dotD1, c0, dc0, c1, dc1, dt, s);
    }

    // Cubic Hermite interior minimum turning point
    float m0 = dotD0 * dt, m1 = dotD1 * dt;
    float a = 6.0 * (d0 - d1) + 3.0 * (m0 + m1);
    float b = -6.0 * (d0 - d1) - 2.0 * (2.0 * m0 + m1);
    float c = m0;
    float u = -c / max(m1 - m0, 1e-6);
    if (abs(a) > 1e-5) {
        float disc = max(b * b - 4.0 * a * c, 0.0);
        u = (-b + sqrt(disc)) / (2.0 * a);
    }

    if (u <= 1e-4 || u >= 0.9999) {
        return _sdfPadeMonotoneSlab(d0, d1, dotD0, dotD1, c0, dc0, c1, dc1, dt, s);
    }

    float u2 = u * u, u3 = u2 * u;
    float h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
    float h10 = u3 - 2.0 * u2 + u;
    float h01 = -2.0 * u3 + 3.0 * u2;
    float h11 = u3 - u2;
    float dMin = min(min(d0, d1), h00 * d0 + h10 * m0 + h01 * d1 + h11 * m1);

    // Hermite colour & colour derivative at turning point u
    vec3 cStar = h00 * c0 + h10 * (dt * dc0) + h01 * c1 + h11 * (dt * dc1);
    float dh00 = 6.0 * u2 - 6.0 * u;
    float dh10 = 3.0 * u2 - 4.0 * u + 1.0;
    float dh01 = -6.0 * u2 + 6.0 * u;
    float dh11 = 3.0 * u2 - 2.0 * u;
    vec3 dcStar = (dh00 * c0 + dh10 * (dt * dc0) + dh01 * c1 + dh11 * (dt * dc1)) / max(dt, 1e-6);

    float dt1 = u * dt;
    float dt2 = (1.0 - u) * dt;

    vec4 res1 = _sdfPadeMonotoneSlab(d0, dMin, dotD0, 0.0, c0, dc0, cStar, dcStar, dt1, s);
    vec4 res2 = _sdfPadeMonotoneSlab(dMin, d1, 0.0, dotD1, cStar, dcStar, c1, dc1, dt2, s);

    vec3 col = res1.rgb + (1.0 - res1.a) * res2.rgb;
    float opacity = 1.0 - (1.0 - res1.a) * (1.0 - res2.a);
    return vec4(col, opacity);
}

#endif // SDF_GEMS_SDF_PADE_AUTO_INT_SLAB_GLSL
