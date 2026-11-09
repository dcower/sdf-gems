// sdfWeylCovariantShell: Zero-bias Weyl-covariant heterogeneous logistic volume rendering & conformal Eikonal calculus.
// Derived by Gemini 4 Argon; implementation & verification by Gemini 3.8 Flash.
// Eliminates parasitic gauge dilation s_dot * d that causes surface distortion in heterogeneous shells.
// See docs/volume/sdfWeylCovariantShell.md.
#ifndef SDF_GEMS_SDF_WEYL_COVARIANT_SHELL_GLSL
#define SDF_GEMS_SDF_WEYL_COVARIANT_SHELL_GLSL

// Numerically robust logistic CDF Phi_s(x) = 1 / (1 + exp(-s x)) clamped against FP32 overflow.
float sdfWeylLogistic(float x, float s) {
    return 1.0 / (1.0 + exp(-clamp(s * x, -35.0, 35.0)));
}

// Pointwise Weyl-covariant extinction density along a ray:
//   sigma_Weyl = max(-s * dotD / (1 + exp(s * d)), 0.0)
// Depends strictly on the covariant derivative D_t d = dotD and is completely free
// of the parasitic gauge dilation s_dot * d that plagues naive heterogeneous volume shells.
float sdfWeylExtinction(float d, float dotD, float s) {
    float u = clamp(s * d, -35.0, 35.0);
    return max((-s * dotD) / (1.0 + exp(u)), 0.0);
}

// Exact zero-bias Weyl-covariant slab opacity between endpoints (d0, s0) and (d1, s1).
// Integrates the optical metric g_ij = s(x)^2 delta_ij via the geometric-mean gauge connection
//   sBar = sqrt(s0 * s1)
// to yield exact unbiased volume opacity without spurious boundary shifts.
float sdfOpacityWeylShell(float d0, float d1, float s0, float s1) {
    float sBar = sqrt(max(s0, 1e-4) * max(s1, 1e-4));
    float c0 = sdfWeylLogistic(d0, sBar);
    float c1 = sdfWeylLogistic(d1, sBar);
    return clamp((c0 - c1) / max(c0, 1e-6), 0.0, 1.0);
}

// Weyl-covariant cubic Hermite turning-point opacity across a ray segment of length dt.
// When grazing a curved silhouette (dotD0 < 0 and dotD1 > 0), finds the O(dt^4) interior turning
// point dMin and geometric log-sharpness sStar = s0^(1 - u*) * s1^(u*), compositing the two half-slabs.
float sdfOpacityWeylHermite(float d0, float d1, float dotD0, float dotD1, float s0, float s1, float dt) {
    if (dotD0 >= 0.0 || dotD1 <= 0.0) {
        return sdfOpacityWeylShell(d0, d1, s0, s1);
    }
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
        return sdfOpacityWeylShell(d0, d1, s0, s1);
    }
    float u2 = u * u, u3 = u2 * u;
    float h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
    float h10 = u3 - 2.0 * u2 + u;
    float h01 = -2.0 * u3 + 3.0 * u2;
    float h11 = u3 - u2;
    float dMin = min(min(d0, d1), h00 * d0 + h10 * m0 + h01 * d1 + h11 * m1);
    float sStar = exp(mix(log(max(s0, 1e-4)), log(max(s1, 1e-4)), u));
    float a1 = sdfOpacityWeylShell(d0, dMin, s0, sStar);
    float a2 = sdfOpacityWeylShell(dMin, d1, sStar, s1);
    return 1.0 - (1.0 - a1) * (1.0 - a2);
}

// Outer fuzzy shell boundary and Lipschitz-safe sphere tracing step bound.
// Given signed distance d, sharpness s, upper bound on ||grad ln s||, and transmittance threshold tau in (0, 0.5):
//   d_shell = d - ln(1/tau - 1) / s
//   t_safe  = d_shell / (1.0 + (ln(1/tau - 1) / s) * gradLogS_norm)
// Returns vec2(d_shell, t_safe).
vec2 sdfWeylShellBound(float d, float s, float gradLogS_norm, float tau) {
    float tauClamped = clamp(tau, 1e-5, 0.4999);
    float offset = log(1.0 / tauClamped - 1.0) / max(s, 1e-4);
    float dShell = d - offset;
    float denom = 1.0 + offset * max(gradLogS_norm, 0.0);
    float tSafe = dShell / denom;
    return vec2(dShell, tSafe);
}

#endif // SDF_GEMS_SDF_WEYL_COVARIANT_SHELL_GLSL
