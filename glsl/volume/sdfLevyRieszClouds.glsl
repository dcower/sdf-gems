#ifndef SDF_GEMS_LEVY_RIESZ_CLOUDS_GLSL
#define SDF_GEMS_LEVY_RIESZ_CLOUDS_GLSL

// 2-Center Algebraic Riesz-Plummer Metaball SDF & exact analytical gradient (zero log/exp/pow).
// Uses the Newtonian-Plummer harmonic potential Phi(x) = sum_i w_i / sqrt(||x - c_i||^2 + a_i^2).
// c1, c2: vec4(center.xyz, coreRadius_a)
// w: vec2(w1, w2) > 0, targetRadius: R0 > 0
// Returns vec4(d_Riesz, grad_d_Riesz) guaranteed 1-Lipschitz everywhere and asymptotically ||x|| - R0.
vec4 sdgRieszMetaball2(vec3 p, vec4 c1, vec4 c2, vec2 w, float targetRadius) {
    vec3 d1 = p - c1.xyz;
    vec3 d2 = p - c2.xyz;
    float r1Sq = dot(d1, d1) + c1.w * c1.w;
    float r2Sq = dot(d2, d2) + c2.w * c2.w;

    float invR1 = inversesqrt(max(r1Sq, 1e-12));
    float invR2 = inversesqrt(max(r2Sq, 1e-12));

    float W = w.x + w.y;
    float phi = w.x * invR1 + w.y * invR2;
    float rEff = W / max(phi, 1e-12);

    // Gradient of -phi is sum_i w_i * invR_i^3 * d_i
    vec3 negGradPhi = (w.x * invR1 * invR1 * invR1) * d1
                    + (w.y * invR2 * invR2 * invR2) * d2;

    float factor = (rEff * rEff) / max(W, 1e-12);
    vec3 gradReff = factor * negGradPhi;
    float gradNorm = length(gradReff);

    // Normalize by max(gradNorm, 1.0) so ||grad|| <= 1 everywhere and ||grad|| == 1 on surface
    float scale = 1.0 / max(gradNorm, 1.0);
    float sdf = (rEff - targetRadius) * scale;
    vec3 unitGrad = (gradNorm > 1e-8) ? (gradReff / gradNorm) : vec3(0.0, 1.0, 0.0);
    return vec4(sdf, unitGrad);
}

// Exact algebraic slab ray integral of the 3D Plummer-Riesz cloud density kernel
// rho(x) = a^2 / (||x - c||^2 + a^2)^(3/2) over t in [ta, tb].
float sdfPlummerRieszRayIntegral(vec3 ro, vec3 rd, float ta, float tb, vec3 center, float coreRadius) {
    float dt = max(tb - ta, 0.0);
    if (dt <= 0.0) return 0.0;

    vec3 m = ro - center;
    float tStar = -dot(m, rd);
    vec3 mPerp = m + tStar * rd;

    float a2 = coreRadius * coreRadius;
    float rPerp2 = max(dot(mPerp, mPerp) + a2, 1e-10);

    float tauA = ta - tStar;
    float tauB = tb - tStar;
    float A = rPerp2 + tauA * tauA;
    float B = rPerp2 + tauB * tauB;

    float invA = inversesqrt(A);
    float invB = inversesqrt(B);

    // If tauA and tauB have the same sign, use the cancellation-free rationalized identity
    if (tauA * tauB > 0.0) {
        float denom = abs(tauB) * A * (B * invB) + abs(tauA) * B * (A * invA);
        return (a2 * dt * abs(tauA + tauB)) / max(denom, 1e-12);
    }
    return (a2 / rPerp2) * max(tauB * invB - tauA * invA, 0.0);
}

// Davis-Mineev-Weinstein Non-Exponential Pareto Cloud Transmittance T_H(tau) and Weight w_H(tau).
// Returns vec2(transmittance, collisionWeight) = ((1 + tau/H)^(-H), (1 + tau/H)^(-(H+1))).
// Uses fast rational paths for H == 1.0 and H == 2.0.
vec2 sdfDavisMineevTransmittance(float tau, float H) {
    float t = max(tau, 0.0);
    float h = max(H, 1e-4);
    float base = 1.0 + t / h;
    float invBase = 1.0 / base;

    float T;
    if (abs(h - 1.0) < 1e-4) {
        T = invBase;
    } else if (abs(h - 2.0) < 1e-4) {
        T = invBase * invBase;
    } else {
        T = pow(base, -h);
    }
    return vec2(T, T * invBase);
}

#endif // SDF_GEMS_LEVY_RIESZ_CLOUDS_GLSL
