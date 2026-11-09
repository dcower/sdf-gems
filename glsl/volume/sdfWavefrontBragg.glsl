// sdfWavefrontBragg: 2x2 Gullstrand-Hamilton wavefront curvature caustics, WKB Maslov index & TPMS Bragg iridescence.
// See docs/volume/sdfWavefrontBragg.md.
#ifndef SDF_GEMS_SDF_WAVEFRONT_BRAGG_GLSL
#define SDF_GEMS_SDF_WAVEFRONT_BRAGG_GLSL

// 2x2 symmetric matrix represented as vec3(K11, K12, K22).

// Extract 2x2 surface shape operator S_surf = E_s^T * H * E_s in the incidence-aligned
// surface tangent basis (u, v_s), given incident ray dIn, unit outward normal n, and 3x3 Hessian H.
vec3 sdfSurfaceShapeOperator(vec3 dIn, vec3 n, mat3 H, out vec3 uAxis, out float cosThetaI) {
    float sgn = (dot(dIn, n) < 0.0) ? 1.0 : -1.0;
    vec3 N = sgn * n;
    cosThetaI = max(-dot(dIn, N), 1e-4);

    vec3 cr = cross(dIn, N);
    float crLen2 = dot(cr, cr);
    if (crLen2 > 1e-12) {
        uAxis = cr * inversesqrt(crLen2);
    } else {
        float s = (N.z >= 0.0) ? 1.0 : -1.0;
        float a = -1.0 / (s + N.z);
        float b = N.x * N.y * a;
        uAxis = vec3(1.0 + s * N.x * N.x * a, s * b, -s * N.x);
    }
    vec3 vSurf = cross(N, uAxis);

    mat3 Hs = sgn * H;
    float s11 = dot(uAxis, Hs * uAxis);
    float s12 = dot(uAxis, Hs * vSurf);
    float s22 = dot(vSurf, Hs * vSurf);
    return vec3(s11, s12, s22);
}

// Exact 2x2 Gullstrand-Hamilton Wavefront Curvature Transfer across a reflective/refractive SDF interface.
// kIn = vec3(K11, K12, K22) in incident frame; eta = eta1 / eta2 (pass eta = 1.0 and isRefraction = false for reflection).
// Returns kOut = vec3(K11_out, K12_out, K22_out).
vec3 sdfGullstrandTransfer(vec3 kIn, vec3 dIn, vec3 n, mat3 H, float eta, bool isRefraction, out vec3 dOut) {
    vec3 uAxis;
    float cosI;
    vec3 sSurf = sdfSurfaceShapeOperator(dIn, n, H, uAxis, cosI);
    vec3 N = (dot(dIn, n) < 0.0) ? n : -n;

    float cosO;
    float gammaNorm;
    float etaRatio;

    if (!isRefraction) {
        dOut = dIn + 2.0 * cosI * N;
        cosO = cosI;
        gammaNorm = 2.0 * cosI;
        etaRatio = 1.0;
    } else {
        float sin2T = eta * eta * max(1.0 - cosI * cosI, 0.0);
        if (sin2T >= 1.0) {
            dOut = dIn + 2.0 * cosI * N;
            cosO = cosI;
            gammaNorm = 2.0 * cosI;
            etaRatio = 1.0;
        } else {
            float cosT = sqrt(max(1.0 - sin2T, 1e-8));
            dOut = eta * dIn + (eta * cosI - cosT) * N;
            cosO = -cosT;
            gammaNorm = eta * cosI - cosT;
            etaRatio = eta;
        }
    }

    float invCosO = 1.0 / (sign(cosO) * max(abs(cosO), 1e-4));
    float kOut11 = etaRatio * kIn.x + gammaNorm * sSurf.x;
    float kOut12 = (etaRatio * cosI * kIn.y + gammaNorm * sSurf.y) * invCosO;
    float kOut22 = (etaRatio * cosI * cosI * kIn.z + gammaNorm * sSurf.z) * (invCosO * invCosO);
    return vec3(kOut11, kOut12, kOut22);
}

// Propagate 2x2 wavefront curvature kOut = vec3(K11, K12, K22) by distance s >= 0.
// Returns vec3(causticIntensity, jacobianDet, maslovIndex) where maslovIndex in {0, 1, 2}.
vec3 sdfWavefrontCausticPropagate(vec3 kOut, float s, float airyEps) {
    float m = 0.5 * (kOut.x + kOut.z);
    float d = 0.5 * (kOut.x - kOut.z);
    float disc = sqrt(max(d * d + kOut.y * kOut.y, 0.0));
    float kappa1 = m + disc;
    float kappa2 = m - disc;

    float f1 = 1.0 + s * kappa1;
    float f2 = 1.0 + s * kappa2;
    float J = f1 * f2;

    float maslov = float(f1 < 0.0) + float(f2 < 0.0);
    float intensity = inversesqrt(J * J + max(airyEps * airyEps, 1e-8));
    return vec3(intensity, J, maslov);
}

// Closed-form O(1) Fourier-Laplace Bragg Iridescence across a Gyroid photonic crystal slab [ta, tb].
// q0: Gyroid lattice wavevector scale (2*pi / period)
// kRGB: Optical wavenumbers 2*pi*n_eff / lambda for (R, G, B), e.g. q0 * vec3(0.55, 0.68, 0.78)
// alpha: Extinction / attenuation coefficient > 0
vec3 sdfGyroidBraggIridescence(vec3 rd, float ta, float tb, float q0, vec3 kRGB, float alpha) {
    float dt = max(tb - ta, 0.0);
    if (dt <= 0.0) return vec3(0.0);

    const vec3 G_DIRS[6] = vec3[6](
        vec3( 1.0,  1.0,  0.0),
        vec3( 1.0, -1.0,  0.0),
        vec3( 0.0,  1.0,  1.0),
        vec3( 0.0,  1.0, -1.0),
        vec3( 1.0,  0.0,  1.0),
        vec3(-1.0,  0.0,  1.0)
    );

    float exp1 = exp(-alpha * dt);
    float exp2 = exp1 * exp1;
    float a2 = max(alpha * alpha, 1e-6);

    vec3 accumRGB = vec3(0.0);
    for (int m = 0; m < 6; ++m) {
        float omegaM = abs(q0 * dot(G_DIRS[m], rd));
        vec3 detune = 2.0 * kRGB - vec3(omegaM);
        vec3 num = vec3(1.0 + exp2) - 2.0 * exp1 * cos(detune * dt);
        vec3 den = vec3(a2) + detune * detune;
        accumRGB += a2 * (num / den);
    }

    return accumRGB / (vec3(1.0) + accumRGB);
}

#endif // SDF_GEMS_SDF_WAVEFRONT_BRAGG_GLSL
