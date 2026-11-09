// sdfChladniOrbitals: Spectral Courant nodal domains, unit-gradient Chladni eigen-SDFs & quantum orbital volume clouds.
// See docs/volume/sdfChladniOrbitals.md.
#ifndef SDF_GEMS_SDF_CHLADNI_ORBITALS_GLSL
#define SDF_GEMS_SDF_CHLADNI_ORBITALS_GLSL

#ifndef SDF_GEMS_HAS_SAFE_SINC
#define SDF_GEMS_HAS_SAFE_SINC
// FP32-safe unnormalized sinc(u) = sin(u) / u
float sdfSafeSinc(float u) {
    float u2 = u * u;
    return (abs(u) > 1e-4) ? (sin(u) / u) : (1.0 - u2 * (1.0 / 6.0 - u2 * (1.0 / 120.0)));
}
#endif

#ifndef SDF_GEMS_HAS_INTEGRATE_COS_SLAB
#define SDF_GEMS_HAS_INTEGRATE_COS_SLAB
// Exact integral of cos(omega * t + phi) over t in [ta, tb]
float sdfIntegrateCosSlab(float omega, float phi, float ta, float tb) {
    float dt = tb - ta;
    float tMid = 0.5 * (ta + tb);
    return dt * sdfSafeSinc(0.5 * omega * dt) * cos(omega * tMid + phi);
}
#endif

#ifndef SDF_GEMS_HAS_INTEGRATE_COS2_SLAB
#define SDF_GEMS_HAS_INTEGRATE_COS2_SLAB
// Exact integral of cos^2(omega * t + phi) over t in [ta, tb] (>= 0)
float sdfIntegrateCos2Slab(float omega, float phi, float ta, float tb) {
    float dt = tb - ta;
    float tMid = 0.5 * (ta + tb);
    return 0.5 * dt * (1.0 + sdfSafeSinc(omega * dt) * cos(2.0 * (omega * tMid + phi)));
}
#endif

// Evaluate 3D Chladni Helmholtz Eigen-SDF and its exact analytical gradient:
// psi(x) = w.x * cos(qx)*cos(qy)*cos(qz)
//        + w.y * sin(qx)*sin(qy)*sin(qz)
//        + w.z * (cos(qx)*sin(qy)*cos(qz) + cos(qy)*sin(qz)*cos(qx) + cos(qz)*sin(qx)*cos(qy))
// Returns vec4(d_Chladni, grad_d_Chladni) with zero 0/0 singularities at nodal intersections.
vec4 sdgChladniEigen(vec3 p, float q, vec3 w) {
    vec3 qp = q * p;
    vec3 c = cos(qp);
    vec3 s = sin(qp);

    float psi1 = c.x * c.y * c.z;
    float psi2 = s.x * s.y * s.z;
    float psi3 = c.x * s.y * c.z + c.y * s.z * c.x + c.z * s.x * c.y;
    float psi = w.x * psi1 + w.y * psi2 + w.z * psi3;

    // First derivatives (divided by q)
    vec3 g1 = vec3(-s.x * c.y * c.z, -c.x * s.y * c.z, -c.x * c.y * s.z);
    vec3 g2 = vec3( c.x * s.y * s.z,  s.x * c.y * s.z,  s.x * s.y * c.z);
    vec3 g3 = vec3(
        c.x * c.y * c.z - c.y * s.x * s.z - c.z * s.x * s.y,
        c.x * c.y * c.z - c.x * s.y * s.z - c.z * s.x * s.y,
        c.x * c.y * c.z - c.x * s.y * s.z - c.y * s.x * s.z
    );
    vec3 g = q * (w.x * g1 + w.y * g2 + w.z * g3);

    // Off-diagonal second derivatives (H_yz, H_xz, H_xy)
    float q2 = q * q;
    vec3 hOff1 = vec3( c.x * s.y * s.z,  s.x * c.y * s.z,  s.x * s.y * c.z);
    vec3 hOff2 = vec3( s.x * c.y * c.z,  c.x * s.y * c.z,  c.x * c.y * s.z);
    vec3 hOff3 = vec3(
        s.x * s.y * s.z - c.x * (c.y * s.z + c.z * s.y),
        s.x * s.y * s.z - c.y * (c.x * s.z + c.z * s.x),
        s.x * s.y * s.z - c.z * (c.x * s.y + c.y * s.x)
    );
    // hOff = vec3(H_yz, H_xz, H_xy)
    vec3 hOff = q2 * (w.x * hOff1 + w.y * hOff2 + w.z * hOff3);

    // Third mixed derivative T = d^3 psi / (dx dy dz)
    float q3 = q2 * q;
    float T1 = -s.x * s.y * s.z;
    float T2 =  c.x * c.y * c.z;
    float T3 =  c.x * s.y * s.z + c.y * s.x * s.z + c.z * s.x * s.y;
    float T = q3 * (w.x * T1 + w.y * T2 + w.z * T3);

    // Total Helmholtz eigenvalue k^2 = 3 * q^2
    float k2 = 3.0 * q2;
    float inv2k2 = 1.0 / max(2.0 * k2, 1e-12);

    // Diagonal of H is -q^2 * psi, so ||H||_F^2 = 3 * q^4 * psi^2 + 2 * ||hOff||^2
    float hFrob2 = 3.0 * q2 * q2 * psi * psi + 2.0 * dot(hOff, hOff);

    // Phase-space energy E(x)
    float E = k2 * psi * psi + dot(g, g) + inv2k2 * hFrob2 + 1e-12;
    float invSqrtE = inversesqrt(E);
    float dChladni = psi * invSqrtE;

    // Matrix-vector product H * g
    vec3 Hg = -q2 * psi * g + vec3(
        hOff.z * g.y + hOff.y * g.z,
        hOff.z * g.x + hOff.x * g.z,
        hOff.y * g.x + hOff.x * g.y
    );

    // Exact gradient of 0.5 * ||H||_F^2 using 3D Harmonic Jet Closure Identity
    vec3 halfGradHFrob2 = (q2 * q2 * psi) * g - (2.0 * q2) * Hg + (2.0 * T) * hOff;

    // Exact gradient of 0.5 * E(x)
    vec3 halfGradE = (k2 * psi) * g + Hg + inv2k2 * halfGradHFrob2;

    vec3 gradD = invSqrtE * (g - (psi / E) * halfGradE);
    return vec4(dChladni, gradD);
}

// Scalar wrapper for sdChladniEigen
float sdChladniEigen(vec3 p, float q, vec3 w) {
    return sdgChladniEigen(p, q, w).x;
}

// Closed-form O(1) 1D ray slab integral of the Born probability cloud
// rho(x) = |psi(x)|^2 for the 4 fundamental cubic Chladni wavevectors (+/-1, +/-1, +/-1)*q.
float sdfChladniBornRayIntegral(vec3 ro, vec3 rd, float ta, float tb, float q) {
    float dt = max(tb - ta, 0.0);
    if (dt <= 0.0) return 0.0;

    const vec3 K_DIRS[4] = vec3[4](
        vec3( 1.0,  1.0,  1.0),
        vec3( 1.0,  1.0, -1.0),
        vec3( 1.0, -1.0,  1.0),
        vec3(-1.0,  1.0,  1.0)
    );

    float omega[4];
    float phi[4];
    for (int i = 0; i < 4; ++i) {
        vec3 kVec = q * K_DIRS[i];
        omega[i] = dot(kVec, rd);
        phi[i] = dot(kVec, ro);
    }

    float sumInt = 0.0;
    for (int i = 0; i < 4; ++i) {
        sumInt += sdfIntegrateCos2Slab(omega[i], phi[i], ta, tb);
        for (int j = i + 1; j < 4; ++j) {
            sumInt += sdfIntegrateCosSlab(omega[i] - omega[j], phi[i] - phi[j], ta, tb);
            sumInt += sdfIntegrateCosSlab(omega[i] + omega[j], phi[i] + phi[j], ta, tb);
        }
    }
    return max(0.0625 * sumInt, 0.0);
}

#endif // SDF_GEMS_SDF_CHLADNI_ORBITALS_GLSL
