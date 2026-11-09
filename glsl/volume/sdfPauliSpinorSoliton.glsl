// sdfPauliSpinorSoliton: Non-Abelian Hermitian Pauli-Spinor SDFs,
// SL(2,C) Photoelastic Stress-Birefringent Volumes, & Hirota determinant(mat4) Soliton Ray Collapse.
// Mathematical framework attributed to Gemini 4 Argon.
// Implementation and verification by Gemini 3.8 Flash (sdf-gems, 2026).
// See docs/volume/sdfPauliSpinorSoliton.md.

#ifndef SDF_GEMS_SDF_PAULI_SPINOR_SOLITON_GLSL
#define SDF_GEMS_SDF_PAULI_SPINOR_SOLITON_GLSL

// -----------------------------------------------------------------------------
// 1. Hermitian Pauli-Spinor SDFs & Exact Von Neumann-Wigner Avoided Crossing
// -----------------------------------------------------------------------------

// Computes the inner and outer avoided-crossing eigenvalues (lambda_-, lambda_+)
// of the 2x2 Hermitian matrix H(x) = d0 * I_2 + b . sigma:
// H = [ d0 + bz     bx - i by ]
//     [ bx + i by   d0 - bz   ]
// packed into vec4 H = vec4(d0, bx, by, bz).
// Returns vec2(lambda_minus, lambda_plus).
vec2 sdfPauliEigenvalues(vec4 H)
{
    float bLen = length(H.yzw);
    return vec2(H.x - bLen, H.x + bLen);
}

// Scalar smooth union via von Neumann-Wigner avoided crossing (coupling > 0).
// Zero exp/log, perfectly symmetric, Lipschitz-bounded.
float sdPauliSpinorUnion(float da, float db, float coupling)
{
    float d0 = 0.5 * (da + db);
    float bz = 0.5 * (da - db);
    return d0 - sqrt(bz * bz + coupling * coupling + 1e-12);
}

// Exact algebraic von Neumann-Wigner avoided crossing smooth union of two SDF+gradient
// fields dga = vec4(da, gxa, gya, gza) and dgb = vec4(db, gxb, gyb, gzb),
// coupled off-diagonally by coupling = vec2(bx, by).
// Returns vec4(lambda_minus, grad_lambda_minus).
// Strictly 1-Lipschitz: ||grad lambda_-|| <= 1 whenever ||ga|| <= 1 and ||gb|| <= 1.
vec4 sdgPauliSpinorUnion(vec4 dga, vec4 dgb, vec2 coupling)
{
    float d0 = 0.5 * (dga.x + dgb.x);
    vec3 grad_d0 = 0.5 * (dga.yzw + dgb.yzw);

    float bz = 0.5 * (dga.x - dgb.x);
    vec3 grad_bz = 0.5 * (dga.yzw - dgb.yzw);

    vec3 b = vec3(coupling.x, coupling.y, bz);
    float bLen = sqrt(dot(b, b) + 1e-12);

    float lambdaMinus = d0 - bLen;
    vec3 gradLambdaMinus = grad_d0 - (bz / bLen) * grad_bz;

    return vec4(lambdaMinus, gradLambdaMinus);
}

// Returns the normalized Bloch spinor vector b_hat = b / max(||b||, 1e-6)
// representing the principal optical / stress anisotropy axis on the Bloch sphere.
vec3 sdfPauliBlochFrame(vec4 H)
{
    vec3 b = H.yzw;
    float len = length(b);
    return b / max(len, 1e-6);
}

// -----------------------------------------------------------------------------
// 2. SL(2,C) Photoelastic Stress-Birefringent Volume Transmission
// -----------------------------------------------------------------------------

// Evaluates transmitted RGB intensity between crossed polarizers (Jones 1941)
// for an optical medium under stressPauli = vec3(bx, by, bz) with isotropic opticalDepth.
// Produces the dark isoclinic cross sin^2(2 phi_az) and rainbow isochromatic
// interference fringes sin^2(pi * ||b_stress|| / lambda).
vec3 sdfPauliBirefringence(float opticalDepth, vec3 stressPauli, vec3 rgbWavelengths)
{
    float bx2 = stressPauli.x * stressPauli.x;
    float bz2 = stressPauli.z * stressPauli.z;
    float denom = bx2 + bz2 + 1e-6;
    float sin2_2phi = (4.0 * bx2 * bz2) / (denom * denom);

    float bMag = length(stressPauli);
    vec3 phase = 3.141592653589793 * bMag / max(rgbWavelengths, vec3(1e-4));
    vec3 sinPhase = sin(phase);
    vec3 isochromatic = sinPhase * sinPhase;

    float atten = exp(-max(opticalDepth, 0.0));
    return atten * sin2_2phi * isochromatic;
}

// -----------------------------------------------------------------------------
// 3. Hirota determinant(mat4) 4-Soliton Tau-Function & 0-Loop Ray Integrals
// -----------------------------------------------------------------------------

// Canonical 4-soliton wavevectors, frequencies, and spectral parameters:
// Solitons 0, 1, 2, 3 in 3D space:
const vec3 SDF_HIROTA_K0 = vec3(1.5, 0.5, 0.0);
const vec3 SDF_HIROTA_K1 = vec3(-1.0, 1.2, 0.2);
const vec3 SDF_HIROTA_K2 = vec3(0.3, -1.4, 0.5);
const vec3 SDF_HIROTA_K3 = vec3(-0.8, -0.6, -1.0);

const vec4 SDF_HIROTA_OMEGA = vec4(1.0, 1.3, 1.6, 2.0);

// Precomputed Cauchy interaction matrix C_ij = 2 * sqrt(kappa_i * kappa_j) / (kappa_i + kappa_j)
// for spectral parameters kappa = vec4(1.0, 1.4, 1.8, 2.2):
// C_00 = 1.0, C_01 = 0.9859296, C_02 = 0.9599657, C_03 = 0.9270248
// C_11 = 1.0, C_12 = 0.9949367, C_13 = 0.9744795
// C_22 = 1.0, C_23 = 0.9949367
// C_33 = 1.0
const mat4 SDF_HIROTA_CAUCHY = mat4(
    1.0000000, 0.9859296, 0.9599657, 0.9270248,
    0.9859296, 1.0000000, 0.9949367, 0.9744795,
    0.9599657, 0.9949367, 1.0000000, 0.9949367,
    0.9270248, 0.9744795, 0.9949367, 1.0000000
);

// Evaluates the 4x4 soliton interaction matrix A(x, t) = I_4 + M(x, t).
mat4 sdfHirotaMatrix(vec3 p, float time)
{
    vec4 eta = vec4(
        dot(SDF_HIROTA_K0, p) - SDF_HIROTA_OMEGA.x * time,
        dot(SDF_HIROTA_K1, p) - SDF_HIROTA_OMEGA.y * time,
        dot(SDF_HIROTA_K2, p) - SDF_HIROTA_OMEGA.z * time,
        dot(SDF_HIROTA_K3, p) - SDF_HIROTA_OMEGA.w * time
    );
    vec4 e = exp(0.5 * clamp(eta, -20.0, 20.0));

    mat4 M;
    M[0] = vec4(SDF_HIROTA_CAUCHY[0][0]*e.x*e.x, SDF_HIROTA_CAUCHY[0][1]*e.x*e.y, SDF_HIROTA_CAUCHY[0][2]*e.x*e.z, SDF_HIROTA_CAUCHY[0][3]*e.x*e.w);
    M[1] = vec4(SDF_HIROTA_CAUCHY[1][0]*e.y*e.x, SDF_HIROTA_CAUCHY[1][1]*e.y*e.y, SDF_HIROTA_CAUCHY[1][2]*e.y*e.z, SDF_HIROTA_CAUCHY[1][3]*e.y*e.w);
    M[2] = vec4(SDF_HIROTA_CAUCHY[2][0]*e.z*e.x, SDF_HIROTA_CAUCHY[2][1]*e.z*e.y, SDF_HIROTA_CAUCHY[2][2]*e.z*e.z, SDF_HIROTA_CAUCHY[2][3]*e.z*e.w);
    M[3] = vec4(SDF_HIROTA_CAUCHY[3][0]*e.w*e.x, SDF_HIROTA_CAUCHY[3][1]*e.w*e.y, SDF_HIROTA_CAUCHY[3][2]*e.w*e.z, SDF_HIROTA_CAUCHY[3][3]*e.w*e.w);

    return mat4(1.0) + M;
}

// Evaluates the Hirota potential Phi(x, t) = 2 d/ds ln det(A)
// and volumetric soliton density rho(x, t) = 2 d^2/ds^2 ln det(A)
// along ray direction rd at point p and time t.
// Returns vec2(Phi, rho).
vec2 sdfHirotaSoliton4(vec3 p, vec3 rd, float time)
{
    vec4 eta = vec4(
        dot(SDF_HIROTA_K0, p) - SDF_HIROTA_OMEGA.x * time,
        dot(SDF_HIROTA_K1, p) - SDF_HIROTA_OMEGA.y * time,
        dot(SDF_HIROTA_K2, p) - SDF_HIROTA_OMEGA.z * time,
        dot(SDF_HIROTA_K3, p) - SDF_HIROTA_OMEGA.w * time
    );
    vec4 e = exp(0.5 * clamp(eta, -20.0, 20.0));

    mat4 M;
    M[0] = vec4(SDF_HIROTA_CAUCHY[0][0]*e.x*e.x, SDF_HIROTA_CAUCHY[0][1]*e.x*e.y, SDF_HIROTA_CAUCHY[0][2]*e.x*e.z, SDF_HIROTA_CAUCHY[0][3]*e.x*e.w);
    M[1] = vec4(SDF_HIROTA_CAUCHY[1][0]*e.y*e.x, SDF_HIROTA_CAUCHY[1][1]*e.y*e.y, SDF_HIROTA_CAUCHY[1][2]*e.y*e.z, SDF_HIROTA_CAUCHY[1][3]*e.y*e.w);
    M[2] = vec4(SDF_HIROTA_CAUCHY[2][0]*e.z*e.x, SDF_HIROTA_CAUCHY[2][1]*e.z*e.y, SDF_HIROTA_CAUCHY[2][2]*e.z*e.z, SDF_HIROTA_CAUCHY[2][3]*e.z*e.w);
    M[3] = vec4(SDF_HIROTA_CAUCHY[3][0]*e.w*e.x, SDF_HIROTA_CAUCHY[3][1]*e.w*e.y, SDF_HIROTA_CAUCHY[3][2]*e.w*e.z, SDF_HIROTA_CAUCHY[3][3]*e.w*e.w);

    mat4 A = mat4(1.0) + M;
    mat4 invA = inverse(A);

    vec4 v = vec4(
        dot(SDF_HIROTA_K0, rd),
        dot(SDF_HIROTA_K1, rd),
        dot(SDF_HIROTA_K2, rd),
        dot(SDF_HIROTA_K3, rd)
    );

    mat4 Aprime;
    Aprime[0] = 0.5 * (v.x + v) * M[0];
    Aprime[1] = 0.5 * (v.y + v) * M[1];
    Aprime[2] = 0.5 * (v.z + v) * M[2];
    Aprime[3] = 0.5 * (v.w + v) * M[3];

    mat4 B = invA * Aprime;
    float phi = 2.0 * (B[0][0] + B[1][1] + B[2][2] + B[3][3]);

    mat4 Adprime;
    vec4 s0 = v.x + v;
    vec4 s1 = v.y + v;
    vec4 s2 = v.z + v;
    vec4 s3 = v.w + v;
    Adprime[0] = 0.25 * (s0 * s0) * M[0];
    Adprime[1] = 0.25 * (s1 * s1) * M[1];
    Adprime[2] = 0.25 * (s2 * s2) * M[2];
    Adprime[3] = 0.25 * (s3 * s3) * M[3];

    mat4 invA_Adprime = invA * Adprime;
    mat4 B2 = B * B;
    float trA_Adprime = invA_Adprime[0][0] + invA_Adprime[1][1] + invA_Adprime[2][2] + invA_Adprime[3][3];
    float trB2 = B2[0][0] + B2[1][1] + B2[2][2] + B2[3][3];
    float rho = 2.0 * (trA_Adprime - trB2);

    return vec2(phi, rho);
}

// Ray integral of the soliton volumetric density rho(s) between ta and tb.
// Collapses to exactly 2 endpoint evaluations with zero raymarching loops:
// integral_{ta}^{tb} rho(ro + s * rd, t) ds = Phi(tb) - Phi(ta)
float sdfHirotaSolitonRayIntegral(vec3 ro, vec3 rd, float ta, float tb, float time)
{
    float phiA = sdfHirotaSoliton4(ro + ta * rd, rd, time).x;
    float phiB = sdfHirotaSoliton4(ro + tb * rd, rd, time).x;
    return phiB - phiA;
}

// Ray integral of the Hirota potential Phi(s) between ta and tb:
// integral_{ta}^{tb} Phi(ro + s * rd, t) ds = 2 * (ln det A(tb) - ln det A(ta))
float sdfHirotaPotentialRayIntegral(vec3 ro, vec3 rd, float ta, float tb, float time)
{
    mat4 Aa = sdfHirotaMatrix(ro + ta * rd, time);
    mat4 Ab = sdfHirotaMatrix(ro + tb * rd, time);
    float detA = max(determinant(Aa), 1e-12);
    float detB = max(determinant(Ab), 1e-12);
    return 2.0 * (log(detB) - log(detA));
}

#endif // SDF_GEMS_SDF_PAULI_SPINOR_SOLITON_GLSL
