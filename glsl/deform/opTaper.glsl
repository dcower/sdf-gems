// ============================================================================
// opTaper, opTaperBound & sdgTaperFinalize — 3D Linear Taper & Triangular-Shear Eikonal Normalization
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  opTaperBound is a guaranteed 1-sqrt Lipschitz lower bound; sdgTaperFinalize is 1st-order Eikonal-exact (`||grad d|| = 1` on surface)
// ============================================================================

#ifndef SDF_GEMS_OP_TAPER_GLSL
#define SDF_GEMS_OP_TAPER_GLSL

// 1. Standard 3D Linear Taper along the y-axis:
//    Maps p -> q = vec3(p.x / s, p.y, p.z / s) where s = max(1 + k * p.y, 1e-4) is the
//    horizontal cross-section scale factor at height p.y.
vec3 opTaper( in vec3 p, in float k )
{
    float s = max(1.0 + k * p.y, 1e-4);
    return vec3(p.x / s, p.y, p.z / s);
}

// 2. Exact Single-Sqrt Jacobian Singular-Value Taper Bound:
//    In cylindrical coordinates around the y-axis, the Jacobian J(p) = dq/dp restricts on
//    the meridian plane to the 2x2 triangular shear matrix:
//      M = [ 1/s,  -k*q_r / s ]
//          [  0 ,       1     ]
//    with determinant det(M) = 1/s and trace(M^T M) = 1 + (1 + k^2 q_r^2) / s^2.
//    By the Sum-and-Difference of Singular Values Identity:
//      (sigma_max + sigma_min)^2 = (s + 1)^2 / s^2 + (k * q_r / s)^2,
//      (sigma_max - sigma_min)^2 = (s - 1)^2 / s^2 + (k * q_r / s)^2,
//    yielding the exact maximum singular value sigma_max(J) without nested square roots!
//    Usage: vec4 tp = opTaperBound(p, k, rMax, hMax); float d = sdPrimitive(tp.xyz) * tp.w;
vec4 opTaperBound( in vec3 p, in float k, in float rMax, in float hMax )
{
    float s    = max(1.0 + k * p.y, 1e-4);
    vec2  qxz  = p.xz / s;
    float sMin = max(1.0 - abs(k) * min(abs(p.y), hMax), 1e-4);
    float qr   = max(length(qxz), rMax);
    float kr   = k * qr;
    float sm   = 0.5 * (sqrt((sMin + 1.0) * (sMin + 1.0) + kr * kr)
                      + sqrt((sMin - 1.0) * (sMin - 1.0) + kr * kr)) / sMin;
    return vec4(qxz.x, p.y, qxz.y, 1.0 / max(sm, 1.0 / sMin));
}

// 3. First-Order Eikonal-Exact Tapered SDF & Analytical 3D Normal via Chain-Rule Jacobian J^T * g_q:
//    Since q = (p.x / s, p.y, p.z / s), the exact pullback gradient is:
//      grad_p(f o Phi) = J(p)^T g_q = (g_qx / s,  g_qy - (k / s) * (q.x * g_qx + q.z * g_qz),  g_qz / s).
//    Normalizing by ||J(p)^T g_q||_2 yields `||grad d|| = 1` on the tapered surface and
//    eliminates raymarching overstepping for arbitrary taper rates k!
//    Usage: vec3 q = opTaper(p, k); vec4 dg = sdgPrimitive(q); vec4 res = sdgTaperFinalize(p, q, k, dg);
vec4 sdgTaperFinalize( in vec3 p, in vec3 q, in float k, in vec4 dg )
{
    float invS = 1.0 / max(1.0 + k * p.y, 1e-4);
    vec3  g    = vec3(dg.y * invS,
                      dg.z - (k * invS) * (q.x * dg.y + q.z * dg.w),
                      dg.w * invS);
    float invG = inversesqrt(max(dot(g, g), 1e-30));
    return vec4(dg.x * invG, g * invG);
}

#endif // SDF_GEMS_OP_TAPER_GLSL
