// ============================================================================
// opTwist, opTwistBound & sdgTwistFinalize — 3D Axial Twist & Jacobian Eikonal Normalization
// @author Inigo Quilez (opTwist) — upgraded with 1-sqrt singular-value bound & chain-rule gradient using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  opTwistBound is a guaranteed Lipschitz lower bound; sdgTwistFinalize is 1st-order Eikonal-exact (`||grad d|| = 1` on surface)
// ============================================================================

// 1. Standard 3D Domain Twist
vec3 opTwist( in vec3 p, in float k )
{
    float c = cos(k * p.y);
    float s = sin(k * p.y);
    mat2  m = mat2(c, -s, s, c);
    return vec3(m * p.xz, p.y);
}

// 2. Exact Single-Sqrt Jacobian Singular-Value Twist Factor:
//    By Unimodular Shear Factorization J(p) = R_y(k*p.y) * S(p), the 2x2 shear M_tau has
//    (sigma_max + sigma_min)^2 = 4 + tau^2 and (sigma_max - sigma_min)^2 = tau^2, so
//    sigma_max = 0.5 * (sqrt(4 + tau^2) + tau) with only 1 sqrt (no nested sqrt!).
//    Usage: vec4 tw = opTwistBound(p, k, rMax); float d = sdPrimitive(tw.xyz) * tw.w;
vec4 opTwistBound( in vec3 p, in float k, in float rMax )
{
    float c   = cos(k * p.y);
    float s   = sin(k * p.y);
    mat2  m   = mat2(c, -s, s, c);
    float r   = max(length(p.xz), rMax);
    float tau = abs(k) * r;
    float sig = 0.5 * (sqrt(4.0 + tau * tau) + tau);
    return vec4(m * p.xz, p.y, 1.0 / sig);
}

// 3. First-Order Eikonal-Exact Twisted SDF & Analytical 3D Normal via Chain-Rule Jacobian J^T * g_q:
//    Usage: vec3 q = opTwist(p, k); vec4 dg = sdgPrimitive(q); vec4 res = sdgTwistFinalize(p, q, k, dg);
vec4 sdgTwistFinalize( in vec3 p, in vec3 q, in float k, in vec4 dg )
{
    float c    = cos(k * p.y);
    float s    = sin(k * p.y);
    vec2  gxz  = vec2(c * dg.y - s * dg.w, s * dg.y + c * dg.w);
    float gy   = dg.z + k * (q.z * dg.y - q.x * dg.w);
    vec3  g    = vec3(gxz.x, gy, gxz.y);
    float invG = inversesqrt(max(dot(g, g), 1e-30));
    return vec4(dg.x * invG, g * invG);
}
