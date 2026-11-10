// ============================================================================
// paUltrametricTree — Non-Archimedean p-adic Cantor Dust, Bruhat-Tits Trees
//                     & Lipschitz-Regularized Monna Wave Pullbacks
//   sdCantorDust1D, sdgCantorDust1D, sdCantorDust3D, sdgCantorDust3D,
//   sdBruhatTitsTree2D, sdBruhatTitsTree3D, sdBruhatTitsCanopy3D,
//   paMonnaTentMap1D, sdMonnaTerrain3D
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  sdCantorDust1D, sdgCantorDust1D, sdCantorDust3D, sdgCantorDust3D:
//           100% Exact Euclidean SDF (interior & exterior) for all alpha in (0, 0.5) and n >= 0.
//         sdBruhatTitsTree2D, sdBruhatTitsTree3D, sdBruhatTitsCanopy3D:
//           Exact union of 2^(n+1)-1 (2D) or (4^(n+1)-1)/3 (3D) tapered round-cone branches
//           iff 2*alpha - alpha^n <= 1 (uniformly for all n iff alpha <= 0.5); exact on exterior,
//           conservative min-union bound inside overlapping joints.
//         sdMonnaTerrain3D:
//           Strict 1-Lipschitz conservative Euclidean SDF bound (normalized by exact Lipschitz slope).
// ============================================================================

#ifndef SDF_GEMS_HAS_PA_ULTRAMETRIC_TREE
#define SDF_GEMS_HAS_PA_ULTRAMETRIC_TREE

// ---------------------------------------------------------------------------
// Internal helper: Exact tapered capsule (round cone) between 3D points a, b
// with radii r1, r2. Written with product radicands for FP32 safety.
// ---------------------------------------------------------------------------
float paRoundCone3D( in vec3 p, in vec3 a, in vec3 b, in float r1, in float r2 )
{
    vec3  ba  = b - a;
    float l2  = dot(ba, ba);
    float rr  = r1 - r2;
    float a2  = l2 - rr * rr;
    float il2 = 1.0 / l2;
    vec3  pa  = p - a;
    float y   = dot(pa, ba);
    float z   = y - l2;
    vec3  w   = pa * l2 - ba * y;
    float x2  = dot(w, w);
    float y2  = y * y * l2;
    float z2  = z * z * l2;
    float k   = sign(rr) * rr * rr * x2;
    if (sign(z) * a2 * z2 > k) return sqrt(x2 + z2) * il2 - r2;
    if (sign(y) * a2 * y2 < k) return sqrt(x2 + y2) * il2 - r1;
    return (sqrt(max(x2 * a2 * il2, 0.0)) + y * rr) * il2 - r1;
}

float paRoundCone2D( in vec2 p, in vec2 a, in vec2 b, in float r1, in float r2 )
{
    return paRoundCone3D(vec3(p, 0.0), vec3(a, 0.0), vec3(b, 0.0), r1, r2);
}

// ---------------------------------------------------------------------------
// 1. Exact 1D & 3D Generalized Cantor Set (p-adic Integer Embedding Z_2 -> R^d)
//    For contraction ratio alpha in (0, 0.5) (classical middle-third: alpha = 1/3),
//    the greedy reflection fold x -> |x| - (1 - alpha) * alpha^k visits a single
//    branch and returns the 100% exact Euclidean SDF of all 2^n intervals (1D)
//    or 8^n cubes (3D) in 1 sqrt.
// ---------------------------------------------------------------------------
vec2 sdgCantorDust1D( in float x, in int n, in float alpha )
{
    float s = 1.0;
    float L = 1.0;
    for (int k = 0; k < 12; k++)
    {
        if (k >= n) break;
        s *= (x >= 0.0) ? 1.0 : -1.0;
        x  = abs(x) - (1.0 - alpha) * L;
        L *= alpha;
    }
    s *= (x >= 0.0) ? 1.0 : -1.0;
    return vec2(abs(x) - L, s);
}

float sdCantorDust1D( in float x, in int n, in float alpha )
{
    return sdgCantorDust1D(x, n, alpha).x;
}

vec4 sdgCantorDust3D( in vec3 p, in int n, in float alpha )
{
    vec2 dx = sdgCantorDust1D(p.x, n, alpha);
    vec2 dy = sdgCantorDust1D(p.y, n, alpha);
    vec2 dz = sdgCantorDust1D(p.z, n, alpha);
    vec3 d  = vec3(dx.x, dy.x, dz.x);
    vec3 s  = vec3(dx.y, dy.y, dz.y);
    vec3 w  = max(d, 0.0);
    float l = length(w);
    float m = max(d.x, max(d.y, d.z));
    vec3  gInside = (d.x >= d.y && d.x >= d.z) ? vec3(s.x, 0.0, 0.0) :
                    ((d.y >= d.z)              ? vec3(0.0, s.y, 0.0) : vec3(0.0, 0.0, s.z));
    vec3  g = (m > 0.0) ? ((l > 1e-20) ? (w / l) * s : gInside) : gInside;
    return vec4(l + min(m, 0.0), g);
}

float sdCantorDust3D( in vec3 p, in int n, in float alpha )
{
    return sdgCantorDust3D(p, n, alpha).x;
}

// ---------------------------------------------------------------------------
// 2. Exact 2D & 3D Bruhat-Tits Tree Skeletons (T_2 and T_2 x T_2)
//    By the Voronoi Bisector Containment Theorem, the single-path greedy fold
//    x -> |x| - (1 - alpha)*L_k yields the exact minimum distance over all
//    2^(n+1)-1 (2D) or (4^(n+1)-1)/3 (3D) tapered round-cone branches iff
//    2*alpha - alpha^n <= 1 (which holds for all n >= 1 iff alpha <= 0.5).
// ---------------------------------------------------------------------------
float sdBruhatTitsTree2D( in vec2 p, in int n, in float alpha, in float r0, in float gammaR )
{
    const float y0 = -0.45, y1 = 0.62, gammaY = 0.64;
    float d  = paRoundCone2D(p, vec2(0.0, -0.78), vec2(0.0, y0), r0 * 1.25, r0);
    float gn = pow(gammaY, float(max(n, 1)));
    float dy = (y1 - y0) * (1.0 - gammaY) / max(1.0 - gn, 1e-6);
    float yk = y0, rk = r0, L = 0.82, stepY = dy;
    float x  = p.x, y = p.y;
    for (int k = 0; k < 10; k++)
    {
        if (k >= n) break;
        float ck  = (1.0 - alpha) * L;
        float yk1 = yk + stepY;
        float rk1 = rk * gammaR;
        x = abs(x);
        d = min(d, paRoundCone2D(vec2(x, y), vec2(0.0, yk), vec2(ck, yk1), rk, rk1));
        x -= ck;
        yk = yk1;
        rk = rk1;
        L *= alpha;
        stepY *= gammaY;
    }
    return d;
}

float sdBruhatTitsTree3D( in vec3 p, in int n, in float alpha, in float r0, in float gammaR )
{
    const float y0 = -0.35, y1 = 0.48, gammaY = 0.62;
    float d  = paRoundCone3D(p, vec3(0.0, -0.64, 0.0), vec3(0.0, y0, 0.0), r0 * 1.3, r0);
    float gn = pow(gammaY, float(max(n, 1)));
    float dy = (y1 - y0) * (1.0 - gammaY) / max(1.0 - gn, 1e-6);
    float yk = y0, rk = r0, L = 0.68, stepY = dy;
    vec3  q  = p;
    for (int k = 0; k < 10; k++)
    {
        if (k >= n) break;
        float ck  = (1.0 - alpha) * L;
        float yk1 = yk + stepY;
        float rk1 = rk * gammaR;
        q.xz = abs(q.xz);
        d = min(d, paRoundCone3D(q, vec3(0.0, yk, 0.0), vec3(ck, yk1, ck), rk, rk1));
        q.xz -= ck;
        yk = yk1;
        rk = rk1;
        L *= alpha;
        stepY *= gammaY;
    }
    return d;
}

// 3D Bruhat-Tits Coral Canopy: 4-regular Bruhat-Tits tree crowned with a 4^n
// Cantor dust slab canopy at the p-adic boundary at infinity.
float sdBruhatTitsCanopy3D( in vec3 p, in int n, in float alpha )
{
    float dTree = sdBruhatTitsTree3D(p, n, alpha, 0.085, 0.66);
    // Crown: 4^n Cantor dust blocks at the top level y1 = 0.48, scaled to L0 = 0.68
    const float L0 = 0.68;
    float dx = sdCantorDust1D(p.x / L0, n, alpha) * L0;
    float dz = sdCantorDust1D(p.z / L0, n, alpha) * L0;
    float dy = abs(p.y - 0.51) - 0.045;
    vec3  d3 = vec3(dx, dy, dz);
    float dCrown = length(max(d3, 0.0)) + min(max(d3.x, max(d3.y, d3.z)), 0.0) - 0.012;
    return min(dTree, dCrown);
}

// ---------------------------------------------------------------------------
// 3. Lipschitz-Regularized Monna Map (p-adic Gray-Code / Tent Wavelet Pullback)
//    Classical Monna digit reversal sum a_k p^k -> sum a_k p^(-k-1) has discontinuous
//    carry jumps (infinite Lipschitz constant). Replacing discontinuous base-p digits
//    with continuous reflected tent digits tau_p(u) = 1 - |2*fract(u) - 1| yields a
//    continuous self-similar p-adic wavelet relief with exact Lipschitz constant.
// ---------------------------------------------------------------------------
float paMonnaTentMap1D( in float x, in int n, in float pBase )
{
    float val   = 0.0;
    float freq  = 1.0;
    float scale = 1.0 / pBase;
    float norm  = 0.0;
    for (int k = 0; k < 10; k++)
    {
        if (k >= n) break;
        float u = fract(x * freq);
        float tent = 1.0 - abs(2.0 * u - 1.0);
        val  += tent * scale;
        norm += scale;
        freq  *= pBase;
        scale /= pBase;
    }
    return (norm > 0.0) ? (val / norm) : 0.0;
}

// Conservative 1-Lipschitz 3D terrain SDF over the 2D regularized Monna surface.
// Each tent mode has derivative |d/dx [scale_k * tent(freq_k * x)]| = 2 * freq_k * scale_k = 2.0,
// so the sum of n modes has exact maximum directional slope L_x = L_z = 2 * n / norm.
float sdMonnaTerrain3D( in vec3 p, in int n, in float pBase, in float amp )
{
    float hx = paMonnaTentMap1D(0.5 * p.x + 0.5, n, pBase);
    float hz = paMonnaTentMap1D(0.5 * p.z + 0.5, n, pBase);
    float h  = amp * (0.5 * (hx + hz) - 0.5);
    float invP = 1.0 / pBase;
    float pNorm = (1.0 - pow(invP, float(max(n, 1)))) / max(1.0 - invP, 1e-6);
    float lipGrad = (pNorm > 0.0) ? (amp * float(n) / (sqrt(2.0) * pNorm)) : 0.0;
    float dSlab = max(abs(p.x) - 0.85, max(abs(p.z) - 0.85, -0.55 - p.y));
    float dSurf = (p.y - h) * inversesqrt(1.0 + lipGrad * lipGrad);
    return max(dSurf, dSlab);
}

#endif // SDF_GEMS_HAS_PA_ULTRAMETRIC_TREE
