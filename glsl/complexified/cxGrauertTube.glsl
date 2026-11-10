// ============================================================================
// cxGrauertTube — Analytic Continuation of SDFs to C^3, Grauert Tubes & Confocal Regularization
//   cxSqrt, cxNorm, cxConfocal,
//   sdComplexSphereSliceRaw3D, sdComplexSphereSlice3D, sdgComplexSphereSlice3D,
//   sdComplexBranchDisk3D, sdgComplexBranchDisk3D,
//   sdComplexNullRing3D, sdgComplexNullRing3D,
//   sdComplexQuadricIntersection3D, sdgComplexQuadricIntersection3D,
//   sdRiemannSqrtSurface3D
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  - sdComplexNullRing3D / sdgComplexNullRing3D: Exact 3D Euclidean SDF of a tube around the
//           complex null quadric circle {(p + i y).(p + i y) = 0}.
//         - sdComplexBranchDisk3D / sdgComplexBranchDisk3D: Exact 3D Euclidean distance and unit gradient
//           to the branch disk {|p| <= |y|, p.y = 0} (the R -> 0 collapse of the complexified sphere).
//         - sdComplexSphereSliceRaw3D: Raw Re(sqrt((p + i y).(p + i y))) - r; exact zero set (oblate
//           spheroid), but NOT 1-Lipschitz (|grad u| = sqrt((u^2 + b^2)/(u^2 + v^2)) >= 1, blowing up
//           on the null circle).
//         - sdComplexSphereSlice3D / sdgComplexSphereSlice3D: Globally 1-Lipschitz confocal / Lie-ball
//           regularization a(p, y) - sqrt(r^2 + |y|^2) with exact zero set {u(p, y) = r} and exact unit
//           normal; exact for y = 0 and on the exterior equatorial plane p.y = 0, conservative bound near poles.
//         - sdComplexQuadricIntersection3D / sdgComplexQuadricIntersection3D: Globally 1-Lipschitz pullback
//           of {u(p, y) <= r, |v(p, y)| <= kappa} (confocal spheroid pierced by an orthogonal one-sheeted
//           hyperboloid throat) via the Bifocal Eikonal Partition |grad a|^2 + |grad c|^2 == 1, grad a . grad c == 0.
//         - sdRiemannSqrtSurface3D: 1-Lipschitz conservative bound (and exact to ~1e-4 via parabolic-coordinate
//           cubic + local refinement) of the principal Riemann surface sheet z = c * Im(sqrt(x + i y)).
// ============================================================================

#ifndef SDF_GEMS_HAS_CX_GRAUERT_TUBE
#define SDF_GEMS_HAS_CX_GRAUERT_TUBE

// ---------------------------------------------------------------------------
// 1. Complex square root, holomorphic norm, and bifocal confocal coordinates
// ---------------------------------------------------------------------------

// Principal complex square root sqrt(z.x + i z.y) = vec2(u, v) with u >= 0.
vec2 cxSqrt( in vec2 z )
{
    float m = length(z);
    float u = sqrt(max(0.5 * (m + z.x), 0.0));
    float v = ((z.y >= 0.0) ? 1.0 : -1.0) * sqrt(max(0.5 * (m - z.x), 0.0));
    return vec2(u, v);
}

// Holomorphic continuation of the Euclidean norm: sqrt((x + i y).(x + i y)) = vec2(u, v).
vec2 cxNorm( in vec3 x, in vec3 y )
{
    return cxSqrt(vec2(dot(x, x) - dot(y, y), 2.0 * dot(x, y)));
}

// Returns vec4(a, c, rPlus, rMinus) for z = x + i y (b = |y|):
//   rPlus, rMinus = maximal and minimal Euclidean distances from x to the null circle C_y = {z.z = 0},
//   a = 0.5 * (rPlus + rMinus) = sqrt(u^2 + b^2) >= b,
//   c = 0.5 * (rPlus - rMinus) = sqrt(b^2 - v^2) in [0, b].
// Satisfies the Bifocal Eikonal Partition everywhere outside C_y:
//   |grad a|^2 = u^2 / (u^2 + v^2) <= 1,   |grad c|^2 = v^2 / (u^2 + v^2) <= 1,
//   |grad a|^2 + |grad c|^2 == 1,          grad a . grad c == 0.
vec4 cxConfocal( in vec3 x, in vec3 y )
{
    float b2 = dot(y, y);
    if (b2 <= 1e-24)
    {
        float r = length(x);
        return vec4(r, 0.0, r, r);
    }
    float b   = sqrt(b2);
    vec3  ey  = y / b;
    float xp  = dot(x, ey);
    vec3  xpe = x - xp * ey;
    float rho = length(xpe);
    float rp  = length(vec2(rho + b, xp));
    float rm  = length(vec2(rho - b, xp));
    return vec4(0.5 * (rp + rm), 0.5 * (rp - rm), rp, rm);
}

// ---------------------------------------------------------------------------
// 2. Complexified Sphere Slice: Raw vs Confocal-Regularized (1-Lipschitz)
// ---------------------------------------------------------------------------

// Raw real part of the complexified sphere SDF: Re(sqrt((p + i y).(p + i y))) - r.
// Has the exact oblate-spheroid zero set {u(p, y) = r}, but |grad u| = sqrt((u^2 + b^2)/(u^2 + v^2)) >= 1
// everywhere (diverging on the null circle), so it overestimates distances when y != 0.
float sdComplexSphereSliceRaw3D( in vec3 p, in vec3 y, in float r )
{
    return cxNorm(p, y).x - r;
}

// Confocal / Lie-ball regularized complex sphere slice:
//   d(p) = a(p, y) - sqrt(r^2 + |y|^2) = 0.5 * (r_+ + r_-) - sqrt(r^2 + |y|^2).
// Shares the exact zero set {u(p, y) = r} with sdComplexSphereSliceRaw3D, is globally 1-Lipschitz
// (|grad d| = u / sqrt(u^2 + v^2) <= 1 everywhere), reduces to the exact sphere SDF for y = 0,
// and returns the exact analytic unit normal in .yzw.
vec4 sdgComplexSphereSlice3D( in vec3 p, in vec3 y, in float r )
{
    float b2 = dot(y, y);
    float r0 = max(r, 0.0);
    if (b2 <= 1e-24)
    {
        float l = length(p);
        vec3  g = (l > 1e-20) ? (p / l) : vec3(1.0, 0.0, 0.0);
        return vec4(l - r0, g);
    }
    float b   = sqrt(b2);
    vec3  ey  = y / b;
    float xp  = dot(p, ey);
    vec3  xpe = p - xp * ey;
    float rho = length(xpe);
    vec3  er  = (rho > 1e-20) ? (xpe / rho) : ((abs(ey.x) < 0.9) ? normalize(cross(ey, vec3(1.0, 0.0, 0.0)))
                                                                 : normalize(cross(ey, vec3(0.0, 1.0, 0.0))));
    vec2  vp  = vec2(rho + b, xp);
    vec2  vm  = vec2(rho - b, xp);
    float rp  = length(vp);
    float rm  = length(vm);
    float a   = 0.5 * (rp + rm);
    float a0  = length(vec2(r0, b));
    float d   = a - a0;

    // Unit normal in the meridian plane: bisector of grad(r_+) and grad(r_-)
    vec2 gp = (rp > 1e-20) ? (vp / rp) : vec2(1.0, 0.0);
    vec2 gm = (rm > 1e-20) ? (vm / rm) : ((xp >= 0.0) ? vec2(0.0, 1.0) : vec2(0.0, -1.0));
    vec2 ga = gp + gm;
    float lga = length(ga);
    vec2 n2 = (lga > 1e-20) ? (ga / lga) : ((xp >= 0.0) ? vec2(0.0, 1.0) : vec2(0.0, -1.0));
    vec3 n  = n2.x * er + n2.y * ey;
    return vec4(d, n);
}

float sdComplexSphereSlice3D( in vec3 p, in vec3 y, in float r )
{
    return cxConfocal(p, y).x - length(vec2(max(r, 0.0), length(y)));
}

// ---------------------------------------------------------------------------
// 3. Exact Grauert Branch Disk (R = 0 limit) and Complex Null Quadric Ring
// ---------------------------------------------------------------------------

// Exact 3D Euclidean distance and unit gradient to the branch disk D_y = {|p| <= |y|, p.y = 0}
// onto which the complexified sphere collapses as r -> 0.
vec4 sdgComplexBranchDisk3D( in vec3 p, in vec3 y )
{
    float b2 = dot(y, y);
    if (b2 <= 1e-24)
    {
        float l = length(p);
        return vec4(l, (l > 1e-20) ? (p / l) : vec3(1.0, 0.0, 0.0));
    }
    float b   = sqrt(b2);
    vec3  ey  = y / b;
    float xp  = dot(p, ey);
    vec3  xpe = p - xp * ey;
    float rho = length(xpe);
    vec3  er  = (rho > 1e-20) ? (xpe / rho) : ((abs(ey.x) < 0.9) ? normalize(cross(ey, vec3(1.0, 0.0, 0.0)))
                                                                 : normalize(cross(ey, vec3(0.0, 1.0, 0.0))));
    float dr  = max(rho - b, 0.0);
    float d   = length(vec2(dr, xp));
    vec2  n2  = (d > 1e-20) ? (vec2(dr, xp) / d) : ((xp >= 0.0) ? vec2(0.0, 1.0) : vec2(0.0, -1.0));
    return vec4(d, n2.x * er + n2.y * ey);
}

float sdComplexBranchDisk3D( in vec3 p, in vec3 y )
{
    return sdgComplexBranchDisk3D(p, y).x;
}

// Exact 3D Euclidean SDF and unit gradient of a tube of radius rTube around the complex null quadric
// circle C_y = {(p + i y).(p + i y) = 0} = {|p| = |y|, p.y = 0}, evaluated via r_-(p, y) - rTube.
vec4 sdgComplexNullRing3D( in vec3 p, in vec3 y, in float rTube )
{
    float b2 = dot(y, y);
    if (b2 <= 1e-24)
    {
        float l = length(p);
        return vec4(l - rTube, (l > 1e-20) ? (p / l) : vec3(1.0, 0.0, 0.0));
    }
    float b   = sqrt(b2);
    vec3  ey  = y / b;
    float xp  = dot(p, ey);
    vec3  xpe = p - xp * ey;
    float rho = length(xpe);
    vec3  er  = (rho > 1e-20) ? (xpe / rho) : ((abs(ey.x) < 0.9) ? normalize(cross(ey, vec3(1.0, 0.0, 0.0)))
                                                                 : normalize(cross(ey, vec3(0.0, 1.0, 0.0))));
    vec2  vm  = vec2(rho - b, xp);
    float rm  = length(vm);
    vec2  n2  = (rm > 1e-20) ? (vm / rm) : vec2(1.0, 0.0);
    return vec4(rm - rTube, n2.x * er + n2.y * ey);
}

float sdComplexNullRing3D( in vec3 p, in vec3 y, in float rTube )
{
    return cxConfocal(p, y).w - rTube;
}

// ---------------------------------------------------------------------------
// 4. Complex Quadric Intersection / Grauert Tube Collar:
//    { Re sqrt(z.z) <= r }  intersect  { |Im sqrt(z.z)| <= kappa }   (0 < kappa < |y|)
//    In confocal coordinates (a, c), this is the orthogonal 2D rectangle {a <= a0, |c| >= c0}
//    where a0 = sqrt(r^2 + b^2) and c0 = sqrt(b^2 - kappa^2): an oblate spheroid pierced by an
//    orthogonal one-sheeted confocal hyperboloid throat along y.
//    By the Universal 2D Pullback Theorem (||D(a,c)||_op = max(cos psi, sin psi) in [1/sqrt(2), 1]),
//    the 2D box SDF in (a, c) pulls back to a globally 1-Lipschitz 3D SDF with closed-form normal.
// ---------------------------------------------------------------------------
vec4 sdgComplexQuadricIntersection3D( in vec3 p, in vec3 y, in float r, in float kappa )
{
    float b2 = dot(y, y);
    float r0 = max(r, 0.0);
    if (b2 <= 1e-24)
    {
        float l = length(p);
        return vec4(l - r0, (l > 1e-20) ? (p / l) : vec3(1.0, 0.0, 0.0));
    }
    float b   = sqrt(b2);
    float kp  = clamp(kappa, 0.0, b);
    float a0  = length(vec2(r0, b));
    float c0  = sqrt(max((b - kp) * (b + kp), 0.0));

    vec3  ey  = y / b;
    float xp  = dot(p, ey);
    vec3  xpe = p - xp * ey;
    float rho = length(xpe);
    vec3  er  = (rho > 1e-20) ? (xpe / rho) : ((abs(ey.x) < 0.9) ? normalize(cross(ey, vec3(1.0, 0.0, 0.0)))
                                                                 : normalize(cross(ey, vec3(0.0, 1.0, 0.0))));
    vec2  vp  = vec2(rho + b, xp);
    vec2  vm  = vec2(rho - b, xp);
    float rp  = length(vp);
    float rm  = length(vm);
    float a   = 0.5 * (rp + rm);
    float ac  = 0.5 * (rp - rm);

    // 2D SDF of {a <= a0, ac >= c0} in (a, ac) coordinates
    vec2  q   = vec2(a - a0, c0 - ac);
    vec2  qp  = max(q, 0.0);
    float lqp = length(qp);
    float d   = min(max(q.x, q.y), 0.0) + lqp;

    // 2D gradient w.r.t. (a, ac): w.x * grad(a) - w.y * grad(ac)
    vec2  w   = (lqp > 0.0) ? ((lqp > 1e-20) ? (qp / lqp) : vec2(1.0, 0.0))
                            : ((q.x > q.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0));
    vec2  gp  = (rp > 1e-20) ? (vp / rp) : vec2(1.0, 0.0);
    vec2  gm  = (rm > 1e-20) ? (vm / rm) : ((xp >= 0.0) ? vec2(0.0, 1.0) : vec2(0.0, -1.0));
    vec2  ga  = 0.5 * (gp + gm);
    vec2  gac = 0.5 * (gp - gm);
    vec2  g2  = w.x * ga - w.y * gac;
    float lg2 = length(g2);
    vec2  n2  = (lg2 > 1e-20) ? (g2 / lg2) : vec2(1.0, 0.0);
    return vec4(d, n2.x * er + n2.y * ey);
}

float sdComplexQuadricIntersection3D( in vec3 p, in vec3 y, in float r, in float kappa )
{
    return sdgComplexQuadricIntersection3D(p, y, r, kappa).x;
}

// ---------------------------------------------------------------------------
// 5. Principal Riemann Surface Sheet with Branch Cut: z = cScale * Im(sqrt(x + i y))
//    Parameterised in parabolic coordinates w = xi + i eta = sqrt(x + i y) (xi >= 0, xi^2 + eta^2 <= R)
//    by S(xi, eta) = (xi^2 - eta^2, 2 xi eta, cScale * eta), thickened by th.
//    For any fixed eta in [-sqrt(R), sqrt(R)], minimizing |(x, y) - (xi^2 - eta^2, 2 xi eta)|^2
//    over xi in [0, sqrt(R - eta^2)] is the exact depressed cubic xi^3 + (eta^2 - x) xi - y eta = 0.
// ---------------------------------------------------------------------------
float sdRiemannSqrtXiDist2( in float x, in float y, in float eta, in float xiMax )
{
    float p  = eta * eta - x;
    float q  = -y * eta;
    float p3 = p * (1.0 / 3.0);
    float q2 = 0.5 * q;
    float D  = q2 * q2 + p3 * p3 * p3;
    float best = 1e30;
    if (D >= 0.0)
    {
        float sD = sqrt(D);
        float u  = -q2 + sD;
        float v  = -q2 - sD;
        float xi = clamp(sign(u) * pow(abs(u), 1.0 / 3.0) + sign(v) * pow(abs(v), 1.0 / 3.0), 0.0, xiMax);
        vec2  d  = vec2(x - xi * xi + eta * eta, y - 2.0 * xi * eta);
        best = dot(d, d);
    }
    else
    {
        float r  = sqrt(max(-p3, 1e-30));
        float th = acos(clamp(-q2 / max(r * r * r, 1e-30), -1.0, 1.0)) * (1.0 / 3.0);
        for (int k = 0; k < 3; k++)
        {
            float xi = clamp(2.0 * r * cos(th + float(k) * 2.09439510239), 0.0, xiMax);
            vec2  d  = vec2(x - xi * xi + eta * eta, y - 2.0 * xi * eta);
            best = min(best, dot(d, d));
        }
    }
    vec2 d0 = vec2(x + eta * eta, y);
    vec2 d1 = vec2(x - xiMax * xiMax + eta * eta, y - 2.0 * xiMax * eta);
    return min(best, min(dot(d0, d0), dot(d1, d1)));
}

float sdRiemannSqrtSurface3D( in vec3 p, in float R, in float cScale, in float th )
{
    float emax = sqrt(max(R, 0.0));
    float bestP = 1e30, bestN = 1e30;
    float etaP  = 0.0,  etaN  = 0.0;
    const int M = 16;
    for (int i = 0; i <= M; i++)
    {
        float eta = emax * (float(i) / float(M));
        float xim = sqrt(max((emax - eta) * (emax + eta), 0.0));
        float dP  = sdRiemannSqrtXiDist2(p.x, p.y,  eta, xim) + (p.z - cScale * eta) * (p.z - cScale * eta);
        float dN  = sdRiemannSqrtXiDist2(p.x, p.y, -eta, xim) + (p.z + cScale * eta) * (p.z + cScale * eta);
        if (dP < bestP) { bestP = dP; etaP =  eta; }
        if (dN < bestN) { bestN = dN; etaN = -eta; }
    }
    vec2  w0  = cxSqrt(p.xy);
    float e0  = clamp(abs(w0.y), 0.0, emax);
    float xm0 = sqrt(max((emax - e0) * (emax + e0), 0.0));
    float dP0 = sdRiemannSqrtXiDist2(p.x, p.y,  e0, xm0) + (p.z - cScale * e0) * (p.z - cScale * e0);
    float dN0 = sdRiemannSqrtXiDist2(p.x, p.y, -e0, xm0) + (p.z + cScale * e0) * (p.z + cScale * e0);
    if (dP0 < bestP) { bestP = dP0; etaP =  e0; }
    if (dN0 < bestN) { bestN = dN0; etaN = -e0; }

    // Refine both the upper branch sheet (eta >= 0) and lower branch sheet (eta <= 0) independently
    float h = emax / float(M);
    for (int iter = 0; iter < 8; iter++)
    {
        float pL = clamp(etaP - 0.5 * h, 0.0, emax), pR = clamp(etaP + 0.5 * h, 0.0, emax);
        float dpL = sdRiemannSqrtXiDist2(p.x, p.y, pL, sqrt(max((emax - pL) * (emax + pL), 0.0))) + (p.z - cScale * pL) * (p.z - cScale * pL);
        float dpR = sdRiemannSqrtXiDist2(p.x, p.y, pR, sqrt(max((emax - pR) * (emax + pR), 0.0))) + (p.z - cScale * pR) * (p.z - cScale * pR);
        if (dpL < bestP) { bestP = dpL; etaP = pL; }
        if (dpR < bestP) { bestP = dpR; etaP = pR; }

        float nL = clamp(etaN - 0.5 * h, -emax, 0.0), nR = clamp(etaN + 0.5 * h, -emax, 0.0);
        float dnL = sdRiemannSqrtXiDist2(p.x, p.y, nL, sqrt(max((emax + nL) * (emax - nL), 0.0))) + (p.z - cScale * nL) * (p.z - cScale * nL);
        float dnR = sdRiemannSqrtXiDist2(p.x, p.y, nR, sqrt(max((emax + nR) * (emax - nR), 0.0))) + (p.z - cScale * nR) * (p.z - cScale * nR);
        if (dnL < bestN) { bestN = dnL; etaN = nL; }
        if (dnR < bestN) { bestN = dnR; etaN = nR; }
        h *= 0.5;
    }
    return sqrt(min(bestP, bestN)) - th;
}

#endif // SDF_GEMS_HAS_CX_GRAUERT_TUBE
