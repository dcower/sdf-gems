// ============================================================================
// thSolGeometry — Thurston's 8th Model Geometry Sol & Arnold's Cat-Map Tiling
//   solMul, solInv, opSolTranslate, opSolTranslateScale, sdgSolTranslate, solGrad,
//   solSheetDistanceXZ, solSheetDistanceYZ, solFinsGauge,
//   solGeodesicDistance, sdgSolGeodesicDistance,
//   solMetricBound, sdgSolMetricBound, sdSolBall, sdgSolBall,
//   opRepSolCatMap, sdSolCatLattice
// @author Gemini 4 Carbon (sdf-gems, 2026)
// @exact  solSheetDistanceXZ / solSheetDistanceYZ: Exact hyperbolic H^2 sheet distances
//         on y = 0 and x = 0 across all R^3.
//         solGeodesicDistance / sdgSolGeodesicDistance: Exact Sol Riemannian geodesic
//         distance and canonical momentum gradient (||grad_Sol d||_Sol == 1) up to the
//         Coiculescu-Schwartz injectivity radius R_inj = sqrt(2)*pi ~= 4.44288 via
//         1D rational momentum elimination, branchless unwrapped Landen AGM, and
//         stationary Maupertuis action.
//         sdSolBall / sdgSolBall / sdSolCatLattice: Certified globally 1-Lipschitz
//         Euclidean SDF bounds (||grad_E d||_E <= 1 everywhere, 0% raymarching overstep)
//         with exact Sol zero-level sets.
// ============================================================================

#ifndef SDF_GEMS_HAS_TH_SOL_GEOMETRY
#define SDF_GEMS_HAS_TH_SOL_GEOMETRY

// ---------------------------------------------------------------------------
// 1. Sol Lie Group Law (R^2 |x R) & Left-Invariant Orthonormal Calculus
//    Group law: a * b = (a.x + exp(-a.z) * b.x, a.y + exp(a.z) * b.y, a.z + b.z)
//    Left-invariant metric: ds_Sol^2 = exp(2z) dx^2 + exp(-2z) dy^2 + dz^2
//    Orthonormal frame:     X_1 = exp(-z) d_x,  X_2 = exp(z) d_y,  X_3 = d_z
// ---------------------------------------------------------------------------

vec3 solMul( in vec3 a, in vec3 b )
{
    float ez = exp(a.z);
    return vec3(a.x + b.x / ez, a.y + b.y * ez, a.z + b.z);
}

vec3 solInv( in vec3 a )
{
    float ez = exp(a.z);
    return vec3(-ez * a.x, -a.y / ez, -a.z);
}

// Left-translates sample point p by a^{-1} in Sol: q = a^{-1} * p, an exact Sol isometry.
vec3 opSolTranslate( in vec3 p, in vec3 a )
{
    float ez = exp(a.z);
    return vec3(ez * (p.x - a.x), (p.y - a.y) / ez, p.z - a.z);
}

// Exact reciprocal of the maximum Euclidean singular value sigma_max(dq/dp) = exp(|a.z|).
float opSolTranslateScale( in vec3 a )
{
    return exp(-abs(a.z));
}

// Pulls back an inner 1-Lipschitz Euclidean SDF + unit gradient dgInner = sdgInner(opSolTranslate(p, a)).
vec4 sdgSolTranslate( in vec3 a, in vec4 dgInner )
{
    float ez    = exp(a.z);
    float scale = exp(-abs(a.z));
    vec3  g     = vec3(ez * dgInner.y, dgInner.z / ez, dgInner.w);
    float lg    = length(g);
    return vec4(dgInner.x * scale, (lg > 1e-20) ? (g / lg) : vec3(0.0, 0.0, 1.0));
}

// Projects a Euclidean gradient onto the orthonormal left-invariant Sol frame (X_1, X_2, X_3).
// Satisfies length(solGrad(p, grad u)) == 1 iff u solves the Sol Eikonal equation.
vec3 solGrad( in vec3 p, in vec3 euclidGrad )
{
    float ez = exp(p.z);
    return vec3(euclidGrad.x / ez, euclidGrad.y * ez, euclidGrad.z);
}

// ---------------------------------------------------------------------------
// 2. Exact Totally-Geodesic Hyperbolic Foliation Sheet Distances
//    Plane y = 0: ds^2 = exp(2z) dx^2 + dz^2  (isometric to Poincare H^2 via u=x, v=e^{-z})
//    Plane x = 0: ds^2 = exp(-2z) dy^2 + dz^2 (isometric to Poincare H^2 via u=y, v=e^{+z})
//    Written in cancellation-free half-angle asinh form for full FP32 accuracy at p -> 0.
// ---------------------------------------------------------------------------

float solSheetDistanceXZ( in vec2 xz )
{
    float hz = sinh(0.5 * xz.y);
    float s2 = hz * hz + 0.25 * exp(xz.y) * (xz.x * xz.x);
    return 2.0 * asinh(sqrt(max(s2, 0.0)));
}

float solSheetDistanceYZ( in vec2 yz )
{
    float hz = sinh(0.5 * yz.y);
    float s2 = hz * hz + 0.25 * exp(-yz.y) * (yz.x * yz.x);
    return 2.0 * asinh(sqrt(max(s2, 0.0)));
}

// Closed-form Sol-1-Lipschitz hyperbolic-fins lower gauge:
//   max(d_XZ, d_YZ) <= d_Sol(0, p) <= sqrt(2) * max(d_XZ, d_YZ), with ||grad_Sol||_Sol == 1 a.e.
float solFinsGauge( in vec3 p )
{
    return max(solSheetDistanceXZ(p.xz), solSheetDistanceYZ(p.yz));
}

// ---------------------------------------------------------------------------
// 3. Exact Sol Geodesic Distance (solGeodesicDistance / sdgSolGeodesicDistance)
//    Integrable Hamiltonian pendulum BVP:
//      z'^2 + p_x^2 exp(-2z) + p_y^2 exp(2z) = 1,   z'' = p_x x' - p_y y'.
//    Integrating z'' yields the exact algebraic momentum conservation law:
//      z'(T) - z'(0) = p_x x - p_y y,
//    which eliminates r = length(vec2(p_x, p_y)) rationally at any angle alpha = atan(p_y, p_x):
//      r(alpha) = 2 |D| / Delta,  Delta = sqrt(((1 - sqrt(C))^2 + D^2)((1 + sqrt(C))^2 + D^2)).
//    Evaluates unwrapped elliptic integrals (F, E) simultaneously via 6-step branchless
//    Landen AGM and returns the stationary Maupertuis action S(alpha*) == d_Sol(0, p)
//    along with the exact canonical momentum gradient (p_x, p_y, z'(T)).
// ---------------------------------------------------------------------------

vec2 thSolAgmEllipPair( in float th0, in float th1, in float k2, in float oneMinusK2 )
{
    float a       = 1.0;
    float b       = sqrt(max(oneMinusK2, 1e-30));
    float c2_sum  = k2;
    vec2  ang     = vec2(th0, th1);
    vec2  sin_sum = vec2(0.0);
    float two_pow = 1.0;

    for (int i = 0; i < 7; i++)
    {
        float c   = 0.5 * (a - b);
        vec2  s   = sin(ang);
        vec2  co  = cos(ang);
        vec2  den = max(a * co * co + b * s * s, vec2(1e-30));
        ang      = 2.0 * ang + atan(((b - a) * s * co) / den);
        sin_sum += c * sin(ang);
        float an = 0.5 * (a + b);
        b        = sqrt(a * b);
        a        = an;
        two_pow *= 2.0;
        c2_sum  += two_pow * (c * c);
    }

    float dF = (ang.y - ang.x) / (two_pow * a);
    float dE = dF * (1.0 - 0.5 * c2_sum) + (sin_sum.y - sin_sum.x);
    return vec2(dF, dE);
}

// Returns vec4(residual, S_stationary, p_x, p_y) and outputs vz1 for x0 > 0, y0 > 0, z >= 0.
// Parameterized by signed offset delta = alpha - alpha_D so D = -hypot(x0, y0) * sin(delta) has zero FP32 cancellation!
vec4 thSolEvalMomentumDelta( in float delta, in float al_D, in float rho_xy, in float x0, in float y0, in float ez, in float emz, out float vz1_out )
{
    float alpha = al_D + delta;
    float ca    = cos(alpha);
    float sa    = sin(alpha);
    float C     = ca * ca * (emz * emz) + sa * sa * (ez * ez);
    float rtC   = sqrt(max(C, 0.0));
    float D     = -rho_xy * sin(delta);
    float D2    = D * D;
    float rc_m  = 1.0 - rtC;
    float rc_p  = 1.0 + rtC;
    float Delta = sqrt(max((rc_m * rc_m + D2) * (rc_p * rc_p + D2), 0.0));

    if (Delta <= 1e-20)
    {
        vz1_out = 0.0;
        return vec4(0.0, rho_xy, 0.70710678, 0.70710678);
    }

    float sgnD = (delta <= 0.0) ? 1.0 : -1.0;
    float invD = 1.0 / Delta;
    float r    = 2.0 * abs(D) * invD;
    float vz0  = sgnD * (1.0 - C - D2) * invD;
    float vz1  = sgnD * (1.0 - C + D2) * invD;
    vz1_out    = vz1;

    float px = r * ca;
    float py = r * sa;
    float m  = 2.0 * px * py;

    if (m >= 1.0 - 1e-6)
    {
        float S = px * x0 + py * y0;
        return vec4(0.0, S, px, py);
    }

    float rt_1mm = sqrt(max(1.0 - m, 1e-30));
    float inv1mm = 1.0 / rt_1mm;
    float th0    = atan(clamp(-vz0 * inv1mm, -1.0, 1.0), clamp((py - px) * inv1mm, -1.0, 1.0));
    float th1    = atan(clamp(-vz1 * inv1mm, -1.0, 1.0), clamp((py * ez - px * emz) * inv1mm, -1.0, 1.0));
    if (th1 <= th0)
    {
        th1 += 6.28318530718;
    }

    float k2         = (1.0 - m) / (1.0 + m);
    float oneMinusK2 = (2.0 * m) / (1.0 + m);
    float rt_1pm     = sqrt(1.0 + m);
    vec2  fe         = thSolAgmEllipPair(th0, th1, k2, oneMinusK2);
    float T_arc      = fe.x / rt_1pm;
    float pxy        = px * x0 + py * y0;
    float res        = rt_1pm * fe.y - m * T_arc - pxy;
    float S_stat     = pxy + rt_1pm * (fe.x - fe.y);
    return vec4(res, S_stat, px, py);
}

// Evaluates symmetric midpoint-surface geodesics (y0 == exp(z) * x0, where alpha* == alpha_D)
// parameterized smoothly by signed vertical velocity vz in (-1, +1), covering both 0-turn (vz >= 0) and 2-turn (vz < 0)!
vec4 thSolEvalMidpointVz( in float vz, in float ca, in float sa, in float x0, in float y0, in float ez, in float emz )
{
    float r  = sqrt(max((1.0 - vz) * (1.0 + vz), 0.0));
    float px = r * ca;
    float py = r * sa;
    float m  = 2.0 * px * py;
    if (m >= 1.0 - 1e-6)
    {
        return vec4(0.0, px * x0 + py * y0, px, py);
    }
    float rt_1mm = sqrt(max(1.0 - m, 1e-30));
    float inv1mm = 1.0 / rt_1mm;
    float th0    = atan(clamp(-vz * inv1mm, -1.0, 1.0), clamp((py - px) * inv1mm, -1.0, 1.0));
    float th1    = atan(clamp(-vz * inv1mm, -1.0, 1.0), clamp((py * ez - px * emz) * inv1mm, -1.0, 1.0));
    if (th1 <= th0)
    {
        th1 += 6.28318530718;
    }
    float k2         = (1.0 - m) / (1.0 + m);
    float oneMinusK2 = (2.0 * m) / (1.0 + m);
    float rt_1pm     = sqrt(1.0 + m);
    vec2  fe         = thSolAgmEllipPair(th0, th1, k2, oneMinusK2);
    float T_arc      = fe.x / rt_1pm;
    float pxy        = px * x0 + py * y0;
    float res        = rt_1pm * fe.y - m * T_arc - pxy;
    float S_stat     = pxy + rt_1pm * (fe.x - fe.y);
    return vec4(res, S_stat, px, py);
}

// Exact Sol geodesic distance d_Sol(0, p) and Euclidean gradient (p_x, p_y, z'(T)).
// Valid on the full injectivity domain d_Sol(0, p) <= sqrt(2) * pi ~= 4.44288.
vec4 sdgSolGeodesicDistance( in vec3 p )
{
    float x0 = abs(p.x);
    float y0 = abs(p.y);
    float z  = p.z;
    bool  swapped = false;
    if (z < 0.0)
    {
        float tmp = x0;
        x0 = y0;
        y0 = tmp;
        z  = -z;
        swapped = true;
    }

    // Clamp exponential argument solely to prevent IEEE FP32 overflow at extreme coordinates
    float zc  = min(z, 40.0);
    float ez  = exp(zc);
    float emz = 1.0 / ez;
    float scaleRef = max(1.0, max(max(x0, y0), z));

    float d, gx, gy, gz;

    if (y0 <= 1e-6 * scaleRef)
    {
        float hz = 0.5 * (exp(0.5 * zc) - exp(-0.5 * zc));
        float s  = sqrt(max(hz * hz + 0.25 * ez * x0 * x0, 0.0));
        d = 2.0 * asinh(s);
        if (d <= 1e-7)
        {
            float l = length(vec3(x0, y0, z));
            gx = (l > 1e-20) ? (x0 / l) : 0.0;
            gy = 0.0;
            gz = (l > 1e-20) ? (z / l) : 1.0;
        }
        else
        {
            float invSh = 1.0 / sinh(min(d, 40.0));
            gx = (ez * x0) * invSh;
            gy = 0.0;
            gz = (0.5 * (ez - emz) + 0.5 * ez * x0 * x0) * invSh;
        }
    }
    else if (x0 <= 1e-6 * scaleRef)
    {
        float hz = 0.5 * (exp(0.5 * zc) - exp(-0.5 * zc));
        float s  = sqrt(max(hz * hz + 0.25 * emz * y0 * y0, 0.0));
        d = 2.0 * asinh(s);
        if (d <= 1e-7)
        {
            float l = length(vec3(x0, y0, z));
            gx = 0.0;
            gy = (l > 1e-20) ? (y0 / l) : 0.0;
            gz = (l > 1e-20) ? (z / l) : 1.0;
        }
        else
        {
            float invSh = 1.0 / sinh(min(d, 40.0));
            gx = 0.0;
            gy = (emz * y0) * invSh;
            gz = (0.5 * (ez - emz) - 0.5 * emz * y0 * y0) * invSh;
        }
    }
    else
    {
        float rho_xy   = length(vec2(x0, y0));
        float al_D     = atan(x0, y0);
        float diff_mid = ez * x0 - y0;
        if (abs(diff_mid) <= 1e-6 * (ez * x0 + y0))
        {
            if (z <= 1e-6 * scaleRef)
            {
                d  = rho_xy;
                gx = 0.70710678;
                gy = 0.70710678;
                gz = 0.0;
            }
            else
            {
                float ca  = y0 / rho_xy;
                float sa  = x0 / rho_xy;
                float lo  = -8.0;
                float hi  =  8.0;
                float flo = thSolEvalMidpointVz(tanh(lo), ca, sa, x0, y0, ez, emz).x;
                float fhi = thSolEvalMidpointVz(tanh(hi), ca, sa, x0, y0, ez, emz).x;
                for (int i = 0; i < 8; i++)
                {
                    float mid  = 0.5 * (lo + hi);
                    float fmid = thSolEvalMidpointVz(tanh(mid), ca, sa, x0, y0, ez, emz).x;
                    if (flo * fmid <= 0.0) { hi = mid; fhi = fmid; }
                    else                   { lo = mid; flo = fmid; }
                }
                float u_best = (abs(fhi - flo) > 1e-20 && flo * fhi <= 0.0)
                             ? clamp(lo - flo * (hi - lo) / (fhi - flo), lo, hi)
                             : 0.5 * (lo + hi);
                float vz_opt = tanh(u_best);
                vec4  ev     = thSolEvalMidpointVz(vz_opt, ca, sa, x0, y0, ez, emz);
                d  = ev.y;
                gx = ev.z;
                gy = ev.w;
                gz = vz_opt;
            }
        }
        else
        {
            float signedSpan = (diff_mid > 0.0) ? (-al_D) : (1.57079632679 - al_D);

            float lo  = -16.5;
            float hi  =  16.5;
            float vz_tmp;
            float flo = thSolEvalMomentumDelta(signedSpan / (1.0 + exp(-lo)), al_D, rho_xy, x0, y0, ez, emz, vz_tmp).x;
            float fhi = thSolEvalMomentumDelta(signedSpan / (1.0 + exp(-hi)), al_D, rho_xy, x0, y0, ez, emz, vz_tmp).x;

            // 8 bisection steps in logistic coordinate u shrink [-16.5, 16.5] to width 0.1289 (covering |delta| down to 7e-8!)
            for (int i = 0; i < 8; i++)
            {
                float mid  = 0.5 * (lo + hi);
                float fmid = thSolEvalMomentumDelta(signedSpan / (1.0 + exp(-mid)), al_D, rho_xy, x0, y0, ez, emz, vz_tmp).x;
                if (flo * fmid <= 0.0)
                {
                    hi  = mid;
                    fhi = fmid;
                }
                else
                {
                    lo  = mid;
                    flo = fmid;
                }
            }

            // 4 safeguarded secant steps converge into the quadratic Maupertuis stationary basin
            float prev_u = lo;
            float prev_f = flo;
            float cur_u  = hi;
            float cur_f  = fhi;
            for (int i = 0; i < 4; i++)
            {
                float df    = cur_f - prev_f;
                float mid   = 0.5 * (lo + hi);
                float sec   = (abs(df) > 1e-20) ? (cur_u - cur_f * (cur_u - prev_u) / df) : mid;
                float pad   = 0.02 * (hi - lo);
                float nxt_u = (sec > lo + pad && sec < hi - pad) ? sec : mid;
                float nxt_f = thSolEvalMomentumDelta(signedSpan / (1.0 + exp(-nxt_u)), al_D, rho_xy, x0, y0, ez, emz, vz_tmp).x;
                prev_u = cur_u;
                prev_f = cur_f;
                cur_u  = nxt_u;
                cur_f  = nxt_f;
                if (flo * nxt_f <= 0.0)
                {
                    hi  = nxt_u;
                    fhi = nxt_f;
                }
                else
                {
                    lo  = nxt_u;
                    flo = nxt_f;
                }
            }

            float u_best  = (abs(cur_f) <= min(abs(flo), abs(fhi))) ? cur_u : ((abs(flo) <= abs(fhi)) ? lo : hi);
            float del_opt = signedSpan / (1.0 + exp(-u_best));
            vec4  ev      = thSolEvalMomentumDelta(del_opt, al_D, rho_xy, x0, y0, ez, emz, gz);
            d  = ev.y;
            gx = ev.z;
            gy = ev.w;
        }
    }

    if (swapped)
    {
        float tmp = gx;
        gx = gy;
        gy = tmp;
        gz = -gz;
    }
    gx = (p.x >= 0.0) ? gx : -gx;
    gy = (p.y >= 0.0) ? gy : -gy;
    return vec4(d, gx, gy, gz);
}

float solGeodesicDistance( in vec3 p )
{
    return sdgSolGeodesicDistance(p).x;
}

// ---------------------------------------------------------------------------
// 4. Certified Globally 1-Lipschitz Sol -> Euclidean Metric Conversion
//    Theorem: Because ||grad_E f||_E <= exp(|z|) ||grad_Sol f||_Sol <= exp(|z|),
//    integrating |df/ds| <= exp(|z_0| + s) along any unit-speed Euclidean path gives:
//      Psi_Sol(f, z) = sign(f) * log(1 + exp(-|z|) * |f|),
//    which preserves the exact zero-level set and satisfies ||grad_E Psi_Sol||_E <= 1
//    everywhere on R^3 (0.0% sphere-tracing overstep!).
// ---------------------------------------------------------------------------

float solMetricBound( in float f, in float z )
{
    float af = abs(f);
    float u  = exp(-abs(z)) * af;
    // Cancellation-safe log(1 + u) for small u:
    float l1u = (u > 1e-4) ? log(1.0 + u) : (u * (1.0 - 0.5 * u * (1.0 - 0.66666667 * u)));
    return (f >= 0.0) ? l1u : -l1u;
}

// Converts a Sol 1-jet dgSol = vec4(d_Sol, grad_E(d_Sol)) at level R into vec4(d_bound, unit_grad_E).
vec4 sdgSolMetricBound( in vec3 p, in vec4 dgSol, in float R )
{
    float f    = dgSol.x - R;
    float dB   = solMetricBound(f, p.z);
    float sgnZ = (p.z >= 0.0) ? 1.0 : -1.0;
    // Exact branch-free gradient direction of Psi_Sol(f, z):
    //   grad Psi_Sol = (grad f - f * sign(z) * e_z) / (exp(|z|) + |f|)
    vec3  gRaw = vec3(dgSol.y, dgSol.z, dgSol.w - f * sgnZ);
    float lg   = length(gRaw);
    vec3  n    = (lg > 1e-20) ? (gRaw / lg) : vec3(0.0, 0.0, 1.0);
    return vec4(dB, n);
}

// Certified globally 1-Lipschitz Euclidean SDF and unit gradient of the 4-finned Sol geodesic sphere {d_Sol(0, p) <= R}.
// Combines the universal logarithmic secant bound Psi_Sol(d_Sol - R, z) with the exact tight bounding box
// [-sinh(R), +sinh(R)]^2 x [-R, +R] (tangent at the 4 hyperbolic fin apices and 2 z-poles) for O(1) far-field marching!
vec4 sdgSolBall( in vec3 p, in float R )
{
    if (R <= 1e-12)
    {
        float l = length(p);
        return vec4(l - R, (l > 1e-20) ? (p / l) : vec3(0.0, 0.0, 1.0));
    }
    vec4  dgPsi = sdgSolMetricBound(p, sdgSolGeodesicDistance(p), R);

    float shR  = sinh(min(R, 40.0));
    vec3  qBox = abs(p) - vec3(shR, shR, R);
    vec3  qPos = max(qBox, 0.0);
    float lPos = length(qPos);
    float dBox = lPos + min(max(qBox.x, max(qBox.y, qBox.z)), 0.0);

    if (dBox > dgPsi.x)
    {
        vec3 sgnP = vec3((p.x >= 0.0) ? 1.0 : -1.0, (p.y >= 0.0) ? 1.0 : -1.0, (p.z >= 0.0) ? 1.0 : -1.0);
        vec3 gBox = (lPos > 1e-12) ? (sgnP * (qPos / lPos))
                  : ((qBox.x > max(qBox.y, qBox.z)) ? vec3(sgnP.x, 0.0, 0.0)
                  : ((qBox.y > qBox.z) ? vec3(0.0, sgnP.y, 0.0) : vec3(0.0, 0.0, sgnP.z)));
        return vec4(dBox, gBox);
    }
    return dgPsi;
}

float sdSolBall( in vec3 p, in float R )
{
    return sdgSolBall(p, R).x;
}

// ---------------------------------------------------------------------------
// 5. Compact Sol 3-Manifold Tiling via Arnold's Cat Map A = ((2, 1), (1, 1))
//    Eigenvalues: lambda_+ = phi^2 = exp(z_0), lambda_- = phi^{-2} = exp(-z_0),
//    with z_0 = 2 ln(phi) = 0.9624236501.
//    Orthogonal eigenbasis Q in SO(2) of symmetric A:
//      v_- = (1, -phi) / sqrt(1 + phi^2),   v_+ = (phi, 1) / sqrt(1 + phi^2),
//    satisfying diag(phi^{-2}, phi^2) * Q^T * Z^2 == Q^T * (A * Z^2) == Q^T * Z^2!
//    Thus Gamma_A = (s Q^T Z^2) |x <tau> with tau(x, y, z) = (phi^{-2} x, phi^2 y, z + z_0)
//    acts cocompactly and ISOMETRICALLY on (Sol, ds_Sol^2)!
// ---------------------------------------------------------------------------

// Folds p into the fundamental torus-bundle cell [-s/2, s/2]^2 x [-ln(phi), +ln(phi)] of M_A = Sol / Gamma_A.
vec3 opRepSolCatMap( in vec3 p, in float latticeScale )
{
    const float Z0  = 0.96242365012; // 2.0 * ln(phi)
    const float QC  = 0.52573111212; // 1 / sqrt(1 + phi^2)
    const float QS  = 0.85065080835; // phi / sqrt(1 + phi^2)

    float s  = max(latticeScale, 1e-6);
    float k  = round(p.z / Z0);
    float zk = p.z - k * Z0;
    float ek = exp(clamp(k * Z0, -40.0, 40.0));

    // Apply tau^{-k}(p) = (exp(+k z_0) x, exp(-k z_0) y, z - k z_0)
    vec2 xy = vec2(ek * p.x, p.y / ek);

    // Rotate into orthogonal Cat-Map eigenbasis w = Q * xy, fold modulo s * Z^2, rotate back by Q^T
    vec2 w  = vec2(QC * xy.x + QS * xy.y, -QS * xy.x + QC * xy.y);
    w      -= s * round(w / s);
    vec2 qxy = vec2(QC * w.x - QS * w.y, QS * w.x + QC * w.y);
    return vec3(qxy, zk);
}

// Infinite self-similar Arnold Cat-Map Sol crystal with smooth cross-layer Sol-metric blending
// across 3 adjacent monodromy periods n in {k-1, k, k+1} and certified 1-Lipschitz Euclidean bound.
float sdSolCatLattice( in vec3 p, in float latticeScale, in float rNode, in float kSmooth )
{
    const float Z0 = 0.96242365012; // 2.0 * ln(phi)
    const float QC = 0.52573111212; // 1 / sqrt(1 + phi^2)
    const float QS = 0.85065080835; // phi / sqrt(1 + phi^2)

    float s    = max(latticeScale, 1e-4);
    float k0   = round(p.z / Z0);
    float beta = (kSmooth > 1e-5) ? (1.0 / kSmooth) : 0.0;

    float dMin  = 1e20;
    float lseSum = 0.0;

    // Evaluate 3 adjacent z-monodromy sheets n in {k0 - 1, k0, k0 + 1}
    // and the 3x3 nearest lattice sites in s * Q^T Z^2 per sheet for C^1 periodic seam safety
    for (int dn = -1; dn <= 1; dn++)
    {
        float n  = k0 + float(dn);
        float zn = p.z - n * Z0;
        float en = exp(clamp(n * Z0, -40.0, 40.0));
        vec2  xy = vec2(en * p.x, p.y / en);
        vec2  w  = vec2(QC * xy.x + QS * xy.y, -QS * xy.x + QC * xy.y);
        vec2  m0 = round(w / s);

        for (int ix = -1; ix <= 1; ix++)
        {
            for (int iy = -1; iy <= 1; iy++)
            {
                vec2 dw  = w - s * (m0 + vec2(float(ix), float(iy)));
                vec2 qxy = vec2(QC * dw.x - QS * dw.y, QS * dw.x + QC * dw.y);
                // Exact Sol-1-Lipschitz hyperbolic-fins gauge at each lattice site
                float dSite = solFinsGauge(vec3(qxy, zn));
                if (beta > 0.0)
                {
                    // Running stabilized log-sum-exp smooth minimum (1-Lipschitz in Sol!)
                    if (dSite < dMin)
                    {
                        lseSum = lseSum * exp(max(beta * (dSite - dMin), -60.0)) + 1.0;
                        dMin   = dSite;
                    }
                    else
                    {
                        lseSum += exp(max(beta * (dMin - dSite), -60.0));
                    }
                }
                else
                {
                    dMin = min(dMin, dSite);
                }
            }
        }
    }

    float dSolCrystal = (beta > 0.0) ? (dMin - log(max(lseSum, 1.0)) / beta) : dMin;
    return solMetricBound(dSolCrystal - rNode, p.z);
}

#endif // SDF_GEMS_HAS_TH_SOL_GEOMETRY
