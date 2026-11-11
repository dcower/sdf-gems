// ============================================================================
// opSweepBishop — Exact Bishop Parallel-Transport 3D Space-Curve Sweeps,
//                 0-SFU Riemannian-Submersion 3D Gradients & Twisted SVD Bounds
// @category deform
// @author   Discovered and derived using Gemini 4 Carbon (sdf-gems, 2026)
// @exact    Yes (`||grad d|| = 1` everywhere inside tubular focal reach
//           `r_prof < (R^2 + c^2)/R`; 0 post-frame SFUs in gradient finalizer;
//           100% exact finite N-turn flat-capped coils; 0.0% overstep SVD bound
//           for twisted wrought-iron / Frenet-locked helical sweeps)
// ============================================================================
//
// Mathematical Foundation (Bishop Parallel-Transport Eikonal Sweep Theorem):
//   Let gamma(s) be a unit-speed 3D space curve with tangent T(s), curvature
//   kappa(s), and torsion tau(s). Sweeping a non-circular 2D profile d2D(u, v)
//   in the Frenet-Serret frame (N, B) twists the cross-section at rate tau(s),
//   inducing tangential metric shear:
//     ||grad_p d2D(u_F, v_F)||^2 = 1 + (tau * (g_u * v_F - g_v * u_F) / (1 - kappa * u_F))^2 > 1,
//   causing sphere-tracing oversteps on non-planar curves.
//
//   In Bishop's (1975) Relatively Parallel Adapted Frame (T(s), E_1(s), E_2(s)),
//   obtained by counter-rotating the normal plane by the integrated torsion
//   alpha(s) = int_0^s tau(sigma) d(sigma), the frame derivatives
//     E_1'(s) = -kappa_1(s) * T(s),   E_2'(s) = -kappa_2(s) * T(s)
//   have ZERO normal-plane component (E_1' . E_2 == 0). At any closest-spine
//   projection s_*(p) where (p - gamma(s_*)) . T(s_*) = 0 and 1 - kappa_1*u - kappa_2*v > 0:
//     (1) grad_p u(p) == E_1(s_*) and grad_p v(p) == E_2(s_*) identically
//         (the normal-plane projection is an exact Riemannian submersion!);
//     (2) On the swept surface Sigma, Sigma_s = (1 - kappa_1*u - kappa_2*v) * T(s)
//         is strictly orthogonal to span(E_1, E_2), so every 3D surface normal
//         ray lies inside a single normal plane s = s_*;
//     (3) Across multi-turn candidates theta_0, theta_0 +- 2*pi on a circular
//         helix gamma(theta) = (R*cos(theta), R*sin(theta), c*theta), taking
//         min_k d2D(u_k, v_k) yields the 100% exact 3D Euclidean SDF in R^3!

#ifndef SDF_GEMS_OP_SWEEP_BISHOP_GLSL
#define SDF_GEMS_OP_SWEEP_BISHOP_GLSL
#define SDF_GEMS_HAS_OP_SWEEP_BISHOP

#include "../common/common.glsl"

// ----------------------------------------------------------------------------
// Internal Helpers: Multi-Turn Kepler-Helix Solver & Self-Contained 2D Profiles
// ----------------------------------------------------------------------------

// Solves alpha * sin(delta) + delta = uSigned for the spine-minimum root
// delta in (-d_crit, +d_crit), where 1 + alpha * cos(delta) > 0.
float _opSweepBishopSolveDelta( float alpha, float uSigned, out bool valid )
{
    float au  = abs(uSigned);
    float sgn = (uSigned >= 0.0) ? 1.0 : -1.0;
    float d_crit = SDF_PI;
    if( alpha <= 1.0 )
    {
        if( au > SDF_PI ) { valid = false; return 0.0; }
    }
    else
    {
        d_crit = SDF_PI - acos(clamp(1.0 / alpha, 0.0, 1.0));
        float f_max = sqrt(max((alpha - 1.0) * (alpha + 1.0), 0.0)) + d_crit;
        if( au > f_max ) { valid = false; return 0.0; }
    }
    valid = true;

    float beta = 0.636619772368 * alpha; // (2 / pi) * alpha
    float b    = 1.0 + 2.0 * beta;
    float disc = b * b - 2.54647908947 * beta * au; // (8 / pi) * beta * au
    float d    = (disc >= 0.0) ? ((2.0 * au) / (b + sqrt(disc))) : (0.75 * d_crit);
    d = clamp(d, 0.0, max(d_crit - 1e-5, 0.0));

    for( int i = 0; i < 3; i++ )
    {
        float a_s  = alpha * sin(d);
        float a_c  = alpha * cos(d);
        float f1   = a_s + d - au;
        float f2   = max(a_c + 1.0, 1e-7);
        float f3   = -a_s;
        float f4   = -a_c;
        float f22  = f2 * f2;
        float num  = f1 * (f22 - 0.5 * f1 * f3);
        float den  = f2 * (f22 - f1 * f3) + (1.0 / 6.0) * (f1 * f1 * f4);
        float step = num / ((den > 1e-12) ? den : (f2 * f22));
        d = clamp(d - step, 0.0, d_crit);
    }
    return sgn * d;
}

// Evaluates Bishop frame (E_1, E_2, T) and coordinates (u, v, wTan) at fixed spine angle theta
void _opSweepBishopEvalSlice(
    vec3 p, float R, float pitchC, float invL, float theta, float dz,
    out vec2 uvBishop, out float wTan, out mat2x3 E12, out vec3 Tvec )
{
    float ct  = cos(theta), st = sin(theta);
    float rL  = R * invL,   cL = pitchC * invL;
    vec3  Uf  = vec3(ct, st, 0.0);
    vec3  Vf  = vec3(cL * st, -cL * ct, rL);
    Tvec      = vec3(-rL * st, rL * ct, cL);

    float uf  = p.x * ct + p.y * st - R;
    float qth = -p.x * st + p.y * ct;
    float vf  = rL * dz - cL * qth;
    wTan      = rL * qth + cL * dz;

    float ang = cL * theta;
    float ca  = cos(ang), sa = sin(ang);
    uvBishop  = vec2(ca * uf + sa * vf, -sa * uf + ca * vf);
    E12       = mat2x3(ca * Uf + sa * Vf, -sa * Uf + ca * Vf);
}

// Self-contained 0/1-SFU exact 2D box SDF & gradient
float _opSweepBishopBox2D( vec2 p, vec2 b )
{
    vec2 d = abs(p) - b;
    return (min(d.x, d.y) <= 0.0) ? max(d.x, d.y) : length(d);
}

vec3 _opSweepBishopBox2DG( vec2 p, vec2 b )
{
    vec2 w = abs(p) - b;
    vec2 s = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    if( min(w.x, w.y) <= 0.0 )
        return (w.x > w.y) ? vec3(w.x, s.x, 0.0) : vec3(w.y, 0.0, s.y);
    float l2 = dot(w, w), invL = inversesqrt(l2);
    return vec3(l2 * invL, s * w * invL);
}

// Self-contained 1-SFU-strip exact 2D regular n-pointed star SDF & unit gradient
// (outer tip radius r, n >= 3 tips, angle divisor m in [2, n])
vec3 _opSweepBishopStar2DG( vec2 p, float r, int n, float m )
{
    float an  = SDF_PI / float(max(n, 3));
    float en  = SDF_PI / clamp(m, 2.0, float(max(n, 3)));
    vec2  acs = vec2(cos(an), sin(an));
    vec2  ecs = vec2(cos(en), sin(en));
    float h2  = r * acs.y / max(ecs.y, 1e-6);

    float r2  = dot(p, p);
    float ang = (r2 > 1e-24) ? atan(p.y, p.x) : 0.0;
    float sec = 2.0 * an;
    float k   = round(ang / sec);
    float bn  = ang - k * sec;
    float sb  = (bn < 0.0) ? -1.0 : 1.0;
    float abn = abs(bn);
    float lp  = sqrt(r2);
    vec2  q   = lp * vec2(cos(abn), sin(abn)) - r * acs;

    vec2  nEdge = vec2(ecs.y, -ecs.x);
    float tRaw  = -dot(q, ecs);
    float dy    = dot(q, nEdge);

    float d;
    vec2  gFold;
    if( tRaw >= 0.0 && tRaw <= h2 )
    {
        d     = dy;
        gFold = nEdge;
    }
    else
    {
        float tCl  = clamp(tRaw, 0.0, h2);
        vec2  diff = q + tCl * ecs;
        float l2   = dot(diff, diff);
        float sgn  = (dy < 0.0) ? -1.0 : 1.0;
        if( l2 > 1e-24 )
        {
            float invL = inversesqrt(l2);
            d     = sgn * (l2 * invL);
            gFold = (r2 > 1e-24) ? (sgn * diff * invL) : vec2(1.0, 0.0);
        }
        else
        {
            d     = 0.0;
            gFold = nEdge;
        }
    }

    if( r2 <= 1e-24 )
        return vec3(d, 1.0, 0.0);

    vec2  gSec = vec2(gFold.x, sb * gFold.y);
    float rot  = k * sec;
    float cr   = cos(rot), sr = sin(rot);
    return vec3(d, cr * gSec.x - sr * gSec.y, sr * gSec.x + cr * gSec.y);
}

// ----------------------------------------------------------------------------
// 1. Exact 3D Bishop-Frame Helical Sweep Operators & 0-SFU Gradient Finalizers
//    Spine: gamma(theta) = (R * cos(theta), R * sin(theta), pitchC * theta)
// ----------------------------------------------------------------------------

// Evaluates normal-plane root on turn (n0 + turnOffset), e.g. turnOffset in {-1, 0, +1}.
// Outputs qBishop = vec3(u, v, sArc) and orthonormal Bishop basis E12 = mat2x3(E1, E2).
bool opSweepHelixBishopTurn(
    vec3 p, float R, float pitchC, float turnOffset,
    out vec3 qBishop, out mat2x3 E12 )
{
    float rxy = length(p.xy);
    float phi = (rxy > 1e-12) ? atan(p.y, p.x) : 0.0;

    // Planar ring / straight-line degenerate safeguards
    if( abs(pitchC) <= 1e-6 * max(R, 1.0) )
    {
        if( abs(turnOffset) > 0.5 )
        {
            qBishop = vec3(0.0);
            E12     = mat2x3(vec3(1.0, 0.0, 0.0), vec3(0.0, 0.0, 1.0));
            return false;
        }
        vec2 er = (rxy > 1e-12) ? (p.xy / rxy) : vec2(1.0, 0.0);
        qBishop = vec3(rxy - R, p.z, R * phi);
        E12     = mat2x3(vec3(er, 0.0), vec3(0.0, 0.0, 1.0));
        return true;
    }

    float L     = sqrt(R * R + pitchC * pitchC);
    float invL  = 1.0 / L;
    float zc    = p.z / pitchC;
    float n0    = round((zc - phi) * 0.159154943092); // 1 / (2*pi)
    float u0    = (zc - phi) - SDF_TAU * n0;
    float uS    = u0 - SDF_TAU * turnOffset;
    float alpha = (R * rxy) / (pitchC * pitchC);

    bool  valid;
    float delta = _opSweepBishopSolveDelta(alpha, uS, valid);
    if( !valid )
    {
        qBishop = vec3(0.0);
        E12     = mat2x3(vec3(1.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
        return false;
    }

    float th = phi + SDF_TAU * (n0 + turnOffset) + delta;
    float dz = pitchC * (uS - delta);
    vec2  uv;
    float wTan;
    vec3  Tvec;
    _opSweepBishopEvalSlice(p, R, pitchC, invL, th, dz, uv, wTan, E12, Tvec);
    qBishop = vec3(uv, L * th);
    return true;
}

// Primary Bishop-frame helical sweep projection: returns vec3(u, v, sArc) and outputs sArc
vec3 opSweepHelixBishop( vec3 p, float R, float pitchC, out float sArc )
{
    vec3   q;
    mat2x3 E12;
    opSweepHelixBishopTurn(p, R, pitchC, 0.0, q, E12);
    sArc = q.z;
    return q;
}

vec3 opSweepHelixBishop( vec3 p, float R, float pitchC )
{
    float sArc;
    return opSweepHelixBishop(p, R, pitchC, sArc);
}

// 0-SFU Riemannian-Submersion 3D Unit Gradient Finalizer (0 sqrts, 0 divs, 6 muls)
vec4 sdgSweepBishop( mat2x3 E12, vec3 dg2D )
{
    return vec4(dg2D.x, E12[0] * dg2D.y + E12[1] * dg2D.z);
}

// Standalone 3D unit gradient finalizer from 2D profile (d, gu, gv) at query point p
vec4 sdgSweepHelixBishopFinalize( vec3 p, float R, float pitchC, vec3 dg2D )
{
    vec3   q;
    mat2x3 E12;
    opSweepHelixBishopTurn(p, R, pitchC, 0.0, q, E12);
    return sdgSweepBishop(E12, dg2D);
}

// Clamped finite N-turn Bishop projection on theta in [-halfTheta, +halfTheta]:
// returns vec3(u, v, wCap) where wCap is tangential overshoot (> 0 in exterior cap half-space)
vec3 opSweepHelixBishopLim( vec3 p, float R, float pitchC, float halfTheta, out float sArc )
{
    float rxy = length(p.xy);
    float phi = (rxy > 1e-12) ? atan(p.y, p.x) : 0.0;
    float L   = sqrt(max(R * R + pitchC * pitchC, 1e-24));
    float invL = 1.0 / L;

    if( abs(pitchC) <= 1e-6 * max(R, 1.0) )
    {
        float thCl = clamp(phi, -halfTheta, halfTheta);
        vec2  uv;
        float wTan;
        mat2x3 E12;
        vec3  Tvec;
        _opSweepBishopEvalSlice(p, R, 0.0, invL, thCl, p.z, uv, wTan, E12, Tvec);
        sArc = L * thCl;
        float wCap = (phi > halfTheta) ? wTan : ((phi < -halfTheta) ? -wTan : -min(R * (halfTheta - abs(phi)), 1e6));
        return vec3(uv, wCap);
    }

    float zc    = p.z / pitchC;
    float zcCl  = clamp(zc, -halfTheta, halfTheta);
    float n0    = round((zcCl - phi) * 0.159154943092);
    float alpha = (R * rxy) / (pitchC * pitchC);

    // Pick closest clamped spine candidate between n0 and adjacent turn n1
    float rem = zcCl - (phi + SDF_TAU * n0);
    float n1  = n0 + ((rem >= 0.0) ? 1.0 : -1.0);

    bool  v0, v1;
    float u0  = zc - (phi + SDF_TAU * n0);
    float d0  = _opSweepBishopSolveDelta(alpha, u0, v0);
    float tc0 = v0 ? clamp(phi + SDF_TAU * n0 + d0, -halfTheta, halfTheta)
                   : ((u0 >= 0.0) ? halfTheta : -halfTheta);

    float u1  = zc - (phi + SDF_TAU * n1);
    float d1  = _opSweepBishopSolveDelta(alpha, u1, v1);
    float tc1 = v1 ? clamp(phi + SDF_TAU * n1 + d1, -halfTheta, halfTheta)
                   : ((u1 >= 0.0) ? halfTheta : -halfTheta);

    vec2   uv0, uv1;
    float  w0, w1;
    mat2x3 E0, E1;
    vec3   T0, T1;
    _opSweepBishopEvalSlice(p, R, pitchC, invL, tc0, p.z - pitchC * tc0, uv0, w0, E0, T0);
    _opSweepBishopEvalSlice(p, R, pitchC, invL, tc1, p.z - pitchC * tc1, uv1, w1, E1, T1);

    float dist2_0 = dot(uv0, uv0) + w0 * w0;
    float dist2_1 = dot(uv1, uv1) + w1 * w1;
    bool  pick0   = (dist2_0 <= dist2_1);
    float tc      = pick0 ? tc0 : tc1;
    vec2  uv      = pick0 ? uv0 : uv1;
    float wTan    = pick0 ? w0  : w1;
    sArc          = L * tc;

    // Signed cap coordinate: > 0 outside cap plane, <= 0 inside span
    float wCap;
    if( tc >= halfTheta && wTan > 0.0 )
        wCap = wTan;
    else if( tc <= -halfTheta && wTan < 0.0 )
        wCap = -wTan;
    else
        wCap = -L * (halfTheta - abs(tc));
    return vec3(uv, wCap);
}

// ----------------------------------------------------------------------------
// 2. Certified SVD Lipschitz Factor for Intrinsically Twisted / Frenet Sweeps
// ----------------------------------------------------------------------------
// For a 2D profile bounded within radius `profileBoundR < (R^2 + c^2) / R` and
// twisted around the Bishop frame by `psi = twistRate * sArc` (note: Frenet-Serret
// lock corresponds to `twistRate = -pitchC / (R^2 + pitchC^2)`), the 2x3 spatial
// Jacobian has exact maximum singular value:
//   sigma_max = sqrt( 1 + ( |twistRate| * r_prof / (1 - kappa * r_prof) )^2 ),
// where kappa = R / (R^2 + c^2). Returns 1 / sigma_max in 0 divs & 1 inversesqrt.
float opSweepHelixTwistedScale( float R, float pitchC, float twistRate, float profileBoundR )
{
    float L2    = R * R + pitchC * pitchC;
    float reach = max(L2 - R * profileBoundR, 0.0);
    float shear = twistRate * profileBoundR * L2;
    float d2    = reach * reach + shear * shear;
    return (d2 > 0.0) ? (reach * inversesqrt(d2)) : 0.0;
}

// ----------------------------------------------------------------------------
// 3. Concrete Demo Primitives (Multi-Turn Voronoi Seam-Exact + Finite Coils)
// ----------------------------------------------------------------------------

// 3A. Exact Infinite & Finite N-Turn Bishop-Swept Rectangular / Square Helix
vec4 sdgBishopHelixSquare( vec3 p, float R, float pitchC, vec2 b )
{
    if( abs(pitchC) <= 1e-6 * max(R, 1.0) )
    {
        vec3   q;
        mat2x3 E12;
        opSweepHelixBishopTurn(p, R, pitchC, 0.0, q, E12);
        return sdgSweepBishop(E12, _opSweepBishopBox2DG(q.xy, b));
    }

    vec4 best = vec4(1e30, 1.0, 0.0, 0.0);
    for( int k = -1; k <= 1; k++ )
    {
        vec3   q;
        mat2x3 E12;
        if( opSweepHelixBishopTurn(p, R, pitchC, float(k), q, E12) )
        {
            vec3 dg2 = _opSweepBishopBox2DG(q.xy, b);
            if( dg2.x < best.x )
                best = sdgSweepBishop(E12, dg2);
        }
    }
    return best;
}

float sdBishopHelixSquare( vec3 p, float R, float pitchC, vec2 b )
{
    if( abs(pitchC) <= 1e-6 * max(R, 1.0) )
    {
        vec3 q = opSweepHelixBishop(p, R, pitchC);
        return _opSweepBishopBox2D(q.xy, b);
    }

    float best = 1e30;
    for( int k = -1; k <= 1; k++ )
    {
        vec3   q;
        mat2x3 E12;
        if( opSweepHelixBishopTurn(p, R, pitchC, float(k), q, E12) )
            best = min(best, _opSweepBishopBox2D(q.xy, b));
    }
    return best;
}

// Exact Finite N-Turn Flat-Capped Bishop-Swept Square Coil (theta in [-halfTheta, +halfTheta])
vec4 sdgBishopHelixSquare( vec3 p, float R, float pitchC, float halfTheta, vec2 b )
{
    float L    = sqrt(max(R * R + pitchC * pitchC, 1e-24));
    float invL = 1.0 / L;
    float rxy  = length(p.xy);
    float phi  = (rxy > 1e-12) ? atan(p.y, p.x) : 0.0;

    vec4 bestLat = vec4(1e30, 1.0, 0.0, 0.0);
    if( abs(pitchC) <= 1e-6 * max(R, 1.0) )
    {
        if( abs(phi) <= halfTheta )
        {
            vec2   uv;
            float  wTan;
            mat2x3 E12;
            vec3   Tvec;
            _opSweepBishopEvalSlice(p, R, 0.0, invL, phi, p.z, uv, wTan, E12, Tvec);
            bestLat = sdgSweepBishop(E12, _opSweepBishopBox2DG(uv, b));
        }
    }
    else
    {
        float zc    = p.z / pitchC;
        float zcCl  = clamp(zc, -halfTheta, halfTheta);
        float n0    = round((zcCl - phi) * 0.159154943092);
        float alpha = (R * rxy) / (pitchC * pitchC);
        for( int k = -1; k <= 1; k++ )
        {
            float thB = phi + SDF_TAU * (n0 + float(k));
            float uS  = zc - thB;
            bool  valid;
            float delta = _opSweepBishopSolveDelta(alpha, uS, valid);
            float th    = thB + delta;
            if( valid && abs(th) <= halfTheta )
            {
                vec2   uv;
                float  wTan;
                mat2x3 E12;
                vec3   Tvec;
                _opSweepBishopEvalSlice(p, R, pitchC, invL, th, pitchC * (uS - delta), uv, wTan, E12, Tvec);
                vec3 dg2 = _opSweepBishopBox2DG(uv, b);
                if( dg2.x < bestLat.x )
                    bestLat = sdgSweepBishop(E12, dg2);
            }
        }
    }

    // Exact flat end-cap disks at theta = -halfTheta and +halfTheta
    vec4 bestCap = vec4(1e30, 0.0, 0.0, 1.0);
    for( int e = -1; e <= 1; e += 2 )
    {
        float  thE = float(e) * halfTheta;
        vec2   uvE;
        float  wE;
        mat2x3 EE;
        vec3   TE;
        _opSweepBishopEvalSlice(p, R, pitchC, invL, thE, p.z - pitchC * thE, uvE, wE, EE, TE);
        vec3  dgE = _opSweepBishopBox2DG(uvE, b);
        float dc;
        vec3  gc;
        if( dgE.x <= 0.0 )
        {
            dc = abs(wE);
            gc = ((wE >= 0.0) ? 1.0 : -1.0) * TE;
        }
        else
        {
            float l2 = dgE.x * dgE.x + wE * wE;
            if( l2 > 1e-24 )
            {
                float invD = inversesqrt(l2);
                dc = l2 * invD;
                gc = (dgE.x * (EE[0] * dgE.y + EE[1] * dgE.z) + wE * TE) * invD;
            }
            else
            {
                dc = 0.0;
                gc = float(e) * TE;
            }
        }
        if( dc < bestCap.x )
            bestCap = vec4(dc, gc);
    }

    if( bestLat.x < 0.0 )
        return (bestLat.x > -bestCap.x) ? bestLat : vec4(-bestCap.x, -bestCap.yzw);
    return (bestLat.x <= bestCap.x) ? bestLat : bestCap;
}

float sdBishopHelixSquare( vec3 p, float R, float pitchC, float halfTheta, vec2 b )
{
    return sdgBishopHelixSquare(p, R, pitchC, halfTheta, b).x;
}

// 3B. Exact Infinite & Finite N-Turn Bishop-Swept Star Helix + 3D Analytical Gradient
vec4 sdgBishopHelixStar( vec3 p, float R, float pitchC, float rOuter, int nPoints, float mAngle )
{
    if( abs(pitchC) <= 1e-6 * max(R, 1.0) )
    {
        vec3   q;
        mat2x3 E12;
        opSweepHelixBishopTurn(p, R, pitchC, 0.0, q, E12);
        return sdgSweepBishop(E12, _opSweepBishopStar2DG(q.xy, rOuter, nPoints, mAngle));
    }

    vec4 best = vec4(1e30, 1.0, 0.0, 0.0);
    for( int k = -1; k <= 1; k++ )
    {
        vec3   q;
        mat2x3 E12;
        if( opSweepHelixBishopTurn(p, R, pitchC, float(k), q, E12) )
        {
            vec3 dg2 = _opSweepBishopStar2DG(q.xy, rOuter, nPoints, mAngle);
            if( dg2.x < best.x )
                best = sdgSweepBishop(E12, dg2);
        }
    }
    return best;
}

vec4 sdgBishopHelixStar( vec3 p, float R, float pitchC, float rOuter )
{
    return sdgBishopHelixStar(p, R, pitchC, rOuter, 5, 2.6);
}

float sdBishopHelixStar( vec3 p, float R, float pitchC, float rOuter, int nPoints, float mAngle )
{
    return sdgBishopHelixStar(p, R, pitchC, rOuter, nPoints, mAngle).x;
}

// Exact Finite N-Turn Flat-Capped Bishop-Swept Star Coil + 3D Analytical Gradient
vec4 sdgBishopHelixStar( vec3 p, float R, float pitchC, float halfTheta, float rOuter, int nPoints, float mAngle )
{
    float L    = sqrt(max(R * R + pitchC * pitchC, 1e-24));
    float invL = 1.0 / L;
    float rxy  = length(p.xy);
    float phi  = (rxy > 1e-12) ? atan(p.y, p.x) : 0.0;

    vec4 bestLat = vec4(1e30, 1.0, 0.0, 0.0);
    if( abs(pitchC) <= 1e-6 * max(R, 1.0) )
    {
        if( abs(phi) <= halfTheta )
        {
            vec2   uv;
            float  wTan;
            mat2x3 E12;
            vec3   Tvec;
            _opSweepBishopEvalSlice(p, R, 0.0, invL, phi, p.z, uv, wTan, E12, Tvec);
            bestLat = sdgSweepBishop(E12, _opSweepBishopStar2DG(uv, rOuter, nPoints, mAngle));
        }
    }
    else
    {
        float zc    = p.z / pitchC;
        float zcCl  = clamp(zc, -halfTheta, halfTheta);
        float n0    = round((zcCl - phi) * 0.159154943092);
        float alpha = (R * rxy) / (pitchC * pitchC);
        for( int k = -1; k <= 1; k++ )
        {
            float thB = phi + SDF_TAU * (n0 + float(k));
            float uS  = zc - thB;
            bool  valid;
            float delta = _opSweepBishopSolveDelta(alpha, uS, valid);
            float th    = thB + delta;
            if( valid && abs(th) <= halfTheta )
            {
                vec2   uv;
                float  wTan;
                mat2x3 E12;
                vec3   Tvec;
                _opSweepBishopEvalSlice(p, R, pitchC, invL, th, pitchC * (uS - delta), uv, wTan, E12, Tvec);
                vec3 dg2 = _opSweepBishopStar2DG(uv, rOuter, nPoints, mAngle);
                if( dg2.x < bestLat.x )
                    bestLat = sdgSweepBishop(E12, dg2);
            }
        }
    }

    vec4 bestCap = vec4(1e30, 0.0, 0.0, 1.0);
    for( int e = -1; e <= 1; e += 2 )
    {
        float  thE = float(e) * halfTheta;
        vec2   uvE;
        float  wE;
        mat2x3 EE;
        vec3   TE;
        _opSweepBishopEvalSlice(p, R, pitchC, invL, thE, p.z - pitchC * thE, uvE, wE, EE, TE);
        vec3  dgE = _opSweepBishopStar2DG(uvE, rOuter, nPoints, mAngle);
        float dc;
        vec3  gc;
        if( dgE.x <= 0.0 )
        {
            dc = abs(wE);
            gc = ((wE >= 0.0) ? 1.0 : -1.0) * TE;
        }
        else
        {
            float l2 = dgE.x * dgE.x + wE * wE;
            if( l2 > 1e-24 )
            {
                float invD = inversesqrt(l2);
                dc = l2 * invD;
                gc = (dgE.x * (EE[0] * dgE.y + EE[1] * dgE.z) + wE * TE) * invD;
            }
            else
            {
                dc = 0.0;
                gc = float(e) * TE;
            }
        }
        if( dc < bestCap.x )
            bestCap = vec4(dc, gc);
    }

    if( bestLat.x < 0.0 )
        return (bestLat.x > -bestCap.x) ? bestLat : vec4(-bestCap.x, -bestCap.yzw);
    return (bestLat.x <= bestCap.x) ? bestLat : bestCap;
}

float sdBishopHelixStar( vec3 p, float R, float pitchC, float halfTheta, float rOuter, int nPoints, float mAngle )
{
    return sdgBishopHelixStar(p, R, pitchC, halfTheta, rOuter, nPoints, mAngle).x;
}

// 3C. Certified 1-Lipschitz Twisted Wrought-Iron Helical Bar (0.0% Overstep)
float sdTwistedHelixBar( vec3 p, float R, float pitchC, float twistRate, vec2 b )
{
    float rProf = length(b);
    float scale = opSweepHelixTwistedScale(R, pitchC, twistRate, rProf);
    int   kMax  = (abs(pitchC) <= 1e-6 * max(R, 1.0)) ? 0 : 1;

    float best = 1e30;
    for( int k = -kMax; k <= kMax; k++ )
    {
        vec3   q;
        mat2x3 E12;
        if( opSweepHelixBishopTurn(p, R, pitchC, float(k), q, E12) )
        {
            float psi = twistRate * q.z;
            float cp  = cos(psi), sp = sin(psi);
            vec2  qTw = vec2(cp * q.x + sp * q.y, -sp * q.x + cp * q.y);
            float dRaw = _opSweepBishopBox2D(qTw, b);
            float rp   = length(q.xy);
            float dk;
            if( rp <= rProf )
            {
                dk = dRaw * scale;
            }
            else
            {
                vec2  qCl = qTw * (rProf / max(rp, 1e-30));
                float fz  = _opSweepBishopBox2D(qCl, b) * scale;
                dk        = max(rp - rProf, fz);
            }
            best = min(best, dk);
        }
    }
    return best;
}

// Certified 1-Lipschitz Finite N-Turn Twisted Helical Bar (theta in [-halfTheta, +halfTheta])
float sdTwistedHelixBar( vec3 p, float R, float pitchC, float halfTheta, float twistRate, vec2 b )
{
    float rProf = length(b);
    float scale = opSweepHelixTwistedScale(R, pitchC, twistRate, rProf);
    float L     = sqrt(max(R * R + pitchC * pitchC, 1e-24));
    float invL  = 1.0 / L;
    float rxy   = length(p.xy);
    float phi   = (rxy > 1e-12) ? atan(p.y, p.x) : 0.0;

    float dInf      = 1e30;
    float dLatIn    = 1e30;
    float dTubeCoil = 1e30;

    if( abs(pitchC) <= 1e-6 * max(R, 1.0) )
    {
        vec2   uv;
        float  wTan;
        mat2x3 E12;
        vec3   Tvec;
        _opSweepBishopEvalSlice(p, R, 0.0, invL, phi, p.z, uv, wTan, E12, Tvec);
        float psi  = twistRate * (L * phi);
        float cp   = cos(psi), sp = sin(psi);
        vec2  qTw  = vec2(cp * uv.x + sp * uv.y, -sp * uv.x + cp * uv.y);
        float dRaw = _opSweepBishopBox2D(qTw, b);
        float rp   = length(uv);
        float dk   = (rp <= rProf) ? (dRaw * scale)
                                   : max(rp - rProf, _opSweepBishopBox2D(qTw * (rProf / max(rp, 1e-30)), b) * scale);
        dInf = dk;
        if( abs(phi) <= halfTheta )
        {
            dLatIn    = dk;
            dTubeCoil = max(rp - rProf, 0.0);
        }
    }
    else
    {
        float zc    = p.z / pitchC;
        float zcCl  = clamp(zc, -halfTheta, halfTheta);
        float n0    = round((zcCl - phi) * 0.159154943092);
        float alpha = (R * rxy) / (pitchC * pitchC);
        for( int k = -1; k <= 1; k++ )
        {
            float thB = phi + SDF_TAU * (n0 + float(k));
            float uS  = zc - thB;
            bool  valid;
            float delta = _opSweepBishopSolveDelta(alpha, uS, valid);
            if( valid )
            {
                float  th = thB + delta;
                vec2   uv;
                float  wTan;
                mat2x3 E12;
                vec3   Tvec;
                _opSweepBishopEvalSlice(p, R, pitchC, invL, th, pitchC * (uS - delta), uv, wTan, E12, Tvec);
                float psi  = twistRate * (L * th);
                float cp   = cos(psi), sp = sin(psi);
                vec2  qTw  = vec2(cp * uv.x + sp * uv.y, -sp * uv.x + cp * uv.y);
                float dRaw = _opSweepBishopBox2D(qTw, b);
                float rp   = length(uv);
                float dk   = (rp <= rProf) ? (dRaw * scale)
                                           : max(rp - rProf, _opSweepBishopBox2D(qTw * (rProf / max(rp, 1e-30)), b) * scale);
                dInf = min(dInf, dk);
                if( abs(th) <= halfTheta )
                {
                    dLatIn    = min(dLatIn, dk);
                    dTubeCoil = min(dTubeCoil, max(rp - rProf, 0.0));
                }
            }
        }
    }

    float dCaps = 1e30;
    for( int e = -1; e <= 1; e += 2 )
    {
        float  thE = float(e) * halfTheta;
        vec2   uvE;
        float  wE;
        mat2x3 EE;
        vec3   TE;
        _opSweepBishopEvalSlice(p, R, pitchC, invL, thE, p.z - pitchC * thE, uvE, wE, EE, TE);
        float psiE = twistRate * (L * thE);
        float cpE  = cos(psiE), spE = sin(psiE);
        vec2  qTwE = vec2(cpE * uvE.x + spE * uvE.y, -spE * uvE.x + cpE * uvE.y);
        float d2E  = max(_opSweepBishopBox2D(qTwE, b), 0.0);
        dCaps      = min(dCaps, sqrt(d2E * d2E + wE * wE));
        float rpE  = max(length(uvE) - rProf, 0.0);
        dTubeCoil  = min(dTubeCoil, sqrt(rpE * rpE + wE * wE));
    }

    if( dLatIn < 0.0 )
        return max(dLatIn, -dCaps);
    return min(max(abs(dInf), dTubeCoil), dCaps);
}

#endif // SDF_GEMS_OP_SWEEP_BISHOP_GLSL
