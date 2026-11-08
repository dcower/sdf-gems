// ============================================================================
// sdf-gems: sdHelicalGear
// Category: primitives3d
// Exact Euclidean 3D SDF of an involute helical gear (infinite face width),
// a finite-face-width version, and the standard (module / pressure angle /
// helix angle) parameterization. Gear axis = +z.
//
// Involute-Helicoid Parallel-Offset Theorem: a helical tooth flank is an
// involute helicoid (the tangent developable of the base helix). With twist
// rate tau = tan(bb) / rb (bb = base helix angle) and untwisted angle
// phi' = phi - tau * z, the function
//
//     f = rb * cos(bb) * (phi' - psi(rho)),   psi(rho) = psib - inv(acos(rb / rho)),
//
// has |grad f| = 1 and straight characteristics lying in the planes tangent
// to the base cylinder: rotating the helicoid about the axis by delta gives
// a PARALLEL surface at normal offset rb * cos(bb) * delta. So f is the exact
// signed distance to the flank wherever its foot point (roll parameter
// t = sqrt(rho^2 - rb^2) / rb - cos(bb)^2 * (phi' - psi)) lies on the flank.
// Closed form: no root finding, no Lipschitz factor.
//
// Orbit-space dihedral fold: the gear is invariant under the screw motion
// (rotate tau*dz, translate dz) and under the half-turn (phi', z) -> (-phi', -z),
// both isometries, and the screw-orbit distance is monotone in |delta phi'|.
// Hence p can be screwed to z = 0 and folded into the chamber phi' in [0, PI/n],
// and the nearest point always lies on the six chamber features:
//   tip cylinder, tip-edge helix, involute flank, straight-helicoid strip
//   (radial flank below the base cylinder), root-edge helix, root cylinder.
// Each helix / strip distance reduces to the Kepler-type equation
// a*sin(x) + x = c (same as sdHelix), solved by 2 Householder steps.
//
// @author Claude Opus 5.5 (sdf-gems, 2026)
// ============================================================================
#include "../common/common.glsl"

#ifndef SDF_GEMS_HAS_HELICAL_GEAR
#define SDF_GEMS_HAS_HELICAL_GEAR

// Solves a*sin(x) + x = |c| on x in [0, min(|c|, PI)] (a >= 0, |c| <= PI) by two
// 3rd-order Householder steps from a parabolic-sine guess; returns sign(c) * x.
float sdHelicalGearKepler( in float a, in float c )
{
    float u  = abs(c);
    float xm = min(u, SDF_PI);
    float be = 0.636619772368 * a;                         // (2/PI) * a
    float b  = 1.0 + 2.0 * be;
    float x  = clamp(2.0 * u / (b + sqrt(max(b * b - 2.54647908947 * be * u, 0.0))), 0.0, xm);
    for( int i = 0; i < 2; i++ )
    {
        float as = a * sin(x), ac = a * cos(x);
        float f1 = as + x - u, f2 = ac + 1.0, f3 = -as, f4 = -ac;
        float f22 = f2 * f2;
        float dn = f2 * (f22 - f1 * f3) + (f1 * f1 * f4) / 6.0;
        x = clamp(x - f1 * (f22 - 0.5 * f1 * f3) / ((abs(dn) > 1e-30) ? dn : 1e-30), 0.0, xm);
    }
    return (c < 0.0) ? -x : x;
}

// Distance^2 from the orbit point (rho, phi' = 0, z = 0) to the helix of radius R
// at untwisted angle c (twist rate tau), and the axial offset dz of its foot point.
vec2 sdHelicalGearHelix( in float rho, in float R, in float c, in float tau )
{
    float w  = sdHelicalGearKepler(rho * R * tau * tau, c);  // angular offset of the foot
    float dz = -rho * R * tau * sin(w);                      // = (w - c) / tau, division-free
    float s  = sin(0.5 * w);
    float dr = rho - R;
    return vec2(dr * dr + 4.0 * rho * R * s * s + dz * dz, dz);
}

// Core: returns vec2(signed distance, axial offset of the nearest point) for an
// infinitely long helical gear (axis z). tbb = tan(base helix angle); tbb = 0
// gives the exact 2D spur-gear SDF of the cross-section.
//   n    : number of teeth          rb : base-cylinder radius
//   ra   : tip radius (ra > rb)     rf : root radius (rf < ra)
//   psib : transverse tooth half-angle at the base cylinder (tooth on +x at z = 0)
vec2 sdHelicalGearCore( in vec3 p, in float n, in float rb, in float ra, in float rf,
                        in float psib, in float tbb )
{
    float tau = tbb / rb;                                // twist rate d(phi)/dz
    float cb2 = 1.0 / (1.0 + tbb * tbb);                 // cos^2(base helix angle)
    float rho = length(p.xy);
    float sec = SDF_TAU / n;
    float ang = atan(p.y, p.x + 1e-30) - tau * p.z;      // screw p to z = 0
    float phs = ang - sec * round(ang / sec);
    float sy  = (phs < 0.0) ? -1.0 : 1.0;                // half-turn (phi', z) -> (-phi', -z)
    float phi = abs(phs);                                // chamber [0, PI/n]

    float rlo  = max(rf, rb);
    float tlo  = sqrt(max((rlo - rb) * (rlo + rb), 0.0)) / rb;
    float ta   = sqrt(max((ra - rb) * (ra + rb), 0.0)) / rb;
    float psil = psib - (tlo - atan(tlo));               // flank angle at its start radius
    float psia = psib - (ta  - atan(ta));                // flank angle at the tip

    // Involute flank (Involute-Helicoid Parallel-Offset Theorem)
    float L  = sqrt(max((rho - rb) * (rho + rb), 0.0));
    float G  = phi - psib + L / rb - atan(L, rb);        // = phi' - psi(rho)
    bool  inside = (rho <= rf) || (rho <= ra && ((rho > rb) ? (G <= 0.0) : (phi <= psib)));

    vec2 best = vec2(1e30, 0.0);
    if( phi <= psia || inside ) best = vec2(abs(rho - ra), 0.0);   // tip cylinder
    if( (phi >= psil || !inside) && abs(rho - rf) < best.x ) best = vec2(abs(rho - rf), 0.0); // root cylinder

    float tq = L / rb - cb2 * G;                         // roll parameter of the foot point
    float sf = rb * sqrt(cb2) * G;
    if( rho > rb && (tq >= tlo || (!inside && tq >= 0.0)) && (tq <= ta || inside) && abs(sf) < best.x )
        best = vec2(abs(sf), sf * tbb * sqrt(cb2));      // dz = f * sin(bb)

    // Convex tip-edge helix (exterior only) / concave root-edge helix (interior only)
    vec2 h = inside ? sdHelicalGearHelix(rho, rf, psil - phi, tau)
                    : sdHelicalGearHelix(rho, ra, psia - phi, tau);
    if( h.x < best.x * best.x ) best = vec2(sqrt(h.x), h.y);

    // Straight-helicoid strip r in [rf, rb] at phi' = psib (radial flank part)
    if( rf < rb )
    {
        float w  = 0.5 * sdHelicalGearKepler(rho * rho * tau * tau, 2.0 * (psib - phi));
        float cw = cos(w), sw = sin(w);
        float e  = rho * cw - clamp(rho * cw, rf, rb);
        float dz = -rho * rho * tau * sw * cw;
        float d2 = e * e + rho * rho * sw * sw + dz * dz;
        if( d2 < best.x * best.x ) best = vec2(sqrt(d2), dz);
    }
    return vec2(inside ? -best.x : best.x, sy * best.y);
}

// Exact SDF of an infinitely long involute helical gear (axis z).
//   tbb : tan(base helix angle); tbb > 0 = right-hand helix.
float sdInvoluteHelicoidGear( in vec3 p, in float n, in float rb, in float ra, in float rf,
                              in float psib, in float tbb )
{
    return sdHelicalGearCore(p, n, rb, ra, rf, psib, tbb).x;
}

// Exact SDF + analytical unit gradient: returns vec4(d, dd/dx, dd/dy, dd/dz).
vec4 sdgInvoluteHelicoidGear( in vec3 p, in float n, in float rb, in float ra, in float rf,
                              in float psib, in float tbb )
{
    float tau = tbb / rb;
    float cb  = inversesqrt(1.0 + tbb * tbb);            // cos(base helix angle)
    float rho = length(p.xy);
    float sec = SDF_TAU / n;
    float ang = atan(p.y, p.x + 1e-30) - tau * p.z;
    float k   = round(ang / sec);
    float phs = ang - sec * k;
    float sy  = (phs < 0.0) ? -1.0 : 1.0;
    float phi = abs(phs);
    vec2  er  = vec2(cos(phi), sin(phi));                // folded point = (rho * er, 0)

    float rlo  = max(rf, rb);
    float tlo  = sqrt(max((rlo - rb) * (rlo + rb), 0.0)) / rb;
    float ta   = sqrt(max((ra - rb) * (ra + rb), 0.0)) / rb;
    float psil = psib - (tlo - atan(tlo));
    float psia = psib - (ta  - atan(ta));

    float L  = sqrt(max((rho - rb) * (rho + rb), 0.0));
    float al = atan(L, rb);
    float G  = phi - psib + L / rb - al;
    bool  inside = (rho <= rf) || (rho <= ra && ((rho > rb) ? (G <= 0.0) : (phi <= psib)));

    // Each candidate keeps (distance, unit direction from its foot point to q)
    vec4 best = vec4(1e30, er, 0.0);
    if( phi <= psia || inside ) best = vec4(abs(rho - ra), er * sign(rho - ra + 1e-30), 0.0);
    if( (phi >= psil || !inside) && abs(rho - rf) < best.x ) best = vec4(abs(rho - rf), er * sign(rho - rf + 1e-30), 0.0);

    float tq = L / rb - cb * cb * G;
    float sf = rb * cb * G;
    if( rho > rb && (tq >= tlo || (!inside && tq >= 0.0)) && (tq <= ta || inside) && abs(sf) < best.x )
    {
        float th = phi - al;                             // base-tangent direction of the foot
        best = vec4(abs(sf), sign(sf + 1e-30) * vec3(-cb * sin(th), cb * cos(th), -tbb * cb));
    }

    // Edge helix (tip if outside, root if inside). Foot vectors are built in the
    // frame of q (cancellation-free), then rotated by phi.
    mat2  rq = mat2(er.x, er.y, -er.y, er.x);
    float R  = inside ? rf : ra;
    float w  = sdHelicalGearKepler(rho * R * tau * tau, (inside ? psil : psia) - phi);
    float sw = sin(w), sh = sin(0.5 * w);
    vec3  v  = vec3(rho - R + 2.0 * R * sh * sh, -R * sw, rho * R * tau * sw);
    float lv = length(v);
    if( lv < best.x ) best = vec4(lv, (lv > 1e-20) ? vec3(rq * v.xy, v.z) / lv : vec3(er, 0.0));

    if( rf < rb )
    {
        w  = 0.5 * sdHelicalGearKepler(rho * rho * tau * tau, 2.0 * (psib - phi));
        sw = sin(w);
        float cw = cos(w);
        float e  = rho * cw - clamp(rho * cw, rf, rb);
        v  = vec3(rho * sw * sw + e * cw, (e - rho * cw) * sw, rho * rho * tau * sw * cw);
        lv = length(v);
        if( lv < best.x ) best = vec4(lv, (lv > 1e-20) ? vec3(rq * v.xy, v.z) / lv : vec3(-er.y, er.x, 0.0));
    }

    // Undo the half-turn and the screw / rotation of the fold
    float s  = inside ? -1.0 : 1.0;
    vec3  g  = s * vec3(best.y, sy * best.z, sy * best.w);
    float a  = sec * k + tau * p.z;
    float ca = cos(a), sa = sin(a);
    return vec4(s * best.x, ca * g.x - sa * g.y, sa * g.x + ca * g.y, g.z);
}

// Finite face width |z| <= h. Exact inside, exact outside wherever the nearest
// point of the infinite gear lies within the slab or the end face is nearer than
// the lateral lower bound; elsewhere (a thin shell around the end-face edges,
// beyond the end planes) a conservative lower bound >= max(d_inf, |z| - h).
float sdInvoluteHelicoidGearCapped( in vec3 p, in float n, in float rb, in float ra, in float rf,
                                    in float psib, in float tbb, in float h )
{
    vec2  di = sdHelicalGearCore(p, n, rb, ra, rf, psib, tbb);
    float az = abs(p.z) - h;
    if( di.x <= 0.0 && az <= 0.0 ) return max(di.x, az);
    // Exact distance to the nearer end face: 2D gear distance in that end plane
    float ze = (p.z >= 0.0) ? h : -h;
    float a  = -tbb / rb * ze;
    float ca = cos(a), sa = sin(a);
    float d2 = sdHelicalGearCore(vec3(ca * p.x - sa * p.y, sa * p.x + ca * p.y, 0.0),
                                 n, rb, ra, rf, psib, 0.0).x;
    float dF = sqrt(max(d2, 0.0) * max(d2, 0.0) + az * az);
    // Lateral distance (only points with |z| <= h) is >= B = max(|d_inf|, |z| - h) and,
    // beyond the end plane, >= the 2D-cone bound: the slice distance to the profile is
    // tau*rho-Lipschitz in z, so it is >= d2 - tau*rho*(h - |z_q|).
    float B  = max(abs(di.x), az);
    if( az > 0.0 && d2 > 0.0 )
    {
        float k = abs(tbb) / rb * length(p.xy);
        B = max(B, (k * d2 <= az) ? dF : (d2 + k * az) * inversesqrt(1.0 + k * k));
    }
    float dL = (abs(p.z + di.y) <= h) ? abs(di.x) : B;
    return min(dF, dL);
}

// Standard full-depth helical gear: n teeth, NORMAL module mn, NORMAL pressure
// angle pan (radians, 20 deg = 0.3490659), helix angle beta at the pitch
// cylinder (radians, > 0 right-hand). Addendum mn, dedendum 1.25 mn. Infinite.
float sdHelicalGear( in vec3 p, in float n, in float mn, in float pan, in float beta )
{
    float cbt = cos(beta);
    float at  = atan(tan(pan), cbt);                     // transverse pressure angle
    float rp  = 0.5 * mn * n / cbt;                      // pitch radius (transverse module)
    return sdInvoluteHelicoidGear(p, n, rp * cos(at), rp + mn, rp - 1.25 * mn,
                                  SDF_PI / (2.0 * n) + tan(at) - at, tan(beta) * cos(at));
}

// Standard helical gear with finite face width |z| <= h.
float sdHelicalGearCapped( in vec3 p, in float n, in float mn, in float pan, in float beta, in float h )
{
    float cbt = cos(beta);
    float at  = atan(tan(pan), cbt);
    float rp  = 0.5 * mn * n / cbt;
    return sdInvoluteHelicoidGearCapped(p, n, rp * cos(at), rp + mn, rp - 1.25 * mn,
                                        SDF_PI / (2.0 * n) + tan(at) - at, tan(beta) * cos(at), h);
}

// Standard helical gear (infinite) with exact analytical gradient vec4(d, grad).
vec4 sdgHelicalGear( in vec3 p, in float n, in float mn, in float pan, in float beta )
{
    float cbt = cos(beta);
    float at  = atan(tan(pan), cbt);
    float rp  = 0.5 * mn * n / cbt;
    return sdgInvoluteHelicoidGear(p, n, rp * cos(at), rp + mn, rp - 1.25 * mn,
                                   SDF_PI / (2.0 * n) + tan(at) - at, tan(beta) * cos(at));
}

#endif // SDF_GEMS_HAS_HELICAL_GEAR
