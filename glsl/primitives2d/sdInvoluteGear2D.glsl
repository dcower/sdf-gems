// ============================================================================
// sdf-gems: sdInvoluteGear2D
// Category: primitives2d
// Exact Euclidean 2D SDF of an involute spur gear, a constant-pitch involute
// spiral (clock spring / scroll), and analytical gradients.
//
// Involute Parallel-Offset Theorem: rotating an involute of the base circle
// (radius rb) by an angle delta yields a PARALLEL curve at normal offset
// rb * delta. Hence the signed Euclidean distance from p = (rho, phi) to the
// involute flank whose polar angle at radius rho is psi(rho) is exactly
//
//     d = rb * (phi - psi(rho)),   psi(rho) = psib - inv(acos(rb / rho)),
//
// with inv(a) = tan(a) - a. No root finding, no iteration, no cubic.
// The gear SDF is completed by an exact Voronoi treatment of the tip arc,
// tip corner, radial flank below the base circle, root corner and root arc
// inside the D_n mirror chamber.
//
// @author Claude Opus 5.5 (sdf-gems, 2026)
// ============================================================================
#include "../common/common.glsl"

#ifndef SDF_GEMS_HAS_INVOLUTE_GEAR_2D
#define SDF_GEMS_HAS_INVOLUTE_GEAR_2D

// Exact SDF of an n-tooth involute gear.
//   n    : number of teeth
//   rb   : base-circle radius (involutes unwind from here)
//   ra   : addendum (tip) radius,   ra > rb
//   rf   : root (dedendum) radius,  rf < ra
//   psib : half-angle of a tooth at the base circle (tooth centered on +x).
//          Below the base circle the flank continues radially at psib.
// Requires a non-pointed tooth (psi(ra) >= 0) and an open root gap (psi(rf) < PI/n).
float sdInvoluteGear2D( in vec2 p, in float n, in float rb, in float ra, in float rf, in float psib )
{
    float rho = length(p);
    float sec = SDF_TAU / n;
    float phi = atan(p.y, p.x + 1e-30);
    phi = abs(phi - sec * round(phi / sec));            // D_n mirror chamber [0, PI/n]
    vec2  q   = rho * vec2(cos(phi), sin(phi));

    // Flank roll limits: t = tan(pressure angle) = sqrt(r^2 - rb^2) / rb
    float rlo  = max(rf, rb);
    float tlo  = sqrt(rlo * rlo - rb * rb) / rb;
    float ta   = sqrt(ra * ra - rb * rb) / rb;
    float psil = psib - (tlo - atan(tlo));              // flank angle at its start radius
    float psia = psib - (ta  - atan(ta));               // flank angle at the tip

    // Tip arc / tip corner
    vec2  ct = vec2(cos(psia), sin(psia));
    float d  = (phi <= psia) ? abs(rho - ra) : length(q - ra * ct);

    // Root arc / root corner, and the radial flank segment [rf, rlo] at angle psil
    vec2  cl = vec2(cos(psil), sin(psil));
    d = min(d, (phi >= psil) ? abs(rho - rf) : length(q - rf * cl));
    vec2  w  = vec2(dot(q, cl), cl.x * q.y - cl.y * q.x);
    d = min(d, length(vec2(w.x - clamp(w.x, rf, rlo), w.y)));

    // Involute flank: signed distance rb * (phi - psi(rho)), valid when the
    // foot point's roll t* = psib - phi + acos(rb / rho) lies on the flank.
    float L   = sqrt(max(rho * rho - rb * rb, 0.0));
    float al  = atan(L, rb);                            // pressure angle at rho
    float sf  = rb * (phi - psib) + L - rb * al;        // = rb * (phi - psi(rho))
    float tf  = psib - phi + al;
    if( rho > rb && tf >= tlo && tf <= ta ) d = min(d, abs(sf));

    bool inside = (rho <= rf) || (rho <= ra && ((rho > rb) ? (sf <= 0.0) : (phi <= psib)));
    return inside ? -d : d;
}

// Standard full-depth spur gear: n teeth, module m (pitch diameter / n),
// pressure angle pa (radians, 20 deg = 0.3490659). Addendum m, dedendum 1.25 m.
float sdSpurGear2D( in vec2 p, in float n, in float m, in float pa )
{
    float rp = 0.5 * m * n;
    return sdInvoluteGear2D(p, n, rp * cos(pa), rp + m, rp - 1.25 * m,
                            SDF_PI / (2.0 * n) + tan(pa) - pa);
}

// Exact SDF + analytical gradient: returns vec3(d, dd/dx, dd/dy).
vec3 sdgInvoluteGear2D( in vec2 p, in float n, in float rb, in float ra, in float rf, in float psib )
{
    float rho = length(p);
    float sec = SDF_TAU / n;
    float ang = atan(p.y, p.x + 1e-30);
    float k   = round(ang / sec);
    float phs = ang - sec * k;
    float sy  = (phs < 0.0) ? -1.0 : 1.0;
    float phi = abs(phs);
    vec2  q   = rho * vec2(cos(phi), sin(phi));

    float rlo  = max(rf, rb);
    float tlo  = sqrt(rlo * rlo - rb * rb) / rb;
    float ta   = sqrt(ra * ra - rb * rb) / rb;
    float psil = psib - (tlo - atan(tlo));
    float psia = psib - (ta  - atan(ta));

    // Each candidate keeps (distance, unit direction from the feature to q).
    // Exactly-on-feature points fall back to the radial direction (never a zero vector).
    vec2  er = (rho > 1e-20) ? q / rho : vec2(1.0, 0.0);
    vec2  ct = vec2(cos(psia), sin(psia));
    vec2  dv = q - ra * ct;
    float dt = length(dv);
    vec3  best = (phi <= psia) ? vec3(abs(rho - ra), (rho >= ra) ? er : -er)
                               : vec3(dt, (dt > 1e-20) ? dv / dt : er);

    vec2  cl = vec2(cos(psil), sin(psil));
    vec2  rv = q - rf * cl;
    float dr = length(rv);
    vec3  c  = (phi >= psil) ? vec3(abs(rho - rf), (rho >= rf) ? er : -er)
                             : vec3(dr, (dr > 1e-20) ? rv / dr : er);
    if( c.x < best.x ) best = c;

    vec2  w  = vec2(dot(q, cl), cl.x * q.y - cl.y * q.x);
    vec2  sv = vec2(w.x - clamp(w.x, rf, rlo), w.y);
    float ds = length(sv);
    if( ds < best.x )
    {
        vec2 u = (ds > 1e-20) ? sv / ds : vec2(0.0, 1.0);
        best = vec3(ds, cl * u.x + vec2(-cl.y, cl.x) * u.y);
    }

    float L  = sqrt(max(rho * rho - rb * rb, 0.0));
    float al = atan(L, rb);
    float sf = rb * (phi - psib) + L - rb * al;
    float tf = psib - phi + al;
    if( rho > rb && tf >= tlo && tf <= ta && abs(sf) < best.x )
    {
        // Normal of the involute through q = tangent of the base circle at angle phi - al
        float th = phi - al;
        best = vec3(abs(sf), vec2(-sin(th), cos(th)) * ((sf >= 0.0) ? 1.0 : -1.0));
    }

    bool  inside = (rho <= rf) || (rho <= ra && ((rho > rb) ? (sf <= 0.0) : (phi <= psib)));
    float s = inside ? -1.0 : 1.0;
    // Undo the mirror and the rotation of the fold
    vec2  g = s * vec2(best.y, best.z * sy);
    float a = sec * k;
    float ca = cos(a), sa = sin(a);
    return vec3(s * best.x, ca * g.x - sa * g.y, sa * g.x + ca * g.y);
}

// Exact SDF of a constant-pitch involute spiral (clock spring, rolled sheet,
// scroll-compressor wrap) of strip half-thickness th: the curve unwinds from
// the base circle rb at angle 0 for 'turns' full turns. Successive turns are
// separated by exactly 2*PI*rb along every normal (an Archimedean spiral is not).
float sdInvoluteSpiral2D( in vec2 p, in float rb, in float turns, in float th )
{
    // Curve: c(t) = rb * (cos t + t sin t, sin t - t cos t), t in [0, 2*PI*turns].
    // Its normal at c(t) is the base-circle tangent at angle t, so the roll of
    // every foot point from p is t = phi + acos(rb / rho) (mod 2*PI).
    float rho = length(p);
    float tm  = SDF_TAU * turns;                         // maximum roll parameter
    vec2  p1  = rb * (vec2(cos(tm), sin(tm)) + tm * vec2(sin(tm), -cos(tm)));
    float d   = min(length(p - vec2(rb, 0.0)), length(p - p1));
    if( rho > rb )
    {
        float L   = sqrt(rho * rho - rb * rb);
        float al  = atan(L, rb);
        float t0  = mod(atan(p.y, p.x + 1e-30) + al, SDF_TAU); // roll of the first candidate turn
        float u   = L / rb - t0;                         // = (string length - unwound arc) / rb
        float km  = floor((tm - t0) / SDF_TAU);          // last turn index with t <= tm
        // rb*|u - 2*PI*k| is convex in k => the clamped nearest turn is the global minimum
        if( km >= 0.0 ) d = min(d, rb * abs(u - SDF_TAU * clamp(round(u / SDF_TAU), 0.0, km)));
    }
    return d - th;
}

#endif // SDF_GEMS_HAS_INVOLUTE_GEAR_2D
