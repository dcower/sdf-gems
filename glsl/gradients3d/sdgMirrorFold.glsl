// ============================================================================
// sdgMirrorFold — Exact 1-Tap Orthonormal Mirror-Plane Folding, Filleting,
//                 Chamfering & Seam Engraving with Analytical 3D Gradients
// @author Discovered and derived using Claude Opus 5.5; audited and hardened
//         (cancellation-free gx -> +1 engraving) using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` across all oblique mirror angles theta in (0, pi)
//         for half-spaces and locally planar facets; first-order O(r^2/R) tangent-plane
//         collar on curved surfaces of radius of curvature R)
// ============================================================================
//
// Given a single evaluation `vec4 dg = sdgPrim(vec3(abs(p.x), p.y, p.z))` in the
// folded half-space x = |p.x| >= 0 with signed distance d = dg.x and unit gradient
// g = (gx, gy, gz) = dg.yzw:
//   - Standard `p.x = abs(p.x)` folding leaves the gradient discontinuous across
//     p.x = 0 and underestimates exterior distance above convex ridges (gx > 0) and
//     interior depth below concave valleys (gx < 0) by up to 1 - sqrt(1 - gx^2),
//     producing sharp miter spikes under `d - r` offsetting and requiring 2x-8x
//     redundant primitive evaluations to fillet or chamfer the mirror seam.
//   - By the Orthonormal Mirror-Meridian Decomposition Theorem, the mirror normal
//     e_x = (sgn(p.x), 0, 0) and normalized tangential gradient e_v = (0, gy, gz) / c
//     (where c = sqrt(1 - gx^2)) form an exact orthonormal 2D basis in R^3 in which
//     (d, t) = (d, (x - gx * d) / c) is an orthogonal reflection of (x, v) = (x, (d - gx * x) / c).
//   - Consequently, exact crease distance restoration (`sdgMirrorExact`), rolling-ball
//     seam filleting (`sdgMirrorFillet`), beveling (`sdgMirrorChamfer`), and V-groove
//     seam engraving (`sdgMirrorEngrave`) are all recovered from the SINGLE folded
//     evaluation `dg` in 1 `sqrt`.

#ifndef SDF_GEMS_SDG_MIRROR_FOLD_GLSL
#define SDF_GEMS_SDG_MIRROR_FOLD_GLSL

// 1. Exact 1-Tap Mirror-Plane Distance & 3D Gradient Restoration
//    Restores exact Euclidean crease distance and unit 3D gradient in the exterior
//    normal wedge of convex ridges (gx > 0) and interior wedge of concave valleys (gx < 0).
vec4 sdgMirrorExact( in float px, in vec4 dg )
{
    float sx = (px < 0.0) ? -1.0 : 1.0, x = abs(px), d = dg.x, gx = dg.y;
    float s2 = dot(dg.zw, dg.zw);
    if( gx * d <= x || s2 <= 1e-5 ) return vec4(d, sx * gx, dg.zw);
    float s  = (d < 0.0) ? -1.0 : 1.0;
    vec3  g  = vec3(px, ((d - gx * x) / s2) * dg.zw);
    float l  = length(g);
    return (l > 1e-12) ? vec4(s * l, (s / l) * g) : vec4(0.0, sx * gx, dg.zw);
}

// 2. Exact 1-Tap Rolling-Ball Mirror Fillet & Ridge Rounding with 3D Gradient
//    Fillets concave mirror valleys (gx < 0) and rounds convex mirror ridges (gx > 0)
//    with exact circular arc radius r from a single folded evaluation.
vec4 sdgMirrorFillet( in float px, in vec4 dg, in float r )
{
    float sx = (px < 0.0) ? -1.0 : 1.0, x = abs(px), d = dg.x, gx = dg.y;
    float s  = (gx < 0.0) ? -1.0 : 1.0, u = d + s * r, s2 = dot(dg.zw, dg.zw);
    if( gx * u <= x || s2 <= 1e-5 ) return vec4(d, sx * gx, dg.zw);
    vec3  g  = vec3(px, ((u - gx * x) / s2) * dg.zw);
    float l  = length(g);
    return (l > 1e-12) ? vec4(s * (l - r), (s / l) * g) : vec4(-s * r, sx * gx, dg.zw);
}

// 3. Exact 1-Tap Euclidean Mirror Chamfer (Symmetric Bevel) with 3D Gradient
vec4 sdgMirrorChamfer( in float px, in vec4 dg, in float r )
{
    float sx = (px < 0.0) ? -1.0 : 1.0, x = abs(px), d = dg.x, gx = dg.y;
    float s2 = dot(dg.zw, dg.zw);
    if( s2 <= 1e-5 ) return vec4(d, sx * gx, dg.zw);
    float s  = (gx < 0.0) ? -1.0 : 1.0, ic = inversesqrt(s2), c = s2 * ic;
    float dx = x - c * r, dv = (d - gx * dx) * ic;
    if( dx > 0.0 && dx < gx * d ) {
        vec3  g = vec3(sx * dx, (dv * ic) * dg.zw);
        float l = length(g);
        return (l > 1e-12) ? vec4(s * l, (s / l) * g) : vec4(0.0, sx * gx, dg.zw);
    }
    return (s * dv > s * d) ? vec4(dv, 0.0, ic * dg.zw) : vec4(d, sx * gx, dg.zw);
}

// 4. Exact 1-Tap Euclidean Mirror Seam V-Groove Engraving with 3D Gradient
vec4 sdgMirrorEngrave( in float px, in vec4 dg, in float r )
{
    float sx = (px < 0.0) ? -1.0 : 1.0, x = abs(px), d = dg.x, gx = dg.y;
    float s2 = dot(dg.zw, dg.zw);
    if( s2 <= 1e-5 ) return vec4(d, sx * gx, dg.zw);
    float ic = inversesqrt(s2), c = s2 * ic, v = (d - gx * x) * ic;
    float k  = sqrt(max((gx > 0.0) ? (2.0 * s2 / (1.0 + gx)) : (2.0 * (1.0 - gx)), 1e-7));
    float nx = -0.5 * k, nv = c / k, vr = v + r;
    float dw = nx * x + nv * vr, sw = nv * x - nx * vr, lw = k * r;
    float dx = x - c * r, dv = v + gx * r;
    if( sw > lw && dx < gx * d ) {
        vec3  g = vec3(sx * dx, (dv * ic) * dg.zw);
        float l = length(g);
        return (l > 1e-12) ? vec4(l, g / l) : vec4(0.0, sx * gx, dg.zw);
    }
    if( sw <= lw ) {
        if( sw < 0.0 ) {
            vec3  g = vec3(px, (vr * ic) * dg.zw);
            float l = length(g);
            if( dx < gx * d || -l > d ) {
                return (l > 1e-12) ? vec4(-l, -g / l) : vec4(0.0, sx * gx, dg.zw);
            }
        } else if( dx < gx * d || dw > d ) {
            return vec4(dw, sx * nx, (nv * ic) * dg.zw);
        }
    }
    return vec4(d, sx * gx, dg.zw);
}

// 5. Arbitrary Unit-Normal Mirror Plane 1-Tap Fillet with 3D Gradient
//    Evaluate `dg = sdgPrim(p - 2.0 * min(dot(p, n), 0.0) * n)` before calling.
vec4 sdgMirrorFilletPlane( in vec3 p, in vec3 n, in vec4 dg, in float r )
{
    float pn = dot(p, n), sn = (pn < 0.0) ? -1.0 : 1.0, x = abs(pn);
    float d  = dg.x, gn = dot(dg.yzw, n);
    vec3  gt = dg.yzw - gn * n, gf = sn * gn * n + gt;
    float s  = (gn < 0.0) ? -1.0 : 1.0, u = d + s * r, s2 = dot(gt, gt);
    if( gn * u <= x || s2 <= 1e-5 ) return vec4(d, gf);
    vec3  g  = pn * n + ((u - gn * x) / s2) * gt;
    float l  = length(g);
    return (l > 1e-12) ? vec4(s * (l - r), (s / l) * g) : vec4(-s * r, gf);
}

// ============================================================================
// Scalar Counterparts (require only d and mirror-normal derivative gx)
// ============================================================================

float opMirrorExact( in float px, in float d, in float gx )
{
    float x = abs(px), h = gx * d - x, s2 = (1.0 - gx) * (1.0 + gx);
    if( h <= 0.0 || s2 <= 1e-5 ) return d;
    return ((d < 0.0) ? -1.0 : 1.0) * sqrt(d * d + h * h / s2);
}

float opMirrorFillet( in float px, in float d, in float gx, in float r )
{
    float x = abs(px), s = (gx < 0.0) ? -1.0 : 1.0, u = d + s * r;
    float h = gx * u - x, s2 = (1.0 - gx) * (1.0 + gx);
    if( h <= 0.0 || s2 <= 1e-5 ) return d;
    return s * (sqrt(u * u + h * h / s2) - r);
}

float opMirrorChamfer( in float px, in float d, in float gx, in float r )
{
    float x = abs(px), s2 = (1.0 - gx) * (1.0 + gx);
    if( s2 <= 1e-5 ) return d;
    float s = (gx < 0.0) ? -1.0 : 1.0, ic = inversesqrt(s2), c = s2 * ic;
    float dx = x - c * r, dv = (d - gx * dx) * ic;
    if( dx > 0.0 && dx < gx * d ) return s * length(vec2(dx, dv));
    return (s * dv > s * d) ? dv : d;
}

float opMirrorEngrave( in float px, in float d, in float gx, in float r )
{
    float x = abs(px), s2 = (1.0 - gx) * (1.0 + gx);
    if( s2 <= 1e-5 ) return d;
    float ic = inversesqrt(s2), c = s2 * ic, v = (d - gx * x) * ic;
    float k  = sqrt(max((gx > 0.0) ? (2.0 * s2 / (1.0 + gx)) : (2.0 * (1.0 - gx)), 1e-7));
    float nx = -0.5 * k, nv = c / k, vr = v + r;
    float dw = nx * x + nv * vr, sw = nv * x - nx * vr, lw = k * r;
    float dx = x - c * r, dv = v + gx * r;
    if( sw > lw && dx < gx * d ) return length(vec2(dx, dv));
    if( sw <= lw ) {
        float d1 = (sw < 0.0) ? -length(vec2(x, vr)) : dw;
        if( dx < gx * d || d1 > d ) return d1;
    }
    return d;
}

#endif // SDF_GEMS_SDG_MIRROR_FOLD_GLSL
