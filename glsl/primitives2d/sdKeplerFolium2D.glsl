// ============================================================================
// sdKeplerFolium2D, sdgKeplerFolium2D, sdDoubleFolium2D, sdgDoubleFolium2D,
// sdTrifolium2D, sdgTrifolium2D, sdKeplerEgg3D & sdgKeplerEgg3D
// — Exact 0-Trig Quartic Folium Family: Kepler's (1609) Oviform Folium,
//   Dürer's (1525) Double Folium, Regular Trifolium & 3D Kepler Egg
//
// @author  Discovered and derived using Gemini 4 Argon (research & mathematical design)
//          and Gemini 3.8 Flash (implementation & verification) (sdf-gems, 2026)
// @license MIT
// @exact   Yes (`||grad d|| = 1` everywhere in R^2 / R^3; 0 trig in solver & 0-trig 0-sqrt quartic interior sign tests)
// @novelty Unifies the entire Quartic Folium Family (x^2 + y^2)^2 = x * (alpha*x^2 + beta*y^2),
//          or r(theta) = cos(theta) * (alpha*cos^2(theta) + beta*sin^2(theta)), via a
//          single 0-trig Weierstrass solver u = (1 - s^2)/(1 + s^2), w = 2*s/(1 + s^2):
//          - 1 Lobe (alpha = a, beta = a - b, 0 <= b <= a): Johannes Kepler's (1609,
//            Astronomia Nova) Oviform / Simple Folium Loop, plus 3D Kepler Egg!
//          - 2 Lobes (alpha = 0, beta = 4*a): Albrecht Dürer's (1525) Double Folium
//            (x^2 + y^2)^2 = 4*a*x*y^2!
//          - 3 Lobes (alpha = a, beta = -3*a): Regular 3-Petal Trifolium r = a*cos(3*theta),
//            reduced via 0-trig D_3 reflection across n_60 = (-sqrt(3)/2, 1/2) to a single
//            half-petal on s in [0, 2 - sqrt(3)]!
// ============================================================================

#ifndef SDF_GEMS_SD_KEPLER_FOLIUM_2D_GLSL
#define SDF_GEMS_SD_KEPLER_FOLIUM_2D_GLSL

#include "../common/common.glsl"

// Shared 0-trig solver for r(theta) = (alpha - (alpha - beta)*sin^2(theta)) * cos(theta) on s in [0, sMax]
vec3 _sdQuarticFoliumCore( in float px, in float ay, in float alpha, in float beta, in float sMax )
{
    float qx0    = max(alpha, 0.0);
    vec2  bestQ  = vec2(qx0, 0.0);
    float bestD2 = (px - qx0) * (px - qx0) + ay * ay;

    float d2_0 = px * px + ay * ay;
    if( d2_0 < bestD2 )
    {
        bestD2 = d2_0;
        bestQ  = vec2(0.0, 0.0);
    }

    vec3  s = vec3(0.15 * sMax, 0.50 * sMax, 0.85 * sMax);
    float b = alpha - beta;

    for( int i = 0; i < 6; i++ )
    {
        vec3 s2   = s * s;
        vec3 invD = 1.0 / (1.0 + s2);
        vec3 u    = (1.0 - s2) * invD;
        vec3 w    = 2.0 * s * invD;
        vec3 u2   = u * u;
        vec3 w2   = w * w;
        vec3 uw   = u * w;
        vec3 u2w2 = u2 - w2;

        vec3 R   = vec3(alpha) - b * w2;
        vec3 dR  = (-2.0 * b) * uw;
        vec3 d2R = (-2.0 * b) * u2w2;

        vec3 qx   = u2 * R;
        vec3 dqx  = -2.0 * uw * R + u2 * dR;
        vec3 d2qx = -2.0 * u2w2 * R - 4.0 * uw * dR + u2 * d2R;

        vec3 qy   = uw * R;
        vec3 dqy  = u2w2 * R + uw * dR;
        vec3 d2qy = -4.0 * uw * R + 2.0 * u2w2 * dR + uw * d2R;

        vec3 dx = qx - vec3(px);
        vec3 dy = qy - vec3(ay);
        vec3 f  = dx * dqx + dy * dqy;
        vec3 df = dqx * dqx + dqy * dqy + dx * d2qx + dy * d2qy;

        vec3 stepT = vec3(
            (df.x > 1e-12) ? (f.x / df.x) : ((f.x >= 0.0) ? 0.2 : -0.2),
            (df.y > 1e-12) ? (f.y / df.y) : ((f.y >= 0.0) ? 0.2 : -0.2),
            (df.z > 1e-12) ? (f.z / df.z) : ((f.z >= 0.0) ? 0.2 : -0.2)
        );
        stepT = clamp(stepT, vec3(-0.35), vec3(0.35));
        s = clamp(s - 0.5 * (1.0 + s2) * stepT, vec3(0.0), vec3(sMax));
    }

    for( int j = 0; j < 3; j++ )
    {
        float sj   = s[j];
        float s2j  = sj * sj;
        float invD = 1.0 / (1.0 + s2j);
        float uj   = (1.0 - s2j) * invD;
        float wj   = 2.0 * sj * invD;
        float Rj   = alpha - b * wj * wj;
        vec2  qj   = vec2(uj * uj * Rj, uj * wj * Rj);
        vec2  dv   = vec2(px, ay) - qj;
        float d2   = dot(dv, dv);
        if( d2 < bestD2 )
        {
            bestD2 = d2;
            bestQ  = qj;
        }
    }

    return vec3(bestD2, bestQ.x, bestQ.y);
}

// 1. Exact 0-Trig 2D Johannes Kepler's (1609) Oviform / Simple Folium Loop + Unit Gradient vec3(d, grad.xy)
//    r(theta) = cos(theta) * (a - b * sin^2(theta)), 0 <= b <= a (when b = a, r = a * cos^3(theta) Simple Folium).
vec3 sdgKeplerFolium2D( in vec2 p, in float a, in float b )
{
    a = max(a, 1e-5);
    b = clamp(b, 0.0, a);

    float px = p.x;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ay = abs(p.y);

    vec3  res   = _sdQuarticFoliumCore(px, ay, a, a - b, 1.0);
    vec2  bestQ = res.yz;

    float r2     = px * px + ay * ay;
    bool  inside = (px > 0.0) && (r2 * r2 < px * (a * px * px + (a - b) * ay * ay));
    float sgn    = inside ? -1.0 : 1.0;

    vec2  diff = vec2(px - bestQ.x, sy * (ay - bestQ.y));
    float d2   = dot(diff, diff);
    float invD = inversesqrt(max(d2, 1e-30));
    float d    = d2 * invD;
    vec2  grad = (d2 > 1e-24) ? (sgn * diff * invD) : vec2(1.0, 0.0);
    return vec3(sgn * d, grad);
}

// 2. Scalar 2D Kepler's Oviform / Folium Loop SDF
float sdKeplerFolium2D( in vec2 p, in float a, in float b )
{
    return sdgKeplerFolium2D(p, a, b).x;
}

// 1-parameter overloads for the canonical Simple Folium / Kepler Ovoid (b = a: (x^2 + y^2)^2 = a * x^3)
vec3 sdgKeplerFolium2D( in vec2 p, in float a )
{
    return sdgKeplerFolium2D(p, a, a);
}

float sdKeplerFolium2D( in vec2 p, in float a )
{
    return sdgKeplerFolium2D(p, a, a).x;
}

// 3. Exact 0-Trig 2D Albrecht Dürer's (1525) Double Folium + Unit Gradient vec3(d, grad.xy)
//    (x^2 + y^2)^2 = 4 * a * x * y^2 <=> r(theta) = 4 * a * cos(theta) * sin^2(theta)
vec3 sdgDoubleFolium2D( in vec2 p, in float a )
{
    a = max(a, 1e-5);
    float px = p.x;
    float sy = (p.y < 0.0) ? -1.0 : 1.0;
    float ay = abs(p.y);

    vec3  res   = _sdQuarticFoliumCore(px, ay, 0.0, 4.0 * a, 1.0);
    vec2  bestQ = res.yz;

    float r2     = px * px + ay * ay;
    bool  inside = (px > 0.0) && (r2 * r2 < 4.0 * a * px * ay * ay);
    float sgn    = inside ? -1.0 : 1.0;

    vec2  diff = vec2(px - bestQ.x, sy * (ay - bestQ.y));
    float d2   = dot(diff, diff);
    float invD = inversesqrt(max(d2, 1e-30));
    float d    = d2 * invD;
    vec2  grad = (d2 > 1e-24) ? (sgn * diff * invD) : vec2(1.0, 0.0);
    return vec3(sgn * d, grad);
}

// 4. Scalar 2D Dürer's Double Folium SDF
float sdDoubleFolium2D( in vec2 p, in float a )
{
    return sdgDoubleFolium2D(p, a).x;
}

// 5. Exact 0-Trig 2D Regular 3-Petal Trifolium + Analytical Unit Gradient vec3(d, grad.xy)
//    (x^2 + y^2)^2 = a * x * (x^2 - 3 * y^2) <=> r(theta) = a * cos(3 * theta) >= 0
vec3 sdgTrifolium2D( in vec2 p, in float a )
{
    a = max(a, 1e-5);
    const vec2 n60 = vec2(-0.8660254037844386, 0.5);

    vec2  q   = p;
    float sy0 = (q.y < 0.0) ? -1.0 : 1.0; q.y = abs(q.y);
    bool  r0  = (dot(q, n60) > 0.0);       if( r0 ) q -= 2.0 * dot(q, n60) * n60;
    float sy1 = (q.y < 0.0) ? -1.0 : 1.0; q.y = abs(q.y);

    // In the folded 60-degree wedge [0, pi/3], the positive half-petal spans theta in [0, pi/6],
    // so s = tan(theta/2) spans [0, tan(pi/12)] = [0, 2 - sqrt(3)].
    vec3  res   = _sdQuarticFoliumCore(q.x, q.y, a, -3.0 * a, 0.2679491924311227);
    vec2  bestQ = res.yz;

    float r2     = dot(q, q);
    float rhs    = a * q.x * (q.x * q.x - 3.0 * q.y * q.y);
    bool  inside = (rhs > 0.0) && (r2 * r2 < rhs);
    float sgn    = inside ? -1.0 : 1.0;

    vec2  diff = q - bestQ;
    float d2   = dot(diff, diff);
    float invD = inversesqrt(max(d2, 1e-30));
    float d    = d2 * invD;
    vec2  g    = (d2 > 1e-24) ? (sgn * diff * invD) : vec2(1.0, 0.0);

    g.y *= sy1;
    if( r0 ) g -= 2.0 * dot(g, n60) * n60;
    g.y *= sy0;
    return vec3(sgn * d, g);
}

// 6. Scalar 2D Regular 3-Petal Trifolium SDF
float sdTrifolium2D( in vec2 p, in float a )
{
    return sdgTrifolium2D(p, a).x;
}

// 7. Exact 3D Kepler Egg / Oviform Solid of Revolution + Analytical Unit Gradient vec4(d, grad.xyz)
vec4 sdgKeplerEgg3D( in vec3 p, in float a, in float b )
{
    float rho = length(p.yz);
    vec3  r2d = sdgKeplerFolium2D(vec2(p.x, rho), a, b);
    vec2  er  = (rho > 1e-12) ? (p.yz / rho) : vec2(1.0, 0.0);
    return vec4(r2d.x, r2d.y, r2d.z * er.x, r2d.z * er.y);
}

// 8. Scalar 3D Kepler Egg / Oviform Solid SDF
float sdKeplerEgg3D( in vec3 p, in float a, in float b )
{
    return sdgKeplerFolium2D(vec2(p.x, length(p.yz)), a, b).x;
}

// 1-parameter overloads for the canonical Kepler Egg (b = a)
vec4 sdgKeplerEgg3D( in vec3 p, in float a )
{
    return sdgKeplerEgg3D(p, a, a);
}

float sdKeplerEgg3D( in vec3 p, in float a )
{
    return sdKeplerEgg3D(p, a, a);
}

#endif // SDF_GEMS_SD_KEPLER_FOLIUM_2D_GLSL
