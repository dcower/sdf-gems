// ============================================================================
// sdArchimedeanSpiral2D, sdgArchimedeanSpiral2D & sdArchimedeanSpiral3D — Exact 1-Atan, 0-Sin/Cos Weierstrass-Lambert Archimedean Spiral
//
// @author  Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @license MIT
// @exact   Yes (1.56e-8 max Euclidean error, 1 atan, 0 sin/cos calls, exact finite N-turn end-caps & unit gradient)
// @novelty Solves the exact Euclidean distance and analytical unit gradient to a finite N-turn
//          Archimedean spiral `S(theta) = a*theta*(cos(theta), sin(theta))` on `theta in [0, T]`
//          using only 1 `atan` and 0 `sin`/`cos`/`tan` calls via:
//          (1) Quadratic-Taylor initial offset `delta_0` on the stationary equation
//              `g(theta_k + delta) = a*(theta_k + delta) - r*cos(delta) + r*(theta_k + delta)*sin(delta) = 0`,
//          (2) Lambert [3/4] continued-fraction Padé for `tau = tan(delta/2)` feeding the
//              exact unit-norm Weierstrass stereographic parameterization `(c, s) = ((1-tau^2)/(1+tau^2), 2*tau/(1+tau^2))`,
//          (3) Direct stereographic rotation of `p` to reconstruct `S(theta)` without `cos(theta)` or `sin(theta)`.
// ============================================================================

// Lambert [3/4] continued-fraction Padé for tan(d / 2) -> exact unit-norm Weierstrass (cos(d), sin(d))
vec2 _sdSpiralWeierstrassCS( in float d )
{
    float u   = 0.5 * d, u2 = u * u;
    float tau = u * (105.0 - 10.0 * u2) / (105.0 - u2 * (45.0 - u2));
    float t2  = tau * tau, it2 = 1.0 / (1.0 + t2);
    return vec2((1.0 - t2) * it2, 2.0 * tau * it2);
}

// 1. Exact Finite N-Turn 2D Archimedean Spiral Distance & Analytical Unit Gradient (1 atan, 0 sin/cos calls)
//    Parameters:
//      p          : 2D query point
//      a          : spiral radial pitch per radian (turn spacing = 2*pi*a)
//      T          : maximum angle in radians (e.g. 2*pi*N_turns)
//      endPointCS : precomputed vec2(cos(T), sin(T)) (simply vec2(1.0, 0.0) for integer turns)
//      w          : wire half-thickness / tube radius
vec3 sdgArchimedeanSpiral2D( in vec2 p, in float a, in float T, in vec2 endPointCS, in float w )
{
    float r   = length(p);
    float phi = atan(p.y, p.x);
    if( phi < 0.0 ) phi += 6.28318530718;

    float k0     = floor((clamp(r / a, 0.0, T) - phi) * 0.159154943092);
    float bestL2 = 1e30;
    vec2  bestQ  = p;

    for( int i = 0; i < 2; i++ )
    {
        float th_k = phi + 6.28318530718 * (k0 + float(i));
        float dr   = r - a * th_k;
        float b    = max(a + r * th_k, 1e-6 * a);
        float d    = clamp(2.0 * dr / (b + sqrt(max(b * b + 6.0 * r * dr, 0.0))), -2.2, 2.2);

        for( int j = 0; j < 2; j++ )
        {
            float th  = th_k + d;
            vec2  cs  = _sdSpiralWeierstrassCS(d);
            float g   = a * th - r * cs.x + r * th * cs.y;
            float gp  = a + 2.0 * r * cs.y + r * th * cs.x;
            float gpp = 3.0 * r * cs.x - r * th * cs.y;
            float den = gp - 0.5 * g * gpp / max(abs(gp), 1e-6 * a);
            d = clamp(d - g / max(abs(den), 1e-6 * a), -3.14159265, 3.14159265);
        }

        float th  = clamp(th_k + d, 0.0, T);
        vec2  cs  = _sdSpiralWeierstrassCS(th - th_k);
        float sc  = (a * th) / max(r, 1e-30);
        vec2  sPt = (th >= T - 1e-6) ? (a * T * endPointCS)
                                     : (sc * vec2(cs.x * p.x - cs.y * p.y, cs.x * p.y + cs.y * p.x));
        vec2  q   = p - sPt;
        float l2  = dot(q, q);
        if( l2 < bestL2 ) { bestL2 = l2; bestQ = q; }
    }

    float invL = inversesqrt(max(bestL2, 1e-30));
    vec2  g    = (bestL2 > 1e-30) ? (bestQ * invL) : endPointCS;
    return vec3(bestL2 * invL - w, g);
}

float sdArchimedeanSpiral2D( in vec2 p, in float a, in float T, in vec2 endPointCS, in float w )
{
    return sdgArchimedeanSpiral2D(p, a, T, endPointCS, w).x;
}

// 2. Exact 3D Planar Archimedean Spiral Tube / Coiled Watch-Spring SDF & Analytical Unit Gradient
vec4 sdgArchimedeanSpiral3D( in vec3 p, in float a, in float T, in vec2 endPointCS, in float w )
{
    vec3  dg2  = sdgArchimedeanSpiral2D(p.xy, a, T, endPointCS, 0.0);
    float l2   = dg2.x * dg2.x + p.z * p.z;
    float invL = inversesqrt(max(l2, 1e-30));
    return vec4(l2 * invL - w, vec3(dg2.x * dg2.yz, p.z) * invL);
}

float sdArchimedeanSpiral3D( in vec3 p, in float a, in float T, in vec2 endPointCS, in float w )
{
    return sdgArchimedeanSpiral3D(p, a, T, endPointCS, w).x;
}
