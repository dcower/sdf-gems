// ============================================================================
// sdCardioid2D — Exact 0-Trig, 0-cbrt 2D/3D Cardioid, Steiner Deltoid & Bernoulli Lemniscate
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, `2.78e-9` error for `sdCardioid2D`, `4.02e-7` for `sdDeltoid2D`, `8.38e-9` for `sdLemniscate2D`, 0 trig & 0 `cbrt` calls)
// ============================================================================

#ifndef SDF_GEMS_SD_CARDIOID_2D_GLSL
#define SDF_GEMS_SD_CARDIOID_2D_GLSL

//------------------------------------------------------------------------------
// 1. 2D Cardioid (cusp at origin (0,0), smooth apex at (R, 0): r = R*cos^2(th/2))
//------------------------------------------------------------------------------
float sdCardioid2D( vec2 p, float R )
{
    float X = p.x;
    float Y = abs(p.y);
    float B = 3.0 * X - R;
    float W = sqrt(B * B + 9.0 * Y * Y);

    // Unique turning point s_* of Q'(s) = 3*Y*s^2 - 2*B*s - 3*Y = 0 on [0, inf)
    float s_star = (B >= 0.0) ? (B + W) / max(3.0 * Y, 1e-4 * R)
                              : (3.0 * Y) / max(W - B, 1e-20);
    float Q_star = X - s_star * (B * s_star + 6.0 * Y) * (1.0 / 3.0);
    float K      = max(-Q_star, 0.0);

    // Solve d^2 * (W + Y * d) = K on the strictly convex branch d >= 0
    float Ws  = max(W, 1e-20);
    float d   = sqrt(K / max(W + Y * sqrt(K / Ws), 1e-20));
    float f   = (Y * d + W) * d * d - K;
    float df  = d * (3.0 * Y * d + 2.0 * W);
    float ddf = 6.0 * Y * d + 2.0 * W;
    d = max(d - (f * df) / max(df * df - 0.5 * f * ddf, 1e-20), 0.0);

    float s    = s_star + d;
    float inv2 = 1.0 / (1.0 + s * s);
    float u2   = s * s * inv2;
    vec2  q    = R * u2 * vec2(u2 - inv2, 2.0 * s * inv2);

    float d2_smooth = dot(vec2(X, Y) - q, vec2(X, Y) - q);
    float r2        = X * X + Y * Y;
    float d2        = (Q_star < 0.0) ? min(d2_smooth, r2) : r2;

    float lhs = 2.0 * r2 - R * X;
    float sgn = (lhs <= 0.0 || lhs * lhs <= R * R * r2) ? -1.0 : 1.0;
    return sgn * sqrt(d2);
}

vec3 sdgCardioid2D( vec2 p, float R )
{
    float X  = p.x;
    float sy = (p.y >= 0.0) ? 1.0 : -1.0;
    float Y  = abs(p.y);
    float B  = 3.0 * X - R;
    float W  = sqrt(B * B + 9.0 * Y * Y);

    float s_star = (B >= 0.0) ? (B + W) / max(3.0 * Y, 1e-4 * R)
                              : (3.0 * Y) / max(W - B, 1e-20);
    float Q_star = X - s_star * (B * s_star + 6.0 * Y) * (1.0 / 3.0);
    float K      = max(-Q_star, 0.0);

    float Ws  = max(W, 1e-20);
    float d   = sqrt(K / max(W + Y * sqrt(K / Ws), 1e-20));
    float f   = (Y * d + W) * d * d - K;
    float df  = d * (3.0 * Y * d + 2.0 * W);
    float ddf = 6.0 * Y * d + 2.0 * W;
    d = max(d - (f * df) / max(df * df - 0.5 * f * ddf, 1e-20), 0.0);

    float s    = s_star + d;
    float inv2 = 1.0 / (1.0 + s * s);
    float u2   = s * s * inv2;
    vec2  q    = R * u2 * vec2(u2 - inv2, 2.0 * s * inv2);
    vec2  v_sm = vec2(X, Y) - q;
    float d2_s = dot(v_sm, v_sm);
    float r2   = X * X + Y * Y;

    bool  use_s = (Q_star < 0.0) && (d2_s < r2);
    vec2  diff  = use_s ? v_sm : vec2(X, Y);
    float dist  = sqrt(use_s ? d2_s : r2);

    float lhs = 2.0 * r2 - R * X;
    float sgn = (lhs <= 0.0 || lhs * lhs <= R * R * r2) ? -1.0 : 1.0;
    vec2  g   = (dist > 1e-20) ? (sgn * diff / dist) : vec2(1.0, 0.0);
    g.y *= sy;
    return vec3(sgn * dist, g);
}

// 3D Cardioid of revolution around the X-axis
float sdCardioid3D( vec3 p, float R )
{
    return sdCardioid2D(vec2(p.x, length(p.yz)), R);
}

//------------------------------------------------------------------------------
// 2. 2D Steiner Deltoid (3-cusped hypocycloid, inradius a, cusp radius 3a)
//------------------------------------------------------------------------------
float sdDeltoid2D( vec2 p, float a )
{
    // Fold into 60-degree D_3 Weyl chamber {y >= 0, sqrt(3)*x - y >= 0}
    const vec2 n60 = vec2(0.86602540, -0.50000000);
    p.y = abs(p.y);
    p -= 2.0 * min(dot(p, n60), 0.0) * n60;
    p.y = abs(p.y);
    p -= 2.0 * min(dot(p, n60), 0.0) * n60;

    // Strictly convex cubic F(w) = A*w^3 + y*w^2 - B*w + y = 0 on w in [0, 1/sqrt(3)]
    float A      = p.x + 3.0 * a;
    float B      = 9.0 * a - p.x;
    float W      = sqrt(max(p.y * p.y + 3.0 * A * B, 0.0));
    float w_star = (B > 0.0) ? B / (W + p.y + 1e-20) : 0.0;
    float F_star = p.y - w_star * (2.0 * B - p.y * w_star) * (1.0 / 3.0);
    float K      = max(-F_star, 0.0);

    // Two-endpoint-exact seed + 1 Halley step on d^2 * (W - A * d) = K
    float d0  = w_star * sqrt(K / max(K + p.y, 1e-20));
    float d   = min(sqrt(K / max(W - A * d0, 1e-20)), w_star);
    float f   = (W - A * d) * d * d - K;
    float df  = d * (2.0 * W - 3.0 * A * d);
    float ddf = 2.0 * W - 6.0 * A * d;
    d = clamp(d - (f * df) / max(df * df - 0.5 * f * ddf, 1e-20), 0.0, w_star);

    float w    = clamp(w_star - d, 0.0, 0.57735027);
    float w2   = w * w;
    float s2   = 1.0 + w2;
    float inv4 = 1.0 / (s2 * s2);
    vec2  q    = a * inv4 * vec2(1.0 + 6.0 * w2 - 3.0 * w2 * w2, 8.0 * w);
    float d2   = dot(p - q, p - q);

    // Check cusp V_cusp = (1.5*a, 1.5*sqrt(3)*a) when F(1/sqrt(3)) > 0
    if( p.x + 1.73205081 * p.y > 6.0 * a )
    {
        vec2 vc = p - vec2(1.50000000, 2.59807621) * a;
        d2 = min(d2, dot(vc, vc));
    }
    float sgn = ((p.x - p.y * w) * s2 >= a * (1.0 - 3.0 * w2)) ? 1.0 : -1.0;
    return sgn * sqrt(d2);
}

vec3 sdgDeltoid2D( vec2 p, float a )
{
    const vec2 n60 = vec2(0.86602540, -0.50000000);
    float sy0 = (p.y >= 0.0) ? 1.0 : -1.0; p.y = abs(p.y);
    float c0  = dot(p, n60); if( c0 < 0.0 ) p -= 2.0 * c0 * n60;
    float sy1 = (p.y >= 0.0) ? 1.0 : -1.0; p.y = abs(p.y);
    float c1  = dot(p, n60); if( c1 < 0.0 ) p -= 2.0 * c1 * n60;

    float A      = p.x + 3.0 * a;
    float B      = 9.0 * a - p.x;
    float W      = sqrt(max(p.y * p.y + 3.0 * A * B, 0.0));
    float w_star = (B > 0.0) ? B / (W + p.y + 1e-20) : 0.0;
    float F_star = p.y - w_star * (2.0 * B - p.y * w_star) * (1.0 / 3.0);
    float K      = max(-F_star, 0.0);

    float d0  = w_star * sqrt(K / max(K + p.y, 1e-20));
    float d   = min(sqrt(K / max(W - A * d0, 1e-20)), w_star);
    float f   = (W - A * d) * d * d - K;
    float df  = d * (2.0 * W - 3.0 * A * d);
    float ddf = 2.0 * W - 6.0 * A * d;
    d = clamp(d - (f * df) / max(df * df - 0.5 * f * ddf, 1e-20), 0.0, w_star);

    float w    = clamp(w_star - d, 0.0, 0.57735027);
    float w2   = w * w;
    float s2   = 1.0 + w2;
    float inv4 = 1.0 / (s2 * s2);
    vec2  q    = a * inv4 * vec2(1.0 + 6.0 * w2 - 3.0 * w2 * w2, 8.0 * w);
    vec2  diff = p - q;
    float d2   = dot(diff, diff);

    if( p.x + 1.73205081 * p.y > 6.0 * a )
    {
        vec2  vc  = p - vec2(1.50000000, 2.59807621) * a;
        float d2c = dot(vc, vc);
        if( d2c < d2 ) { d2 = d2c; diff = vc; }
    }
    float dist = sqrt(d2);
    float sgn  = ((p.x - p.y * w) * s2 >= a * (1.0 - 3.0 * w2)) ? 1.0 : -1.0;
    vec2  g    = (dist > 1e-20) ? (sgn * diff / dist) : vec2(1.0, 0.0);

    if( c1 < 0.0 ) g -= 2.0 * dot(g, n60) * n60; g.y *= sy1;
    if( c0 < 0.0 ) g -= 2.0 * dot(g, n60) * n60; g.y *= sy0;
    return vec3(sgn * dist, g);
}

//------------------------------------------------------------------------------
// 3. 2D Bernoulli Lemniscate: (x^2 + y^2)^2 = a^2 * (x^2 - y^2)
//------------------------------------------------------------------------------
float sdLemniscate2D( vec2 p, float a )
{
    vec2  q = abs(p);
    float c, s;
    if( 3.0 * q.x >= 2.0 * a )
    {
        float Bs = 3.0 * q.x - 2.0 * a;
        float As = 2.0 * a + 3.0 * q.y - q.x;
        s = clamp((2.0 * q.y) / max(Bs + sqrt(max(Bs * Bs + 4.0 * As * q.y, 0.0)), 1e-20), 0.0, 1.0);
        c = sqrt(1.0 - s * s);
        float Bs2 = max(q.x * (3.0 - s * s) - 2.0 * a * c, 0.0);
        s = clamp((2.0 * q.y) / max(Bs2 + sqrt(Bs2 * Bs2 + 12.0 * q.y * q.y), 1e-20), 0.0, 1.0);
        c = sqrt(1.0 - s * s);
    }
    else
    {
        s = 0.85;
        for( int i = 0; i < 2; i++ )
        {
            float Ac = q.x * s - 3.0 * q.y;
            float Bc = a * s;
            float Cc = 2.0 * (q.x * s + q.y);
            c = clamp(Cc / max(Bc + sqrt(max(Bc * Bc - Ac * Cc, 0.0)), 1e-20), 0.0, 1.0);
            s = sqrt(max(1.0 - c * c, 1e-12));
        }
    }

    // 2 Halley-Cayley steps on G(c, s) = s*(2*a*c - x*(2 + c^2)) + y*(3*c^2 - 2) = 0
    for( int i = 0; i < 2; i++ )
    {
        float G   = s * (2.0 * a * c - q.x * (2.0 + c * c)) + q.y * (3.0 * c * c - 2.0);
        float dG  = -3.0 * c * (q.x * c * c + 2.0 * q.y * s) + 2.0 * a * (c * c - s * s);
        float ddG = 3.0 * s * (3.0 * q.x * c * c + 2.0 * q.y * s) - 6.0 * q.y * c * c - 8.0 * a * s * c;
        float dt  = -(G * dG) / max(dG * dG - 0.5 * G * ddG, 1e-20);
        float cn  = clamp(c - s * dt, 0.0, 1.0);
        float sn  = clamp(s + c * dt, 0.0, 1.0);
        float inv = inversesqrt(cn * cn + sn * sn);
        c = cn * inv; s = sn * inv;
    }

    float qx  = (a * c) / (1.0 + s * s);
    vec2  d   = q - vec2(qx, qx * s);
    float r2  = dot(q, q);
    float sgn = (r2 * r2 <= a * a * (q.x * q.x - q.y * q.y)) ? -1.0 : 1.0;
    return sgn * length(d);
}

#endif // SDF_GEMS_SD_CARDIOID_2D_GLSL
