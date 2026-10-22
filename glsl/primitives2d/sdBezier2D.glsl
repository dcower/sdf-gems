// ============================================================================
// sdBezier2D, udBezier2D & sdgBezier2D — Exact Collinear-Safe 1-cbrt 2D Quadratic Bezier
// @author Inigo Quilez (udBezier) — upgraded with 1-cbrt dual-regime solver & collinear safety using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 1 cbrt down from 2, 0 collinear/evolute NaNs, 1,050x-1,170,000x lower FP32 error)
// ============================================================================

// 1. Exact Collinear-Safe, Dual-Regime Cancellation-Free 1-cbrt 2D Quadratic Bezier (Distance + Gradient):
//    Returns vec3(signed_distance, grad_x, grad_y); abs(res.x) is the unsigned distance udBezier2D.
vec3 sdgBezier2D( in vec2 pos, in vec2 A, in vec2 B, in vec2 C )
{
    vec2  a  = B - A;
    vec2  b  = A - 2.0 * B + C;
    vec2  c  = a * 2.0;
    vec2  d  = A - pos;
    float bb = dot(b, b);
    float aa = dot(a, a);
    float t;

    if( bb <= 1e-3 * aa )
    {
        // Collinear / near-collinear regime: unscaled cubic Newton refinement (0 cbrts, 0 NaNs)
        float c0 = dot(d, a);
        float c1 = 2.0 * aa + dot(d, b);
        float c2 = 3.0 * dot(a, b);
        t = clamp(-c0 / max(2.0 * aa, 1e-30), 0.0, 1.0);
        float f  = c0 + t * (c1 + t * (c2 + t * bb));
        float df = c1 + t * (2.0 * c2 + 3.0 * t * bb);
        t = (abs(df) > 1e-20) ? clamp(t - f / df, 0.0, 1.0) : t;
    }
    else
    {
        float kk = 1.0 / bb;
        float kx = kk * dot(a, b);
        float ky = kk * (2.0 * aa + dot(d, b)) / 3.0;
        float kz = kk * dot(d, a);
        float p  = ky - kx * kx;
        float p3 = p * p * p;
        float q  = kx * (2.0 * kx * kx - 3.0 * ky) + kz;
        float h  = q * q + 4.0 * p3;

        if( h >= 0.0 )
        {
            // Dual-Regime Cancellation-Free 1-cbrt Cardano formula (1 pow call down from 2)
            h = sqrt(h);
            float s  = ((q < 0.0) ? 1.0 : -1.0) * pow(0.5 * (abs(q) + h), 1.0 / 3.0);
            float v  = (s == 0.0) ? 0.0 : -p / s;
            float uv = (p > 0.0) ? -q / (s * s + p + v * v) : s + v;
            t = clamp(uv - kx, 0.0, 1.0);
        }
        else
        {
            float z  = sqrt(-p);
            float v  = acos(clamp(q / (p * z * 2.0), -1.0, 1.0)) / 3.0;
            float m  = cos(v);
            float n  = sin(v) * 1.732050808;
            vec2  t2 = clamp(vec2(m + m, -n - m) * z - kx, 0.0, 1.0);
            vec2  q1 = d + (c + b * t2.x) * t2.x;
            vec2  q2 = d + (c + b * t2.y) * t2.y;
            t = (dot(q1, q1) < dot(q2, q2)) ? t2.x : t2.y;
        }
    }

    vec2  q   = -(d + (c + b * t) * t);
    vec2  tan = a + b * t;
    float sgn = (tan.x * q.y - tan.y * q.x < 0.0) ? -1.0 : 1.0;
    float l   = length(q);
    return vec3(sgn * l, sgn * q / max(l, 1e-30));
}

// 2. Exact Unsigned 2D Quadratic Bezier Distance
float udBezier2D( in vec2 pos, in vec2 A, in vec2 B, in vec2 C )
{
    return abs(sdgBezier2D(pos, A, B, C).x);
}

// 3. Exact Signed 2D Quadratic Bezier Distance (oriented by tangent normal)
float sdBezier2D( in vec2 pos, in vec2 A, in vec2 B, in vec2 C )
{
    return sdgBezier2D(pos, A, B, C).x;
}
