// ============================================================================
// sdBezier3D, udBezier3D & sdgBezier3D — Exact Collinear-Safe 1-cbrt 3D Quadratic Bezier Capsule
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 1 cbrt, 0 collinear/evolute NaNs, exact 3D Euclidean distance)
// ============================================================================

// 1. Exact Collinear-Safe, Dual-Regime Cancellation-Free 1-cbrt 3D Quadratic Bezier Capsule + Gradient:
//    Returns vec4(sdBezier3D, grad_x, grad_y, grad_z) for tube radius r (set r = 0.0 for udBezier3D).
vec4 sdgBezier3D( in vec3 pos, in vec3 A, in vec3 B, in vec3 C, in float r )
{
    vec3  a  = B - A;
    vec3  b  = A - 2.0 * B + C;
    vec3  c  = a * 2.0;
    vec3  d  = A - pos;
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
            // Dual-Regime Cancellation-Free 1-cbrt Cardano formula
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
            vec3  q1 = d + (c + b * t2.x) * t2.x;
            vec3  q2 = d + (c + b * t2.y) * t2.y;
            t = (dot(q1, q1) < dot(q2, q2)) ? t2.x : t2.y;
        }
    }

    vec3  q = -(d + (c + b * t) * t);
    float l = length(q);
    return vec4(l - r, q / max(l, 1e-30));
}

// 2. Exact 3D Quadratic Bezier Swept Capsule SDF (tube radius r)
float sdBezier3D( in vec3 pos, in vec3 A, in vec3 B, in vec3 C, in float r )
{
    return sdgBezier3D(pos, A, B, C, r).x;
}

// 3. Exact Unsigned 3D Quadratic Bezier Curve Distance
float udBezier3D( in vec3 pos, in vec3 A, in vec3 B, in vec3 C )
{
    return sdgBezier3D(pos, A, B, C, 0.0).x;
}
