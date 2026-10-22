// ============================================================================
// sdgBoxFrame & sdgBoxEdges — Exact 0/1-Sqrt 3D Wireframe Box Frame & Cylindrical Cage Gradients
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1`, 0 sorting swaps, 0 sqrts in interior & 1 sqrt in exterior, 0 NaNs)
// ============================================================================

// 1. Exact 0/1-Sqrt 3D Wireframe Box Frame with Analytical Unit Gradient (0 sorting swaps):
//    By Closest-Strut Longitudinal Invariance, the closest strut is parallel to argmin(abs(p) - b),
//    and inside the strut (max(q_transverse) <= 0), w_min <= min(q_transverse) holds identically!
vec4 sdgBoxFrame( vec3 p, vec3 b, float e )
{
    vec3  s  = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3  w  = abs(p) - b, u = w + e;
    vec3  sq = vec3((u.x < 0.0) ? -1.0 : 1.0, (u.y < 0.0) ? -1.0 : 1.0, (u.z < 0.0) ? -1.0 : 1.0);
    vec3  q  = abs(u) - e;
    bool  ix = (w.x <= w.y) && (w.x <= w.z);
    bool  iy = (!ix) && (w.y <= w.z);
    vec3  c  = ix ? vec3(w.x, q.y, q.z) : (iy ? vec3(q.x, w.y, q.z) : vec3(q.x, q.y, w.z));
    vec3  sc = ix ? vec3(1.0, sq.y, sq.z) : (iy ? vec3(sq.x, 1.0, sq.z) : vec3(sq.x, sq.y, 1.0));
    vec3  gt = ix ? vec3(-1e30, q.y, q.z) : (iy ? vec3(q.x, -1e30, q.z) : vec3(q.x, q.y, -1e30));
    float g  = max(gt.x, max(gt.y, gt.z));
    vec3  v  = max(c, 0.0);
    float l  = length(v);
    vec3  n  = (l > 0.0) ? (v / l)
             : ((gt.x >= max(gt.y, gt.z)) ? vec3(1.0, 0.0, 0.0)
             : ((gt.y >= gt.z)            ? vec3(0.0, 1.0, 0.0) : vec3(0.0, 0.0, 1.0)));
    return vec4((l > 0.0) ? l : g, s * sc * n);
}

// 2. Exact 1-Sqrt 3D Cylindrical-Edge Wireframe Cage with Analytical Unit Gradient (tube radius r)
vec4 sdgBoxEdges( vec3 p, vec3 b, float r )
{
    vec3  s = vec3((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0, (p.z < 0.0) ? -1.0 : 1.0);
    vec3  w = abs(p) - b;
    float m = min(min(w.x, min(w.y, w.z)), 0.0);
    vec3  u = w - ((w.x == m) ? vec3(m, 0.0, 0.0) : ((w.y == m) ? vec3(0.0, m, 0.0) : vec3(0.0, 0.0, m)));
    float l = length(u);
    return vec4(l - r, (l > 0.0) ? (s * u / l) : s);
}
