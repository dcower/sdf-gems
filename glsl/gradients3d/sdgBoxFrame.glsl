// ============================================================================
// sdgBoxFrame & sdgBoxEdges — Exact 1-Sqrt 3D Wireframe Box Frame & Cylindrical Cage Gradients
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1`, 0 sorting swaps, 1 sqrt, 0 NaNs; interior corner cubes corrected, see docs)
// ============================================================================

// 1. Exact 1-Sqrt 3D Wireframe Box Frame with Analytical Unit Gradient (0 sorting swaps):
//    By Closest-Strut Longitudinal Invariance, the closest strut is parallel to argmin(abs(p) - b).
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
    vec3  v  = max(c, 0.0);
    float l  = length(v);
    if( l > 0.0 ) return vec4(l, s * sc * (v / l));

    // Inside the solid (1 sqrt). The nearest non-solid point is either across an outer face (depth -max w)
    // or in the hollow, which needs TWO coordinates below b - 2e: depth = sqrt of the two smallest t_i^2 with
    // t = max(w + 2e, 0). Inside the 8 corner cubes the inner planes w_j = -2e are NOT surface, so the former
    // axis-aligned answer max(q_j, q_k) under-estimated the depth there and its gradient pointed at an interior
    // point (found by the lipschitzcache experiment; fixed with Claude Fable 5.2, sdf-gems 2026).
    // An axis with b_i <= 2e has no hollow side: its struts span the whole box.
    vec3  t  = vec3((b.x > 2.0 * e) ? max(w.x + 2.0 * e, 0.0) : 1e15,
                    (b.y > 2.0 * e) ? max(w.y + 2.0 * e, 0.0) : 1e15,
                    (b.z > 2.0 * e) ? max(w.z + 2.0 * e, 0.0) : 1e15);
    bool  ex = (t.x >= t.y) && (t.x >= t.z);                 // drop the largest t: keep the two smallest
    bool  ey = (!ex) && (t.y >= t.z);
    vec3  tt = ex ? vec3(0.0, t.y, t.z) : (ey ? vec3(t.x, 0.0, t.z) : vec3(t.x, t.y, 0.0));
    float h  = length(tt);
    float wm = max(w.x, max(w.y, w.z));                      // -wm = depth below the nearest outer face
    vec3  no = (w.x >= max(w.y, w.z)) ? vec3(1.0, 0.0, 0.0)
             : ((w.y >= w.z)          ? vec3(0.0, 1.0, 0.0) : vec3(0.0, 0.0, 1.0));
    if( h < -wm )
    {
        if( h > 0.0 ) return vec4(-h, s * (-tt / h));
        // Exactly on an inner face (FP32 surface shell): axis normal of the active transverse plane.
        vec3 ni = (gt.x >= max(gt.y, gt.z)) ? vec3(1.0, 0.0, 0.0)
                : ((gt.y >= gt.z)           ? vec3(0.0, 1.0, 0.0) : vec3(0.0, 0.0, 1.0));
        return vec4(0.0, s * sc * ni);
    }
    return vec4(wm, s * no);
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
