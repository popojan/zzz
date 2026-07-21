/*  zzz live loop — BUFFER B tab (trace cache).  Set its iChannel0 = Buffer A.
 *
 *  Computes the detector r(x) ONCE per screen column (row 0 only), so the Image
 *  pass samples it with a texelFetch instead of a per-pixel cos sum. The zero
 *  count scales with the window so spikes reach the threshold at any depth
 *  (matching Buffer A's prime decision). Also stashes the scroll window in .gba.
 */
const int NTRACE = 8000;  // trace zero-count cap (matches Buffer A's NLIM)

// n = zeros to sum, chosen by the caller to resolve the current window.
float RdetVis(float x, int n){
    float lx=log(x), s=0.0;
    for(int k=0;k<NTRACE;k++){ if(k>=n) break; s += pcos(texelFetch(iChannel0, idx2(k),0).r * lx); }
    float psi = 1.0 + 1.0/(x - x*x*x) - (4.0/sqrt(x))*s;
    return psi/envf(x, n);
}

void mainImage(out vec4 O, in vec2 F){
    if(int(F.y) != 0){ O = vec4(0.0); return; }          // only row 0 holds the trace
    int pc; float xLo, xHi; primeWindow(iChannel0, iTime, pc, xLo, xHi);
    int nz = int(texelFetch(iChannel0, idx2(0),0).b + 0.5);
    // enough zeros to resolve THIS window (reach ~0.5*n; want x well inside it), so the
    // trace's spikes reach the threshold wherever the window sits -- matching the prime
    // DECISION in Buffer A. Scales with depth, so you only pay for what you're viewing.
    int n = int(2.8*xHi) + 128; if(n>nz) n=nz; if(n<8) n=8;
    float x = xLo + (xHi-xLo) * (F.x / iResolution.x);   // this column's number-line position
    // .r = the trace r(x); .gba = the (pixel-independent) window, so the Image pass
    // reads pc/xLo/xHi from one texel instead of re-running primeWindow per pixel.
    O = vec4(RdetVis(x, n), float(pc), xLo, xHi);
}
