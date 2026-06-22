/*  zzz live loop — IMAGE tab.   Published: https://www.shadertoy.com/view/s3BSzD
 *    iChannel0 = Buffer A   (loop state: zeros, prime-indicator, control)
 *    iChannel1 = Buffer B   (per-column trace cache r(x), stored in row 0)
 *  The heavy detector sum now lives in Buffer B; here we just sample it.
 */
bool isPrimeLive(int x){ return texelFetch(iChannel0, idx2(x),0).g > 0.5; }
float traceAt(int px){            // cached r(x) at screen column px (Buffer B, row 0)
    px = clamp(px, 0, int(iResolution.x)-1);
    return texelFetch(iChannel1, ivec2(px,0), 0).r;
}

// ---- 8< ---- GLSL Number Printing - @P_Malin ---- 8< ----
// Creative Commons CC0 1.0 Universal (CC-0)  https://www.shadertoy.com/view/4sBSWW
float DigitBin( const int x ){
    return x==0?480599.0:x==1?139810.0:x==2?476951.0:x==3?476999.0:x==4?350020.0:x==5?464711.0:x==6?464727.0:x==7?476228.0:x==8?481111.0:x==9?481095.0:0.0;
}
float PrintValue( vec2 vStringCoords, float fValue, float fMaxDigits, float fDecimalPlaces ){
    if ((vStringCoords.y < 0.0) || (vStringCoords.y >= 1.0)) return 0.0;
    bool bNeg = ( fValue < 0.0 ); fValue = abs(fValue);
    float fLog10Value = log2(abs(fValue)) / log2(10.0);
    float fBiggestIndex = max(floor(fLog10Value), 0.0);
    float fDigitIndex = fMaxDigits - floor(vStringCoords.x);
    float fCharBin = 0.0;
    if(fDigitIndex > (-fDecimalPlaces - 1.01)) {
        if(fDigitIndex > fBiggestIndex) { if((bNeg) && (fDigitIndex < (fBiggestIndex+1.5))) fCharBin = 1792.0; }
        else {
            if(fDigitIndex == -1.0) { if(fDecimalPlaces > 0.0) fCharBin = 2.0; }
            else {
                float v = fValue;
                if(fDigitIndex < 0.0) { v = fract( fValue ); fDigitIndex += 1.0; }
                fCharBin = DigitBin(int(floor(mod((abs(v / pow(10.0, fDigitIndex))), 10.0))));
            }
        }
    }
    return floor(mod((fCharBin / pow(2.0, floor(fract(vStringCoords.x) * 4.0) + (floor(vStringCoords.y * 5.0) * 4.0))), 2.0));
}
// ---- 8< ---------------------------------------------------------- 8< ----

// point-to-segment distance (pixels); round joins/caps fall out of the vertex distance
float segDist(vec2 p, vec2 a, vec2 b){
    vec2 pa=p-a, ba=b-a;
    float h=clamp(dot(pa,ba)/dot(ba,ba), 0.0, 1.0);
    return length(pa - ba*h);
}

void mainImage(out vec4 O, in vec2 F){
    vec2 R=iResolution.xy, uv=F/R;
    // window cached by Buffer B (.g=pc, .b=xLo, .a=xHi) — one fetch instead of a
    // per-pixel primeWindow scan; also makes the markers match the trace exactly.
    vec4 win = texelFetch(iChannel1, ivec2(0,0), 0);
    int pc = int(win.g+0.5); float pcf = float(pc);
    float xLo = win.b, xHi = win.a;
    float x = xLo + (xHi-xLo)*uv.x;                       // this column's number-line position

    float baseY=0.10, scaleY=0.46, yThr=baseY+0.30*scaleY;
    int px = int(F.x);
    float yC = baseY + clamp(traceAt(px),  -0.3,1.7)*scaleY;        // uv units (fill-under)
    float yLp=(baseY+clamp(traceAt(px-1),-0.3,1.7)*scaleY)*R.y;     // neighbour heights, pixels
    float yMp= yC*R.y;
    float yRp=(baseY+clamp(traceAt(px+1),-0.3,1.7)*scaleY)*R.y;
    float reveal = smoothstep(pcf+3.0, pcf-0.5, x);     // trace exists up to the current prime

    vec3 col = mix(vec3(0.02,0.03,0.06), vec3(0.04,0.05,0.11), uv.y);
    col += vec3(0.40,0.12,0.10)*reveal*smoothstep(0.0035,0.0,abs(uv.y-yThr));
    col += vec3(0.05,0.22,0.38)*reveal*step(uv.y,yC)*smoothstep(yC,baseY,uv.y)*0.5;
    // round stroke: min distance (px) to the two adjacent polyline segments. The shared
    // apex is a vertex, so its distance field is radial => peaks round off, joins included.
    vec2 Pp = vec2(0.0, F.y);
    float dpx = min(segDist(Pp, vec2(-1.0,yLp), vec2(0.0,yMp)),
                    segDist(Pp, vec2( 0.0,yMp), vec2(1.0,yRp)));
    float line = 1.0 - smoothstep(2.0, 3.6, dpx);       // ~2.8px stroke, soft round edge
    float glow = exp(-dpx*dpx/240.0);                   // soft round halo
    float band = 1.0 - 0.80*smoothstep(0.60, 0.80, uv.y);   // dim trace where the number sits
    col += vec3(0.25,0.85,1.0)*reveal*(line + 0.22*glow)*band;

    int ni=int(floor(x+0.5));
    if(ni>=2 && float(ni)<=pcf+0.5 && isPrimeLive(ni))   // cheap tests gate the texelFetch
        col += vec3(0.1,0.95,0.35)*reveal*smoothstep(0.6,0.0,abs(x-float(ni)))*0.30;

    col += vec3(0.6,1.0,0.75)*smoothstep(1.3,0.0,abs(x-float(pc)))*(0.55+0.45*sin(iTime*9.0));

    float val=float(pc), digits=max(1.0,ceil(log(val+0.5)/log(10.0)));
    vec2 fs=vec2(34.0,62.0)/520.0*R.y;
    vec2 org=vec2(R.x*0.5 - fs.x*0.5*digits, R.y*0.72);
    col = mix(col, vec3(1.0,1.0,0.95), PrintValue((F-org)/fs, val, digits, 0.0));

    O=vec4(col,1.0);
}
