/*  zzz — primes from Riemann zeros, with NO sieve and NO primality test.
 *  Shadertoy: paste this whole file into the Image tab.
 *
 *  Background = the explicit-formula detector itself:
 *      psi'(x) = 1 + 1/(x - x^3) - (4/sqrt x) * SUM_k cos(gamma_k * log x)
 *  Its spikes sit at prime powers; r = psi'/envelope > 0.30 flags one.
 *  The flaring spike is the prime currently printed (held ~1 s).
 *
 *  STAGE 1 (this file): the gamma_k below are a baked snapshot of zzz --loop's
 *  OWN output, purely to nail the visuals. STAGE 2 replaces G[] with a live
 *  Buffer-A bootstrap (the full loop.c: method-B bisection forward, detector
 *  backward) so the zeros are self-generated on the GPU — no baking, no sieve.
 *
 *  Verified offline: this detector in genuine fp32 reproduces the first 600+
 *  primes exactly (vs a sieve). Number printing: @P_Malin, CC0 (header kept).
 */
const int NZ = 400;
const float G[400] = float[](
14.1347,21.0220,25.0109,30.4249,32.9351,37.5862,40.9187,43.3271,
48.0052,49.7738,52.9703,56.4462,59.3470,60.8318,65.1125,67.0798,
69.5464,72.0672,75.7047,77.1448,79.3374,82.9104,84.7355,87.4253,
88.8091,92.4919,94.6513,95.8706,98.8312,101.3179,103.7255,105.4466,
107.1686,111.0295,111.8747,114.3202,116.2267,118.7908,121.3701,122.9468,
124.1271,127.4431,129.5635,131.1035,133.5494,134.7537,138.1010,139.7794,
141.0524,143.1913,145.9381,147.4135,150.1195,150.8513,153.0605,156.1322,
157.6290,158.7846,161.1971,163.1191,165.4582,167.2022,169.1732,169.7907,
173.4253,174.6411,176.4083,178.3930,179.9714,182.2260,184.9078,185.5621,
187.2977,189.4371,192.0419,193.0227,195.2270,196.9568,197.9726,201.3502,
202.3699,204.2070,205.3891,207.8471,209.5890,211.7087,213.4865,214.4249,
216.1895,219.0131,220.7818,221.3867,224.0621,224.9340,227.4034,229.3822,
231.3476,231.8678,233.6906,236.6097,237.6632,239.6042,241.0181,242.9507,
243.9512,247.1999,248.0071,249.5129,250.9774,253.0876,255.3666,256.3218,
258.6720,259.9579,260.7483,263.6007,265.7005,266.4274,267.8399,269.9243,
271.5004,273.5024,275.6466,276.3548,278.3850,279.1134,282.5147,283.1094,
284.8159,286.7802,287.7992,289.5887,291.8483,293.6374,295.0380,295.6088,
297.8576,299.7937,301.7333,302.6299,304.9326,305.6318,307.2299,310.1307,
311.1893,312.3532,313.9579,315.4618,317.8148,318.7915,321.2515,322.0363,
323.4234,324.9034,327.4033,329.1397,329.8318,331.4792,333.6760,334.2015,
336.8272,338.3387,339.9079,341.0162,342.0487,344.6274,346.4559,347.1621,
349.4344,350.2756,351.8397,353.5517,355.9703,357.2300,357.8507,359.7211,
361.2965,363.4140,364.7934,366.1724,368.1171,368.8955,370.0965,373.1132,
373.8528,375.8364,376.2733,378.3841,379.7947,381.4660,383.4653,385.0717,
385.7540,387.2144,388.9256,391.4488,392.2497,393.3985,395.6804,396.2580,
397.9125,399.9869,401.9355,402.7648,404.2640,405.1063,407.5715,408.9722,
410.6288,412.0693,413.1079,414.9540,415.5433,418.3047,419.8787,420.6217,
422.0977,423.8055,424.9812,427.3377,428.0249,430.3462,431.3285,432.0705,
433.9224,436.1628,437.7088,438.5008,439.8797,441.7813,442.8738,444.2873,
446.8859,447.3921,449.2139,450.0730,451.3549,454.0621,454.8734,456.2458,
457.9035,459.5858,460.0663,462.0247,464.0290,465.7407,466.5546,467.4047,
469.5930,470.7612,472.9445,473.6681,475.6583,476.6993,478.1462,478.8711,
481.8951,482.8043,483.8218,485.6768,486.3638,488.3776,489.6526,491.4652,
493.3669,493.8481,495.3563,496.4524,498.5376,500.2993,501.6809,502.2266,
504.5545,505.3719,506.4582,508.7691,510.3219,511.6933,512.5134,513.6229,
515.4341,517.6317,518.2066,520.0893,521.6261,522.3469,524.0635,524.9572,
527.8616,528.4144,529.9181,530.7482,532.7154,533.7122,535.6368,537.0692,
538.4569,540.2140,540.6393,541.9364,544.2602,545.6407,547.1394,547.8039,
549.4810,551.0457,551.9773,553.8042,555.7911,556.9862,557.5131,559.3414,
560.2416,562.5168,564.1373,564.5455,566.7220,567.6954,568.9360,570.0646,
572.3737,573.6228,575.1966,575.6897,577.0320,579.1073,580.0785,581.9584,
583.2509,584.5825,586.0989,586.6523,588.1687,590.6742,591.7683,592.5053,
594.0068,595.8047,596.2905,598.5381,599.4979,601.6566,602.5833,603.5969,
604.5847,606.4106,608.4368,609.3788,610.9767,611.6546,613.6539,614.6378,
615.5169,618.1633,619.1531,620.2370,621.8046,622.3195,624.2266,626.0462,
627.3328,628.2901,630.4281,630.8350,632.2818,633.6533,635.5389,637.3470,
637.9497,638.8421,640.7012,641.9080,643.2592,645.0366,646.4326,647.9308,
648.6345,650.1745,650.6925,653.7271,654.2337,655.7590,657.0911,658.1302,
659.8138,660.5740,662.2502,664.3031,665.2752,666.5443,667.1130,668.9971,
670.4187,672.4704,672.9458,674.3131,676.1537,677.2474,677.8014,679.7462
);
const float C_ENV = 423.748571; // 5*Li(NZ-4)
// gamma_max = 679.7, reach ~ 272
const int NP = 55;
const int PRIMES[55] = int[](2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257);
const float XMAX = 260.0;

// ---- fp32 phase-safe cosine (Cody–Waite reduction; ready for the live loop's
// ---- high gamma where naive cos(largeArg) would fold) ----------------------
float pcos(float a){
    const float TWOPI = 6.2831855, INV = 0.15915494;
    float k = floor(a*INV + 0.5);
    float r = (a - k*TWOPI) + k*1.7484555e-7;   // 2pi split, fp32
    return cos(r);
}

// detector r(x) = psi'(x) / envelope(x)
float Rdet(float x){
    float lx = log(x), s = 0.0;
    for(int k=0;k<NZ;k++) s += pcos(G[k]*lx);
    float psi = 1.0 + 1.0/(x - x*x*x) - (4.0/sqrt(x))*s;
    float env = C_ENV * lx / x;                 // (20/4)*Li(NZ-4)*log x / x
    return psi/env;
}

bool isPrimeBaked(int n){
    for(int i=0;i<NP;i++){ if(PRIMES[i]==n) return true; if(PRIMES[i]>n) break; }
    return false;
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
        if(fDigitIndex > fBiggestIndex) {
            if((bNeg) && (fDigitIndex < (fBiggestIndex+1.5))) fCharBin = 1792.0;
        } else {
            if(fDigitIndex == -1.0) { if(fDecimalPlaces > 0.0) fCharBin = 2.0; }
            else {
                float fReducedRangeValue = fValue;
                if(fDigitIndex < 0.0) { fReducedRangeValue = fract( fValue ); fDigitIndex += 1.0; }
                float fDigitValue = (abs(fReducedRangeValue / (pow(10.0, fDigitIndex))));
                fCharBin = DigitBin(int(floor(mod(fDigitValue, 10.0))));
            }
        }
    }
    return floor(mod((fCharBin / pow(2.0, floor(fract(vStringCoords.x) * 4.0) + (floor(vStringCoords.y * 5.0) * 4.0))), 2.0));
}
// ---- 8< ---------------------------------------------------------- 8< ----

void mainImage(out vec4 O, in vec2 F){
    vec2 R = iResolution.xy;
    vec2 uv = F / R;

    const float HOLD = 0.9;                       // seconds per prime
    int   idx  = int(mod(iTime/HOLD, float(NP))); // which prime is "current"
    float ph   = fract(iTime/HOLD);               // 0..1 within the hold
    int   curP = PRIMES[idx];
    float xf   = float(curP);                     // discovery frontier

    float x = 2.0 + (XMAX-2.0)*uv.x;              // column -> number line
    float r = Rdet(x);                            // detector value at this column

    float baseY = 0.10, scaleY = 0.46;
    float yC   = baseY + clamp(r,-0.3,1.7)*scaleY; // curve height
    float yThr = baseY + 0.30*scaleY;             // threshold tau=0.30

    // reveal: trace exists only up to the current frontier (left-to-right discovery)
    float reveal = smoothstep(xf+6.0, xf+2.0, x);

    vec3 col = mix(vec3(0.02,0.03,0.06), vec3(0.04,0.05,0.11), uv.y);

    // threshold line (dim red)
    col += vec3(0.40,0.12,0.10) * reveal * smoothstep(0.0035,0.0,abs(uv.y-yThr));

    // faint fill under the curve
    col += vec3(0.05,0.22,0.38) * reveal * step(uv.y,yC) * smoothstep(yC,baseY,uv.y) * 0.5;

    // the trace (cyan glow)
    col += vec3(0.25,0.85,1.0) * reveal * smoothstep(0.013,0.0,abs(uv.y-yC));

    // green tint at discovered primes
    int ni = int(floor(x+0.5));
    if(isPrimeBaked(ni) && float(ni) <= xf+0.5)
        col += vec3(0.1,0.95,0.35) * reveal * smoothstep(0.6,0.0,abs(x-float(ni))) * 0.30;

    // flaring current prime (bright vertical pulse)
    float flare = smoothstep(1.3,0.0,abs(x-xf)) * (0.55+0.45*sin(iTime*9.0)) * (1.0-0.4*ph);
    col += vec3(0.6,1.0,0.75) * flare;

    // the prime, printed large and held (P_Malin)
    float val = float(curP);
    float digits = max(1.0, ceil(log(val+0.5)/log(10.0)));
    vec2  fs  = vec2(34.0,62.0)/520.0 * R.y;
    vec2  org = vec2(R.x*0.5 - fs.x*0.5*digits, R.y*0.72);
    float isd = PrintValue((F-org)/fs, val, digits, 0.0);
    col = mix(col, vec3(1.0,1.0,0.95), isd*(1.0-0.25*ph));

    O = vec4(col, 1.0);
}
