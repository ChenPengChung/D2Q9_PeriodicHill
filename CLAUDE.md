# D2Q9 Periodic Hill å°ˆæ??‹ç™¼è¨˜é?

## å°ˆæ?æ¦‚è¿°
ä½¿ç”¨ Lattice Boltzmann Method (LBM) æ¨¡æ“¬ Periodic Hill æµå ´
- æ¨¡å?ï¼šD2Q9 + MRT ç¢°æ?ç®—å?
- ?·è«¾?¸ï?Re = 50
- ç¶²æ ¼ï¼šNY=128, NZ=64 (??buffer ??NY6=135, NZ6=70)
- ?å??»ç¶²?¼ï?Z ?¹å?ä½¿ç”¨ tanh ?†ä?

---

## 2026-02-05 ±Ò°Ê/­«±Ò¡]Àô¹ÒÅÜ¼Æ³Æ§Ñ¡^
±`¨£¨SÅª¨ìªì©lÀÉªº­ì¦]
1. `RESTART_VTK` ¨S¦³³]©w
2. `OUTPUT_DIR` ¨S¦³«ü¨ì¦³ `latest_vtk.txt` ªº¸ê®Æ§¨

¥¿½T±Ò°Ê¤è¦¡¡]PowerShell¡^
```powershell
# ¥Î latest_vtk.txt
$env:OUTPUT_DIR="output_test"
$env:RESTART_VTK="LATEST"
.\main.exe

# ª½±µ«ü©wÀÉ®×
$env:RESTART_VTK="output/flow_036000.vtk"
.\main.exe
```

§Ö³tÀË¬d¥Ø«eÀô¹ÒÅÜ¼Æ
```powershell
$env:OUTPUT_DIR
$env:RESTART_VTK
```


## 2025-01-14 ?²åº¦æª¢æŸ¥

### å·²å??æ¨¡çµ?
| æª”æ? | ?€??| èªªæ? |
|------|------|------|
| `variables.h` | ??| ?©ç??ƒæ•¸?‡ç¶²?¼é?ç½?|
| `model.h` | ??| å±±å¡å¹¾ä??½æ•¸ HillFunction() |
| `MRT_Matrix.h` | ??| MRT ?©é™£å®šç¾©ï¼?025-01-15 å·²ä¿®å¾©ï? |
| `MRT_Process.h` | ? ï? | m_vector å·¨é?ï¼Œé?è¼¯å?è£œå? |
| `initializationTool.h` | ??| ?å??»ç¶²?¼ã€Lagrange ?’å€?|
| `InterpolationHillISLBM.h` | ??| 7é»æ??¼å·¨?†ï?ç¶“æª¢?¥ç„¡èªæ??¯èª¤ï¼?|

### å·²ä¿®å¾©å?é¡?(2025-01-15)

1. **MRT_Matrix.h** - äºŒç¶­????„å?ç¼ºå??—è??†é?ç¬?   - `Matrix` å·¨é?ï¼šå???`{...}` å¾Œè?ä¸Šé€—è?
   - `Matrix_Inverse` å·¨é?ï¼šå???`{...}` å¾Œè?ä¸Šé€—è?

2. **InterpolationHillISLBM.h** - ç¶“æª¢?¥ç¬¬50è¡Œè?æ³•æ­£ç¢ºï??¡é?ä¿®æ”¹

### å¾…å¯¦ä½œé???
- [ ] ä¸»ç?å¼?main.cu
- [ ] MRT ç¢°æ??‹ç?å®Œæ•´?è¼¯
- [ ] Streaming æ­¥é?
- [ ] ?Šç?æ¢ä»¶?•ç? (?±æ??§é??Œã€å???
- [ ] ?å??–æ???- [ ] çµæ?è¼¸å‡º

---

## ?‹ç™¼ç­†è?

- `y_h` ??main.cu ä¸­ç??¨å?è®Šæ•¸ï¼Œinitialization.h ?…å«å¾Œå¯?´æ¥è®€å¯?- ?…å«?†å??€?¨å…¨?Ÿå?ç¾©ä?å¾?
### å¾…å?ç¾©å…¨?Ÿè???(main.cu)

ä»¥ä?è®Šæ•¸??`initialization.h` ä¸­ä½¿?¨ï??€??`main.cu` ä¸­å®£?Šï?

```cpp
double rho[NY6 * NZ6];           // å¯†åº¦??double v[NY6 * NZ6];             // Y?¹å??Ÿåº¦
double w[NY6 * NZ6];             // Z?¹å??Ÿåº¦
double f[9][NY6 * NZ6];          // ?†ä??½æ•¸ (D2Q9)
double Force[2];                 // å¤–å???(Fy, Fz)
double xi_h[NZ6];                // ?¡å?æ¬¡å?Zåº§æ?
```

**æ³¨æ?**ï¼š`#include "initialization.h"` å¿…é??¾åœ¨?™ä??¨å?è®Šæ•¸å®??ä¹‹å?

---

## 2025-01-14 BFL ?Šç?æ¢ä»¶?¤æ–·å¼?
### å±±ä?å¹¾ä?çµæ?
- å·¦å?ä¸˜ï?`y ??[0, 54/28]` ??`[0, 1.9286]`
- ?³å?ä¸˜ï?`y ??[LY - 54/28, LY]` ??`[7.07, 9.0]`
- ?±æ??§é??Œé€?¥å¾Œå½¢?å??´å±±ä¸?
### ?°å??¤æ–·?½æ•¸ (initializationTool.h)

```c
//5. ?¤æ–·?¯å¦??+Y ?¹å??Šç?è¨ˆç?é»?(?å?å·¦å?ä¸?
#define HillHalfWidth (54.0/28.0)

#define IsBoundary_LeftHill_PlusY(y, z) \
(                                        \
    (y) >= 0.0 &&                        \
    (y) + minSize <= HillHalfWidth &&    \
    (z) > HillFunction(y) &&             \
    (z) <= HillFunction((y) + minSize)   \
)

//6. ?¤æ–·?¯å¦??-Y ?¹å??Šç?è¨ˆç?é»?(?å??³å?ä¸?
#define IsBoundary_RightHill_MinusY(y, z) \
(                                          \
    (y) <= LY &&                           \
    (y) - minSize >= (LY - HillHalfWidth) && \
    (z) > HillFunction(y) &&               \
    (z) <= HillFunction((y) - minSize)     \
)
```

### ?è¼¯èªªæ?
- å·¦å?ä¸˜ï?è¨ˆç?é»åœ¨æµé??€ï¼Œå? +Y ?¹å?ç§»å? minSize ?ƒç¢°?°æ›²??- ?³å?ä¸˜ï?è¨ˆç?é»åœ¨æµé??€ï¼Œå? -Y ?¹å?ç§»å? minSize ?ƒç¢°?°æ›²??
---

## 2026-01-16 å®Œæ•´ç¨‹å?ç¢¼å¯©??
### ç¨‹å?ç¢¼ç?æ§‹ç¸½è¦?
å°ˆæ??…å« 7 ?‹æ??­æ?ï¼?1. `variables.h` - ?©ç??ƒæ•¸å®šç¾©
2. `model.h` - å±±å¡å¹¾ä?æ¨¡å?
3. `MRT_Matrix.h` - MRT è½‰æ??©é™£
4. `MRT_Process.h` - MRT ç¢°æ??•ç?
5. `initializationTool.h` - ?å??–å·¥?·å‡½??6. `InterpolationHillISLBM.h` - ISLBM ?’å€¼å·¨??7. `initialization.h` - ?…è?æª?
### ?„æ¨¡çµ„è©³ç´°å???
#### ??variables.h (æ­?¸¸)
- å®šç¾©?©ç??ƒæ•¸ï¼šRe=50, LY=9.0, LZ=3.036
- ç¶²æ ¼è¨­å?ï¼šNY=128, NZ=64

- Buffer ç¶²æ ¼ï¼šNY6=135, NZ6=70
- CFL=0.6, minSize è¨ˆç?æ­?¢º
- tau=0.6833, omega_2 ??omega_7 å®šç¾©æ­?¢º
- Uref = Re*niu

#### ??model.h (æ­?¸¸)
- HillFunction() å¯¦ä?å®Œæ•´
- å·¦å?ä¸˜ï?6 æ®µä?æ¬¡å??…å? (y ??[0, 54/28])
- ?³å?ä¸˜ï?6 æ®µä?æ¬¡å??…å? (y ??[LY-54/28, LY])
- ?±æ??§é??Œè??†æ­£ç¢?
#### ? ï? MRT_Matrix.h (?€æª¢æŸ¥)
- **?é? 1**ï¼šç¬¬ 10 è¡Œå?ç¾?`s1 = omega_1` ä½?`omega_1` ?ªåœ¨ variables.h ä¸­å?ç¾?  - ?¹æ?è¨»è§£?‰è©²??free parameterï¼Œå»ºè­°æ”¹??`s1 = 1.0`
- **?é? 2**ï¼šç¬¬ 9??1??4??6 è¡Œç¼ºå°‘å???- Matrix M[9][9] å®šç¾©æ­?¢ºï¼ˆå·²??2025-01-15 ä¿®å¾©ï¼?- Matrix_Inverse M_I[9][9] å®šç¾©æ­?¢ºï¼ˆå·²??2025-01-15 ä¿®å¾©ï¼?
#### ??MRT_Process.h (?´é??¯èª¤)
- **?é? 1**ï¼šé™£?—ç´¢å¼•è???  - ?„è?ä½¿ç”¨ `M[i][9]` ä½?D2Q9 ?ªæ? 9 ?‹æ–¹??(ç´¢å? 0-8)
  - æ­?¢º?‰ç‚º `M[i][8]`
- **?é? 2**ï¼šå?è¡Œç¼ºå°‘å???- **?é? 3**ï¼šè??¸å‘½?éŒ¯èª?  - ä½¿ç”¨ F9_in ä½?D2Q9 ?ªæ? F0-F8
  - ?‰æ”¹??F8_in
- **?é? 4**ï¼šm_vector å·¨é??è¼¯ä¸å???  - ?…è?ç®?momentï¼Œæœªå¯¦ä?ç¢°æ?é¬†å??Œé€†è???
#### ??initializationTool.h (?Ÿèƒ½å®Œæ•´ï¼Œæ?å°éŒ¯èª?
- **?é?**ï¼šç¬¬ 6 è¡?`HillHalfWidth` å®šç¾©?¯èª¤
  - ?¾åœ¨ï¼š`#define HillHalfWidth 9.0*(54.0/28.0)` ??17.36
  - æ­?¢ºï¼š`#define HillHalfWidth (54.0/28.0)` ??1.93

**å·²å¯¦ä½œå???*ï¼?- ??HillFunction_Inverse_Left/Right (äºŒå??œå?æ³?
- ??tanhFunction å·¨é? (?å??»ç¶²??
- ??GetNonuniParameter() (è¨ˆç?ä¼¸ç¸®?ƒæ•¸ aï¼Œä½¿?©ç?ç©ºé?è¨ˆç?é»ç??€å°è??¢ç??¼æ™¶?¼å¤§å°?minSize)
  - **?®ç?**ï¼šèª¿??tanh ?å??»ç¶²?¼ç?ä¼¸ç¸®?ƒæ•¸ aï¼Œç¢ºä¿æ?å°ç¶²?¼é?è·?= minSize (?¶æ ¼å¤§å?)
  - **?¹æ?**ï¼šä??†æ??œå? a ??(0.1, 1.0)
  - **?œéµ?¤æ–·**ï¼?    ```c
    x_temp[0] = tanhFunction(total, minSize, a_mid, 0, (NZ6-7));
    x_temp[1] = tanhFunction(total, minSize, a_mid, 1, (NZ6-7));
    dx = x_temp[1] - x_temp[0];  // ç¬?0 ?‡ç¬¬ 1 ?‹è?ç®—é??„é?è·?    if( dx - minSize >= 0.0 )    // dx >= minSize ??a ?¯ä»¥?´å¤§
        a_temp[0] = a_mid;
    else
        a_temp[1] = a_mid;
    // ?¶æ?æ¢ä»¶ï¼š|dx - minSize| < 1e-14
    ```
  - **?©ç??ç¾©**ï¼šè?å£é¢?„è?(k=0)?„ç¶²?¼é?è·ç²¾ç¢ºç??¼æ™¶?¼å¤§å°ï?ç¢ºä? ISLBM ?’å€¼ç²¾åº?- ??Lagrange_6th() (6 ?æ???
- ??GetParameter_6th() (?è?ç®—æ??¼æ???
- ??IsLeftHill_Boundary_yPlus() (å·¦ä? +Y ?Šç??¤æ–·)
- ??IsRightHill_Boundary_yMinus() (?³ä? -Y ?Šç??¤æ–·)
- ??IsLeftHill_Boundary_Diagonal45() (å·¦ä? 45Â° ?Šç??¤æ–·)
- ??IsRightHill_Boundary_Diagonal135() (?³ä? 135Â° ?Šç??¤æ–·)
- ??IsLeftHill_Boundary_DiagonalMinus45() (å·¦ä? -45Â° ?Šç??¤æ–·)
- ??IsRightHill_Boundary_DiagonalMinus135() (?³ä? -135Â° ?Šç??¤æ–·)
- ??Left_q_yPlus() (å·¦ä? +Y BFL è·é›¢)
- ??Right_q_yMinus() (?³ä? -Y BFL è·é›¢)
- ??Left_q_Diagonal45() (å·¦ä? 45Â° BFL è·é›¢)
- ??Right_q_Diagonal135() (?³ä? 135Â° BFL è·é›¢)
- ??Left_q_DiagonalMinus45() (å·¦ä? -45Â° BFL è·é›¢)
- ??Right_q_DiagonalMinus135() (?³ä? -135Â° BFL è·é›¢)

#### ??InterpolationHillISLBM.h (æ­?¸¸)
- å®šç¾© 7 é»æ??¼å·¨??Intrpl7()
- F0 ??F8 ??2D ?’å€¼å·¨??- ?¯æ´?å??»ç¶²?¼ç? Lagrange ?’å€?- èªæ?æª¢æŸ¥?¡èª¤ï¼?025-01-15 å·²ç¢ºèªï?

#### ??initialization.h (æ­?¸¸)
- ç°¡å–®?…è?æª”ï??…å« initializationTool.h

### ?€ä¿®å¾©?„éŒ¯èª¤æ???
#### é«˜å„ª?ˆç? (?»æ–·?§éŒ¯èª?
1. **MRT_Process.h ???è¶Šç?** - å¿…é?ç«‹å³ä¿®å¾©

---

## 2025-01-20 ?’å€¼è??å??»ç¶²?¼å¿«?–èª¿?´ï?AI ?”åŠ©ï¼?
### ?°å?/è®Šæ›´
- ?°å??¨å? `nonuni_a`ï¼š`globalVariables.h` å®??ï¼Œ`main.cpp` å®šç¾©ï¼Œ`GenerateMesh_Z` è¨ˆç?ä¸€æ¬¡å??±ç”¨ï¼Œé¿?å?æ¬¡å‘¼??`GetNonuniParameter()`??- `GenerateMesh_Z` ?§éƒ¨ï¼šä½¿??`nonuni_a` ?Ÿæ? `xi_h`?`z_global`ï¼Œé?è¼¯ä?è®Šï??…ç§»?¤é?è¤‡æ? a??- ?°å? `BuildXiWeights`ï¼ˆ`initialization.h`ï¼‰ï?å°å?ä¸€ `(j,k)` ä¸€æ¬¡æ€§è?ç®—ä??—ç? `RelationXi`/`GetParameter_6th2`ï¼Œå¯«?¥å???`XiPara*_h[?][index_xi + row*NZ6]`ï¼›å…§??y ?±æ??•ç?ï¼Œä½¿?¨é??ˆç? `nonuni_a`ï¼Œé¿??8 ?é?ç®—ã€?- `GetIntrplParameter_Xi`ï¼šå…«?‹é€Ÿåº¦?¹å??¹å‘¼??`BuildXiWeights`ï¼Œå? `GetXiParameter` ä¿ç?ä½†ä??ä½¿?¨ã€?- ä¿®æ­£ `GetParameter_6th2` ç¬¬ä??—ä½¿??`pos_z` ?„éŒ¯èª¤ï?çµ±ä??¨å?ä¸€çµ?`pos_z2`??
#### ä¸»è?ç¨‹å?ç¢¼ç?æ®?
`globalVariables.h` / `main.cpp`ï¼šæ–°å¢å…±?¨ç??å??»å???`nonuni_a`??```cpp
// globalVariables.h
extern double nonuni_a;

// main.cpp
double nonuni_a = 0.0;
```

`initialization.h`ï¼šåœ¨?¢ç¶²?¼æ?è¨ˆç?ä¸€æ¬?aï¼Œä?å¾Œéƒ½?¨å?ä¸€?¸å€¼ã€?```cpp
void GenerateMesh_Z() {
    ...
    nonuni_a = GetNonuniParameter();      // ?ªç?ä¸€æ¬?    const double a = nonuni_a;
    xi_h[k] = tanhFunction(LXi, minSize, a, (k-3), (NZ6-7)) - minSize/2.0;
    ...
}
```



`initialization.h`ï¼š`BuildXiWeights` å°‡å?ä¸€ `(j,k)` ?„ä??—æ??ä?æ¬¡ç?å®Œï?ä¾?8 ?‹æ–¹?‘å…±?¨ã€?```cpp
inline void BuildXiWeights(double* XiPara_h[7], double pos_z, double pos_y,
                           int index_xi, int j, int k) {
    int jj = j;
    if (jj < 3) jj += NY6 - 7;
    if (jj > NY6 - 4) jj -= (NY6 - 7);
    double L = LZ - HillFunction(pos_y) - minSize;
    double pos_xi = pos_z - (HillFunction(pos_y) + minSize/2.0);
    double j_cont = Inverse_tanh_index(pos_xi, L, minSize, nonuni_a, (NZ6-7));

    auto fillRow = [&](int rowOffset, int storeRow) {
        int rowIdx = jj + rowOffset;
        if (rowIdx < 0) rowIdx += NY6;
        if (rowIdx >= NY6) rowIdx -= NY6;
        double H = HillFunction(y_global[rowIdx]);
        double Lrow = LZ - H - minSize;
        double pos_z_row = tanhFunction(Lrow, minSize, nonuni_a, j_cont, (NZ6-7)) - minSize/2.0;
        double rel[7];
        RelationXi(j_cont, Lrow, minSize, nonuni_a, (NZ6-7), rel);
        GetParameter_6th2(XiPara_h, pos_z_row, rel, storeRow, index_xi);
    };
    fillRow(-3,0); fillRow(-2,1); fillRow(-1,2); fillRow(0,3);
    fillRow(1,4);  fillRow(2,5);  fillRow(3,6);
}
```

`GetIntrplParameter_Xi`ï¼šå…«?‹æ–¹?‘å…±??`BuildXiWeights`ï¼Œé¿?é?ç®—ã€?```cpp
int idx = j * NZ6 + k;
BuildXiWeights(XiParaF1_h, z_global[idx],         y_global[j]-minSize, idx, j, k);
BuildXiWeights(XiParaF2_h, z_global[idx]-minSize, y_global[j],         idx, j, k);
BuildXiWeights(XiParaF3_h, z_global[idx],         y_global[j]+minSize, idx, j, k);
BuildXiWeights(XiParaF4_h, z_global[idx]+minSize, y_global[j],         idx, j, k);
BuildXiWeights(XiParaF5_h, z_global[idx]-minSize, y_global[j]-minSize, idx, j, k);
BuildXiWeights(XiParaF6_h, z_global[idx]-minSize, y_global[j]+minSize, idx, j, k);
BuildXiWeights(XiParaF7_h, z_global[idx]+minSize, y_global[j]+minSize, idx, j, k);
BuildXiWeights(XiParaF8_h, z_global[idx]+minSize, y_global[j]-minSize, idx, j, k);
```

`GetParameter_6th2` ?¼å«ä¿®æ­£ï¼š`(j+1)` ???ä½¿ç”¨?Œä?çµ?`pos_z2`??```cpp
RelationXi(j_cont, LT, minSize, a, (NZ6-7), RelationXi_4);
GetParameter_6th2(XiPara_h, pos_z2, RelationXi_4, 4, index_xi);
```

### ?ªæ”¹?•ç??è¼¯
- ç¶²æ ¼?Ÿæ??¬å??HillFunction?Y ?‘æ??¼ã€BFL ?å??–ã€MRT/æ¼”å?ä¸»æ?ç¨‹å??ªæ”¹??- ?¢æ???`GetXiParameter` ä»åœ¨æª”æ?ä¸­ï??¯åˆª?¯ç?ï¼Œä?å½±éŸ¿?¾è?æµç???2. **MRT_Process.h ç¼ºå??†è?** - ç·¨è­¯?¯èª¤
3. **initializationTool.h HillHalfWidth ?¯èª¤** - å½±éŸ¿?Šç??¤æ–·

#### ä¸­å„ª?ˆç?
4. **MRT_Matrix.h omega_1 ?ªå?ç¾?* - ?¯èƒ½å°è‡´ç·¨è­¯?¯èª¤
5. **MRT_Matrix.h ç¼ºå??†è?** - ç·¨è­¯?¯èª¤

#### å¾…å¯¦ä½œå???- [ ] MRT ç¢°æ??‹ç?å®Œæ•´?è¼¯ (é¬†å? + ?†è???
- [ ] Streaming æ­¥é?
- [ ] ?Šç?æ¢ä»¶?•ç??è¼¯
- [ ] ä¸»ç?å¼?main.cu
- [ ] ?å??–æ??´è¨­å®?- [ ] çµæ?è¼¸å‡º?Ÿèƒ½

### D2Q9 æ¨¡å??Ÿåº¦?¹å?

```
?¹å?ç·¨è??‡é€Ÿåº¦?‘é?ï¼?F0: (0, 0)
F1: (1, 0)   F2: (0, 1)
F3: (-1, 0)  F4: (0, -1)
F5: (1, 1)   F6: (-1, 1)
F7: (-1, -1) F8: (1, -1)
```

---

## 2026-01-16 ?¯èª¤ä¿®å¾©å®Œæ?

### å·²ä¿®å¾©ç?é«˜å„ª?ˆç??¯èª¤

#### 1. initializationTool.h:20 - HillHalfWidth å®šç¾©?¯èª¤ ??- **?é?**ï¼š`#define HillHalfWidth 9.0*(54.0/28.0)` ?¼ç‚º 17.36
- **ä¿®å¾©**ï¼šæ”¹??`#define HillHalfWidth (54.0/28.0)` ??1.9286
- **å½±éŸ¿**ï¼šä¿®æ­????BFL ?Šç??¤æ–·?½æ•¸?„è?ç®—ç???
#### 2. MRT_Process.h - ???è¶Šç??¯èª¤ ??- **?é?**ï¼šæ??‰è?ä½¿ç”¨ `M[i][9]` ??`F9_in`ï¼Œä? D2Q9 æ¨¡å??ªæ? 0-8 ç´¢å?
- **ä¿®å¾©**ï¼šå? `M[i][9]` ?¹ç‚º `M[i][8]`ï¼Œ`F9_in` ?¹ç‚º `F8_in`
- **å½±éŸ¿**ï¼šé¿?è??¶é?è¶Šç?å­˜å?

#### 3. MRT_Process.h - ç¼ºå??†è? ??- **?é?**ï¼šm0-m8 ?„è?å®šç¾©å¾Œç¼ºå°‘å???- **ä¿®å¾©**ï¼šæ?è¡Œç?å°¾å?ä¸Šå???`;`
- **å½±éŸ¿**ï¼šä¿®æ­?·¨?†è?æ³•éŒ¯èª?
### å·²ä¿®å¾©ç?ä¸­å„ª?ˆç??¯èª¤

#### 4. MRT_Matrix.h:10 - omega_1 ?ªå?ç¾???- **?é?**ï¼š`s1 = omega_1` ä½?omega_1 ?ªåœ¨ variables.h å®šç¾©
- **ä¿®å¾©**ï¼šæ”¹??`s1 = omega_2` (?¹æ? variables.h:36 å®šç¾©)
- **å½±éŸ¿**ï¼šä½¿?¨æ­£ç¢ºç?é¬†å??ƒæ•¸

#### 5. MRT_Matrix.h - ç¼ºå??†è? ??- **?é?**ï¼šRelaxation å·¨é??„è?ç¼ºå??†è?
- **ä¿®å¾©**ï¼šæ??‰è?çµå°¾è£œä??†è? `;`
- **å½±éŸ¿**ï¼šä¿®æ­?·¨?†è?æ³•éŒ¯èª?
### ä¿®å¾©ç¸½ç?

?€?‰é˜»?·æ€§ç·¨è­¯éŒ¯èª¤å·²ä¿®å¾©ï¼Œç?å¼ç¢¼?®å??€?‹ï?
- ??èªæ??¯èª¤ï¼šå…¨?¨ä¿®å¾?- ?????è¶Šç?ï¼šå·²ä¿®æ­£
- ???ªå?ç¾©è??¸ï?å·²ä¿®æ­?- ? ï? ?Ÿèƒ½å®Œæ•´?§ï?MRT ç¢°æ??è¼¯ä»å?å¯¦ä?

---

## 2026-01-16 å¤–å??…å¯¦ä½œè?è«?
### ?é??†æ?ï¼šHalf-Way Correction Force

?®å?ç¨‹å??„å??¢æ•£??*ç¼ºå? Half-Way (Guo) ä¿®æ­£**ï¼Œå??¨ä»¥ä¸‹å?é¡Œï?

#### 1. ç¼ºå? Half-Way ä¿®æ­£??**?¾ç?**ï¼šå??…åª?¨å¸¸?¸ä??¸ï?ä¸å« `(1 - s/2)` ? å?ï¼Œä?ä¸å«?Ÿåº¦ `u` ?¸é???
**å½±éŸ¿**ï¼?- ?›å??Ÿåº¦?„å½±?¿åª?‰ä??ç²¾åº?- ? ?é›¢ 1 ?–æ??Ÿè?å¤§æ??ƒå??¥å?å·?- ?¡æ?æ­?¢º?•æ?æµé?å°å??›ç??¿æ?

#### 2. ?›é??†é?ä¸å???**?¾ç?**ï¼šåªå°å« Y ?†é??„é€Ÿåº¦?¹å?? å…¥å¤–å?
- D2Q9 æ¨¡å?ä¸­ï??ªæ? F1(+Y), F3(-Y), F5(+Y,+Z), F6(-Y,+Z), F7(-Y,-Z), F8(+Y,-Z) ? å…¥??- F2(+Z), F4(-Z) æ²’æ?? å…¥?›é?

**?Ÿå?**ï¼šç?å¼å? `Force[0]` ?¶æ??®ä?è»¸å??„é??›ï?æ²?Y ?¹å?ï¼‰ï??¶ä??†é???0

### Half-Way (Guo) ?›æ¨¡?‹ä¿®æ­?–¹æ¡?
#### ä¿®æ­£ 1ï¼šå?å¢é? (Force Distribution)

?¨ç¢°?å?ï¼Œå?æ¯å€‹æ–¹?‘ç??›é?ä½¿ç”¨ Guo ?¬å?ï¼?
```c
Fi = wi * (1 - sÎ±/2) * [ (ei - u)/cs^2 + (eiÂ·u) ei / cs^4 ] Â· F
```

**?ƒæ•¸èªªæ?**ï¼?- `wi`ï¼šD2Q9 æ¬Šé?ä¿‚æ•¸
  - w0 = 4/9 (?œæ­¢?¹å?)
  - w1,2,3,4 = 1/9 (?´å?ï¼šÂ±Y, Â±Z)
  - w5,6,7,8 = 1/36 (?œå?ï¼šå?è§’ç?)
- `sÎ±`ï¼šå???moment ?„é?å¼›ç?
  - BGK: ä½¿ç”¨ `(1 - 1/(2?))`
  - MRT: ä½¿ç”¨å°æ??©ç? `sÎ±` ??- `ei`ï¼šç¬¬ i ?‹é€Ÿåº¦?¹å?
- `u`ï¼šå?è§€?Ÿåº¦
- `F`ï¼šé?ç©å??‘é?
- `cs = 1/??`ï¼šæ™¶?¼è²??
#### ä¿®æ­£ 2ï¼šé€Ÿåº¦ä¿®æ­£ (Velocity Correction)

è¨ˆç?å®è??Ÿåº¦?‚ä½¿?¨å?æ­¥ä¿®æ­??

```c
u = (Î£ ei fi + 0.5 * F * dt) / rho
```

?–å??©æ­¥ï¼?```c
u_raw = (Î£ ei fi) / rho
u = u_raw + 0.5 * F * dt / rho
```

?™æ¨£ momentum ?´æ–°?‡å?å¢é?ä¸€?´ï?ä¿è?äºŒé?ç²¾åº¦??
### D2Q9 æ¨¡å??„å??…å???
#### ?Ÿåº¦?¹å??‡æ???
```
?¹å?  ?Ÿåº¦?‘é? (ey, ez)  æ¬Šé? wi    ?¯å¦?«Y?†é?
F0    (0, 0)             4/9       ??F1    (1, 0)             1/9       ??(+Y)
F2    (0, 1)             1/9       ??(ç´?Z)
F3    (-1, 0)            1/9       ??(-Y)
F4    (0, -1)            1/9       ??(ç´?Z)
F5    (1, 1)             1/36      ??(+Y+Z)
F6    (-1, 1)            1/36      ??(-Y+Z)
F7    (-1, -1)           1/36      ??(-Y-Z)
F8    (1, -1)            1/36      ??(+Y-Z)
```

#### ?®å?å¯¦ä??é?

**?®è»¸?›å?è¨?*ï¼?- ?‡è¨­?ªæ? Fy ??0ï¼ŒFz = 0
- ? æ­¤ F2, F4 (ç´?Z ?¹å?) ??`ei Â· F = 0`ï¼Œä?? å???- ?™æ˜¯ç°¡å??‡è¨­ï¼Œé©?¨æ–¼ç´?Y ?¹å?é©…å?æµ?
**?šç”¨å¤–å?å¯¦ä?**ï¼?- å¦‚é??¯æ´ä»»æ??¹å?å¤–å? (Fy, Fz)ï¼Œé?å°?*?€??9 ?‹æ–¹??*??Guo ?¬å??†é?
- ?…å« `(1 - sÎ±/2)` ? å??Œé€Ÿåº¦?¸é???`u`

### å¾…å¯¦ä½œé???
- [ ] å¯¦ä? Guo ?›æ¨¡?‹ç??›å??è?ç®?- [ ] ? å…¥ `(1 - sÎ±/2)` ä¿®æ­£? å?
- [ ] å¯¦ä??Ÿåº¦?¸é???`(eiÂ·u) ei / cs^4`
- [ ] ä¿®æ­£å®è??Ÿåº¦è¨ˆç?ï¼ˆå?æ­¥ä¿®æ­??
- [ ] ?¯æ´?‘é?å¤–å? `(Fy, Fz)` ?Œé??®è»¸
- [ ] å°æ???9 ?‹æ–¹?‘æ­£ç¢ºå??å???
### ?©ç??ç¾©

**?ºä??€è¦?Half-Way ä¿®æ­£**ï¼?1. **äºŒé?ç²¾åº¦**ï¼šç¢ºä¿é›¢??Navier-Stokes ?¹ç??„æ??“ç²¾åº¦é???O(dtÂ²)
2. **æ­?¢º?•é?å®ˆæ?**ï¼šå?å¢é??‡é€Ÿåº¦?´æ–°?¨å??‚é?æ­¥ä???3. **æ¸›å??¼å??ˆæ?**ï¼šé¿?å??…å??¥ç??¸å€¼å?å·?
**?ƒè€ƒæ???*ï¼?- Guo, Z., Zheng, C., & Shi, B. (2002). Discrete lattice effects on the forcing term in the lattice Boltzmann method. *Physical Review E*, 65(4), 046308.

---

## 2026-01-17 Xi ?’å€¼æ??ç?è¨˜ï?GetIntrplParameter_Xi / GetXiParameter / GetParameterXiï¼?
### ?ºä?éº¼é?è¦å??‚è¼¸??`pos_y` ??`pos_z`ï¼?
è£œå?ï¼š`GetIntrplParameter_Xi()` ?¬èº«??`void` ?¡å??¸ï??Ÿæ­£?€è¦?`pos_z`/`pos_y` ?„æ˜¯å®ƒåœ¨è¿´å?ä¸­å‘¼?«ç? `GetXiParameter(...)`ï¼ˆCPU ?ˆå???`GetParameterXi(...)`ï¼‰ï??¨ä???`(y,z)` è½‰æ? `pos_xi` ä¸¦ç”¢?Ÿæ??¼æ??ã€?
Periodic Hill ?„å?å£é?åº¦ç‚º `HillFunction(y)`ï¼ˆéš¨ `y` ?¹è?ï¼‰ï?? æ­¤å¾ç‰©?†åº§æ¨?`(y, z)` è½‰æ??¡å?æ¬¡å???`xi` ?‚ï??€è¦å??‚çŸ¥?“ï?
- ?™å€‹ä?ç½®ç?åº•å?é«˜åº¦ `HillFunction(pos_y)`
- ?™å€‹ä?ç½®ç??‚ç›´åº§æ? `pos_z`

?œéµ?¨æ–¼ï¼šå??¨é€šé?é«˜åº¦?‡é›¢å£è??¢éƒ½è·?`y` ?‰é?ï¼Œæ?ä»¥å?ä¸€??`pos_z` ?¨ä??Œç? `pos_y` ?ƒå??°ä??Œç? `pos_xi`??
`pos_xi` ?„æ ¸å¿ƒé?ä¿‚å?ï¼ˆæ?å¿µä?å°±æ˜¯ä½ ç?è¨˜è£¡??`pos_z - Hill - 0.5*minSize`ï¼‰ï?

```cpp
L      = LZ - HillFunction(pos_y) - minSize;
pos_xi = LXi * (pos_z - (HillFunction(pos_y) + minSize/2.0)) / L;
```

ä¹Ÿå?æ­¤åœ¨ä½¿ç”¨ä¸Šï???`y+` ??`y-` ?ƒå??‰åˆ°ä¸å???`HillFunction(yÂ±minSize)`ï¼Œæ?ä»¥å³ä½?`pos_z` ä¸€æ¨??`pos_xi` ä»æ?ä¸å?ï¼ˆæ??é?è¦å??‹ç?ï¼‰ã€?
### ?©å€‹æ?æ¨™ï?`XiPara` / `Pos_xi`ï¼‰ç??ç¾©

ä»?`GetParameterXi(...)` / `GetXiParameter(...)` ?™é??½å?ä¾†èªªï¼Œå¸¸è¦‹æ??©å€‹ã€Œæ?æ¨™é??å??¸ï?
- `XiPara`ï¼ˆæ? GPU ?ˆç? `double *XiPara*_h[7]`ï¼‰ï??¨ä?**å­˜æ??¼æ???*??th order ?ƒæ? 7 ?‹æ??ï??€ä»¥ç¬¬ä¸€ç¶­æ˜¯ 7ï¼›ç¬¬äºŒç¶­?¯ã€Œè?å­˜åˆ°?ªå€‹ç¶²?¼é??ç?ç´¢å???- `Pos_xi`ï¼ˆä?å¦?`xi_h`ï¼‰ï?`xi` ?¹å???*åº§æ????**ï¼Œæ?ä¾›çµ¦ `GetParameter_6th(...)` ?»å?ä½?stencil ä¸¦è?ç®?7 é»æ??ã€?
### ?ºä?éº¼ã€ŒçŸ©??¬¬äºŒå€‹å??¸ã€æ???`NYD6*NZ6`ï¼?
`z_h` ?¨ç?å¼è£¡?¯ä»¥ `(j,k)` ??2D ?†ä?ä½¿ç”¨ï¼Œä??šå¸¸?¨ä?ç¶­é€??è¨˜æ†¶é«”å??¾ä¸¦ flattenï¼?- `index = j*NZ6 + k`
- ç¸½é???= `NYD6 * NZ6`

??`xi` ?’å€¼æ??æ???`y` ?¹è?ï¼ˆå???`HillFunction(y)` ?ƒå½±??`L` ??`pos_xi`ï¼‰ï??€ä»¥æ???`(j,k)` ?½é?è¦å??ªä?çµ?7 é»æ??ï?? æ­¤æ¯å€‹æ??é™£?—ç??·åº¦å°±æ???`NYD6*NZ6`??
### ?œéµç¨‹å?ç¢¼ï?å¾?`periodic hill_6thIBLBM` ?œå‡ºï¼?
æª”æ?ï¼š`periodic hill_6thIBLBM/50hill_4GPU/initialization.h`

```cpp
void GetXiParameter(
    double *XiPara_h[7],    double pos_z,       double pos_y,
    double *Pos_xi,         int IdxToStore,     int k  )
{
    double L = LZ - HillFunction(pos_y) - minSize;
    double pos_xi = LXi * (pos_z - (HillFunction(pos_y)+minSize/2.0)) / L;

    if( k >= 3 && k <= 6 ){
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, 3 );
    } else if ( k >= NZ6-7 && k <= NZ6-4 ) {
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, NZ6-10 );
    } else {
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, k-3 );
    }
}

void GetIntrplParameter_Xi() {
    for( int j = 3; j < NYD6-3; j++ ){
    for( int k = 3; k < NZ6-3;  k++ ){
        GetXiParameter( XiParaF3_h,  z_h[j*NZ6+k],         y_h[j]-minSize, xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF4_h,  z_h[j*NZ6+k],         y_h[j]+minSize, xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF5_h,  z_h[j*NZ6+k]-minSize, y_h[j],         xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF6_h,  z_h[j*NZ6+k]+minSize, y_h[j],         xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF15_h, z_h[j*NZ6+k]-minSize, y_h[j]-minSize, xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF16_h, z_h[j*NZ6+k]-minSize, y_h[j]+minSize, xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF17_h, z_h[j*NZ6+k]+minSize, y_h[j]-minSize, xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF18_h, z_h[j*NZ6+k]+minSize, y_h[j]+minSize, xi_h, j*NZ6+k, k );
    }}
}
```

ï¼ˆå??§æ?å¿µï?æª”æ?ï¼š`initialization.h`

```cpp
void GetParameterXi(double** XiPara , double pos_y , double pos_z , double* Pos_xi , double now , double start){
    double L = LZ - HillFunction( pos_y ) - minSize;
    double pos_xi = (LXi / L) * (pos_z - HillFunction( pos_y ) - minSize/2.0);
    ...
}
```

### ä½ ä¿®?¹ç? `GetParameterXi(...)` ?¯å¦ç­‰åƒ¹ï¼Ÿï??¿å?é¡¯å??¡å?æ¬¡å?ï¼?
ä½ æ??ºç??¹æ??¸å??¯ï??ˆç??¸å?ä½ç½®æ¯”ä? `rate`ï¼ˆç???`xi/LXi`ï¼‰ï??æ??Œä??‹æ?ä¾‹æ??æ??‹å??ƒå? `y_ref = y_global[now_y]` ?„å¯¦é«”åº§æ¨?`pos_z_mapped`ï¼Œæ?å¾Œç›´?¥åœ¨ `Pos_z` ä¸Šç? 6th æ¬Šé???
?¨ç›®?å?æ¡ˆç? `tanhFunction(...)` ç¶²æ ¼å®šç¾©ä¸‹ï??™æ¨£?šåœ¨?¸å­¸ä¸Šæ˜¯**ç­‰åƒ¹**?„ï?æ¬Šé??ƒä??´ï?å·®ç•°?ªæ?ä¾†è‡ªæµ®é?èª¤å·®ï¼‰ï?? ç‚ºï¼?- ?±å?ç¾©å¯å¾?`z(y,k) - Hill(y) - minSize/2` ??`xi(k)` å°å?ä¸€çµ?`k` ??*ç·šæ€§æ?ä¾‹ç¸®??*ï¼ˆ`tanhFunction` å°é•·åº¦å???`L` ?¯ç??§ç?ï¼‰ã€?- 6th Lagrange æ¬Šé??¨åº§æ¨™å?ä»¿å?è®Šæ?ï¼ˆç??§ç¸®??+ å¹³ç§»ï¼‰ä?ä¸è?ï¼›ç”¨ `xi_h` ?¶åº§æ¨™ç?æ¬Šé?ï¼Œæ??¹ç”¨å°æ? `y_ref` ??`z` åº§æ????ç®—æ??ï??¬è³ªä¸Šç??¹ã€?
ä½ ç??¬ä¸­?™å…©è¡Œç??©ç??ç¾©ä¹Ÿå¯ä»¥é€™æ¨£?‹ï?
- `rate = (pos_z - Hill(pos_y) - minSize/2) / L(pos_y)` ?¶å¯¦å°±æ˜¯ `rate = pos_xi / LXi`ï¼ˆåª?¯ä?ä¸å?é¡¯å?å¯«å‡º `pos_xi`ï¼‰ã€?- `pos_z_mapped = rate * L(y_ref) + Hill(y_ref) + minSize/2` ?¯æ??Œä???`rate` ? å??å??ƒå???`z`??
å¯¦ä?ä¸Šæ³¨?å…©é»ï?
- ä½ ç?å¼ç?æ®µä¸­ `double pos_z = ...;` ?ƒå??ƒæ•¸ `pos_z` ?å?ï¼Œæ?ç·¨è­¯å¤±æ?ï¼›è??¹å???`pos_z_mapped`ï¼ˆæ?é¡ä¼¼ï¼‰ã€?- `Pos_z` å¿…é???**`y_ref = y_global[now_y]` ?????*??`z` åº§æ????ï¼ˆå?æ¨?? flatten ç´¢å?ï¼‰ï??Œä? `now_y/now_z` ?„ç´¢å¼•è?è·?`Pos_z` ?„å„²å­˜æ–¹å¼ä??´ï??æ?ç­‰åƒ¹??
### ?œæ–¼?²é¢ä¸‹ã€?D/3D ?é?ç½®é€??æ¬Šé??æ˜¯?¦æ?å¼•å…¥é¡å?èª¤å·®ï¼?
ä½ æ??ºç??´è¦º?¯å??†ç?ï¼šå??œç”¨?Œå¹³?¢ç›´è§’åº§æ¨™ã€ç??³æ??»ç?ï¼Œ`y` ä¸å???`Hill(y)` ä¸å?ï¼Œç¢ºå¯¦æ?å°è‡´?Œä??‹ç‰©??`z` å°æ??°ä??Œç??¡å?æ¬¡åº§æ¨™ï??‹èµ·ä¾†å??¯ã€Œå??‹é??±ç”¨?Œä?çµ?`xi` æ¬Šé??æ??¢ç?èª¤å·®??
ä½†é€™ä»½ç¨‹å?ï¼ˆ`GetXiParameter` + `F3/F4..._Intrpl7` ?„å¯«æ³•ï??¶å¯¦?¯åœ¨??*?²ç?åº§æ? (y, xi) ?„å¼µ?ç??’å€?*ï¼Œè€Œä??¯åœ¨ (y, z) ?„ç›´è§’ç¶²?¼ä??’å€¼ï?
- ç¶²æ ¼è³‡æ?é»æ˜¯ `(y_j, xi_k)` ??tensor-productï¼›å??‰åˆ°?©ç?åº§æ??æ?è®Šæ? `(y_j, z(y_j,xi_k))` ?„ã€Œæ›²?¢ç¶²?¼ã€ã€?- ? æ­¤?¨å? 2D/3D ?’å€¼æ?ï¼Œå…§å±¤å??¨å?ä¸€??`xi_dep`ï¼ˆç”±?®æ?é»?`(pos_y,pos_z)` ç®—å‡ºï¼‰åœ¨æ¯å€?`y_j` ä¸Šæ? `f(y_j, xi_dep)`ï¼Œå?å±¤å?æ²?`y` ?’åˆ° `y_dep`ï¼Œé€™æ˜¯æ¨™æ???separable/interpolation-in-computational-space ä½œæ???- ?›å¥è©±èªªï¼Œç?å¼?*ä¸æ˜¯**?‡è¨­ä¸å? `y` ?‚ã€Œå?ä¸€??`z`?è??±ç”¨?Œä?å¥—æ??ï??Œæ˜¯??`xi` ?Šä???`y` ?„å??´æ–¹?‘å?æ­???–å?ï¼Œåœ¨?Œä???`xi` ä¸Šæ??¼ã€?
ä»ç„¶?ƒæ?èª¤å·®ï¼Œä?ä¸»è??¯ã€Œæ??¼æˆª?·èª¤å·®ã€ï?è·Ÿç¶²?¼è§£?åº¦/?æ•¸/?½æ•¸?‰æ?åº¦æ??œï?ï¼Œä??¯å??ºæ??è¢«?¯èª¤?±ç”¨??
å¦‚æ?ä½ æƒ³?šæ›´?Œç›´è§’ç‰©?†ç©º?“ã€ç??šæ?ï¼ˆå›ºå®šç‰©??`z` ?¨ä???`y` ä¸Šå??¼ï?ï¼Œé‚£å°±é?è¦å?æ¯å€?`y` stencil é»å??ªç? `xi(y_j, z_dep)`ï¼Œè?æ¯ä???`y_j` ?½ç”¨**ä¸å???`xi` æ¬Šé?**?šå…§å±¤æ??¼ï????è®Šæ??å¼µ?ç?ï¼ˆæ??´æ?è²´ç?ï¼‰æ??¼æ?ç¨‹ï??æœ¬?‡å¯¦ä½œè??œåº¦?½æ??é¡¯ä¸Šå???
### `xi_h[0..2]`ï¼ˆä»¥?Šå°¾ç«?bufferï¼‰æœª?å??–æ˜¯?¦æ??‰å?é¡Œï?

?¨å?ä½œè€…ç? `periodic hill_6thIBLBM/50hill_4GPU` ä¸­ï?`xi_h` ?„å?å§‹å??»æ??ªå??Œé? buffer?å??“ï?
- `periodic hill_6thIBLBM/50hill_4GPU/initialization.h:130`ï½`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:136`
  - `for( int k = bfr; k < NZ6-bfr; k++ )`ï¼Œå…¶ä¸?`bfr = 3`
  - ?€ä»¥åª?ƒç? `k = 3 .. NZ6-4`
  - `k = 0,1,2` ??`k = NZ6-3, NZ6-2, NZ6-1` ?½ä??ƒè¢«è³¦å€¼ï?å±¬æ–¼ buffer/ghost layersï¼?
å°æ?ç¨‹å?ç¢¼ï?`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:110` ?§ç? `GenerateMesh_Z()`ï¼‰ï?

```cpp
int bfr = 3;
...
for( int k = bfr; k < NZ6-bfr; k++ ){
    xi_h[k] = tanhFunction( LXi, minSize, a, (k-3), (NZ6-7) ) - minSize/2.0;
}
```

**å°æ?è§?œ¬èº«é€šå¸¸ä¸æ??‰å?é¡?*ï¼Œå?? æ˜¯å¾Œç??’å€¼æ??ç??é?ç½®è? kernel è¨ˆç?ä¹Ÿå?æ¨?·³??buffer ç¯„å?ï¼?- `GetIntrplParameter_Xi()` ?ªå? `k = 3 .. NZ6-4` ?å?è¨ˆç?æ¬Šé?  
  - `periodic hill_6thIBLBM/50hill_4GPU/initialization.h:231`ï½`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:241`
- CUDA kernel ?´æ¥?’é™¤ `k <= 2` ??`k >= NZ6-3` ??cellï¼ˆä??ƒç”¨??buffer ä½ç½®??`idx_xi = j*NZ6 + k`ï¼? 
  - `periodic hill_6thIBLBM/50hill_4GPU/evolution.h:80`ï¼ˆ`if( ... k <= 2 || k >= NZ6-3 ) return;`ï¼?
å°æ?ç¨‹å?ç¢¼ï??é?ç½®åªç®?interiorï¼‰ï?

```cpp
for( int j = 3; j < NYD6-3; j++ ){
for( int k = 3; k < NZ6-3;  k++ ){
    GetXiParameter( XiParaF3_h, z_h[j*NZ6+k], y_h[j]-minSize, xi_h, j*NZ6+k, k );
    ...
}}
```

å°æ?ç¨‹å?ç¢¼ï?kernel ?´æ¥è·³é? bufferï¼‰ï?

```cpp
if( i <= 2 || i >= NX6-3 || k <= 2 || k >= NZ6-3 ) return;
```

?Œä?å°±ç??¯åœ¨? è?ä¸‹å???`k = 3..6`ï¼Œæ??ä??ƒç”¨?ºå? stencil èµ·é? `n=3`ï¼Œå?æ­¤å??¨ä??€è¦ç¢°??`xi_h[0..2]`ï¼?- `periodic hill_6thIBLBM/50hill_4GPU/initialization.h:192`ï½`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:198`

å°æ?ç¨‹å?ç¢¼ï?? è??Šç??ºå???`n=3`ï¼Œæ?ä»¥ç”¨?°ç???`xi_h[3..9]`ï¼‰ï?

```cpp
if( k >= 3 && k <= 6 ){
    GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, 3 );
}
```

**ä½†æ?ä¸€?‹ã€Œå®¹?“èª¤?ƒã€ç??¯ä???*ï¼šç?å¼æ??Šæ•´??`xi_h[0..NZ6-1]` ?½è¼¸???·è?ï¼Œå??«æœª?å??–ç? buffer ?¼ï?
- è¼¸å‡º `meshXi.DAT` ?ƒæ? `xi_h[0..NZ6-1]` ?¨éƒ¨å¯«å‡ºä¾? 
  - `periodic hill_6thIBLBM/50hill_4GPU/initialization.h:170`ï½`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:175`
- `cudaMemcpy(xi_d, xi_h, NZ6*sizeof(double), ...)` ä¹Ÿæ??Šæœª?å??–ç? buffer ä¸€èµ·æ‹·??GPU  
  - `periodic hill_6thIBLBM/50hill_4GPU/initialization.h:178`

å°æ?ç¨‹å?ç¢¼ï?è¼¸å‡º/?·è??…å« bufferï¼‰ï?

```cpp
for( int k = 0; k < NZ6; k++ ){
    fprintf( meshXi, "%.15lf\n", xi_h[k] );
}
...
CHECK_CUDA( cudaMemcpy(xi_d, xi_h, NZ6*sizeof(double), cudaMemcpyHostToDevice) );
```

? æ­¤ï¼?- **æ¨¡æ“¬è¨ˆç???*ï¼šé€šå¸¸å®‰å…¨ï¼ˆå???buffer ä¸æ?è¢«ç”¨?°ï???- **å¾Œè????¤éŒ¯??*ï¼š`meshXi.DAT` ?„å? 3 ?‹è??€å¾?3 ?‹å€¼å¯?½æ˜¯?¨æ?/?ƒåœ¾?¼ï?å®¹æ?? æ?èª¤è§£??
å»ºè­°ï¼ˆé?å¿…è?ï¼Œä??½é¿?è¸©?·ï?ï¼?- ?¨è¼¸??`meshXi.DAT` ?‚åªè¼¸å‡º `k=3..NZ6-4`ï¼›æ?
- ?ç¢º??`xi_h[0..2]`?`xi_h[NZ6-3..NZ6-1]` è¨­æ??‰æ?ç¾©ç??¼ï?ä¾‹å?è¤‡è£½?Šç??æ?è¨­æ? `NAN` è®“é™¤?¯ä??¼ç??ºä?ï¼‰ã€?
### ?ºä?éº?`XiParameter*_d` ä¸æ?è¶Šç??è€?buffer ?€ä¹Ÿä??ƒè¢«?¨åˆ°ï¼?
?–ç„¶ `XiParaF*_d[i]` ?„é?ç½®å¤§å°æ˜¯ `NYD6*NZ6`ï¼ˆæ??‹é?ä¸€çµ„æ??ï?ï¼Œä?ç¨‹å??¨ã€Œç”¢?Ÿã€è??Œä½¿?¨ã€å…©?‹é?æ®µéƒ½??buffer ?’é™¤?‰ï?
- **?¢ç??æ®µ**ï¼šåªå¡?`k = 3 .. NZ6-4`ï¼ˆ`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:231`ï½`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:241`ï¼?- **ä½¿ç”¨?æ®µ**ï¼škernel ?´æ¥ `return` ??`k <= 2` ??`k >= NZ6-3`ï¼ˆ`periodic hill_6thIBLBM/50hill_4GPU/evolution.h:80`ï¼?
?€ä»¥å?ä½ æ??°ç??€? è?å£é¢?„æ??ˆé?ï¼ˆä?å¦?`k=3,4,5`ï¼Œå???`idx_xi = NZ6*j + 3/4/5`ï¼‰ï?
- **ä¸æ?è¶Šç?**ï¼š`idx_xi` ä»åœ¨ `0 .. NYD6*NZ6-1` ç¯„å??§ã€?- **?ƒè¢«ä½¿ç”¨**ï¼šå???`k=3,4,5` ?½ä???`k<=2` ?„æ??¤ç??å…§??- **ä¸ä?è³?`xi_h[0..2]`**ï¼šæ???stencil ?ƒå? `n=3` ?‹å??–ï??¿å?ç¢°åˆ° buffer ??`xi_h[0..2]`ï¼‰ã€?
---

## 2026-01-18 stream_collide CUDA Kernel ?¼å«æ©Ÿåˆ¶

### stream_collide ä¸æ˜¯ for è¿´å?ï¼Œè€Œæ˜¯ CUDA kernel

`stream_collide` ??`stream_collide_Buffer` ??**CUDA `__global__` kernel**ï¼Œä?ä½¿ç”¨?³çµ± CPU ??for è¿´å??‚å??‘ç??Œè¿´?ˆã€æ˜¯?é? **CUDA grid/block ?ç½®**ä¾†å¯¦?¾å¹³è¡Œè?ç®—ã€?
### Grid/Block ?ç½®

æª”æ?ï¼š`periodic hill_6thIBLBM/50hill_4GPU/evolution.h:785`ï½`793`

```cpp
// stream_collide ?¨é€™ç?ï¼ˆä¸»è¦è?ç®—å??Ÿï?
dim3 griddim(  NX6/NT+1, NYD6, NZ6);
dim3 blockdim( NT, 1, 1);

// stream_collide_Buffer ?¨é€™ç?ï¼ˆY?¹å??Šç? buffer ?€ï¼?dim3 griddimBuf(NX6/NT+1, 1, NZ6);
dim3 blockdimBuf(NT, 4, 1);
```

### Kernel ?¼å«ä½ç½®

æª”æ?ï¼š`periodic hill_6thIBLBM/50hill_4GPU/evolution.h`

```cpp
// evolution.h:795 - ?•ç? j=3 ??buffer ?€
stream_collide_Buffer<<<griddimBuf, blockdimBuf, 0, stream1>>>(
    f_old[0], f_old[1], ..., f_old[18],
    f_new[0], f_new[1], ..., f_new[18],
    XPara0_d[0..6], XPara2_d[0..6],
    YPara0_d[0..6], YPara2_d[0..6],
    XiParaF3_d[0..6], XiParaF4_d[0..6], ...
    u, v, w, rho_d, Force_d, 3, ...  // æ³¨æ??€å¾Œå???3 è¡¨ç¤º j=3
);

// evolution.h:823 - ?•ç? j=NYD6-7 ??buffer ?€
stream_collide_Buffer<<<griddimBuf, blockdimBuf, 0, stream1>>>(
    ...
    u, v, w, rho_d, Force_d, NYD6-7, ...  // æ³¨æ??€å¾Œå???NYD6-7
);

// evolution.h:859 - ?•ç?ä¸»è?è¨ˆç??€?Ÿï?interiorï¼?stream_collide<<<griddim, blockdim, 0, stream0>>>(
    f_old[0], f_old[1], ..., f_old[18],
    f_new[0], f_new[1], ..., f_new[18],
    XPara0_d[0..6], XPara2_d[0..6],
    YPara0_d[0..6], YPara2_d[0..6],
    XiParaF3_d[0..6], XiParaF4_d[0..6], ...
    u, v, w, rho_d, Force_d, ...
);
```

### Kernel ?§éƒ¨?„ç´¢å¼•è?ç®?
??kernel ?§éƒ¨ï¼Œæ???`blockIdx`/`threadIdx` è¨ˆç??ºå??‰ç? `(i, j, k)` ç´¢å?ï¼?
æª”æ?ï¼š`periodic hill_6thIBLBM/50hill_4GPU/evolution.h:80` ?„è?

```cpp
int i = blockIdx.x * blockDim.x + threadIdx.x;  // X ?¹å?
int j = blockIdx.y;                              // Y ?¹å?
int k = blockIdx.z;                              // Z ?¹å?

// è·³é? buffer ?€??if( i <= 2 || i >= NX6-3 || k <= 2 || k >= NZ6-3 ) return;
```

### CUDA kernel ç­‰æ???CPU ä¸‰å±¤ for è¿´å?

ä¸Šè¿° CUDA ?ç½®ç­‰å???CPU ä¸Šç?ä¸‰å±¤ for è¿´å?ï¼?
```cpp
// CPU ç­‰æ?å¯«æ?ï¼ˆå?ä¾›ç?è§??å¯¦é???GPU å¹³è??·è?ï¼?for(int i = 3; i < NX6-3; i++)      // X ?¹å?ï¼ˆç”± griddim.x * blockdim.x å±•é?ï¼?for(int j = 0; j < NYD6; j++)        // Y ?¹å?ï¼ˆç”± griddim.y å±•é?ï¼?for(int k = 3; k < NZ6-3; k++)       // Z ?¹å?ï¼ˆç”± griddim.z å±•é?ï¼?{
    // stream + collide è¨ˆç??è¼¯
}
```

### ?ºä?éº¼é?è¦å??‹è???Buffer ?€ï¼?
`stream_collide_Buffer` å°ˆé??•ç? Y ?¹å??Šç?ï¼ˆ`j=3` ??`j=NYD6-7`ï¼‰ï??Ÿå?ï¼?1. ?™ä?ä½ç½®?€è¦ç‰¹æ®Šç??±æ??§é??Œæ?ä»¶è???2. ?†é??•ç??¯ä»¥ä½¿ç”¨ä¸å???streamï¼ˆ`stream1`ï¼‰ï??‡ä¸»è¨ˆç?ï¼ˆ`stream0`ï¼‰å¹³è¡ŒåŸ·è¡?3. `blockdimBuf` ??`y=4` è¡¨ç¤ºä¸€æ¬¡è???4 ?—ï?`j=3,4,5,6` ??`j=NYD6-7, NYD6-6, NYD6-5, NYD6-4`ï¼?
### ?·è??†å?

1. `stream_collide_Buffer`ï¼ˆj=3ï¼‰åœ¨ `stream1` ä¸Šå???2. `stream_collide_Buffer`ï¼ˆj=NYD6-7ï¼‰åœ¨ `stream1` ä¸Šæ???3. `AccumulateUbulk` ??`stream1` ä¸Šè?ç®—å¹³?‡é€Ÿåº¦
4. `stream_collide`ï¼ˆä¸»è¨ˆç??€ï¼‰åœ¨ `stream0` ä¸Šå???
?±æ–¼ä½¿ç”¨ä¸å???CUDA streamï¼ŒBuffer ?•ç??‡ä¸»è¨ˆç??¯ä»¥?¨å??ç??·è???
---

## 2026-01-18 BFL ?Šç?æ¢ä»¶?å??–ç?è¨?
### BFL ?Šç??•ç??„æ ¸å¿ƒæ?å¿?
**?å??¹å?ï¼ˆé›¢?‹å??¢ç??¹å?ï¼‰æ??€è¦å??Šç??•ç?**ï¼Œä??¯å…¥å°„æ–¹?‘ã€?
?¶ç?å­æ??°å??¢å??ƒå?å½ˆï?
- **?¥å??¹å?**ï¼šç?å­?*?²å…¥**å£é¢?„æ–¹??- **?å??¹å?**ï¼šç?å­?*?¢é?**å£é¢?„æ–¹????**?™å€‹æ??€è¦?BFL ?•ç?**

### ?Šç??¤æ–·?½æ•¸?„è¨­è¨ˆé?è¼?
?Šç??¤æ–·?½æ•¸æª¢æ¸¬?„æ˜¯ï¼šã€Œå??™å€‹è?ç®—é?å¾€**?¥å??¹å?**ç§»å?ï¼Œæ?ä¸æ??åˆ°å£é¢ï¼Ÿã€?- å¦‚æ??ƒæ??????™å€‹é???*?å??¹å?**?„é??Œè?ç®—é?
- è©?*?å??¹å?**?„å?ä½ˆå‡½?¸é?è¦ç”¨ BFL ?´æ–°

### D2Q9 ?¹å??‡é??Œè??†å??‰é?ä¿?
| ?Šç??¤æ–·?½æ•¸ | ?¥å??¹å?ï¼ˆæ??†ï? | ?å??¹å?ï¼ˆé›¢?‹ï? | ?€è¦æ›´??| ?¨èª°?„å€¼ä??’å€?|
|-------------|-----------------|-----------------|---------|---------------|
| `IsLeftHill_Boundary_yPlus` | F3 (-Y) ?å·¦ä¸?| **F1 (+Y)** | F1 | F3 |
| `IsRightHill_Boundary_yMinus` | F1 (+Y) ?å³ä¸?| **F3 (-Y)** | F3 | F1 |
| `IsLeftHill_Boundary_Diagonal45` | F7 (-Y,-Z) ?å·¦ä¸?| **F5 (+Y,+Z)** | F5 | F7 |
| `IsRightHill_Boundary_Diagonal135` | F8 (+Y,-Z) ?å³ä¸?| **F6 (-Y,+Z)** | F6 | F8 |

### BFLInitialization ?½æ•¸çµæ?

æª”æ?ï¼š`initialization.h:129`

```cpp
void BFLInitialization() {
    for(int j = 3; j < NY6-3; j++){
        for(int k = 3; k < NZ6-3; k++){

            // F1 ?„é??Œè??†ï?F3 ?¥å??å·¦ä¸???F1 ?å??¢é?
            if(IsLeftHill_Boundary_yPlus(y_global[j], z_global[j*NZ6+k])){
                double q1 = Left_q_yPlus(y_global[j], z_global[j*NZ6+k]);
                double delta1 = minSize * (1.0 - 2.0*q1);
                // BFL ?å?é»åœ¨ +Y ?¹å?: y + delta, z ä¸è?
                GetParameter_6th(YBFLParaF3_h, y_global[j]+delta1, ...);  // ??F3 ?’å€?                GetXiParameter(XiBFLParaF3_h, z_global[j*NZ6+k], y_global[j]+delta1, ...);
                Q1_h[j*NZ6+k] = q1;
            }

            // F3 ?„é??Œè??†ï?F1 ?¥å??å³ä¸???F3 ?å??¢é?
            if(IsRightHill_Boundary_yMinus(y_global[j], z_global[j*NZ6+k])){
                double q3 = Right_q_yMinus(y_global[j], z_global[j*NZ6+k]);
                double delta3 = minSize * (1.0 - 2.0*q3);
                // BFL ?å?é»åœ¨ -Y ?¹å?: y - delta, z ä¸è?
                GetParameter_6th(YBFLParaF1_h, y_global[j]-delta3, ...);  // ??F1 ?’å€?                GetXiParameter(XiBFLParaF1_h, z_global[j*NZ6+k], y_global[j]-delta3, ...);
                Q3_h[j*NZ6+k] = q3;
            }

            // F5 ?„é??Œè??†ï?F7 ?¥å??å·¦ä¸???F5 ?å??¢é?
            if(IsLeftHill_Boundary_Diagonal45(y_global[j], z_global[j*NZ6+k])){
                double q5 = Left_q_Diagonal45(y_global[j], z_global[j*NZ6+k]);
                double delta5 = minSize * (1.0 - 2.0*q5);
                // BFL ?å?é»åœ¨ (+Y,+Z) ?¹å?: y + delta, z + delta
                GetParameter_6th(YBFLParaF7_h, y_global[j]+delta5, ...);  // ??F7 ?’å€?                GetXiParameter(XiBFLParaF7_h, z_global[j*NZ6+k]+delta5, y_global[j]+delta5, ...);
                Q5_h[j*NZ6+k] = q5;
            }

            // F6 ?„é??Œè??†ï?F8 ?¥å??å³ä¸???F6 ?å??¢é?
            if(IsRightHill_Boundary_Diagonal135(y_global[j], z_global[j*NZ6+k])){
                double q6 = Right_q_Diagonal135(y_global[j], z_global[j*NZ6+k]);
                double delta6 = minSize * (1.0 - 2.0*q6);
                // BFL ?å?é»åœ¨ (-Y,+Z) ?¹å?: y - delta, z + delta
                GetParameter_6th(YBFLParaF8_h, y_global[j]-delta6, ...);  // ??F8 ?’å€?                GetXiParameter(XiBFLParaF8_h, z_global[j*NZ6+k]+delta6, y_global[j]-delta6, ...);
                Q6_h[j*NZ6+k] = q6;
            }
        }
    }
}
```

### delta ?„è?ç®?
**delta** ??BFL ?å?é»ç›¸å°æ–¼è¨ˆç?é»ç??ç§»?ï?

```cpp
delta = minSize * (1.0 - 2.0*q)
```

?¶ä¸­ `q` ?¯è?ç®—é??°å??¢ç??¡å?æ¬¡è??¢ï?`0 < q < 1`ï¼‰ã€?
- ??`q = 0`ï¼ˆè?ç®—é??¨å??¢ä?ï¼‰ï?`delta = minSize`
- ??`q = 0.5`ï¼ˆè?ç®—é??¨ä¸­?“ï?ï¼š`delta = 0`
- ??`q = 1`ï¼ˆè?ç®—é??¢å??¢ä??¼ï?ï¼š`delta = -minSize`

å°æ–¼?œå?ï¼?5Â°/135Â°ï¼‰ï???`q` ?¯ç”¨?™ä»½ç¨‹å??„å?ç¾©ï?`q = |?y|/minSize`ï¼Œç??Œæ²¿ link ?„ç„¡? æ¬¡è·é›¢ï¼‰ï???`delta` ä»ç”¨ `delta = minSize*(1-2q)`ï¼›æ??‘ç??Œå¯¦?›ä?ç§»é•·åº¦ã€æ??ƒæ˜¯ `sqrt(2.0)*delta`ï¼ˆä?ä½ åœ¨åº§æ?ä¸Šå??„æ˜¯ `yÂ±delta, zÂ±delta`ï¼‰ã€?
### æ¬Šé?????½å?è¦å?

`YBFLParaF3_h` ??`XiBFLParaF3_h` ?„å‘½?æ?ç¾©ï?
- **F3** ä»?¡¨?™å€‹æ??é™£?—æ˜¯?¨ä?**?’å€?F3 ?†ä??½æ•¸**??- ?’å€¼ç??œæ??¨ä?**?´æ–° F1**ï¼ˆF3 ?„å??¹å?ï¼?
?Œç?ï¼?- `YBFLParaF1_h` ???’å€?F1ï¼Œæ›´??F3
- `YBFLParaF7_h` ???’å€?F7ï¼Œæ›´??F5
- `YBFLParaF8_h` ???’å€?F8ï¼Œæ›´??F6

---

## ?è?ç¨‹å?ç¢¼ç´¢å¼?
### initializationTool.h - ?å??–å·¥?·å‡½??
| ç·¨è? | ?½æ•¸?ç¨± | è¡Œè? | ?Ÿèƒ½èªªæ? |
|------|---------|------|----------|
| 1 | `HillFunction_Inverse_Left()` | 29-43 | å·¦å?ä¸˜å??½æ•¸ï¼ˆä??†æ?ï¼?|
| 2 | `HillFunction_Inverse_Right()` | 55-68 | ?³å?ä¸˜å??½æ•¸ï¼ˆä??†æ?ï¼?|
| 3 | `tanhFunction` | 82-85 | ?™æ›²æ­???å??»ç¶²?¼åº§æ¨™è??›å·¨??|
| 4 | `GetNonuniParameter()` | 96-117 | è¨ˆç??å??»ç¶²?¼ä¼¸ç¸®å???a |
| 5 | `Lagrange_6th()` | 135-138 | ?­é? Lagrange ?’å€¼åŸºåº•å‡½??|
| 6 | `GetParameter_6th()` | 151-159 | ?¢ç??­é? Lagrange ?’å€¼é??ç½®æ¬Šé? |
| 7 | `IsLeftHill_Boundary_yPlus()` | 174-185 | ?¤æ–· F1 (+Y) ?Šç?è¨ˆç?é»?|
| 8 | `IsRightHill_Boundary_yMinus()` | 198-209 | ?¤æ–· F3 (-Y) ?Šç?è¨ˆç?é»?|
| 9 | `IsLeftHill_Boundary_Diagonal45()` | 222-237 | ?¤æ–· F5 (+Y,+Z) ?Šç?è¨ˆç?é»?|
| 10 | `IsRightHill_Boundary_Diagonal135()` | 250-265 | ?¤æ–· F6 (-Y,+Z) ?Šç?è¨ˆç?é»?|
| 11 | `IsLeftHill_Boundary_DiagonalMinus45()` | 282-297 | ?¤æ–· F8 ?Šç?ï¼ˆæ???1?‚ç”¨ï¼?|
| 12 | `IsRightHill_Boundary_DiagonalMinus135()` | 315-330 | ?¤æ–· F7 ?Šç?ï¼ˆæ???1?‚ç”¨ï¼?|
| 13 | `Left_q_yPlus()` | 342-348 | è¨ˆç?å·¦ä? +Y ?¹å? q ??|
| 14 | `Right_q_yMinus()` | 360-366 | è¨ˆç??³ä? -Y ?¹å? q ??|
| 15 | `Left_q_Diagonal45()` | 378-404 | è¨ˆç?å·¦ä? 45Â° ?¹å? q ??|
| 16 | `Right_q_Diagonal135()` | 416-441 | è¨ˆç??³ä? 135Â° ?¹å? q ??|
| 17 | `Left_q_DiagonalMinus45()` | 453-478 | è¨ˆç?å·¦ä? -45Â° ?¹å? q ??|
| 18 | `Right_q_DiagonalMinus135()` | 490-515 | è¨ˆç??³ä? -135Â° ?¹å? q ??|

### initialization.h - ?å??–ä¸»?½æ•¸

| ?½æ•¸?ç¨± | è¡Œè? | ?Ÿèƒ½èªªæ? |
|---------|------|----------|
| `InitialUsingDftFunc()` | 7-31 | ?å??–å?ä½ˆå‡½?¸è?å¤–å???|
| `GenerateMesh_Y()` | 34-47 | å»ºç? Y ?¹å??‡å‹»ç¶²æ ¼ |
| `GenerateMesh_Z()` | 49-74 | å»ºç? Z ?¹å??å??»ç¶²?¼ï??«å±±ä¸˜ï? |
| `GetXiParameter()` | 75-90 | è¨ˆç? Xi ?¹å??’å€¼æ???|
| `GetIntrplParameter_Y()` | 92-97 | ?é?ç½?Y ?¹å??’å€¼æ???|
| `GetIntrplParameter_Xi()` | 100-127 | ?é?ç½?Xi ?¹å??’å€¼æ??ï??€??F1-F8ï¼?|
| `BFLInitialization()` | 129-183 | BFL ?Šç?æ¢ä»¶?å??–ï?æ¬Šé???q ?¼ï? |

### ?œéµè®Šæ•¸å°ç…§è¡?
| è®Šæ•¸?ç¨± | ?¨é€?|
|---------|------|
| `y_global[NY6]` | Y ?¹å??©ç?åº§æ???? |
| `z_global[NY6*NZ6]` | Z ?¹å??©ç?åº§æ????ï¼?D flattenï¼?|
| `xi_h[NZ6]` | ?¡å?æ¬¡å? Z åº§æ?ï¼ˆä??«å±±ä¸˜å½±?¿ï? |
| `XiParaF*_h[7]` | F* ?¹å???Xi ?’å€¼æ??ï?7é»ï? |
| `YBFLParaF*_h[7]` | BFL ?¨ç? Y ?¹å??’å€¼æ???|
| `XiBFLParaF*_h[7]` | BFL ?¨ç? Xi ?¹å??’å€¼æ???|
| `Q*_h[]` | ?„æ–¹?‘ç? BFL q ?¼ï?è¨ˆç?é»åˆ°å£é¢è·é›¢ï¼‰|

---

## 2026-01-18 ?ªå?å§‹å?æ¸…å–®ï¼ˆä¸»ç¨‹å??°å¯«?‚è?çµ±æ•´ï¼?
ä»¥ä??¯ã€Œåœ¨?®å?ç¨‹å?æµç?ä¸­ï?ä¸ä?å®šæ?è¢«å¯«?¥ã€ç?è³‡æ??€?Ÿï??¥å?çºŒç?å¼æ?è®€?°å??‘ï??€?¨é?ç½??å??–é?æ®µå? `memset`?å¡« 0?æ?å¡?`NAN`/sentinel??
### ç¶²æ ¼/åº§æ?

- `xi_h` ??buffer ?€ï¼š`initialization.h:58` ?ªå¯« `k=3..NZ6-4` ??`xi_h[0..2]` ??`xi_h[NZ6-3..NZ6-1]` ?ªå?å§‹å???- `z_global` ??buffer ?€ï¼š`initialization.h:67` ?ªå¯« `k=3..NZ6-4`ï¼Œå¦å¤–åªè£?`k=2` ??`k=NZ6-3`ï¼ˆ`initialization.h:71-72`ï¼‰â? `k=0,1,NZ6-2,NZ6-1` ?ªå?å§‹å?ï¼ˆæ???`j` ?½ä?æ¨????- `Force[1]`ï¼š`InitialUsingDftFunc()` ?ªè³¦??`Force[0]`ï¼ˆ`initialization.h:30`ï¼‰â? ?¥ä???`Force` ?¶æ? `(Fy,Fz)`ï¼Œå? `Force[1]` ?€è¦æ?ç¢ºå?å§‹å???
### ?’å€¼æ??ï???BFLï¼?
- `YPara0_h[*][i]`?`YPara2_h[*][i]`ï¼š`GetIntrplParameter_Y()` ?ªç? `i=3..NY6-4`ï¼ˆ`initialization.h:92-97`ï¼‰â? buffer `i<=2` ??`i>=NY6-3` ?ªå?å§‹å???- `XiParaF1_h..XiParaF8_h`ï¼š`GetIntrplParameter_Xi()` ?ªç? `j=3..NYD6-4`?`k=3..NZ6-4`ï¼ˆ`initialization.h:100-127`ï¼‰â? ä»»ä??¨é€™å€‹ç??å???`(j,k)` æ¬Šé??ªå?å§‹å???
### BFL æ¬Šé???q ??
- `Q1_h/Q3_h/Q5_h/Q6_h`ï¼šåª?¨å???`Is*Boundary*()` ?ç??‚æ?å¯«å…¥ï¼ˆ`initialization.h:140-180`ï¼‰â? ?é??Œé?ä½ç½®?ªå?å§‹å???- `YBFLParaF3_h/YBFLParaF1_h/YBFLParaF7_h/YBFLParaF8_h` ??`XiBFLParaF3_h/XiBFLParaF1_h/XiBFLParaF7_h/XiBFLParaF8_h`ï¼šå?ä¸Šï??ªåœ¨?Šç?é»å¡«æ¬Šé? ???é??Œé?ä½ç½®?ªå?å§‹å???
### ?®å???repo ?§ã€Œæ‰¾ä¸åˆ°å®???ç?ç¬¦è?ï¼ˆåªå­˜åœ¨?¼ç?è¨˜ï?å°šæœª?¨å¯ç·¨è­¯?„ç?å¼æ??ºç¾ï¼?
`initialization.h` ?´æ¥ä½¿ç”¨ä½†åœ¨ç¨‹å?ç¢¼è£¡å°šæœª?œå??°å®£?Šï?ä½ èªª?Œå?è¨­å·²å®???ç???‰¹ï¼‰ï?
- æµå ´/?†ä?ï¼š`rho[]`, `v[]`, `w[]`, `f[9][]`, `Force[]`
- ç¶²æ ¼åº§æ?ï¼š`y_global[]`, `z_global[]`, `xi_h[]`, ä»¥å?ä½ åœ¨æ¬Šé??Ÿæ??‚ç”¨??`y_h[]`, `z_h[]`
- ä¸€?¬æ??¼æ??ï?`YPara0_h`, `YPara2_h`, `XiParaF1_h..XiParaF8_h`
- BFLï¼š`Q1_h`, `Q3_h`, `Q5_h`, `Q6_h`, `YBFLParaF3_h/YBFLParaF1_h/YBFLParaF7_h/YBFLParaF8_h`, `XiBFLParaF3_h/XiBFLParaF1_h/XiBFLParaF7_h/XiBFLParaF8_h`

### D2Q9 ?Ÿåº¦?¹å?å®šç¾©

```
     F6(-1,+1)  F2(0,+1)  F5(+1,+1)
              \    |    /
               \   |   /
     F3(-1,0) ?â€?F0 ?•â? F1(+1,0)
               /   |   \
              /    |    \
     F7(-1,-1)  F4(0,-1)  F8(+1,-1)
```

| ?¹å? | ?Ÿåº¦?‘é? (ey, ez) | æ¬Šé? |
|------|------------------|------|
| F0 | (0, 0) | 4/9 |
| F1 | (+1, 0) | 1/9 |
| F2 | (0, +1) | 1/9 |
| F3 | (-1, 0) | 1/9 |
| F4 | (0, -1) | 1/9 |
| F5 | (+1, +1) | 1/36 |
| F6 | (-1, +1) | 1/36 |
| F7 | (-1, -1) | 1/36 |
| F8 | (+1, -1) | 1/36 |

---

## 2026-01-18 MRT ?©é™£å·¨é?ä½¿ç”¨?¹å?

### å·¨é?å®šç¾©ä½ç½®

æª”æ?ï¼š`MRT_Matrix.h`

### ä¸‰å€‹ä¸»è¦å·¨??
| å·¨é??ç¨± | å±•é?å¾Œå®£??| ?¨é€?|
|---------|-----------|------|
| `Matrix` | `double M[9][9] = {...}` | ?†ä??½æ•¸ ???©ç©º?“è??›çŸ©??|
| `Matrix_Inverse` | `double M_I[9][9] = {...}` | ?©ç©º?????†ä??½æ•¸?†è??›çŸ©??|
| `Relaxation` | `double s0=0; ... double s8=omega_7;` | é¬†å??ƒæ•¸ï¼ˆå?è§’çŸ©???ç´ ï? |

### evolution.h ä¸­ç?ä½¿ç”¨?¹å?

æª”æ?ï¼š`evolution.h:64-68`

```cpp
// MRT ?©é™£?‡é?å¼›å???(å·¨é?å±•é?å¾Œæ?å®?? M[9][9], M_I[9][9], s0~s8)
Matrix;
Matrix_Inverse;
Relaxation;
```

**æ³¨æ?**ï¼?- å·¨é?å·²å??«è??¸å®£?Šï?**ä¸é?è¦?*?å?å®?? `s0~s8`
- å±•é?å¾Œå³?¯ç›´?¥ä½¿??`M[i][j]`?`M_I[i][j]`?`s0~s8`

### ?¯èª¤ç¤ºç?ï¼ˆæ?å°è‡´?è?å®??ï¼?
```cpp
// ???¯èª¤ï¼šs0~s8 ?ƒè¢«?è?å®??
double s0, s1, s2, s3, s4, s5, s6, s7, s8;  // ?‹å?å®??
Relaxation;  // å·¨é??§å?å®??ä¸€æ¬???ç·¨è­¯?¯èª¤
```

### æ­?¢ºç¤ºç?

```cpp
// ??æ­?¢ºï¼šç›´?¥ä½¿?¨å·¨??Matrix;           // ??double M[9][9] = {...};
Matrix_Inverse;   // ??double M_I[9][9] = {...};
Relaxation;       // ??double s0=0; double s1=omega_2; ... double s8=omega_7;

// ?¾åœ¨?¯ä»¥?´æ¥ä½¿ç”¨ M, M_I, s0~s8
```

### å¸¸è??¯èª¤ï¼šå·¨?†å?ç¨±æ‹¼??
| ?¯èª¤å¯«æ? | æ­?¢ºå¯«æ? |
|---------|---------|
| `Inverse_Matrix;` | `Matrix_Inverse;` |
| `Relaxtion;` | `Relaxation;` |

### é¬†å??ƒæ•¸?©ç??ç¾©

```
s0 = 0        // å¯†åº¦å®ˆæ?ï¼ˆä?é¬†å?ï¼?s1 = omega_2  // ?½é?é¬†å?
s2 = 1.0      // free parameter
s3 = 0        // Y?•é?å®ˆæ?ï¼ˆä?é¬†å?ï¼?s4 = 1.0      // free parameter
s5 = 0        // Z?•é?å®ˆæ?ï¼ˆä?é¬†å?ï¼?s6 = 1.0      // free parameter
s7 = omega_7  // ?ªå??‰å?é¬†å?ï¼ˆè?é»æ»¯ä¿‚æ•¸?¸é?ï¼?s8 = omega_7  // ?ªå??‰å?é¬†å?ï¼ˆè?é»æ»¯ä¿‚æ•¸?¸é?ï¼?```

?¶ä¸­ `omega_2` ??`omega_7` å®šç¾©??`variables.h:36-37`ï¼?
```cpp
#define omega_2  1/(niu/9.0 + 0.5)
#define omega_7  1/(niu/3.0 + 0.5)
```

### MRT ç¢°æ?æµç?ï¼ˆæ­??MRT_Process.hï¼?
```cpp
// 1. å®???©é™£?‡å???Matrix;
Matrix_Inverse;
Relaxation;

// 2. è¨ˆç???m = M * f ï¼ˆä½¿??m_vector å·¨é?ï¼?m_vector;  // ??m0, m1, ..., m8

// 3. è¨ˆç?å¹³è¡¡?‹çŸ© meq = M * feq ï¼ˆä½¿??meq å·¨é?ï¼?meq;  // ??meq0, meq1, ..., meq8

// 4. ç¢°æ?ï¼šf = f - M^(-1) * S * (m - meq) ï¼ˆä½¿??collision å·¨é?ï¼?collision;  // ??F0_in, F1_in, ..., F8_in ?´æ–°
```

### CUDA Constant Memory ?ˆæœ¬ï¼ˆå¯?¸ï?

å¦‚æ?è¦åœ¨ CUDA kernel ä¸­ä½¿?¨ï??¯ä»¥?Ÿç”¨ constant memory ?ˆæœ¬ä»¥æ??‡æ??½ï?

```cpp
// ??main.cu ?‹é ­
#define USE_CUDA_CONSTANT
#include "MRT_Matrix.h"

// ?å???initRelaxationToGPU();  // å°‡é?å¼›å??¸å‚³??GPU

// ??kernel ä¸­ç›´?¥ä½¿??__global__ void kernel() {
    // ä½¿ç”¨ d_M[i][j], d_M_I[i][j], d_S[i]
    // å¾?constant memory è®€?–ï?ä¸é?è¦æ?æ¬¡é??°å»ºç«‹é™£??}
```

---

## 2026-01-19 evolution.h Equilibrium æª¢æŸ¥?‡ä¿®å¾?
### Equilibrium ?†ä??½æ•¸é©—è?çµæ?

æª”æ?ï¼š`evolution.h:153-165`

#### ??æ­?¢º?„éƒ¨??
**1. å¯†åº¦è¨ˆç? (Line 153)**
```cpp
double rho_s = F0_in + F1_in + F2_in + F3_in + F4_in + F5_in + F6_in + F7_in + F8_in;
```
ç¬¦å? D2Q9 å¯†åº¦å®šç¾©ï¼š`? = Î£ fi`

**2. ?Ÿåº¦è¨ˆç? (Lines 154-155)**
```cpp
double v1 = (F1_in + F5_in + F8_in - (F3_in + F6_in + F7_in)) / rho_s;  // Y?¹å??Ÿåº¦
double w1 = (F2_in + F5_in + F6_in - (F4_in + F7_in + F8_in)) / rho_s;  // Z?¹å??Ÿåº¦
```
ç¬¦å? D2Q9 ?•é?å®šç¾©ï¼?- `?uy = F1 + F5 + F8 - F3 - F6 - F7`ï¼ˆå« +Y ?†é??„æ–¹??- ??-Y ?†é??„æ–¹?‘ï?
- `?uz = F2 + F5 + F6 - F4 - F7 - F8`ï¼ˆå« +Z ?†é??„æ–¹??- ??-Z ?†é??„æ–¹?‘ï?

**3. å¹³è¡¡?‹å?ä½ˆå‡½?¸å…¬å¼?(Lines 157-165)**

?€??9 ?‹æ–¹?‘ç?æ¬Šé??‡å…¬å¼éƒ½æ­?¢ºï¼Œç¬¦?ˆæ?æº?D2Q9 å¹³è¡¡?‹ï?
```
f_eq[i] = w[i] * ? * (1 + 3(eÂ·u) + 4.5(eÂ·u)Â² - 1.5uÂ²)
```

| ?¹å? | æ¬Šé? | ?Ÿåº¦?‘é? | ?¬å?é©—è? |
|-----|------|---------|---------|
| F0 | 4/9 | (0,0) | ??`1 - 1.5*udot` |
| F1 | 1/9 | (+1,0) | ??`1 + 3v + 4.5vÂ² - 1.5*udot` |
| F2 | 1/9 | (0,+1) | ??`1 + 3w + 4.5wÂ² - 1.5*udot` |
| F3 | 1/9 | (-1,0) | ??`1 - 3v + 4.5vÂ² - 1.5*udot` |
| F4 | 1/9 | (0,-1) | ??`1 - 3w + 4.5wÂ² - 1.5*udot` |
| F5 | 1/36 | (+1,+1) | ??`1 + 3(v+w) + 4.5(v+w)Â² - 1.5*udot` |
| F6 | 1/36 | (-1,+1) | ??`1 + 3(-v+w) + 4.5(-v+w)Â² - 1.5*udot` |
| F7 | 1/36 | (-1,-1) | ??`1 + 3(-v-w) + 4.5(-v-w)Â² - 1.5*udot` |
| F8 | 1/36 | (+1,-1) | ??`1 + 3(v-w) + 4.5(v-w)Â² - 1.5*udot` |

### å·²ä¿®å¾©ç??¯èª¤

#### 1. Line 8: ä¸­æ??’è?èªæ??¯èª¤ ??- **?é?**ï¼š`#include "variables.h"ï¼š`ï¼ˆä¸­?‡å…¨å½¢å??Ÿï?
- **ä¿®å¾©**ï¼šåˆª?¤å?????`#include "variables.h"`

#### 2. Line 56: è¿´å?ç¯„å??¯èª¤ ??- **?é?**ï¼š`for(int j = 3 ; j <= NZ6-4 ; j++)`ï¼ˆj ?‰è©²??Y ?¹å?ï¼Œç”¨ NY6ï¼?- **ä¿®å¾©**ï¼š`for(int j = 3 ; j < NY6-3 ; j++)`

#### 3. Line 57: è¿´å?ç¯„å?ä¿®æ­£ ??- **?é?**ï¼š`for(int k = 3 ; k <= NZ6-4 ; k++)`
- **ä¿®å¾©**ï¼š`for(int k = 3 ; k < NZ6-3 ; k++)`

#### 4. Lines 108, 121, 134, 147: BFL q>0.5 ?‹ç?å­å„ª?ˆé?åºéŒ¯èª???
**?é?**ï¼šæ•´?¸é™¤æ³•å??´éŒ¯èª¤ç???```cpp
// ?¯èª¤ï¼?/2*q1 = (1/2)*q1 = 0*q1 = 0ï¼ˆæ•´?¸é™¤æ³•ï?
F1_in = (1/2*q1)*f3_old[idx_xi] + ((2*q1-1)/2*q1)*f1_old[idx_xi];
```

**ä¿®å¾©**ï¼šä½¿?¨æµ®é»æ•¸?¤æ?
```cpp
// æ­?¢ºï¼?.0/(2.0*q1) å¾—åˆ°æ­?¢º?„å???F1_in = (1.0/(2.0*q1))*f3_old[idx_xi] + ((2.0*q1-1.0)/(2.0*q1))*f1_old[idx_xi];
```

**BFL ç·šæ€§æ??¼å…¬å¼èªª??*ï¼???`q > 0.5` ?‚ï?ä½¿ç”¨ç·šæ€§æ??¼ï?
```
F_?å? = (1/(2q)) * F_?¥å?(å£é¢é»? + ((2q-1)/(2q)) * F_?å?(è¨ˆç?é»?
```

### ä¿®å¾©å¾Œç?ç¨‹å?ç¢¼ç???
| è¡Œè? | ä¿®å¾©?§å®¹ |
|-----|---------|
| 8 | ?ªé™¤ä¸­æ??’è? |
| 56 | `j < NY6-3` |
| 57 | `k < NZ6-3` |
| 108 | `1.0/(2.0*q1)` ??`(2.0*q1-1.0)/(2.0*q1)` |
| 121 | `1.0/(2.0*q3)` ??`(2.0*q3-1.0)/(2.0*q3)` |
| 134 | `1.0/(2.0*q5)` ??`(2.0*q5-1.0)/(2.0*q5)` |
| 147 | `1.0/(2.0*q6)` ??`(2.0*q6-1.0)/(2.0*q6)` |

### Equilibrium ??MRT_Process.h ?„é?ä¿?
`MRT_Process.h` ä¸­ç? `Equilibrium` å·¨é???`evolution.h` ä¸­æ?å¯«ç?å¹³è¡¡?‹è?ç®—å??½ç›¸?Œï?

**MRT_Process.h å·¨é??ˆæœ¬**ï¼?```cpp
#define Equilibrium \
    double udot = uy*uy + uz*uz; \
    double F0_eq = (4.0/9.0)  * rho_local * (1.0 - 1.5*udot); \
    ...
```
- è¼¸å…¥ï¼š`rho_local`, `uy`, `uz`
- è¼¸å‡ºï¼š`F0_eq ~ F8_eq`

**evolution.h ?‹å¯«?ˆæœ¬**ï¼ˆLines 153-165ï¼‰ï?
```cpp
double rho_s = F0_in + ... + F8_in;
double v1 = (...) / rho_s;
double w1 = (...) / rho_s;
double udot = v1*v1 + w1*w1;
const double F0_eq = (4./9.) * rho_s * (1.0-1.5*udot);
...
```
- ?´æ¥å¾å?ä½ˆå‡½?¸è?ç®—å?åº¦è??Ÿåº¦
- ?¶å?è¨ˆç?å¹³è¡¡??
?©ç¨®?¹å??¸å­¸ä¸Šç??¹ï?`evolution.h` ?„ç??¬æ›´?©å???stream-collide æµç?ä¸­ä½¿?¨ï?? ç‚ºå®ƒç›´?¥å? post-streaming ?„å?ä½ˆå‡½?¸é?å§‹è?ç®—ã€?
---

## 2026-01-19 å¤–å??…å??‹ä¿®æ­????
### ?Ÿèƒ½?®ç?

?¨é€±æ??§é??•æ?ä¸­ï??å?å¤–å??¯æ ¹??Poiseuille æµè§£?è§£ä¼°ç??„ï?ä½†ç”±?¼ï?
1. Periodic Hill ?„å¹¾ä½•å½¢?€? æ?é¡å?å£“å??å¤±
2. ?å??»ç¶²?¼ç??¸å€¼èª¤å·?3. MRT ç¢°æ??„æ•¸?¼é?æ»?
å¯¦é?æµé€Ÿå¯?½å??¢ç›®æ¨™å€?`Uref`??*?•æ?å¤–å?èª¿æ•´**ä½¿ç”¨æ¯”ä??§åˆ¶?¨è?å¯¦é?æµé€Ÿè¶¨è¿‘ç›®æ¨™æ??Ÿã€?
### ?§åˆ¶?¨å…¬å¼?
```
F_new = F_old + Î² ? (U_target - U_actual) ? U_ref / L
```

| ?ƒæ•¸ | ?ç¾© | ?¸å€?|
|------|------|------|
| `Î²` | ?§åˆ¶å¢ç? | `max(0.001, 3/Re)` |
| `U_target` | ?®æ?æµé€?| `Uref` |
| `U_actual` | å¯¦é?å¹³å?æµé€?| `Ub_avg` |
| `L` | ?¹å¾µ?·åº¦ | `LZ` |

**?§åˆ¶?è¼¯**ï¼?- ??`U_actual < U_target` ??èª¤å·®?ºæ­£ ??å¢å?å¤–å? ??? é€Ÿæ?é«?- ??`U_actual > U_target` ??èª¤å·®?ºè? ??æ¸›å?å¤–å? ??æ¸›é€Ÿæ?é«?
### ?°å??½æ•¸

#### 1. `AccumulateUbulk()` (evolution.h:235-242)

```cpp
void AccumulateUbulk(double* v_field, double* Ub_sum_ptr) {
    for(int j = 3; j < NY6-3; j++) {
        for(int k = 3; k < NZ6-3; k++) {
            int idx = j * NZ6 + k;
            *Ub_sum_ptr += v_field[idx];
        }
    }
}
```

- **?Ÿèƒ½**ï¼šç´¯ç©?Y ?¹å?å¹³å??Ÿåº¦
- **?¼å«?‚æ?**ï¼šæ??‹æ??“æ­¥
- **è¼¸å…¥**ï¼š`v_field` - Y ?¹å??Ÿåº¦??- **è¼¸å‡º**ï¼š`Ub_sum_ptr` - ç´¯ç??„é€Ÿåº¦ç¸½å?ï¼ˆé€é??‡æ??´æ–°ï¼?
#### 2. `ModifyForcingTerm()` (evolution.h:261-278)

```cpp
void ModifyForcingTerm(double* Force, double* Ub_sum_ptr, int NDTFRC) {
    // 1. è¨ˆç??‚é??‡ç©º?“å¹³?‡é€Ÿåº¦
    int num_cells = (NY6 - 6) * (NZ6 - 6);
    double Ub_avg = (*Ub_sum_ptr) / (double)(num_cells * NDTFRC);

    // 2. è¨ˆç??§åˆ¶å¢ç?
    double beta = fmax(0.001, 3.0 / (double)Re);

    // 3. èª¿æ•´å¤–å?
    Force[0] = Force[0] + beta * (Uref - Ub_avg) * Uref / LZ;

    // 4. è¼¸å‡º??§è³‡è?
    printf("Force Update: Ub_avg = %.6f, Uref = %.6f, Force = %.5e\n",
           Ub_avg, Uref, Force[0]);

    // 5. ?ç½®ç´¯å???    *Ub_sum_ptr = 0.0;
}
```

- **?Ÿèƒ½**ï¼šä½¿?¨æ?ä¾‹æ§?¶å™¨èª¿æ•´å¤–å?
- **?¼å«?‚æ?**ï¼šæ? `NDTFRC` æ­¥ï?å»ºè­° 10000 æ­¥ï?
- **è¼¸å…¥**ï¼?  - `Force` - å¤–å????
  - `Ub_sum_ptr` - ç´¯ç??„é€Ÿåº¦ç¸½å?
  - `NDTFRC` - ç´¯ç??„æ??“æ­¥??- **è¼¸å‡º**ï¼?  - ?´æ–°å¾Œç? `Force[0]`
  - ?ç½® `*Ub_sum_ptr = 0`

### ä¸»ç?å¼ä½¿?¨ç?ä¾?
```cpp
// main.cpp
#include "evolution.h"

int main() {
    // ... ?å???...

    double Ub_sum = 0.0;           // ç´¯ç??„å¹³?‡é€Ÿåº¦
    int force_update_count = 0;    // ç´¯ç??„æ??“æ­¥??    const int NDTFRC = 10000;      // æ¯å?å°‘æ­¥ä¿®æ­£ä¸€æ¬?
    for(int t = 0; t < loop; t++) {
        // 1. Stream + Collide
        stream_collide(...);

        // 2. ?±æ??§é???        periodicSW(...);

        // 3. ç´¯ç?å¹³å??Ÿåº¦
        AccumulateUbulk(v, &Ub_sum);
        force_update_count++;

        // 4. æ¯?NDTFRC æ­¥ä¿®æ­????        if(force_update_count >= NDTFRC) {
            ModifyForcingTerm(Force, &Ub_sum, NDTFRC);
            force_update_count = 0;
        }

        // 5. äº¤æ? f_old ??f_new ?‡æ?
        std::swap(f0_old, f0_new);
        std::swap(f1_old, f1_new);
        // ... ?¶ä??¹å? ...
    }
    return 0;
}
```

### ?ƒæ•¸å»ºè­°

| ?ƒæ•¸ | å»ºè­°??| èªªæ? |
|------|--------|------|
| `NDTFRC` | 10000 | å¤ªå??ƒé??ªï?å¤ªå¤§?ƒæ”¶?‚æ…¢ |
| `Î²` ä¸‹é? | 0.001 | ?¿å?é«?Re ?‚æ§?¶å??Šé?å°?|
| `Î²` ä¸Šé? | 3/Re | ä½?Re ?‚è?å¤§ï?é«?Re ?‚è?å°?|

### ?‡å?å§?GPU ?ˆæœ¬?„å·®??
| ?…ç›® | ?Ÿå? GPU ?ˆæœ¬ | ?®å? CPU ?ˆæœ¬ |
|------|--------------|--------------|
| ?Ÿåº¦ç´¯ç? | GPU kernel + cudaMemcpy | CPU ?´æ¥è¿´å? |
| å¤?GPU ?Œæ­¥ | MPI_Reduce + MPI_Bcast | ä¸é?è¦?|
| è¨˜æ†¶é«?| Ub_avg_d[NX6*NZ6] | Ub_sum (ç´”é?) |

---

## 2026-01-19 ?±æ??§é??Œæ?ä»?periodicSW()

### ?Ÿèƒ½èªªæ?

`periodicSW()` å¯¦ç¾ Y ?¹å?ï¼ˆStream-Wiseï¼Œä¸»æµå ´?¹å?ï¼‰ç??±æ??§é??Œæ?ä»¶ã€?
### ç¶²æ ¼?ç½®

```
Y ?¹å?ç´¢å? (j):
  0   1   2  |  3   4   5  ...  NY6-6  NY6-5  NY6-4  |  NY6-3  NY6-2  NY6-1
  ?buffer?? |  ?â??€?€?€?€?€?€?€?€?€ è¨ˆç??€???€?€?€?€?€?€?€?€?€?€?€?€?€?€?? |  ?â??€?€buffer?€?€?€??```

### è¤‡è£½?è¼¯

#### å·¦å´ Buffer å¡«å?ï¼ˆå??³é?è¤‡è£½ï¼?```cpp
int idx_right = (i+NY6-6)*NZ6 + k;   // ä¾†æ?ï¼šj = NY6-6, NY6-5, NY6-4
int buffer_left = i*NZ6 + k;         // ?®æ?ï¼šj = 0, 1, 2
```

#### ?³å´ Buffer å¡«å?ï¼ˆå?å·¦é?è¤‡è£½ï¼?```cpp
int idx_left = (i+3)*NZ6 + k;        // ä¾†æ?ï¼šj = 3, 4, 5
int buffer_right = (i+NY6-3)*NZ6+k;  // ?®æ?ï¼šj = NY6-3, NY6-2, NY6-1
```

### è¤‡è£½?§å®¹

- ?†ä??½æ•¸ï¼š`f0_new ~ f8_new`
- å®è??ï?`v`, `w`, `rho_d`

### å·²ä¿®æ­??é¡?
å°?`k` ?„ç??å? `0 ~ NZ6` ?¹ç‚º `3 ~ NZ6-3`ï¼Œåªè¤‡è£½ Z ?¹å??„æ??ˆè?ç®—å??Ÿã€?
---

## 2026-01-20 æµå ´æ¢ç??é?è¨ºæ–·

### ?é??è¿°

æ¨¡æ“¬çµæ??ºç¾**?œå?æ¢ç??–æ?**ï¼?- æ¢ç??¹å?ï¼šç? 45Â° ?œå?
- æ¢ç??±æ?ï¼šç? 7-8 ?¼ï???7 é»?stencil ?¸é?ï¼?- ?€?‹ï??¸å€¼ç©©å®šï?æ²’æ??¼æ•£ï¼Œä??é¡¯?¯ç³»çµ±æ€§èª¤å·?
### æµå ´?–å??¹å¾µ?†æ?

```
è§€å¯Ÿåˆ°?„ç¾è±¡ï?
1. æ¢ç?æ²?45Â° ?¹å? ???—ç¤º Y ??Z ?¹å??‰æ?ç¨®è€¦å??é?
2. æ¢ç??±æ?ç´?7-8 ?????¯èƒ½??6 ?æ??¼ç? 7 é»?stencil ?‰é?
3. ?¸å€¼ç©©å®???ä¸æ˜¯?¼æ•£?é?ï¼Œæ˜¯ç³»çµ±?§ç?æ¬Šé?/ç´¢å??¯èª¤
```

### å·²æ??¤ç??é?

ç¶“é?è©³ç´°æª¢æŸ¥ï¼Œä»¥ä¸‹é???*ç¢ºè??¡èª¤**ï¼?
#### ??å¹³è¡¡?‹å?ä½ˆå‡½?¸è?ç®?(evolution.h:202-214)
- å¯†åº¦è¨ˆç?æ­?¢ºï¼š`rho_s = Î£ Fi_in`
- ?Ÿåº¦è¨ˆç?æ­?¢ºï¼š`v1 = (F1+F5+F8-F3-F6-F7)/rho_s`
- å¹³è¡¡?‹å…¬å¼æ­£ç¢ºï?ç¬¦å?æ¨™æ? D2Q9 ?¬å?

#### ??MRT ç¢°æ?å·¨é? (MRT_Process.h)
- m_vector å·¨é?ï¼šæ­£ç¢ºè?ç®?`m = M * f`
- meq å·¨é?ï¼šæ­£ç¢ºè?ç®?`meq = M * feq`
- collision å·¨é?ï¼šæ­£ç¢ºè?ç®?`f = f - M_I * S * (m - meq)`

#### ???±æ??Šç?æ¢ä»¶ (evolution.h:237-276)
- å·?buffer `(0,1,2)` ??å¾?`(NY6-6, NY6-5, NY6-4)` è¤‡è£½ ??- ??buffer `(NY6-3, NY6-2, NY6-1)` ??å¾?`(3, 4, 5)` è¤‡è£½ ??
#### ???å???(initialization.h:10-34)
- ?¨å ´?å??–ç‚º `rho=1, v=0, w=0` ??- ?†ä??½æ•¸?¨å¹³è¡¡æ??å?????
#### ??Lagrange ?’å€¼æ??è?ç®?(initializationTool.h:154-162)
- `GetParameter_6th()` è¨ˆç?æ­?¢º??Lagrange ?ºå??½æ•¸

### ?¯èƒ½?„å?é¡Œä?æº?
#### ?”´ ?‘é? 1ï¼š`GetXiParameter` ??stencil èµ·é??¤æ–·ä½¿ç”¨?®æ?é»ç? k

**æª”æ?**ï¼š`initialization.h:79-93`

```cpp
void GetXiParameter(
    double** XiPara_h,    double pos_z,       double pos_y,
    double *Pos_xi,         int IdxToStore,     int k  )
{
    ...
    if( k >= 3 && k <= 6 ){
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, 3 );
    } else if ( k >= NZ6-7 && k <= NZ6-4 ) {
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, NZ6-10 );
    } else {
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, k-3 );
    }
}
```

**?é?**ï¼šå‚³?¥ç? `k` ??*?®æ?é»?*?„ç´¢å¼•ï?ä½†å???streaming ä¾†èªªï¼Œæ?è©²ç”¨**ä¾†æ?é»?*?„ä?ç½®ä?æ±ºå? stencil èµ·é???
ä¾‹å? F2 (+Z)ï¼?- ?®æ?é»?`(j, k)`
- ä¾†æ?é»?`(j, k-1)`ï¼ˆZ ?¹å?å¾€ä¸‹ä??¼ï?
- ç¨‹å???`k`ï¼ˆç›®æ¨™é?ï¼‰åˆ¤??stencil èµ·é?
- ?‰è©²?¨ä?æºé??„ä?ç½®ä??¤æ–·

**æ½›åœ¨å½±éŸ¿**ï¼šç•¶ä¾†æ?é»å??®æ?é»è·¨è¶?stencil ?Šç??‚ï?å¦?`k=7` ?®æ?ä½†ä?æºé??‰è©²??`start=3`ï¼‰ï??ƒä½¿?¨éŒ¯èª¤ç? stencil èµ·é???
#### ?”´ ?‘é? 2ï¼š`F1_Intrpl7` / `F3_Intrpl7` ?„æ??¼é?åºè?æ¬Šé??¯èƒ½ä¸åŒ¹??
**æª”æ?**ï¼š`interpolationHillISLBM.h:14-36`

`F1_Intrpl7` ??`F3_Intrpl7` ?„æ??¼é?åºï?
```
å¤–å±¤: Y ?¹å??’å€?(y_0..y_6)
?§å±¤: Xi ?¹å??’å€?(xi_0..xi_6)
```

`Y_XI_Intrpl7` (?¨æ–¼ F5-F8) ?„æ??¼é?åºï?
```
å¤–å±¤: Xi ?¹å??’å€?(xi_0..xi_6)
?§å±¤: Y ?¹å??’å€?(y_0..y_6)
```

**?©è€…é?åºç›¸??*ï¼é??¶ç?è«–ä? 2D å¼µé?ç©æ??¼ç??œè??†å??¡é?ï¼Œä?å¦‚æ?æ¬Šé?è¨ˆç??¹å??‡ä½¿?¨æ–¹å¼ä?ä¸€?´ï??¯èƒ½?¢ç??é???
#### ?”´ ?‘é? 3ï¼šCLAUDE.md ?ƒè€ƒç??¬è??¾æ?ç¨‹å???Xi æ¬Šé?è¨ˆç?ä¸ä???
**CLAUDE.md:438-441 ?„å??ƒç???*ï¼?```cpp
GetXiParameter( XiParaF3_h,  z_h[j*NZ6+k],         y_h[j]-minSize, xi_h, j*NZ6+k, k );
GetXiParameter( XiParaF4_h,  z_h[j*NZ6+k],         y_h[j]+minSize, xi_h, j*NZ6+k, k );
GetXiParameter( XiParaF5_h,  z_h[j*NZ6+k]-minSize, y_h[j],         xi_h, j*NZ6+k, k );
```

**?¾æ?ç¨‹å? (initialization.h:115-129)**ï¼?```cpp
// F1: pos_y = y_global[j]-minSize
GetXiParameter( XiParaF1_h,  z_global[j*NZ6+k],         y_global[j]-minSize, xi_h, j*NZ6+k, k );
// F5: pos_y = y_global[j]-minSize  ???‡å??ƒç??¬ä??Œï??ƒè€ƒç??¬æ˜¯ y_h[j]
GetXiParameter( XiParaF5_h,  z_global[j*NZ6+k]-minSize, y_global[j]-minSize, xi_h, j*NZ6+k, k );
```

**å·®ç•°**ï¼?- ?ƒè€ƒç???F5 ??`y_h[j]`ï¼ˆä??ç§»ï¼?- ?¾æ?ç¨‹å? F5 ??`y_global[j]-minSize`ï¼ˆæ??ç§»ï¼?
**?©ç??ç¾©**ï¼šF5 (+Y,+Z) ?„ä?æºé???`(y-?, z-?)`ï¼Œæ?ä»?`pos_y = y-?` ?¯å??†ç??‚ä??™è??ƒè€ƒç??¬ä??Œï??€è¦ç¢ºèªå“ª?‹æ˜¯æ­?¢º?„ã€?
### å»ºè­°?„è¨º?·æ­¥é©?
#### æ­¥é? 1ï¼šè¼¸?ºæ??ç¸½?Œé?è­?
?¨æ??‹ç‰¹å®?`(j, k)` ä½ç½®ï¼Œè¼¸?ºæ???7 ?‹æ??ä¸¦æª¢æŸ¥ç¸½å??¯å¦??1ï¼?
```cpp
// ? åœ¨ GetIntrplParameter_Xi() çµå°¾
if(j == 64 && k == 35) {
    double sum = 0.0;
    for(int i = 0; i < 7; i++) {
        sum += XiParaF1_h[i][j*NZ6+k];
        printf("XiParaF1_h[%d] = %e\n", i, XiParaF1_h[i][j*NZ6+k]);
    }
    printf("Sum of weights = %f (should be 1.0)\n", sum);
}
```

#### æ­¥é? 2ï¼šè¼¸??pos_xi ç¯„å?é©—è?

ç¢ºè? `pos_xi` ?¨å??†ç??å…§ï¼? ??LXiï¼‰ï?

```cpp
// ? åœ¨ GetXiParameter() ??if(pos_xi < 0 || pos_xi > LXi) {
    printf("WARNING: pos_xi = %f out of range [0, %f] at k=%d\n", pos_xi, LXi, k);
}
```

#### æ­¥é? 3ï¼šç°¡?–æ¸¬è©?- ?¹ç”¨?‡å‹»ç¶²æ ¼

?«æ???Z ?¹å??¹æ??‡å‹»ç¶²æ ¼ï¼Œç?æ¢ç??¯å¦æ¶ˆå¤±ï¼?
```cpp
// ??GenerateMesh_Z() ä¸­ï??«æ??¹ç”¨ï¼?for( int k = 3; k < NZ6-3; k++ ){
    z_global[j*NZ6+k] = HillFunction(y_global[j]) + (LZ - HillFunction(y_global[j])) * (k-3) / (NZ6-7);
}
```

#### æ­¥é? 4ï¼šå–®?¹å?æ¸¬è©¦

?«æ??ªåŸ·è¡?F0ï¼ˆé?æ­¢ï???F2/F4ï¼ˆç? Z ?¹å?ï¼‰ï??œé? Y ?¹å??Œæ??‘æ–¹?‘ï?

```cpp
// ??stream_collide ä¸­æš«?‚è¨»è§??ï¼?// F1_Intrpl7(...);
// F3_Intrpl7(...);
// Y_XI_Intrpl7(f5_old, ...);
// Y_XI_Intrpl7(f6_old, ...);
// Y_XI_Intrpl7(f7_old, ...);
// Y_XI_Intrpl7(f8_old, ...);

// ?¹ç‚º?´æ¥è¤‡è£½ï¼?F1_in = f1_old[idx_xi];
F3_in = f3_old[idx_xi];
F5_in = f5_old[idx_xi];
F6_in = f6_old[idx_xi];
F7_in = f7_old[idx_xi];
F8_in = f8_old[idx_xi];
```

å¦‚æ?æ¢ç?æ¶ˆå¤±ï¼Œå?é¡Œåœ¨ Y ?¹å??–æ??‘æ??¼ï?å¦‚æ?ä»æ?æ¢ç?ï¼Œå?é¡Œåœ¨ Z ?¹å??’å€¼ã€?
### å¾…ä¿®å¾©é??®æ???
| ?ªå?ç´?| ?…ç›® | ?€??|
|-------|------|------|
| é«?| ç¢ºè? `GetXiParameter` ??stencil èµ·é??¤æ–·?è¼¯ | ??å¾…é?è­?|
| é«?| æ¯”å? CLAUDE.md ?ƒè€ƒç??¬è??¾æ?ç¨‹å???Xi æ¬Šé?è¨ˆç? | ??å¾…ç¢ºèª?|
| ä¸?| é©—è? `F1_Intrpl7` / `F3_Intrpl7` ?„æ??¼é?åºæ­£ç¢ºæ€?| ??å¾…é?è­?|
| ä½?| æª¢æŸ¥ `cell_z` è¨ˆç???`GetXiParameter` ?„ä??´æ€?| ??å¾…ç¢ºèª?|

### ?¸é?æª”æ?ç´¢å?

| æª”æ? | ?œéµè¡Œè? | ?§å®¹ |
|------|---------|------|
| `interpolationHillISLBM.h` | 14-24 | `F1_Intrpl7` å·¨é? |
| `interpolationHillISLBM.h` | 26-36 | `F3_Intrpl7` å·¨é? |
| `interpolationHillISLBM.h` | 46-56 | `Y_XI_Intrpl7` å·¨é? |
| `initialization.h` | 79-93 | `GetXiParameter` ?½æ•¸ |
| `initialization.h` | 104-130 | `GetIntrplParameter_Xi` ?½æ•¸ |
| `evolution.h` | 123-131 | ?’å€¼å·¨?†å‘¼??|
| `initializationTool.h` | 154-162 | `GetParameter_6th` ?½æ•¸ |

### è£œå?ï¼šD2Q9 Streaming ?¹å??‡ä?æºé?å°æ?

| ?Ÿåº¦?¹å? | ?Ÿåº¦?‘é? | Streaming ä¾†æ?é»?| ?™è¨» |
|---------|---------|-----------------|------|
| F0 | (0, 0) | ?Ÿåœ° | ä¸é??’å€?|
| F1 | (+1, 0) | (y-?, z) | Y ?¹å??ç§» |
| F2 | (0, +1) | (y, z-?) | Z ?¹å??ç§» |
| F3 | (-1, 0) | (y+?, z) | Y ?¹å??ç§» |
| F4 | (0, -1) | (y, z+?) | Z ?¹å??ç§» |
| F5 | (+1, +1) | (y-?, z-?) | ?œå??ç§» |
| F6 | (-1, +1) | (y+?, z-?) | ?œå??ç§» |
| F7 | (-1, -1) | (y+?, z+?) | ?œå??ç§» |
| F8 | (+1, -1) | (y-?, z+?) | ?œå??ç§» |

## 2026-01-20 ç¨‹å?ç¢¼æ”¹??

```cpp
//çµ¦æ?ä¸€?‹ç·¨?Ÿï??¢ç?è©²Y?¼æ?å°æ??„ä??‹ç„¡? æ¬¡?–åº§æ¨?void RelationXi(double nonintindex, double L , double MinSize , double a , int N , double* RelationXi){
    int i = (int)floor(nonintindex);
    RelationXi[0] = tanhFunction( L , MinSize , a, i-3 , N) - MinSize/2.0;
    RelationXi[1] = tanhFunction( L , MinSize , a, i-2 , N) - MinSize/2.0;
    RelationXi[2] = tanhFunction( L , MinSize , a, i-1 , N) - MinSize/2.0;
    RelationXi[3] = tanhFunction( L , MinSize , a, i , N) - MinSize/2.0;
    RelationXi[4] = tanhFunction( L , MinSize , a, i+1 , N) - MinSize/2.0;
    RelationXi[5] = tanhFunction( L , MinSize , a, i+2 , N) - MinSize/2.0;
    RelationXi[6] = tanhFunction( L , MinSize , a, i+3 , N) - MinSize/2.0;
}



void GetParameter_6th(
    double *Para_h[7],      double Position,
    double *Pos,            int i,              int n  )
{
    Para_h[0][i] = Lagrange_6th(Position, Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+5], Pos[n+6]);
    Para_h[1][i] = Lagrange_6th(Position, Pos[n+1], Pos[n],   Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+5], Pos[n+6]);
    Para_h[2][i] = Lagrange_6th(Position, Pos[n+2], Pos[n],   Pos[n+1], Pos[n+3], Pos[n+4], Pos[n+5], Pos[n+6]);
    Para_h[3][i] = Lagrange_6th(Position, Pos[n+3], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+4], Pos[n+5], Pos[n+6]);
    Para_h[4][i] = Lagrange_6th(Position, Pos[n+4], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+5], Pos[n+6]);
    Para_h[5][i] = Lagrange_6th(Position, Pos[n+5], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+6]);
    Para_h[6][i] = Lagrange_6th(Position, Pos[n+6], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+5]);
}
void GetParameter_6th2(double** XiPara , double pos_z ,  double* RelationXi , int r , int index_xi){
    XiPara[0][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[0],  RelationXi[1],  RelationXi[2] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[1][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[1],  RelationXi[0],  RelationXi[2] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[2][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[2],  RelationXi[0],  RelationXi[1] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[3][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[3],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[4][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[4],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[5], RelationXi[6]); 
    XiPara[5][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[5],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[4], RelationXi[6]); 
    XiPara[6][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[6],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[4], RelationXi[5]);    
}//pos_xi?ºæ?ç®—é?å¾Œç??¡å?æ¬¡å?Zåº§æ? 




void GetXiParameter(double* XiPara_h[7], double pos_z, double pos_y, int index_xi, int j, int k ) 
{   
    //è¶Šç??²å?
    if(j<3) j = j + NY6-7 ; 
    if(j>NY6-4) j = j - (NY6-7) ;
    //Z?¹å?ä¿‚æ•¸ç·¨è?ç³»çµ±(ç¬¬ä???  = (j+relation)*NZ6 + k  ; relation:1~7
    double L = LZ - HillFunction(pos_y) - minSize;
    //æ¯ä??‹yä½ç½®?Œè?æ¯ä??‹è?ç®—é??½ä?ä¸€æ¨?
    double pos_xi = LXi * (pos_z - (HillFunction(pos_y)+minSize/2.0)) / L;
    double a = GetNonuniParameter();
     //æ±‚è§£?®æ?é»ç?ç·¨è?(? å??å??»ç³»çµ?
    double j_cont = Inverse_tanh_index( pos_xi , L , minSize , GetNonuniParameter() , (NZ6-7) );
    //Z?¹å??å“¡èµ·å?ç·¨è?
    int cell_z1 = k-3;
    if( k <= 6 ) cell_z1 = 3;
    if( k >= NZ6-7 ) cell_z1 = NZ6-10; 
    //çµ¦æ?ä¸€?‹double çµ¦å‡º?¸å??‰ä??‹å…§?’ç??¡å?æ¬¡å?åº§æ?
    double RelationXi_0[7] ; //(j-3)
    double LT = LZ - HillFunction(y_global[j-3]) - minSize;
    double pos_z =  tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_0);//å¯«å…¥ç¬¬ä?å¥—å??´æ–¹?‘ç„¡? æ¬¡?–Z???
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_0 , 0 , index_xi); //XiPara_h[0][index_xi + 0*NZ6] ~ XiPara_h[6][index_xi + 0*NZ6]
    
    double RelationXi_1[7] ; //(j-2)
    LT = LZ - HillFunction(y_global[j-2]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_1);//å¯«å…¥ç¬¬ä?å¥—å??´æ–¹?‘ç„¡? æ¬¡?–Z???
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_1 , 1 , index_xi); //XiPara_h[0][index_xi + 1*NZ6] ~ XiPara_h[6][index_xi + 1*NZ6]
    
    double RelationXi_2[7] ; //(j-1)
    LT = LZ - HillFunction(y_global[j-1]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_2);//å¯«å…¥ç¬¬ä?å¥—å??´æ–¹?‘ç„¡? æ¬¡?–Z???
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_2 , 2 , index_xi); //XiPara_h[0][index_xi + 2*NZ6] ~ XiPara_h[6][index_xi + 2*NZ6]
    
    double RelationXi_3[7] ; //(j)
    LT = LZ - HillFunction(y_global[j]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_3);//å¯«å…¥ç¬¬å?å¥—å??´æ–¹?‘ç„¡? æ¬¡?–Z???
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_3 , 3 , index_xi); //XiPara_h[0][index_xi + 3*NZ6] ~ XiPara_h[6][index_xi + 3*NZ6]

    double RelationXi_4[7] ; //(j+1)
    LT = LZ - HillFunction(y_global[j+1]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_4);//å¯«å…¥ç¬¬ä?å¥—å??´æ–¹?‘ç„¡? æ¬¡?–Z???
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_4 , 4 , index_xi); //XiPara_h[0][index_xi + 4*NZ6] ~ XiPara_h[6][index_xi + 4*NZ6]

    double RelationXi_5[7] ; //(j+2)
    LT = LZ - HillFunction(y_global[j+2]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_5);//å¯«å…¥ç¬¬å…­å¥—å??´æ–¹?‘ç„¡? æ¬¡?–Z???
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_5 , 5 , index_xi); //XiPara_h[0][index_xi + 5*NZ6] ~ XiPara_h[6][index_xi + 5*NZ6]
   
    double RelationXi_6[7] ; //(j+3)
    LT = LZ - HillFunction(y_global[j+3]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_6);//å¯«å…¥ç¬¬ä?å¥—å??´æ–¹?‘ç„¡? æ¬¡?–Z???
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_6 , 6 , index_xi); //XiPara_h[0][index_xi + 6*NZ6] ~ XiPara_h[6][index_xi + 6*NZ6]
}
```
---

## 2026-01-23 Xi æ¬Šé?è¦†å¯«?é?ï¼ˆä??¥è?è«–ï?

### ?é?ä¾†æ?ï¼ˆé??µç?å¼æ®µï¼?
1) `initializationTool.h` / `GetParameter_6th2(...)`

```cpp
XiPara[2][index_xi + r*NZ6] = Lagrange_6th(...);
```

2) `initialization.h` / `GetXiParameter(...)`

```cpp
GetParameter_6th2(XiPara_h, pos_z2, RelationXi_0, 0, index_xi);
GetParameter_6th2(XiPara_h, pos_z2, RelationXi_1, 1, index_xi);
...
GetParameter_6th2(XiPara_h, pos_z2, RelationXi_6, 6, index_xi);
```

3) `initialization.h` / `GetIntrplParameter_Xi()` ??`BFLInitialization()`

```cpp
GetXiParameter(XiParaF1_h, ...);
...
GetXiParameter(XiParaF8_h, ...);
```

### è¦†å¯«?¯ä?éº¼æ???- ä¸å?ç©ºé?é»å¯«?°å?ä¸€?‹æ??ä?ç½®ï?å¾Œå¯«è¦†è??ˆå¯«ï¼ˆä??¯è??Œï???
### ?ºä?éº¼æ?è¦†å¯«ï¼ˆä?å­ï?
? ç‚º `index_xi = j*NZ6 + k`ï¼Œå¯¦?›å¯«?¥ä?ç½®æ˜¯ï¼?
```
write_idx = index_xi + r*NZ6 = (j + r)*NZ6 + k
```

?€ä»¥æ??¼ç?ï¼?
```
(j=3, k=3, r=1) -> write_idx = (4,3)
(j=4, k=3, r=0) -> write_idx = (4,3)
```

? æ­¤ `(j=3, r=1)` ?„æ??æ?è¢?`(j=4, r=0)` è¦†è???
### çµè?
- ?®å???`index_xi + r*NZ6` ?ªæ˜¯??r ? å???j ä½ç§»ï¼Œç??Œæ? 7 å±¤æ??ç??¨å?ä¸€å¼?2D ?¼å?ä¸Šï??ƒç³»çµ±æ€§è?å¯«ã€?

## 2026-01-23 ?ºé¿?è??¶é??†ç‚¸ï¼ŒY?¹å??ªç”¨ä¸‰é??¼å??é?
1.?ˆèª¿??RelaxationXi ??? 

---

## 2026-01-24 ?Šç??•ç?ä¿®æ­£?‡æ•¸?¼ç©©å®šæ€§å·¥ä½?
### ?é?è¨ºæ–·?ç?

#### 1. ?å??é?ï¼šç?å¼ç™¼??(t=10 ?‹å??ºç¾ nan)
- **?‡ç?**ï¼š`rho = -54151.9` (t=10)ï¼Œè??Ÿç™¼??- **?¤éŒ¯è¼¸å‡º**ï¼šåœ¨ k=152 (ä¸Šé???NZ6-4) ?¼ç¾ F4 ?ºè???(-1.38293)
- **?¹æœ¬?Ÿå?**ï¼šF4 (-Z?¹å?) ?„ä?æºä?ç½®æ˜¯ z+?ï¼Œç•¶ k=152 ?‚æ?è¶…å‡ºè¨ˆç???
#### 2. cellZ_F* è¨ˆç??é?
- **?é?**ï¼š`GetXiParameter` ä¸?`cell_z1 = k-3` ä½¿ç”¨?„æ˜¯?®æ?é»?kï¼Œè€Œé?ä¾†æ?é»?- **ä¿®æ­£**ï¼šæ”¹?ºåŸº??`j_cont`ï¼ˆæ?å°„å??„é€??åº§æ?ï¼‰è?ç®?```cpp
// ?Šç?
int cell_z1 = k-3;

// ?°ç?  
int k_source = static_cast<int>(std::round(j_cont)) + 3;
int cell_z1 = k_source - 3;
if(cell_z1 < 0) cell_z1 = 0;
if(cell_z1 > NZ6-7) cell_z1 = NZ6-7;
```

#### 3. ?Šç??€?Ÿç??’å€¼å??¨å?é¡?- **?é?**ï¼?-point stencil ?¨é??Œé?è¿‘æ?å­˜å? buffer ?€?Ÿï?å°è‡´ Lagrange å¤–æ¨?¢ç?æ¥µç«¯æ¬Šé?
- **å½±éŸ¿ç¯„å?**ï¼?  - ä¸‹é???(k ??5)ï¼šF2, F5, F6 ?„ä?æºå¯?½è???  - ä¸Šé???(k ??NZ6-6=150)ï¼šF4, F7, F8 ?„ä?æºå¯?½è???  - F1, F3 ??Xi stencil ä¹Ÿæ??—å½±??
### ?€çµ‚è§£æ±ºæ–¹æ¡ˆï??Šç??€?Ÿä½¿?¨ç°¡??streaming

```cpp
// evolution.h ä¸­ç??Šç??•ç??è¼¯
if( k <= 5 ) {
    // ä¸‹é??Œé?è¿‘ï?ä½¿ç”¨ç°¡å–® streaming ??bounce-back
    F1_in = f1_old[(j-1)*NZ6 + k];
    F3_in = f3_old[(j+1)*NZ6 + k];
    F2_in = f4_old[idx_xi];  // bounce-back
    F5_in = f7_old[idx_xi];  // bounce-back
    F6_in = f8_old[idx_xi];  // bounce-back
    F4_Intrpl7(...);  // æ­?¸¸?’å€¼ï?? é›¢ä¸‹é??Œï?
    F7_in = f7_old[(j+1)*NZ6 + k+1];
    F8_in = f8_old[(j-1)*NZ6 + k+1];
    
} else if( k >= NZ6-6 ) {
    // ä¸Šé??Œé?è¿‘ï?ä½¿ç”¨ç°¡å–® streaming ??bounce-back
    F1_in = f1_old[(j-1)*NZ6 + k];
    F3_in = f3_old[(j+1)*NZ6 + k];
    F2_in = f2_old[j*NZ6 + k-1];      // ç°¡å–® streaming
    F5_in = f5_old[(j-1)*NZ6 + k-1];  // ç°¡å–® streaming
    F6_in = f6_old[(j+1)*NZ6 + k-1];  // ç°¡å–® streaming
    F4_in = f2_old[idx_xi];  // bounce-back
    F7_in = f5_old[idx_xi];  // bounce-back
    F8_in = f6_old[idx_xi];  // bounce-back
    
} else {
    // ?§éƒ¨?€?Ÿï?æ­?¸¸ ISLBM ?’å€?    F1_Intrpl3(...); F3_Intrpl3(...);
    F2_Intrpl7(...); F4_Intrpl7(...);
    Y_XI_Intrpl3(...);  // F5~F8
}
```

### ä¿®æ”¹?„æ?æ¡ˆæ???
| æª”æ? | ä¿®æ”¹?§å®¹ |
|------|----------|
| `initialization.h` | `GetXiParameter` ä¸?cell_z1 è¨ˆç??¹ç”¨ j_cont |
| `evolution.h` | ?°å??Šç??€?Ÿæ?ä»¶åˆ¤?·ï?ä½¿ç”¨ç°¡å? streaming ?¿ä»£?’å€?|

### æ¸¬è©¦çµæ?
- ??ç¨‹å?ç©©å??‹è???t=1000+
- ??å¯†åº¦?¶æ???~0.988
- ? ï? ?ºç¾ **Checkerboard instability**ï¼ˆæ??¤æ ¼?¯ç›ªï¼?
### Checkerboard ?é??†æ?

**è§€å¯?*ï¼?- Y ?¹å?ä¸ƒé??’å€???æ¢ç?è¼ƒç?
- Y ?¹å?ä¸‰é??’å€???æ¢ç?è¼ƒå?

**?¯èƒ½?Ÿå?**ï¼?1. å¥‡å¶è§?€?(odd-even decoupling)
2. ?å??»ç¶²?¼ç??¸å€¼è‰²??3. Lagrange ?’å€¼ç??¯ç›ª?¹æ€?
**å¾…å?è©¦è§£æ±ºæ–¹æ¡?*ï¼?- [ ] Z ?¹å?ä¹Ÿé??ºä??æ???- [ ] èª¿æ•´ MRT é¬†å??ƒæ•¸
- [ ] ? å…¥?¸å€¼æ¿¾æ³¢å™¨

### ä¸‹æ¬¡å·¥ä??¹å?

1. **Z ?¹å??é?æ¸¬è©¦**
   - ?°å? `F2_Intrpl3` ??`F4_Intrpl3` å·¨é?
   - ä¿®æ”¹ `evolution.h` ä½¿ç”¨ä¸‰é??ˆæœ¬
   
2. **è©•ä¼°?é??ˆæ?**
   - ??checkerboard æ¶ˆå¤±ä¸”ç??œå?????ç¹¼ç?ä½¿ç”¨
   - ?¥ç²¾åº¦ä?è¶????ƒæ…®?¶ä??‘åˆ¶?¹æ?

### ?®å??ç½®?˜è?

| ?…ç›® | è¨­å???|
|------|--------|
| Y ?¹å??’å€?| 3 ??Lagrange |
| Z ?¹å??’å€?| 7 ??Lagrange |
| ä¸‹é??Œè???(k??) | ç°¡å? streaming + bounce-back |
| ä¸Šé??Œè???(k??50) | ç°¡å? streaming + bounce-back |
| ?§éƒ¨?€??| ISLBM ?’å€?|
| Re | 1 |
| CFL | 0.2 |
| ç¶²æ ¼ | NY=200, NZ=150 (??buffer: NY6=207, NZ6=156) |

---

## 2026-01-28 ?å¤§çªç ´ï¼šå???Streaming ?Šç???Re=500 ç©©å??‹è?

### ä»Šæ—¥?¸å??æ?

| ?…ç›® | ?€??| èªªæ? |
|------|------|------|
| **Re=500 ç©©å??‹è?** | ??| ?å?è·‘åˆ° t=39000+ ?¡ç™¼??|
| **?•æ? Streaming ?Šç?** | ??| æ¼¸é€²å??´å¤§è§??å±¤è¨­è¨ˆå???|
| **è³ªé??šé?èª¤å·®** | ??| ?¹ç”¨ div(?u) è¨ˆç?å±€?¨è³ª?å???|

---

### 1. ?å¤§?¼ç¾ï¼šStreaming Layer ??Interpolation ?Šç??„ç¨ç«‹æ€?
#### ?œéµæ´å?

ä¹‹å??¯èª¤èªç‚º `streaming_lower` å¿…é? >= `interpolation_lower`ï¼Œä?å¯¦é?ä¸?*?©è€…æ˜¯?¨ç???*ï¼?
| ?ƒæ•¸ | ä½œç”¨ | è¨ˆç??‚æ? | ?¯å¦?•æ?èª¿æ•´ |
|------|------|----------|--------------|
| `interpolation_lower/upper` | æ±ºå? XiPara[] ?å??„æ??¼æ??é??‹ï?3é»æ?7é»ï?| **?å??–æ?**?è?ç®?| ?¦ï??€?ç?æ¬Šé?ï¼‰|
| `streaming_lower/upper` | æ±ºå??¯å¦è·³é??’å€¼æ”¹??streaming | **æ¯å€‹æ??“æ­¥**?³æ??¤æ–· | ???¯ä»¥?•æ?èª¿æ•´ |

#### ?è¼¯èªªæ?

```cpp
// interpolation_lower = 25 ?å‘³?—ï?
// k < 25ï¼šXiPara å­˜ã€Œä?é»æ??¼ã€æ???// k >= 25ï¼šXiPara å­˜ã€Œä?é»æ??¼ã€æ???
// streaming_lower ?å‘³?—ï?
// k <= streaming_lowerï¼šè·³?æ??¼ï???streamingï¼ˆæ?ä¿å?ï¼?// k > streaming_lowerï¼šè???XiPara[] ?šæ???```

#### ?è?ä¿®æ­£

**ä¹‹å??„éŒ¯èª¤è???*ï¼š`streaming_lower >= interpolation_lower`ï¼ˆå??ˆï?

**æ­?¢º?†è§£**ï¼š`streaming_lower` ?¯ä»¥ < `interpolation_lower`ï¼?- ?™æ¨£ä¸‰é??’å€¼å?å°±èƒ½ä½œç‚º**ç©©å?ç·©è??€**?¼æ®ä½œç”¨
- æ¼¸é€²é??¾é?åºï?streaming ??ä¸‰é??’å€???ä¸ƒé??’å€?
---

### 2. ?•æ? Streaming ?Šç?è¨­è?

#### è¨­è??Ÿç?

```
?‚é??²å? ??
t=0:      |===streaming (k??0)===|---ä¸ƒé??’å€?--|===streaming===|
                                                 
t=50000:  |==streaming (k??0)==|--ä¸‰é?--|--ä¸ƒé?--|--ä¸‰é?--|==streaming==|

t=100000: |=s(k??0)=|--ä¸‰é?--|-------ä¸ƒé??’å€?(å®Œæ•´)-------|--ä¸‰é?--|=s=|
                     ??                                     ??                ç·©è??€?‹æ”¾ï¼?                           ç·©è??€?‹æ”¾ï¼?```

#### ?€çµ‚è§£?å±¤çµæ? (t >= 100000)

```
k=0    k=10       k=25                    k=236      k=255  k=262
|======|==========|========================|==========|======|
stream   ä¸‰é??’å€?       ä¸ƒé??’å€?         ä¸‰é??’å€?   stream
 (10å±?  (15å±¤ç·©è¡?       (ä¸»è§£?å?)         (19å±¤ç·©è¡?  (7å±?
```

#### ?ƒæ•¸è¨­å?ï¼ˆvariables.hï¼?
```cpp
//=== ?•æ? Streaming ?Šç??ƒæ•¸ï¼ˆæ¼¸?²å??´å¤§è§??å±¤ï?===//
// ?ºå??„æ??¼é??Œï?æ±ºå?æ¬Šé?é¡å?ï¼Œå?å§‹å??‚è?ç®—ï?
#define     interpolation_lower  (25)                // ä¸ƒé??§æ?ä¸‹ç?
#define     interpolation_upper  (NZ6-26)            // ä¸ƒé??§æ?ä¸Šç?

// ?å??¼ï?ä¿å?ï¼Œæ›´å¤§ç? streaming ?€?Ÿï?
#define     streaming_lower_init     (50)            // ?å?ä¸‹ç? (k <= 50 ??streaming)
#define     streaming_upper_init     (NZ6-51)        // ?å?ä¸Šç? (k >= NZ6-51 ??streaming)

// ?®æ??¼ï?æ¯?interpolation ?´æ??²ï??‹æ”¾ä¸‰é??’å€¼ç·©è¡å?ï¼?#define     streaming_lower_target   (10)            // ?®æ?ä¸‹ç?
#define     streaming_upper_target   (NZ6-7)         // ?®æ?ä¸Šç?

// ?æ¸¡?‚é?è¨­å?
#define     ramp_start_time          (0)             // ?‹å?æ¼¸é€²ç??‚é?æ­?#define     ramp_end_time            (100000)        // å®Œæ?æ¼¸é€²ç??‚é?æ­?
// ?¨å?è®Šæ•¸å®??ï¼ˆåœ¨ main.cpp ä¸­å?ç¾©ï?
extern int streaming_lower;  // ?•æ?ä¸‹ç?
extern int streaming_upper;  // ?•æ?ä¸Šç?
```

#### ?´æ–°?½æ•¸å¯¦ä?ï¼ˆmain.cppï¼?
```cpp
//-----------------------------------------------------------------------------
// 2.12 ?•æ? Streaming ?Šç?ï¼ˆæ¼¸?²å??´å¤§è§??å±¤ï?
//-----------------------------------------------------------------------------
int streaming_lower = streaming_lower_init;  // ?•æ?ä¸‹ç?ï¼Œå?å§‹ç‚ºä¿å???int streaming_upper = streaming_upper_init;  // ?•æ?ä¸Šç?ï¼Œå?å§‹ç‚ºä¿å???
// ä½¿ç”¨ tanh å¹³æ??æ¸¡?´æ–° streaming ?Šç?
void UpdateStreamingBounds(int t) {
    if (t >= ramp_end_time) {
        // ?æ¸¡å®Œæ?ï¼Œä½¿?¨ç›®æ¨™å€?        streaming_lower = streaming_lower_target;
        streaming_upper = streaming_upper_target;
    } else if (t <= ramp_start_time) {
        // å°šæœª?‹å?ï¼Œä½¿?¨å?å§‹å€?        streaming_lower = streaming_lower_init;
        streaming_upper = streaming_upper_init;
    } else {
        // ?æ¸¡?Ÿï?ä½¿ç”¨ tanh å¹³æ??æ¸¡
        double progress = (double)(t - ramp_start_time) / (ramp_end_time - ramp_start_time);
        // tanh å¹³æ?ï¼šå? [0,1] ? å???[0,1]ï¼Œä?ä¸­é??æ¸¡?´å¹³æ»?        double smooth_ratio = 0.5 * (1.0 + tanh(6.0 * (progress - 0.5)));
        
        // è¨ˆç??¶å??Šç???        streaming_lower = streaming_lower_init - 
            (int)(smooth_ratio * (streaming_lower_init - streaming_lower_target));
        streaming_upper = streaming_upper_init + 
            (int)(smooth_ratio * (streaming_upper_target - streaming_upper_init));
    }
}
```

#### ?‚é?è¿´å?ä¸­ç??¼å«ï¼ˆmain.cppï¼?
```cpp
for(t = 0; t < loop; t++) {
    // ?´æ–°?•æ? streaming ?Šç?ï¼ˆæ¼¸?²å??´å¤§è§??å±¤ï?
    UpdateStreamingBounds(t);
    
    // ... ?¶é??‚é?è¿´å??è¼¯ ...
}
```

#### è¼¸å‡º??§ï¼ˆmain.cppï¼?
```cpp
// è¼¸å‡º?¶å? streaming ?Šç?ï¼ˆé¡¯ç¤ºè§£?å±¤?´å¤§?²åº¦ï¼?cout << "[Streaming Bounds] lower=" << streaming_lower 
     << " upper=" << streaming_upper 
     << " (target: " << streaming_lower_target << "/" << streaming_upper_target << ")" << endl;
```

---

### 3. è³ªé??šé?èª¤å·®è¨ˆç??¹é€²ï?evolution.hï¼?
#### ?Šç?ï¼šå?åº¦å??¹æ?

```cpp
// ?Šç? - ç°¡å–®å¯†åº¦??double ComputeMaxLocalMassError(double* rho_d) {
    double max_error = 0.0;
    for(int j = 4; j < NY6-4; j++) {
        for(int k = streaming_lower+1; k < streaming_upper-1; k++) {
            double local_mass = 0.0;
            for(int dj = -1; dj <= 1; dj++) {
                for(int dk = -1; dk <= 1; dk++) {
                    local_mass += rho_d[(j + dj) * NZ6 + (k + dk)];
                }
            }
            double error = std::fabs(local_mass - 9.0) / 9.0;
            if(error > max_error) max_error = error;
        }
    }
    return max_error;
}
```

#### ?°ç?ï¼šè³ª?é€šé???º¦?¹æ?

```cpp
//==========================================
//4.è¨ˆç??€?Ÿè³ª?å???- è³ªé??šé?æ·¨æ??¥æ??ºå·®
// å°æ??‹æ ¼é»è?ç®? |æµå‡ºè³ªé??šé? - æµå…¥è³ªé??šé?|
// è³ªé??šé? = ?*v (Y?¹å?) + ?*w (Z?¹å?)
// ä½¿ç”¨ä¸­å?å·®å?: div(?u) ??[(?v)_{j+1} - (?v)_{j-1}]/2 + [(?w)_{k+1} - (?w)_{k-1}]/2
//==========================================
double ComputeMassFluxError(double* rho_d, double* v_field, double* w_field) {
    double total_error = 0.0;
    
    for(int j = 4; j < NY6-4; j++) {
        for(int k = 4; k < NZ6-4; k++) {
            int idx = j * NZ6 + k;
            int idx_jp1 = (j+1) * NZ6 + k;  // j+1
            int idx_jm1 = (j-1) * NZ6 + k;  // j-1
            int idx_kp1 = j * NZ6 + (k+1);  // k+1
            int idx_km1 = j * NZ6 + (k-1);  // k-1
            
            // Y?¹å?è³ªé??šé?å·?(?v)_{j+1} - (?v)_{j-1}
            double flux_y = (rho_d[idx_jp1] * v_field[idx_jp1]) - (rho_d[idx_jm1] * v_field[idx_jm1]);
            
            // Z?¹å?è³ªé??šé?å·?(?w)_{k+1} - (?w)_{k-1}
            double flux_z = (rho_d[idx_kp1] * w_field[idx_kp1]) - (rho_d[idx_km1] * w_field[idx_km1]);
            
            // è³ªé?å®ˆæ?èª¤å·® = |div(?u)| = |???v)/?‚y + ???w)/?‚z|
            double div_rho_u = std::fabs(flux_y + flux_z);
            total_error += div_rho_u;
        }
    }
    
    // ?å‚³ç¸½è³ª?é€šé?èª¤å·®ï¼ˆé™¤ä»¥æ ¼é»æ•¸å¾—åˆ°å¹³å??¼ï?
    int num_cells = (NY6 - 8) * (NZ6 - 8);
    return total_error / (double)num_cells;
}
```

---

### 4. Re=500 ç©©å??‹è?è¨˜é?

#### æ¨¡æ“¬?ƒæ•¸

| ?ƒæ•¸ | ??| èªªæ? |
|------|-----|------|
| Re | 500 | ?·è«¾??|
| tau | 0.6833 | é¬†å??‚é? |
| CFL | 0.8 | CFL ??|
| NY | 512 | Y ?¹å?ç¶²æ ¼??|
| NZ | 256 | Z ?¹å?ç¶²æ ¼??|
| NY6 | 519 | ??buffer ??Y ç¶²æ ¼??|
| NZ6 | 262 | ??buffer ??Z ç¶²æ ¼??|
| Ma_theoretical | 0.336666 | ?†è? Mach ?¸ï?è¶…é? Ma_max=0.3ï¼‰|

#### ?‹è?è¼¸å‡ºï¼ˆt=39000 ?‚ï?

```
Time=39000 ; Average Density=1 ; Density Correction=... ; Mass Flux Err=...
[Streaming Bounds] lower=25 upper=236 (target: 10/255)
[t=39000] Mach stats: max=0.3 (j=...,k=...), avg=0.17...
```

#### æµå ´?†æ?çµæ?

```
============================================================
Global Velocity Statistics (t=34000~39000)
============================================================
    Time     Uy_max     Uy_min    Uy_mean     Uz_max     Uz_min
------------------------------------------------------------
   34000    0.17320   -0.01115    0.10711    0.05660   -0.04676
   35000    0.17320   -0.01402    0.10701    0.05144   -0.04878
   36000    0.17320   -0.01578    0.10697    0.04710   -0.05055
   37000    0.17320   -0.01486    0.10728    0.04853   -0.05171
   38000    0.17320   -0.01069    0.10776    0.05344   -0.05275
   39000    0.17320   -0.01112    0.10809    0.05523   -0.05398

============================================================
Recirculation Zone Analysis
============================================================
k= 12: Backflow region Y=[0.65, 7.77], 406 points
k= 25: Backflow region Y=[0.70, 3.87], 181 points
k= 42: Backflow region Y=[1.21, 2.14], 54 points
k= 64: No backflow detected
```

#### ?œéµè§€å¯?
1. **Uy_max = 0.17320** è¢?Ma_max ?åˆ¶ä½?2. **?æ??€?ç¢ºå­˜åœ¨**ï¼šé?è¿‘å??¨æ?å¤§ç??å?æµï?ç¬¦å? Periodic Hill ?©ç??¹å¾µï¼?3. **å¯†åº¦ç©©å???~1.0**ï¼Œè³ª?å??†è‰¯å¥?
---

### 5. å¾?Edit6 ?å??ˆæœ¬?°ä?å¤©ç?æ¼”é€²æ­·ç¨?
#### Edit6 ?å??ˆæœ¬?¹é?
- ?ªé©?‰æ€§å…§?’æ??ä½¿?¨å…¨?Ÿè??¸è³¦??- `streaming_lower/upper` ?ºé???`#define` å¸¸æ•¸
- ?ºå??Šç?ï¼Œç„¡æ³•å??‹èª¿??
#### ä»Šæ—¥?¹é€?
| ?¹å??…ç›® | ?Šç? | ?°ç? |
|----------|------|------|
| streaming ?Šç? | `#define streaming_lower (25)` ?œæ? | `extern int streaming_lower` ?•æ??¨å?è®Šæ•¸ |
| ?Šç??´æ–° | ??| `UpdateStreamingBounds(t)` ?¨æ???tanh å¹³æ??æ¸¡ |
| è³ªé?èª¤å·®è¨ˆç? | å¯†åº¦?Œæ–¹æ³?| è³ªé??šé???º¦ div(?u) |
| è§??å±¤ç???| ?ºå? | æ¼¸é€²å??´å¤§ï¼ˆstreaming ??ä¸‰é? ??ä¸ƒé?ï¼‰|

---

### 6. ä¿®æ”¹?„æ?æ¡ˆå??´æ???
#### variables.h ä¿®æ”¹

```cpp
// === ?Ÿæœ¬ (?œæ?) ===
#define     streaming_lower      (25)
#define     streaming_upper      (NZ6-26)

// === ?¹ç‚º (?•æ?) ===
#define     interpolation_lower  (25)
#define     interpolation_upper  (NZ6-26)

//=== ?•æ? Streaming ?Šç??ƒæ•¸ï¼ˆæ¼¸?²å??´å¤§è§??å±¤ï?===//
#define     streaming_lower_init     (50)
#define     streaming_upper_init     (NZ6-51)
#define     streaming_lower_target   (10)
#define     streaming_upper_target   (NZ6-7)
#define     ramp_start_time          (0)
#define     ramp_end_time            (100000)

extern int streaming_lower;
extern int streaming_upper;
```

#### main.cpp ä¿®æ”¹

```cpp
// === ?°å??¨å?è®Šæ•¸?‡æ›´?°å‡½??===
int streaming_lower = streaming_lower_init;
int streaming_upper = streaming_upper_init;

void UpdateStreamingBounds(int t) {
    if (t >= ramp_end_time) {
        streaming_lower = streaming_lower_target;
        streaming_upper = streaming_upper_target;
    } else if (t <= ramp_start_time) {
        streaming_lower = streaming_lower_init;
        streaming_upper = streaming_upper_init;
    } else {
        double progress = (double)(t - ramp_start_time) / (ramp_end_time - ramp_start_time);
        double smooth_ratio = 0.5 * (1.0 + tanh(6.0 * (progress - 0.5)));
        streaming_lower = streaming_lower_init - 
            (int)(smooth_ratio * (streaming_lower_init - streaming_lower_target));
        streaming_upper = streaming_upper_init + 
            (int)(smooth_ratio * (streaming_upper_target - streaming_upper_init));
    }
}

// === ?‚é?è¿´å?ä¸­å‘¼??===
for(t = 0; t < loop; t++) {
    UpdateStreamingBounds(t);  // ?°å?
    // ...
}

// === è¼¸å‡º??§ ===
cout << "[Streaming Bounds] lower=" << streaming_lower 
     << " upper=" << streaming_upper 
     << " (target: " << streaming_lower_target << "/" << streaming_upper_target << ")" << endl;
```

#### evolution.h ä¿®æ”¹

```cpp
// === ?¹ç‚ºè³ªé??šé?èª¤å·® ===
double ComputeMassFluxError(double* rho_d, double* v_field, double* w_field) {
    double total_error = 0.0;
    for(int j = 4; j < NY6-4; j++) {
        for(int k = 4; k < NZ6-4; k++) {
            int idx = j * NZ6 + k;
            int idx_jp1 = (j+1) * NZ6 + k;
            int idx_jm1 = (j-1) * NZ6 + k;
            int idx_kp1 = j * NZ6 + (k+1);
            int idx_km1 = j * NZ6 + (k-1);
            double flux_y = (rho_d[idx_jp1] * v_field[idx_jp1]) - (rho_d[idx_jm1] * v_field[idx_jm1]);
            double flux_z = (rho_d[idx_kp1] * w_field[idx_kp1]) - (rho_d[idx_km1] * w_field[idx_km1]);
            double div_rho_u = std::fabs(flux_y + flux_z);
            total_error += div_rho_u;
        }
    }
    int num_cells = (NY6 - 8) * (NZ6 - 8);
    return total_error / (double)num_cells;
}

// === main.cpp ?¼å«?¹ç‚º ===
" ; Mass Flux Err=" << ComputeMassFluxError(rho, v, w) << endl;
```

---

### 7. ä¸‹ä?æ­¥å·¥ä½œæ–¹??
1. **æ¸¬è©¦æ¼¸é€²å??Šç??ˆæ?**
   - è§€å¯?streaming bounds ?¨æ??“è??–ç?è¼¸å‡º
   - ç¢ºè??¨é?æ¸¡æ??“æ¨¡?¬ç©©å®?
2. **?‘æˆ°?´é??·è«¾??*
   - ??Re=500 ç©©å?å¾Œï??—è©¦ Re=700, Re=1000
   - ?¯èƒ½?€è¦èª¿??ramp_end_time ?–å?å§‹é??Œå€?
3. **è©•ä¼°ä¸‰é??’å€¼ç·©è¡å??ˆæ?**
   - æ¯”è????¡ç·©è¡å??„ç©©å®šæ€§å·®??   - ?†æ?æµå ´ç²¾åº¦å½±éŸ¿

---

## 2026-01-28 ä¿®æ­£ï¼šå??æ®µ?æ¸¡è§?±º t=53000 å´©æ½°?é?

### ?é??è¿°

æ¨¡æ“¬??t??3000 ?‚é€??å´©æ½°?©æ¬¡?‚å??ç™¼?¾ï?

```
t=53000 ?‚ï?
- progress = 53000/100000 = 0.53
- smooth_ratio ??0.5 * (1 + tanh(6*(0.53-0.5))) ??0.59
- streaming_lower = 50 - 0.59*(50-10) ??26
```

æ­¤æ? `streaming_lower ??26` æ­?¥½?¥è? `interpolation_lower = 25`ï¼Œé€™æ˜¯ä¸€??*?¨ç??æ¸¡é»?*ï¼?
### ?é??¹å?

?Ÿè¨­è¨ˆä½¿?¨å–®?æ®µ?æ¸¡ï¼?- streaming 50 ??10ï¼ˆè·¨è¶Šä?ä¸‰é??’å€¼å??Œä?é»æ??¼å??„é??Œï?
- ??streaming_lower ?¥è? interpolation_lower ?‚ï?ç³»çµ±å¾?streaming ?´æ¥è·³åˆ°ä½¿ç”¨ä¸ƒé??’å€?- ?™å€‹ç?è®Šå¯?½å??´æ•¸?¼ä?ç©©å?

### è§?±º?¹æ?ï¼šå??æ®µ?æ¸¡

#### è¨­è??Ÿç?

```
=== ç¬¬ä??æ®µ (t=0 ~ 100000) ===
streaming 50 ??25ï¼ˆé??¾ä?é»æ??¼å?ï¼?   
t=0:      |====streaming (k??0)====|---ä¸ƒé??’å€?--|====streaming====|
t=100000: |==streaming (k??5)==|-------ä¸ƒé??’å€?------|==streaming==|
                               ??                      ??                          interpolation_lower    interpolation_upper

=== ç¬¬ä??æ®µ (t=100000 ~ 200000) ===
streaming 25 ??10ï¼ˆé??¾ä?é»æ??¼ç·©è¡å?ï¼?
t=100000: |==streaming (k??5)==|-------ä¸ƒé??’å€?------|==streaming==|
t=200000: |=s(k??0)=|--ä¸‰é?--|-------ä¸ƒé??’å€?------|--ä¸‰é?--|=s=|
              ??   ??                                ??     ??           target  interpolation_lower        interpolation_upper
```

#### ?ƒæ•¸ä¿®æ”¹ï¼ˆvariables.hï¼?
```cpp
//=== ?•æ? Streaming ?Šç??ƒæ•¸ï¼ˆå??æ®µæ¼¸é€²å??´å¤§è§??å±¤ï?===//
// ?å??¼ï?ä¿å?ï¼Œæ›´å¤§ç? streaming ?€?Ÿï?
#define     streaming_lower_init     (50)            // ?å?ä¸‹ç? (k <= 50 ??streaming)
#define     streaming_upper_init     (NZ6-51)        // ?å?ä¸Šç? (k >= NZ6-51 ??streaming)

// === ç¬¬ä??æ®µï¼šé??¾ä?é»æ??¼å? (streaming ??interpolation_lower) ===
#define     streaming_lower_phase1   (interpolation_lower)  // ç¬¬ä??æ®µ?®æ?: 25
#define     streaming_upper_phase1   (interpolation_upper)  // ç¬¬ä??æ®µ?®æ?: NZ6-26
#define     phase1_start_time        (0)             // ç¬¬ä??æ®µ?‹å?
#define     phase1_end_time          (100000)        // ç¬¬ä??æ®µçµæ?

// === ç¬¬ä??æ®µï¼šé??¾ä?é»æ??¼ç·©è¡å? (interpolation_lower ??target) ===
#define     streaming_lower_target   (10)            // ?€çµ‚ç›®æ¨™ä???#define     streaming_upper_target   (NZ6-7)         // ?€çµ‚ç›®æ¨™ä???#define     phase2_start_time        (100000)        // ç¬¬ä??æ®µ?‹å?
#define     phase2_end_time          (200000)        // ç¬¬ä??æ®µçµæ?
```

#### ?´æ–°?½æ•¸å¯¦ä?ï¼ˆmain.cppï¼?
```cpp
// ä½¿ç”¨ tanh å¹³æ??æ¸¡?´æ–° streaming ?Šç?ï¼ˆå??æ®µ?ˆæœ¬ï¼?// ç¬¬ä??æ®µ (t=0~100000)ï¼šstreaming 50??5ï¼ˆé??¾ä?é»æ??¼å?ï¼?// ç¬¬ä??æ®µ (t=100000~200000)ï¼šstreaming 25??0ï¼ˆé??¾ä?é»æ??¼ç·©è¡å?ï¼?void UpdateStreamingBounds(int t) {
    if (t >= phase2_end_time) {
        // ç¬¬ä??æ®µå®Œæ?ï¼Œä½¿?¨æ?çµ‚ç›®æ¨™å€?        streaming_lower = streaming_lower_target;
        streaming_upper = streaming_upper_target;
    } else if (t >= phase2_start_time) {
        // ç¬¬ä??æ®µï¼šå? interpolation ?Šç??æ¸¡?°æ?çµ‚ç›®æ¨™ï??‹æ”¾ä¸‰é??’å€¼å?ï¼?        double progress = (double)(t - phase2_start_time) / (phase2_end_time - phase2_start_time);
        double smooth_ratio = 0.5 * (1.0 + tanh(6.0 * (progress - 0.5)));
        
        streaming_lower = streaming_lower_phase1 - 
            (int)(smooth_ratio * (streaming_lower_phase1 - streaming_lower_target));
        streaming_upper = streaming_upper_phase1 + 
            (int)(smooth_ratio * (streaming_upper_target - streaming_upper_phase1));
    } else if (t >= phase1_start_time) {
        // ç¬¬ä??æ®µï¼šå??å??¼é?æ¸¡åˆ° interpolation ?Šç?ï¼ˆé??¾ä?é»æ??¼å?ï¼?        double progress = (double)(t - phase1_start_time) / (phase1_end_time - phase1_start_time);
        double smooth_ratio = 0.5 * (1.0 + tanh(6.0 * (progress - 0.5)));
        
        streaming_lower = streaming_lower_init - 
            (int)(smooth_ratio * (streaming_lower_init - streaming_lower_phase1));
        streaming_upper = streaming_upper_init + 
            (int)(smooth_ratio * (streaming_upper_phase1 - streaming_upper_init));
    } else {
        // å°šæœª?‹å?ï¼Œä½¿?¨å?å§‹å€?        streaming_lower = streaming_lower_init;
        streaming_upper = streaming_upper_init;
    }
}
```

### è¨­è??ªé?

1. **?¿å?è·¨è??¨ç?é»?*ï¼šç¬¬ä¸€?æ®µçµæ???interpolation_lowerï¼Œç¬¬äºŒé?æ®µæ?ç¹¼ç??‘ä?
2. **?…å?ç©©å??‚é?**ï¼šæ??æ®µ??100000 æ­¥ï?è®“ç³»çµ±å??†é©??3. **ä¿ç? tanh å¹³æ??æ¸¡**ï¼šé¿?é?æ¢¯è·³è®Šå¸¶ä¾†ç??¸å€¼è???4. **æ¼¸é€²å??‹æ”¾**ï¼šå??‹æ”¾ä¸ƒé??’å€¼å?ï¼ˆç²¾åº¦é?ï¼‰ï?ç©©å?å¾Œå??‹æ”¾ä¸‰é?ç·©è??€

### ?æ??‚é?ç·?
| ?‚é? | streaming_lower | streaming_upper | ?‹æ”¾?€??|
|------|-----------------|-----------------|----------|
| t=0 | 50 | 211 | streaming only |
| t=50000 | 30~35 | 225~230 | ä¸ƒé??’å€¼å??¨å??‹æ”¾ |
| t=100000 | 25 | 236 | ä¸ƒé??’å€¼å?å®Œå…¨?‹æ”¾ |
| t=150000 | 16~20 | 245~250 | ä¸‰é?ç·©è??€?¨å??‹æ”¾ |
| t=200000 | 10 | 255 | å®Œå…¨?‹æ”¾ï¼ˆæ?çµ‚ç??‹ï? |

---

### Git Commit è¨˜é?

```
commit e5f3c1b
Author: ChenPengChung
Date:   2026-01-28

feat: ?å¤§?¼ç¾ - ?•æ? Streaming ?Šç?æ¼¸é€²å??´å¤§è§??å±?
## ?¸å??¼ç¾ï¼šStreaming Layer ??Interpolation ?Šç??„ç¨ç«‹æ€?- interpolation_lower/upper: æ±ºå? XiPara[] ?å?ä»€éº¼é??‹ç??’å€¼æ??ï??å??–æ?ï¼?- streaming_lower/upper: æ±ºå??¯å¦è·³é??’å€¼æ”¹??streamingï¼ˆæ?æ­¥å³?‚åˆ¤?·ï??¯å??‹èª¿?´ï?
- ?è?ä¿®æ­£ï¼šstreaming_lower ?¯ä»¥ < interpolation_lowerï¼Œè?ä¸‰é??’å€¼å?ä½œç‚ºç©©å?ç·©è?

## æ¼¸é€²å?è¨­è?
- ä½¿ç”¨ tanh å¹³æ??æ¸¡?½æ•¸?¿å??æ¢¯è·³è?
- ?å??? streaming_lower=50, streaming_upper=NZ6-51
- ?®æ??? streaming_lower=10, streaming_upper=NZ6-7
- ?æ¸¡?‚é?: t=0 ~ t=100000

## ?¶ä??¹é€?- è³ªé??šé?èª¤å·®è¨ˆç??¹ç”¨ div(?u) ?¬å?

3 files changed, 81 insertions(+), 18 deletions(-)
```




?®å??„å?é¡Œæ˜¯ï¼Œç•¶?‘CFLä¸‹é?ï¼Œæœ¬ä¾†é??Œä?å¢œé??ªï??½æ??ºå?é¡Œï??¸è??¼ï¼£ï¼¦ï¼¬ï¼?.8 ?¦å?ï¼Œç•¶?‚é?æ­¥æ¨?²æ?ï¼ŒMa?ƒæ?é«˜æ?å¾Œå¹³?‡è???.3?‚ï??´å€‹ç??œé?å§‹å‡º?¾é??ªï?? æ­¤?‹å?çµæ?ä¸å¯ä¿¡ï??‰ä?éº¼æ–¹æ³•å¯ä»¥æ??¶é¦¬èµ«æ•¸?›è?ï¼Œé›£?“ç›®?é??Ÿæ°£?¯å”¯ä¸€è§?–¹ï¼Œä??¯é??Ÿå™¨?ƒä?æ­?¢º?°å½±?³å??›ä¿®æ­??ä¹Ÿä¸¦?å¥½?¹æ?ï¼Œæœ¨?åª?šä?å®¶æººå¾€??56->512 ä½†æ˜¯?¹æœ¬é¦¬èµ«?¸é?é«˜å??´ç??œä??¯ä¿¡?„å?é¡Œä??¶ç„¡æ³•è§£æ±ºThiukHard (ä¸è??¹code?ˆè?è«?


```
