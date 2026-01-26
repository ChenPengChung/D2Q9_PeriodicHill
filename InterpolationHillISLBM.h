#ifndef INTERPOLATIONHILLISLBM_FILE
#define INTERPOLATIONHILLISLBM_FILE
//定義二階插值公式
#define Intrpl3(f1, a1, f2, a2, f3, a3) \
(       \
    f1*a1 + f2*a2 + f3*a3   \
)
//定義六階插值公式
#define Intrpl7(f1, a1, f2, a2, f3, a3, f4, a4, f5, a5, f6, a6, f7, a7)     \
(       \
    f1*a1 + f2*a2 + f3*a3 + f4*a4 + f5*a5 + f6*a6 + f7*a7   \
)
//定義F0_in的計算公式
//在此，F0_in為物理空間計算點的碰撞前插值後一般態分佈函數
#define F0_Intrpl7(f,j,k) \
    idx = j * nface + k ; \
    F0_in = f[idx] ; /*輸入指標變數使用指標變數相對應之一維連續陣列*/

//如下定義對應D2Q9模型的F1 (+Y方向)
// f: 分佈函數陣列
// idx_f: f 的基底索引 = j_c * nface + k_c (stencil 起點)
// idx_xi: Xi 權重的基底索引 = j * nface + k (目標點)
// Xi權重存取: xi_w[idx_xi + row*NY6*nface], w=0~6 權重編號, row=0~6 y方向線編號
//2026.01.25 將 j_c , k_c 改成 cell_z內插成起始陣列 
// 注意: cell_z 儲存的是 k_start (z方向起始索引)，需要加上 y 列偏移
// N=0 對應 j-1 列, N=1 對應 j 列, N=2 對應 j+1 列
#define F1_Intrpl3(f,j,k,cell_z,idx_y,idx_xi,y_0,y_1,y_2,xi_0,xi_1,xi_2,xi_3,xi_4,xi_5,xi_6) \
do { /*降階版本*/ /*先插 Xi 方向，再插 Y 方向*/\
    const int jm1 = (j == 3) ? (NY6 - 4) : (j - 1); /*處理週期邊界*/ \
    const int jp1 = (j == NY6 - 4) ? 3 : (j + 1); \
     \
    const int base0 = jm1*nface + cell_z[idx_xi+0*NY6*NZ6]; \
    const int base1 = j*nface + cell_z[idx_xi+1*NY6*NZ6]; \
    const int base2 = jp1*nface + cell_z[idx_xi+2*NY6*NZ6]; \
    F1_in = Intrpl3( \
        Intrpl7( f[base0],xi_0[idx_xi+0*NY6*nface],f[base0+1],xi_1[idx_xi+0*NY6*nface],f[base0+2],xi_2[idx_xi+0*NY6*nface],f[base0+3],xi_3[idx_xi+0*NY6*nface],f[base0+4],xi_4[idx_xi+0*NY6*nface],f[base0+5],xi_5[idx_xi+0*NY6*nface],f[base0+6],xi_6[idx_xi+0*NY6*nface] ), y_0[idx_y], \
        Intrpl7( f[base1],xi_0[idx_xi+1*NY6*nface],f[base1+1],xi_1[idx_xi+1*NY6*nface],f[base1+2],xi_2[idx_xi+1*NY6*nface],f[base1+3],xi_3[idx_xi+1*NY6*nface],f[base1+4],xi_4[idx_xi+1*NY6*nface],f[base1+5],xi_5[idx_xi+1*NY6*nface],f[base1+6],xi_6[idx_xi+1*NY6*nface] ), y_1[idx_y], \
        Intrpl7( f[base2],xi_0[idx_xi+2*NY6*nface],f[base2+1],xi_1[idx_xi+2*NY6*nface],f[base2+2],xi_2[idx_xi+2*NY6*nface],f[base2+3],xi_3[idx_xi+2*NY6*nface],f[base2+4],xi_4[idx_xi+2*NY6*nface],f[base2+5],xi_5[idx_xi+2*NY6*nface],f[base2+6],xi_6[idx_xi+2*NY6*nface] ), y_2[idx_y] \
    ); \
} while(0)

#define F1_Intrpl7(f,j,k,j_c,k_c,idx_y,idx_xi,y_0,y_1,y_2,y_3,y_4,y_5,y_6,xi_0,xi_1,xi_2,xi_3,xi_4,xi_5,xi_6) \
do { \
    const int idx_f = j_c * nface + k_c; \
    F1_in = Intrpl7( \
        Intrpl7( f[idx_f],             xi_0[idx_xi+0*NY6*nface], f[idx_f+1],             xi_1[idx_xi+0*NY6*nface], f[idx_f+2],             xi_2[idx_xi+0*NY6*nface], f[idx_f+3],             xi_3[idx_xi+0*NY6*nface], f[idx_f+4],             xi_4[idx_xi+0*NY6*nface], f[idx_f+5],             xi_5[idx_xi+0*NY6*nface], f[idx_f+6],             xi_6[idx_xi+0*NY6*nface] ), y_0[idx_y], \
        Intrpl7( f[idx_f+1*nface],     xi_0[idx_xi+1*NY6*nface], f[idx_f+1*nface+1],     xi_1[idx_xi+1*NY6*nface], f[idx_f+1*nface+2],     xi_2[idx_xi+1*NY6*nface], f[idx_f+1*nface+3],     xi_3[idx_xi+1*NY6*nface], f[idx_f+1*nface+4],     xi_4[idx_xi+1*NY6*nface], f[idx_f+1*nface+5],     xi_5[idx_xi+1*NY6*nface], f[idx_f+1*nface+6],     xi_6[idx_xi+1*NY6*nface] ), y_1[idx_y], \
        Intrpl7( f[idx_f+2*nface],     xi_0[idx_xi+2*NY6*nface], f[idx_f+2*nface+1],     xi_1[idx_xi+2*NY6*nface], f[idx_f+2*nface+2],     xi_2[idx_xi+2*NY6*nface], f[idx_f+2*nface+3],     xi_3[idx_xi+2*NY6*nface], f[idx_f+2*nface+4],     xi_4[idx_xi+2*NY6*nface], f[idx_f+2*nface+5],     xi_5[idx_xi+2*NY6*nface], f[idx_f+2*nface+6],     xi_6[idx_xi+2*NY6*nface] ), y_2[idx_y], \
        Intrpl7( f[idx_f+3*nface],     xi_0[idx_xi+3*NY6*nface], f[idx_f+3*nface+1],     xi_1[idx_xi+3*NY6*nface], f[idx_f+3*nface+2],     xi_2[idx_xi+3*NY6*nface], f[idx_f+3*nface+3],     xi_3[idx_xi+3*NY6*nface], f[idx_f+3*nface+4],     xi_4[idx_xi+3*NY6*nface], f[idx_f+3*nface+5],     xi_5[idx_xi+3*NY6*nface], f[idx_f+3*nface+6],     xi_6[idx_xi+3*NY6*nface] ), y_3[idx_y], \
        Intrpl7( f[idx_f+4*nface],     xi_0[idx_xi+4*NY6*nface], f[idx_f+4*nface+1],     xi_1[idx_xi+4*NY6*nface], f[idx_f+4*nface+2],     xi_2[idx_xi+4*NY6*nface], f[idx_f+4*nface+3],     xi_3[idx_xi+4*NY6*nface], f[idx_f+4*nface+4],     xi_4[idx_xi+4*NY6*nface], f[idx_f+4*nface+5],     xi_5[idx_xi+4*NY6*nface], f[idx_f+4*nface+6],     xi_6[idx_xi+4*NY6*nface] ), y_4[idx_y], \
        Intrpl7( f[idx_f+5*nface],     xi_0[idx_xi+5*NY6*nface], f[idx_f+5*nface+1],     xi_1[idx_xi+5*NY6*nface], f[idx_f+5*nface+2],     xi_2[idx_xi+5*NY6*nface], f[idx_f+5*nface+3],     xi_3[idx_xi+5*NY6*nface], f[idx_f+5*nface+4],     xi_4[idx_xi+5*NY6*nface], f[idx_f+5*nface+5],     xi_5[idx_xi+5*NY6*nface], f[idx_f+5*nface+6],     xi_6[idx_xi+5*NY6*nface] ), y_5[idx_y], \
        Intrpl7( f[idx_f+6*nface],     xi_0[idx_xi+6*NY6*nface], f[idx_f+6*nface+1],     xi_1[idx_xi+6*NY6*nface], f[idx_f+6*nface+2],     xi_2[idx_xi+6*NY6*nface], f[idx_f+6*nface+3],     xi_3[idx_xi+6*NY6*nface], f[idx_f+6*nface+4],     xi_4[idx_xi+6*NY6*nface], f[idx_f+6*nface+5],     xi_5[idx_xi+6*NY6*nface], f[idx_f+6*nface+6],     xi_6[idx_xi+6*NY6*nface] ), y_6[idx_y]  \
    ); \
} while(0)

//如下定義對應D2Q9模型的F3 (-Y方向)
// 注意: cell_z 儲存的是 k_start，需要加上 y 列偏移
// F3 來自 j+1 方向，所以 N=0 對應 j 列, N=1 對應 j+1 列, N=2 對應 j+2 列
#define F3_Intrpl3(f,j,k,cell_z,idx_y,idx_xi,y_0,y_1,y_2,xi_0,xi_1,xi_2,xi_3,xi_4,xi_5,xi_6) \
do { /*降階版本*/\
    const int jm1 = (j == 3) ? (NY6 - 4) : (j - 1); /*處理週期邊界*/ \
    const int jp1 = (j == NY6 - 4) ? 3 : (j + 1); \
     \
    const int base0 = jm1*nface + cell_z[idx_xi+0*NY6*NZ6]; \
    const int base1 = j*nface + cell_z[idx_xi+1*NY6*NZ6]; \
    const int base2 = jp1*nface + cell_z[idx_xi+2*NY6*NZ6]; \
    F3_in = Intrpl3( \
        Intrpl7( f[base0],xi_0[idx_xi+0*NY6*nface],f[base0+1],xi_1[idx_xi+0*NY6*nface],f[base0+2],xi_2[idx_xi+0*NY6*nface],f[base0+3],xi_3[idx_xi+0*NY6*nface],f[base0+4],xi_4[idx_xi+0*NY6*nface],f[base0+5],xi_5[idx_xi+0*NY6*nface],f[base0+6],xi_6[idx_xi+0*NY6*nface] ), y_0[idx_y], \
        Intrpl7( f[base1],xi_0[idx_xi+1*NY6*nface],f[base1+1],xi_1[idx_xi+1*NY6*nface],f[base1+2],xi_2[idx_xi+1*NY6*nface],f[base1+3],xi_3[idx_xi+1*NY6*nface],f[base1+4],xi_4[idx_xi+1*NY6*nface],f[base1+5],xi_5[idx_xi+1*NY6*nface],f[base1+6],xi_6[idx_xi+1*NY6*nface] ), y_1[idx_y], \
        Intrpl7( f[base2],xi_0[idx_xi+2*NY6*nface],f[base2+1],xi_1[idx_xi+2*NY6*nface],f[base2+2],xi_2[idx_xi+2*NY6*nface],f[base2+3],xi_3[idx_xi+2*NY6*nface],f[base2+4],xi_4[idx_xi+2*NY6*nface],f[base2+5],xi_5[idx_xi+2*NY6*nface],f[base2+6],xi_6[idx_xi+2*NY6*nface] ), y_2[idx_y] \
    ); \
} while(0)
#define F3_Intrpl7(f,j,k,j_c,k_c,idx_y,idx_xi,y_0,y_1,y_2,y_3,y_4,y_5,y_6,xi_0,xi_1,xi_2,xi_3,xi_4,xi_5,xi_6) \
do { \
    const int idx_f = j_c * nface + k_c; \
    F3_in = Intrpl7( \
        Intrpl7( f[idx_f],             xi_0[idx_xi+0*NY6*nface], f[idx_f+1],             xi_1[idx_xi+0*NY6*nface], f[idx_f+2],             xi_2[idx_xi+0*NY6*nface], f[idx_f+3],             xi_3[idx_xi+0*NY6*nface], f[idx_f+4],             xi_4[idx_xi+0*NY6*nface], f[idx_f+5],             xi_5[idx_xi+0*NY6*nface], f[idx_f+6],             xi_6[idx_xi+0*NY6*nface] ), y_0[idx_y], \
        Intrpl7( f[idx_f+1*nface],     xi_0[idx_xi+1*NY6*nface], f[idx_f+1*nface+1],     xi_1[idx_xi+1*NY6*nface], f[idx_f+1*nface+2],     xi_2[idx_xi+1*NY6*nface], f[idx_f+1*nface+3],     xi_3[idx_xi+1*NY6*nface], f[idx_f+1*nface+4],     xi_4[idx_xi+1*NY6*nface], f[idx_f+1*nface+5],     xi_5[idx_xi+1*NY6*nface], f[idx_f+1*nface+6],     xi_6[idx_xi+1*NY6*nface] ), y_1[idx_y], \
        Intrpl7( f[idx_f+2*nface],     xi_0[idx_xi+2*NY6*nface], f[idx_f+2*nface+1],     xi_1[idx_xi+2*NY6*nface], f[idx_f+2*nface+2],     xi_2[idx_xi+2*NY6*nface], f[idx_f+2*nface+3],     xi_3[idx_xi+2*NY6*nface], f[idx_f+2*nface+4],     xi_4[idx_xi+2*NY6*nface], f[idx_f+2*nface+5],     xi_5[idx_xi+2*NY6*nface], f[idx_f+2*nface+6],     xi_6[idx_xi+2*NY6*nface] ), y_2[idx_y], \
        Intrpl7( f[idx_f+3*nface],     xi_0[idx_xi+3*NY6*nface], f[idx_f+3*nface+1],     xi_1[idx_xi+3*NY6*nface], f[idx_f+3*nface+2],     xi_2[idx_xi+3*NY6*nface], f[idx_f+3*nface+3],     xi_3[idx_xi+3*NY6*nface], f[idx_f+3*nface+4],     xi_4[idx_xi+3*NY6*nface], f[idx_f+3*nface+5],     xi_5[idx_xi+3*NY6*nface], f[idx_f+3*nface+6],     xi_6[idx_xi+3*NY6*nface] ), y_3[idx_y], \
        Intrpl7( f[idx_f+4*nface],     xi_0[idx_xi+4*NY6*nface], f[idx_f+4*nface+1],     xi_1[idx_xi+4*NY6*nface], f[idx_f+4*nface+2],     xi_2[idx_xi+4*NY6*nface], f[idx_f+4*nface+3],     xi_3[idx_xi+4*NY6*nface], f[idx_f+4*nface+4],     xi_4[idx_xi+4*NY6*nface], f[idx_f+4*nface+5],     xi_5[idx_xi+4*NY6*nface], f[idx_f+4*nface+6],     xi_6[idx_xi+4*NY6*nface] ), y_4[idx_y], \
        Intrpl7( f[idx_f+5*nface],     xi_0[idx_xi+5*NY6*nface], f[idx_f+5*nface+1],     xi_1[idx_xi+5*NY6*nface], f[idx_f+5*nface+2],     xi_2[idx_xi+5*NY6*nface], f[idx_f+5*nface+3],     xi_3[idx_xi+5*NY6*nface], f[idx_f+5*nface+4],     xi_4[idx_xi+5*NY6*nface], f[idx_f+5*nface+5],     xi_5[idx_xi+5*NY6*nface], f[idx_f+5*nface+6],     xi_6[idx_xi+5*NY6*nface] ), y_5[idx_y], \
        Intrpl7( f[idx_f+6*nface],     xi_0[idx_xi+6*NY6*nface], f[idx_f+6*nface+1],     xi_1[idx_xi+6*NY6*nface], f[idx_f+6*nface+2],     xi_2[idx_xi+6*NY6*nface], f[idx_f+6*nface+3],     xi_3[idx_xi+6*NY6*nface], f[idx_f+6*nface+4],     xi_4[idx_xi+6*NY6*nface], f[idx_f+6*nface+5],     xi_5[idx_xi+6*NY6*nface], f[idx_f+6*nface+6],     xi_6[idx_xi+6*NY6*nface] ), y_6[idx_y]  \
    ); \
} while(0)

//如下定義為D2Q9模型的F2 (+Z方向) - 只需要 Z 方向插值
// F2 來自 j 列（y 不變），所以只用 N=1 對應的 cell_z
#define F2_Intrpl7(f, j, k, cell_z , idx_y, idx_xi, xi_0, xi_1, xi_2, xi_3, xi_4, xi_5, xi_6)    \
do { \
    const int base1 = j*nface + cell_z[idx_xi+1*NY6*NZ6]; \
    F2_in = Intrpl7( f[base1],xi_0[idx_xi+1*NY6*nface],f[base1+1],xi_1[idx_xi+1*NY6*nface],f[base1+2],xi_2[idx_xi+1*NY6*nface],f[base1+3],xi_3[idx_xi+1*NY6*nface],f[base1+4],xi_4[idx_xi+1*NY6*nface],f[base1+5],xi_5[idx_xi+1*NY6*nface],f[base1+6],xi_6[idx_xi+1*NY6*nface] ); \
} while(0)


//如下定義為D2Q9模型的F4 (-Z方向) - 只需要 Z 方向插值
// F4 來自 j 列（y 不變），所以只用 N=1 對應的 cell_z
#define F4_Intrpl7(f, j, k,  cell_z , idx_y, idx_xi, xi_0, xi_1, xi_2, xi_3, xi_4, xi_5, xi_6)    \
do { \
    const int base1 = j*nface + cell_z[idx_xi+1*NY6*NZ6]; \
    F4_in = Intrpl7( f[base1],xi_0[idx_xi+1*NY6*nface],f[base1+1],xi_1[idx_xi+1*NY6*nface],f[base1+2],xi_2[idx_xi+1*NY6*nface],f[base1+3],xi_3[idx_xi+1*NY6*nface],f[base1+4],xi_4[idx_xi+1*NY6*nface],f[base1+5],xi_5[idx_xi+1*NY6*nface],f[base1+6],xi_6[idx_xi+1*NY6*nface] ); \
} while(0)



//註解:權重根據起點判斷選擇性為 0
//如下定義為D2Q9模型的F5,F6,F7,F8 (斜向) - 需要 Y 和 Z 
// 注意: cell_z 儲存的是 k_start，需要加上 y 列偏移
#define Y_XI_Intrpl3(f, F_in, j, k, cell_z , idx_y, idx_xi, y_0, y_1, y_2, xi_0, xi_1, xi_2, xi_3, xi_4, xi_5, xi_6)    \
do { /*降階版本*/\
    const int jm1 = (j == 3) ? (NY6 - 4) : (j - 1); /*處理週期邊界*/ \
    const int jp1 = (j == NY6 - 4) ? 3 : (j + 1); \
     \
    const int base0 = jm1*nface + cell_z[idx_xi+0*NY6*NZ6]; /*在哪一個y層+起點座標*/\
    const int base1 = j*nface + cell_z[idx_xi+1*NY6*NZ6]; \
    const int base2 = jp1*nface + cell_z[idx_xi+2*NY6*NZ6]; \
    F_in = Intrpl3( \
        Intrpl7( f[base0],xi_0[idx_xi+0*NY6*nface],f[base0+1],xi_1[idx_xi+0*NY6*nface],f[base0+2],xi_2[idx_xi+0*NY6*nface],f[base0+3],xi_3[idx_xi+0*NY6*nface],f[base0+4],xi_4[idx_xi+0*NY6*nface],f[base0+5],xi_5[idx_xi+0*NY6*nface],f[base0+6],xi_6[idx_xi+0*NY6*nface] ), y_0[idx_y], \
        Intrpl7( f[base1],xi_0[idx_xi+1*NY6*nface],f[base1+1],xi_1[idx_xi+1*NY6*nface],f[base1+2],xi_2[idx_xi+1*NY6*nface],f[base1+3],xi_3[idx_xi+1*NY6*nface],f[base1+4],xi_4[idx_xi+1*NY6*nface],f[base1+5],xi_5[idx_xi+1*NY6*nface],f[base1+6],xi_6[idx_xi+1*NY6*nface] ), y_1[idx_y], \
        Intrpl7( f[base2],xi_0[idx_xi+2*NY6*nface],f[base2+1],xi_1[idx_xi+2*NY6*nface],f[base2+2],xi_2[idx_xi+2*NY6*nface],f[base2+3],xi_3[idx_xi+2*NY6*nface],f[base2+4],xi_4[idx_xi+2*NY6*nface],f[base2+5],xi_5[idx_xi+2*NY6*nface],f[base2+6],xi_6[idx_xi+2*NY6*nface] ), y_2[idx_y] \
    ); \
} while(0)

#define Y_XI_Intrpl7(f, F_in, j, k, j_c, k_c, idx_y, idx_xi, y_0, y_1, y_2, y_3, y_4, y_5, y_6, xi_0, xi_1, xi_2, xi_3, xi_4, xi_5, xi_6)    \
do { \
    const int idx_f = j_c * nface + k_c; \
    F_in = Intrpl7( \
        Intrpl7( f[idx_f],             xi_0[idx_xi+0*NY6*nface], f[idx_f+1],             xi_1[idx_xi+0*NY6*nface], f[idx_f+2],             xi_2[idx_xi+0*NY6*nface], f[idx_f+3],             xi_3[idx_xi+0*NY6*nface], f[idx_f+4],             xi_4[idx_xi+0*NY6*nface], f[idx_f+5],             xi_5[idx_xi+0*NY6*nface], f[idx_f+6],             xi_6[idx_xi+0*NY6*nface] ), y_0[idx_y], \
        Intrpl7( f[idx_f+1*nface],     xi_0[idx_xi+1*NY6*nface], f[idx_f+1*nface+1],     xi_1[idx_xi+1*NY6*nface], f[idx_f+1*nface+2],     xi_2[idx_xi+1*NY6*nface], f[idx_f+1*nface+3],     xi_3[idx_xi+1*NY6*nface], f[idx_f+1*nface+4],     xi_4[idx_xi+1*NY6*nface], f[idx_f+1*nface+5],     xi_5[idx_xi+1*NY6*nface], f[idx_f+1*nface+6],     xi_6[idx_xi+1*NY6*nface] ), y_1[idx_y], \
        Intrpl7( f[idx_f+2*nface],     xi_0[idx_xi+2*NY6*nface], f[idx_f+2*nface+1],     xi_1[idx_xi+2*NY6*nface], f[idx_f+2*nface+2],     xi_2[idx_xi+2*NY6*nface], f[idx_f+2*nface+3],     xi_3[idx_xi+2*NY6*nface], f[idx_f+2*nface+4],     xi_4[idx_xi+2*NY6*nface], f[idx_f+2*nface+5],     xi_5[idx_xi+2*NY6*nface], f[idx_f+2*nface+6],     xi_6[idx_xi+2*NY6*nface] ), y_2[idx_y], \
        Intrpl7( f[idx_f+3*nface],     xi_0[idx_xi+3*NY6*nface], f[idx_f+3*nface+1],     xi_1[idx_xi+3*NY6*nface], f[idx_f+3*nface+2],     xi_2[idx_xi+3*NY6*nface], f[idx_f+3*nface+3],     xi_3[idx_xi+3*NY6*nface], f[idx_f+3*nface+4],     xi_4[idx_xi+3*NY6*nface], f[idx_f+3*nface+5],     xi_5[idx_xi+3*NY6*nface], f[idx_f+3*nface+6],     xi_6[idx_xi+3*NY6*nface] ), y_3[idx_y], \
        Intrpl7( f[idx_f+4*nface],     xi_0[idx_xi+4*NY6*nface], f[idx_f+4*nface+1],     xi_1[idx_xi+4*NY6*nface], f[idx_f+4*nface+2],     xi_2[idx_xi+4*NY6*nface], f[idx_f+4*nface+3],     xi_3[idx_xi+4*NY6*nface], f[idx_f+4*nface+4],     xi_4[idx_xi+4*NY6*nface], f[idx_f+4*nface+5],     xi_5[idx_xi+4*NY6*nface], f[idx_f+4*nface+6],     xi_6[idx_xi+4*NY6*nface] ), y_4[idx_y], \
        Intrpl7( f[idx_f+5*nface],     xi_0[idx_xi+5*NY6*nface], f[idx_f+5*nface+1],     xi_1[idx_xi+5*NY6*nface], f[idx_f+5*nface+2],     xi_2[idx_xi+5*NY6*nface], f[idx_f+5*nface+3],     xi_3[idx_xi+5*NY6*nface], f[idx_f+5*nface+4],     xi_4[idx_xi+5*NY6*nface], f[idx_f+5*nface+5],     xi_5[idx_xi+5*NY6*nface], f[idx_f+5*nface+6],     xi_6[idx_xi+5*NY6*nface] ), y_5[idx_y], \
        Intrpl7( f[idx_f+6*nface],     xi_0[idx_xi+6*NY6*nface], f[idx_f+6*nface+1],     xi_1[idx_xi+6*NY6*nface], f[idx_f+6*nface+2],     xi_2[idx_xi+6*NY6*nface], f[idx_f+6*nface+3],     xi_3[idx_xi+6*NY6*nface], f[idx_f+6*nface+4],     xi_4[idx_xi+6*NY6*nface], f[idx_f+6*nface+5],     xi_5[idx_xi+6*NY6*nface], f[idx_f+6*nface+6],     xi_6[idx_xi+6*NY6*nface] ), y_6[idx_y]  \
    ); \
} while(0)



#endif
