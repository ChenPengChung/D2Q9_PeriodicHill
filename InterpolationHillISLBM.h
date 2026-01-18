#ifndef INTERPOLATIONHILLISLBM_FILE
#define INTERPOLATIONHILLISLBM_FILE
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
//如下定義對應D3Q19模型的F3 
#define F1_Intrpl7(f,j,k,j_c,k_c,idx_y,idx_xi,y_0,y_1,y_2,y_3,y_4,y_5,y_6,xi_0,xi_1,xi_2,xi_3,xi_4,xi_5,xi_6)\
    idx = j_c * nface + k_c ;/*權重參數獨立坐標系先z方向內插再往y方向內插*/\
    F1_in = Intrpl7(    \
        Intrpl7( f[idx],         xi_0[idx_xi], f[idx+1],         xi_1[idx_xi], f[idx+2],         xi_2[idx_xi], f[idx+3],         xi_3[idx_xi], f[idx+4],         xi_4[idx_xi], f[idx+5],         xi_5[idx_xi], f[idx+6],         xi_6[idx_xi] ), y_0[idx_y],   \
        Intrpl7( f[idx+nface],   xi_0[idx_xi], f[idx+nface+1],   xi_1[idx_xi], f[idx+2+nface],   xi_2[idx_xi], f[idx+3+nface],   xi_3[idx_xi], f[idx+4+nface],   xi_4[idx_xi], f[idx+5+nface],   xi_5[idx_xi], f[idx+6+nface],   xi_6[idx_xi] ), y_1[idx_y],   \
        Intrpl7( f[idx+2*nface], xi_0[idx_xi], f[idx+2*nface+1], xi_1[idx_xi], f[idx+2+2*nface], xi_2[idx_xi], f[idx+3+2*nface], xi_3[idx_xi], f[idx+4+2*nface], xi_4[idx_xi], f[idx+5+2*nface], xi_5[idx_xi], f[idx+6+2*nface], xi_6[idx_xi] ), y_2[idx_y],   \
        Intrpl7( f[idx+3*nface], xi_0[idx_xi], f[idx+3*nface+1], xi_1[idx_xi], f[idx+2+3*nface], xi_2[idx_xi], f[idx+3+3*nface], xi_3[idx_xi], f[idx+4+3*nface], xi_4[idx_xi], f[idx+5+3*nface], xi_5[idx_xi], f[idx+6+3*nface], xi_6[idx_xi] ), y_3[idx_y],   \
        Intrpl7( f[idx+4*nface], xi_0[idx_xi], f[idx+1+4*nface], xi_1[idx_xi], f[idx+2+4*nface], xi_2[idx_xi], f[idx+3+4*nface], xi_3[idx_xi], f[idx+4+4*nface], xi_4[idx_xi], f[idx+5+4*nface], xi_5[idx_xi], f[idx+6+4*nface], xi_6[idx_xi] ), y_4[idx_y],   \
        Intrpl7( f[idx+5*nface], xi_0[idx_xi], f[idx+1+5*nface], xi_1[idx_xi], f[idx+2+5*nface], xi_2[idx_xi], f[idx+3+5*nface], xi_3[idx_xi], f[idx+4+5*nface], xi_4[idx_xi], f[idx+5+5*nface], xi_5[idx_xi], f[idx+6+5*nface], xi_6[idx_xi] ), y_5[idx_y],   \
        Intrpl7( f[idx+6*nface], xi_0[idx_xi], f[idx+1+6*nface], xi_1[idx_xi], f[idx+2+6*nface], xi_2[idx_xi], f[idx+3+6*nface], xi_3[idx_xi], f[idx+4+6*nface], xi_4[idx_xi], f[idx+5+6*nface], xi_5[idx_xi], f[idx+6+6*nface], xi_6[idx_xi] ), y_6[idx_y]    \
    );
//如下定義對應D3Q19模型的F4
#define F3_Intrpl7(f,j,k,j_c,k_c,idx_y,idx_xi,y_0,y_1,y_2,y_3,y_4,y_5,y_6,xi_0,xi_1,xi_2,xi_3,xi_4,xi_5,xi_6)\
    idx = j_c * nface + k_c ;/*權重參數獨立坐標系先z方向內插再往y方向內插*/\
    F3_in = Intrpl7(    \
        Intrpl7( f[idx],         xi_0[idx_xi], f[idx+1],         xi_1[idx_xi], f[idx+2],         xi_2[idx_xi], f[idx+3],         xi_3[idx_xi], f[idx+4],         xi_4[idx_xi], f[idx+5],         xi_5[idx_xi], f[idx+6],         xi_6[idx_xi] ), y_0[idx_y],   \
        Intrpl7( f[idx+nface],   xi_0[idx_xi], f[idx+nface+1],   xi_1[idx_xi], f[idx+2+nface],   xi_2[idx_xi], f[idx+3+nface],   xi_3[idx_xi], f[idx+4+nface],   xi_4[idx_xi], f[idx+5+nface],   xi_5[idx_xi], f[idx+6+nface],   xi_6[idx_xi] ), y_1[idx_y],   \
        Intrpl7( f[idx+2*nface], xi_0[idx_xi], f[idx+2*nface+1], xi_1[idx_xi], f[idx+2+2*nface], xi_2[idx_xi], f[idx+3+2*nface], xi_3[idx_xi], f[idx+4+2*nface], xi_4[idx_xi], f[idx+5+2*nface], xi_5[idx_xi], f[idx+6+2*nface], xi_6[idx_xi] ), y_2[idx_y],   \
        Intrpl7( f[idx+3*nface], xi_0[idx_xi], f[idx+3*nface+1], xi_1[idx_xi], f[idx+2+3*nface], xi_2[idx_xi], f[idx+3+3*nface], xi_3[idx_xi], f[idx+4+3*nface], xi_4[idx_xi], f[idx+5+3*nface], xi_5[idx_xi], f[idx+6+3*nface], xi_6[idx_xi] ), y_3[idx_y],   \
        Intrpl7( f[idx+4*nface], xi_0[idx_xi], f[idx+1+4*nface], xi_1[idx_xi], f[idx+2+4*nface], xi_2[idx_xi], f[idx+3+4*nface], xi_3[idx_xi], f[idx+4+4*nface], xi_4[idx_xi], f[idx+5+4*nface], xi_5[idx_xi], f[idx+6+4*nface], xi_6[idx_xi] ), y_4[idx_y],   \
        Intrpl7( f[idx+5*nface], xi_0[idx_xi], f[idx+1+5*nface], xi_1[idx_xi], f[idx+2+5*nface], xi_2[idx_xi], f[idx+3+5*nface], xi_3[idx_xi], f[idx+4+5*nface], xi_4[idx_xi], f[idx+5+5*nface], xi_5[idx_xi], f[idx+6+5*nface], xi_6[idx_xi] ), y_5[idx_y],   \
        Intrpl7( f[idx+6*nface], xi_0[idx_xi], f[idx+1+6*nface], xi_1[idx_xi], f[idx+2+6*nface], xi_2[idx_xi], f[idx+3+6*nface], xi_3[idx_xi], f[idx+4+6*nface], xi_4[idx_xi], f[idx+5+6*nface], xi_5[idx_xi], f[idx+6+6*nface], xi_6[idx_xi] ), y_6[idx_y]    \
    );
//如下定義為對應D3Q19模型的F5
#define F2_Intrpl7(f, j, k, j_c, k_c, idx_y, idx_xi, xi_0, xi_1, xi_2, xi_3, xi_4, xi_5, xi_6)    \
    idx = j*nface + k_c ;    \
    F2_in = Intrpl7( f[idx], xi_0[idx_xi], f[idx+1], xi_1[idx_xi], f[idx+2], xi_2[idx_xi], f[idx+3], xi_3[idx_xi], f[idx+4], xi_4[idx_xi], f[idx+5], xi_5[idx_xi], f[idx+6], xi_6[idx_xi] );
//如下定義為對應D3Q19模型的F6
#define F4_Intrpl7(f, j, k, j_c, k_c, idx_y, idx_xi, xi_0, xi_1, xi_2, xi_3, xi_4, xi_5, xi_6)    \
    idx = j*nface + k_c ;    \
    F4_in = Intrpl7( f[idx], xi_0[idx_xi], f[idx+1], xi_1[idx_xi], f[idx+2], xi_2[idx_xi], f[idx+3], xi_3[idx_xi], f[idx+4], xi_4[idx_xi], f[idx+5], xi_5[idx_xi], f[idx+6], xi_6[idx_xi] );
//如下定義為對應D3Q19模型的F15,16,17,18對於D2Q9模型則對應到F5,F6,F7,F8
#define Y_XI_Intrpl7(f, F_in, i, j, k, i_c, j_c, k_c, idx_x, idx_y, idx_xi, y_0, y_1, y_2, y_3, y_4, y_5, y_6, xi_0, xi_1, xi_2, xi_3, xi_4, xi_5, xi_6)    \
    idx = j_c*nface + k_c ;  \
    F_in = Intrpl7(    /*//idx_x , idx_y , idx_xi作為權重陣列的第二編號*/\
        Intrpl7( f[idx],   y_0[idx_y], f[idx+nface],   y_1[idx_y], f[idx+2*nface],   y_2[idx_y], f[idx+3*nface],   y_3[idx_y], f[idx+4*nface],   y_4[idx_y], f[idx+5*nface],   y_5[idx_y], f[idx+6*nface],   y_6[idx_y] ), xi_0[idx_xi],   \
        Intrpl7( f[idx+1], y_0[idx_y], f[idx+nface+1], y_1[idx_y], f[idx+2*nface+1], y_2[idx_y], f[idx+3*nface+1], y_3[idx_y], f[idx+4*nface+1], y_4[idx_y], f[idx+5*nface+1], y_5[idx_y], f[idx+6*nface+1], y_6[idx_y] ), xi_1[idx_xi],   \
        Intrpl7( f[idx+2], y_0[idx_y], f[idx+nface+2], y_1[idx_y], f[idx+2*nface+2], y_2[idx_y], f[idx+3*nface+2], y_3[idx_y], f[idx+4*nface+2], y_4[idx_y], f[idx+5*nface+2], y_5[idx_y], f[idx+6*nface+2], y_6[idx_y] ), xi_2[idx_xi],   \
        Intrpl7( f[idx+3], y_0[idx_y], f[idx+nface+3], y_1[idx_y], f[idx+2*nface+3], y_2[idx_y], f[idx+3*nface+3], y_3[idx_y], f[idx+4*nface+3], y_4[idx_y], f[idx+5*nface+3], y_5[idx_y], f[idx+6*nface+3], y_6[idx_y] ), xi_3[idx_xi],   \
        Intrpl7( f[idx+4], y_0[idx_y], f[idx+nface+4], y_1[idx_y], f[idx+2*nface+4], y_2[idx_y], f[idx+3*nface+4], y_3[idx_y], f[idx+4*nface+4], y_4[idx_y], f[idx+5*nface+4], y_5[idx_y], f[idx+6*nface+4], y_6[idx_y] ), xi_4[idx_xi],   \
        Intrpl7( f[idx+5], y_0[idx_y], f[idx+nface+5], y_1[idx_y], f[idx+2*nface+5], y_2[idx_y], f[idx+3*nface+5], y_3[idx_y], f[idx+4*nface+5], y_4[idx_y], f[idx+5*nface+5], y_5[idx_y], f[idx+6*nface+5], y_6[idx_y] ), xi_5[idx_xi],   \
        Intrpl7( f[idx+6], y_0[idx_y], f[idx+nface+6], y_1[idx_y], f[idx+2*nface+6], y_2[idx_y], f[idx+3*nface+6], y_3[idx_y], f[idx+4*nface+6], y_4[idx_y], f[idx+5*nface+6], y_5[idx_y], f[idx+6*nface+6], y_6[idx_y] ), xi_6[idx_xi]    \
    );
#endif