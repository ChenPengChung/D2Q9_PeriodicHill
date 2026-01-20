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
#define F1_Intrpl7(f,j,k,j_c,k_c,idx_y,idx_xi,y_0,y_1,y_2,y_3,y_4,y_5,y_6,xi_0,xi_1,xi_2,xi_3,xi_4,xi_5,xi_6) \
    do { \
        const int idx_r0 = (j_c + 0) * nface + k_c; \
        const int idx_r1 = (j_c + 1) * nface + k_c; \
        const int idx_r2 = (j_c + 2) * nface + k_c; \
        const int idx_r3 = (j_c + 3) * nface + k_c; \
        const int idx_r4 = (j_c + 4) * nface + k_c; \
        const int idx_r5 = (j_c + 5) * nface + k_c; \
        const int idx_r6 = (j_c + 6) * nface + k_c; \
        F1_in = Intrpl7( \
            Intrpl7( f[idx_r0],         xi_0[idx_r0],         f[idx_r0+1],         xi_1[idx_r0],         f[idx_r0+2],         xi_2[idx_r0],         f[idx_r0+3],         xi_3[idx_r0],         f[idx_r0+4],         xi_4[idx_r0],         f[idx_r0+5],         xi_5[idx_r0],         f[idx_r0+6],         xi_6[idx_r0] ), y_0[idx_y], \
            Intrpl7( f[idx_r1],         xi_0[idx_r1],         f[idx_r1+1],         xi_1[idx_r1],         f[idx_r1+2],         xi_2[idx_r1],         f[idx_r1+3],         xi_3[idx_r1],         f[idx_r1+4],         xi_4[idx_r1],         f[idx_r1+5],         xi_5[idx_r1],         f[idx_r1+6],         xi_6[idx_r1] ), y_1[idx_y], \
            Intrpl7( f[idx_r2],         xi_0[idx_r2],         f[idx_r2+1],         xi_1[idx_r2],         f[idx_r2+2],         xi_2[idx_r2],         f[idx_r2+3],         xi_3[idx_r2],         f[idx_r2+4],         xi_4[idx_r2],         f[idx_r2+5],         xi_5[idx_r2],         f[idx_r2+6],         xi_6[idx_r2] ), y_2[idx_y], \
            Intrpl7( f[idx_r3],         xi_0[idx_r3],         f[idx_r3+1],         xi_1[idx_r3],         f[idx_r3+2],         xi_2[idx_r3],         f[idx_r3+3],         xi_3[idx_r3],         f[idx_r3+4],         xi_4[idx_r3],         f[idx_r3+5],         xi_5[idx_r3],         f[idx_r3+6],         xi_6[idx_r3] ), y_3[idx_y], \
            Intrpl7( f[idx_r4],         xi_0[idx_r4],         f[idx_r4+1],         xi_1[idx_r4],         f[idx_r4+2],         xi_2[idx_r4],         f[idx_r4+3],         xi_3[idx_r4],         f[idx_r4+4],         xi_4[idx_r4],         f[idx_r4+5],         xi_5[idx_r4],         f[idx_r4+6],         xi_6[idx_r4] ), y_4[idx_y], \
            Intrpl7( f[idx_r5],         xi_0[idx_r5],         f[idx_r5+1],         xi_1[idx_r5],         f[idx_r5+2],         xi_2[idx_r5],         f[idx_r5+3],         xi_3[idx_r5],         f[idx_r5+4],         xi_4[idx_r5],         f[idx_r5+5],         xi_5[idx_r5],         f[idx_r5+6],         xi_6[idx_r5] ), y_5[idx_y], \
            Intrpl7( f[idx_r6],         xi_0[idx_r6],         f[idx_r6+1],         xi_1[idx_r6],         f[idx_r6+2],         xi_2[idx_r6],         f[idx_r6+3],         xi_3[idx_r6],         f[idx_r6+4],         xi_4[idx_r6],         f[idx_r6+5],         xi_5[idx_r6],         f[idx_r6+6],         xi_6[idx_r6] ), y_6[idx_y]  \
        ); \
    } while(0);
//如下定義對應D3Q19模型的F4
#define F3_Intrpl7(f,j,k,j_c,k_c,idx_y,idx_xi,y_0,y_1,y_2,y_3,y_4,y_5,y_6,xi_0,xi_1,xi_2,xi_3,xi_4,xi_5,xi_6) \
    do { \
        const int idx_r0 = (j_c + 0) * nface + k_c; \
        const int idx_r1 = (j_c + 1) * nface + k_c; \
        const int idx_r2 = (j_c + 2) * nface + k_c; \
        const int idx_r3 = (j_c + 3) * nface + k_c; \
        const int idx_r4 = (j_c + 4) * nface + k_c; \
        const int idx_r5 = (j_c + 5) * nface + k_c; \
        const int idx_r6 = (j_c + 6) * nface + k_c; \
        F3_in = Intrpl7( \
            Intrpl7( f[idx_r0],         xi_0[idx_r0],         f[idx_r0+1],         xi_1[idx_r0],         f[idx_r0+2],         xi_2[idx_r0],         f[idx_r0+3],         xi_3[idx_r0],         f[idx_r0+4],         xi_4[idx_r0],         f[idx_r0+5],         xi_5[idx_r0],         f[idx_r0+6],         xi_6[idx_r0] ), y_0[idx_y], \
            Intrpl7( f[idx_r1],         xi_0[idx_r1],         f[idx_r1+1],         xi_1[idx_r1],         f[idx_r1+2],         xi_2[idx_r1],         f[idx_r1+3],         xi_3[idx_r1],         f[idx_r1+4],         xi_4[idx_r1],         f[idx_r1+5],         xi_5[idx_r1],         f[idx_r1+6],         xi_6[idx_r1] ), y_1[idx_y], \
            Intrpl7( f[idx_r2],         xi_0[idx_r2],         f[idx_r2+1],         xi_1[idx_r2],         f[idx_r2+2],         xi_2[idx_r2],         f[idx_r2+3],         xi_3[idx_r2],         f[idx_r2+4],         xi_4[idx_r2],         f[idx_r2+5],         xi_5[idx_r2],         f[idx_r2+6],         xi_6[idx_r2] ), y_2[idx_y], \
            Intrpl7( f[idx_r3],         xi_0[idx_r3],         f[idx_r3+1],         xi_1[idx_r3],         f[idx_r3+2],         xi_2[idx_r3],         f[idx_r3+3],         xi_3[idx_r3],         f[idx_r3+4],         xi_4[idx_r3],         f[idx_r3+5],         xi_5[idx_r3],         f[idx_r3+6],         xi_6[idx_r3] ), y_3[idx_y], \
            Intrpl7( f[idx_r4],         xi_0[idx_r4],         f[idx_r4+1],         xi_1[idx_r4],         f[idx_r4+2],         xi_2[idx_r4],         f[idx_r4+3],         xi_3[idx_r4],         f[idx_r4+4],         xi_4[idx_r4],         f[idx_r4+5],         xi_5[idx_r4],         f[idx_r4+6],         xi_6[idx_r4] ), y_4[idx_y], \
            Intrpl7( f[idx_r5],         xi_0[idx_r5],         f[idx_r5+1],         xi_1[idx_r5],         f[idx_r5+2],         xi_2[idx_r5],         f[idx_r5+3],         xi_3[idx_r5],         f[idx_r5+4],         xi_4[idx_r5],         f[idx_r5+5],         xi_5[idx_r5],         f[idx_r5+6],         xi_6[idx_r5] ), y_5[idx_y], \
            Intrpl7( f[idx_r6],         xi_0[idx_r6],         f[idx_r6+1],         xi_1[idx_r6],         f[idx_r6+2],         xi_2[idx_r6],         f[idx_r6+3],         xi_3[idx_r6],         f[idx_r6+4],         xi_4[idx_r6],         f[idx_r6+5],         xi_5[idx_r6],         f[idx_r6+6],         xi_6[idx_r6] ), y_6[idx_y]  \
        ); \
    } while(0);
//如下定義為對應D3Q19模型的F5
#define F2_Intrpl7(f, j, k, j_c, k_c, idx_y, idx_xi, xi_0, xi_1, xi_2, xi_3, xi_4, xi_5, xi_6)    \
    idx = j*nface + k_c ;    \
    F2_in = Intrpl7( f[idx], xi_0[idx_xi], f[idx+1], xi_1[idx_xi], f[idx+2], xi_2[idx_xi], f[idx+3], xi_3[idx_xi], f[idx+4], xi_4[idx_xi], f[idx+5], xi_5[idx_xi], f[idx+6], xi_6[idx_xi] );
//如下定義為對應D3Q19模型的F6
#define F4_Intrpl7(f, j, k, j_c, k_c, idx_y, idx_xi, xi_0, xi_1, xi_2, xi_3, xi_4, xi_5, xi_6)    \
    idx = j*nface + k_c ;    \
    F4_in = Intrpl7( f[idx], xi_0[idx_xi], f[idx+1], xi_1[idx_xi], f[idx+2], xi_2[idx_xi], f[idx+3], xi_3[idx_xi], f[idx+4], xi_4[idx_xi], f[idx+5], xi_5[idx_xi], f[idx+6], xi_6[idx_xi] );
//如下定義為對應D3Q19模型的F15,16,17,18對於D2Q9模型則對應到F5,F6,F7,F8
#define Y_XI_Intrpl7(f, F_in, j, k,  j_c, k_c, idx_y, idx_xi, y_0, y_1, y_2, y_3, y_4, y_5, y_6, xi_0, xi_1, xi_2, xi_3, xi_4, xi_5, xi_6)    \
    do { \
        const int idx_r0 = (j_c + 0) * nface + k_c; \
        const int idx_r1 = (j_c + 1) * nface + k_c; \
        const int idx_r2 = (j_c + 2) * nface + k_c; \
        const int idx_r3 = (j_c + 3) * nface + k_c; \
        const int idx_r4 = (j_c + 4) * nface + k_c; \
        const int idx_r5 = (j_c + 5) * nface + k_c; \
        const int idx_r6 = (j_c + 6) * nface + k_c; \
        F_in = Intrpl7( \
            Intrpl7( f[idx_r0],         xi_0[idx_r0],         f[idx_r0+1],         xi_1[idx_r0],         f[idx_r0+2],         xi_2[idx_r0],         f[idx_r0+3],         xi_3[idx_r0],         f[idx_r0+4],         xi_4[idx_r0],         f[idx_r0+5],         xi_5[idx_r0],         f[idx_r0+6],         xi_6[idx_r0] ), y_0[idx_y], \
            Intrpl7( f[idx_r1],         xi_0[idx_r1],         f[idx_r1+1],         xi_1[idx_r1],         f[idx_r1+2],         xi_2[idx_r1],         f[idx_r1+3],         xi_3[idx_r1],         f[idx_r1+4],         xi_4[idx_r1],         f[idx_r1+5],         xi_5[idx_r1],         f[idx_r1+6],         xi_6[idx_r1] ), y_1[idx_y], \
            Intrpl7( f[idx_r2],         xi_0[idx_r2],         f[idx_r2+1],         xi_1[idx_r2],         f[idx_r2+2],         xi_2[idx_r2],         f[idx_r2+3],         xi_3[idx_r2],         f[idx_r2+4],         xi_4[idx_r2],         f[idx_r2+5],         xi_5[idx_r2],         f[idx_r2+6],         xi_6[idx_r2] ), y_2[idx_y], \
            Intrpl7( f[idx_r3],         xi_0[idx_r3],         f[idx_r3+1],         xi_1[idx_r3],         f[idx_r3+2],         xi_2[idx_r3],         f[idx_r3+3],         xi_3[idx_r3],         f[idx_r3+4],         xi_4[idx_r3],         f[idx_r3+5],         xi_5[idx_r3],         f[idx_r3+6],         xi_6[idx_r3] ), y_3[idx_y], \
            Intrpl7( f[idx_r4],         xi_0[idx_r4],         f[idx_r4+1],         xi_1[idx_r4],         f[idx_r4+2],         xi_2[idx_r4],         f[idx_r4+3],         xi_3[idx_r4],         f[idx_r4+4],         xi_4[idx_r4],         f[idx_r4+5],         xi_5[idx_r4],         f[idx_r4+6],         xi_6[idx_r4] ), y_4[idx_y], \
            Intrpl7( f[idx_r5],         xi_0[idx_r5],         f[idx_r5+1],         xi_1[idx_r5],         f[idx_r5+2],         xi_2[idx_r5],         f[idx_r5+3],         xi_3[idx_r5],         f[idx_r5+4],         xi_4[idx_r5],         f[idx_r5+5],         xi_5[idx_r5],         f[idx_r5+6],         xi_6[idx_r5] ), y_5[idx_y], \
            Intrpl7( f[idx_r6],         xi_0[idx_r6],         f[idx_r6+1],         xi_1[idx_r6],         f[idx_r6+2],         xi_2[idx_r6],         f[idx_r6+3],         xi_3[idx_r6],         f[idx_r6+4],         xi_4[idx_r6],         f[idx_r6+5],         xi_5[idx_r6],         f[idx_r6+6],         xi_6[idx_r6] ), y_6[idx_y]  \
        ); \
    } while(0);
#endif
