#ifndef INTERPOLATIONHILLISLBM_3PT_FILE
#define INTERPOLATIONHILLISLBM_3PT_FILE
//=============================================================================
// 2階（3點）Lagrange 插值版本
// 優點：減少 Runge 振盪，計算更快，記憶體需求更少
// 適用於 CFL < 0.5 的情況
//=============================================================================

// 3點線性插值公式
#define Intrpl3(f1, a1, f2, a2, f3, a3) \
    (f1*a1 + f2*a2 + f3*a3)

// F0: 靜止粒子，不需要插值
#define F0_Intrpl3(f, j, k) \
    idx = j * nface + k; \
    F0_in = f[idx];

// F1 (+Y方向): 2D插值，先Xi後Y
#define F1_Intrpl3(f, j, k, j_c, k_c, idx_y, idx_xi, y_0, y_1, y_2, xi_0, xi_1, xi_2) \
do { /*降階版本*/\
    const int idx_f = j_c * nface + k_c; /*j_c = j-1 ; k_c = cell_z */\
    F1_in = Intrpl3( \
        Intrpl3( f[idx_f],             xi_0[idx_xi+0*NY6*nface], f[idx_f+1],             xi_1[idx_xi+0*NY6*nface], f[idx_f+2],             xi_2[idx_xi+0*NY6*nface]), y_0[idx_y], \
        Intrpl3( f[idx_f+1*nface],     xi_0[idx_xi+1*NY6*nface], f[idx_f+1*nface+1],     xi_1[idx_xi+1*NY6*nface], f[idx_f+1*nface+2],     xi_2[idx_xi+1*NY6*nface]), y_1[idx_y], \
        Intrpl3( f[idx_f+2*nface],     xi_0[idx_xi+2*NY6*nface], f[idx_f+2*nface+1],     xi_1[idx_xi+2*NY6*nface], f[idx_f+2*nface+2],     xi_2[idx_xi+2*NY6*nface]), y_2[idx_y]  \
    ); \
} while(0)

// F3 (-Y方向): 2D插值，先Xi後Y
#define F3_Intrpl3(f, j, k, j_c, k_c, idx_y, idx_xi, y_0, y_1, y_2, xi_0, xi_1, xi_2) \
do { /*降階版本*/\
    const int idx_f = j_c * nface + k_c; /*j_c = j-1 ; k_c = cell_z */\
    F3_in = Intrpl3( \
        Intrpl3( f[idx_f],             xi_0[idx_xi+0*NY6*nface], f[idx_f+1],             xi_1[idx_xi+0*NY6*nface], f[idx_f+2],             xi_2[idx_xi+0*NY6*nface]), y_0[idx_y], \
        Intrpl3( f[idx_f+1*nface],     xi_0[idx_xi+1*NY6*nface], f[idx_f+1*nface+1],     xi_1[idx_xi+1*NY6*nface], f[idx_f+1*nface+2],     xi_2[idx_xi+1*NY6*nface]), y_1[idx_y], \
        Intrpl3( f[idx_f+2*nface],     xi_0[idx_xi+2*NY6*nface], f[idx_f+2*nface+1],     xi_1[idx_xi+2*NY6*nface], f[idx_f+2*nface+2],     xi_2[idx_xi+2*NY6*nface]), y_2[idx_y]  \
    ); \
} while(0)

// F2 (+Z方向): 1D Xi插值
#define F2_Intrpl3(f, j, k, j_c, k_c, idx_y, idx_xi, xi_0, xi_1, xi_2) \
do { \
    const int idx_f = j * nface + k_c; \
    F2_in = Intrpl3( f[idx_f], xi_0[idx_xi+1*NY6*nface], f[idx_f+1], xi_1[idx_xi+1*NY6*nface], f[idx_f+2], xi_2[idx_xi+1*NY6*nface]); \
} while(0)
// F4 (-Z方向): 1D Xi插值
#define F4_Intrpl3(f, j, k, j_c, k_c, idx_y, idx_xi, xi_0, xi_1, xi_2) \
do { \
    const int idx_f = j * nface + k_c; \
    F4_in = Intrpl3( f[idx_f], xi_0[idx_xi+1*NY6*nface], f[idx_f+1], xi_1[idx_xi+1*NY6*nface], f[idx_f+2], xi_2[idx_xi+1*NY6*nface]); \
} while(0)
// F5,F6,F7,F8 (對角方向): 2D插值，先Xi後Y
#define Y_XI_Intrpl3(f, F_in, j, k, j_c, k_c, idx_y, idx_xi, y_0, y_1, y_2, xi_0, xi_1, xi_2) \
do { /*降階版本*/\
    const int idx_f = j_c * nface + k_c; /*j_c = j-1 ; k_c = cell_z */\
    F_in = Intrpl3( \
        Intrpl3( f[idx_f],             xi_0[idx_xi+0*NY6*nface], f[idx_f+1],             xi_1[idx_xi+0*NY6*nface], f[idx_f+2],             xi_2[idx_xi+0*NY6*nface]), y_0[idx_y], \
        Intrpl3( f[idx_f+1*nface],     xi_0[idx_xi+1*NY6*nface], f[idx_f+1*nface+1],     xi_1[idx_xi+1*NY6*nface], f[idx_f+1*nface+2],     xi_2[idx_xi+1*NY6*nface]), y_1[idx_y], \
        Intrpl3( f[idx_f+2*nface],     xi_0[idx_xi+2*NY6*nface], f[idx_f+2*nface+1],     xi_1[idx_xi+2*NY6*nface], f[idx_f+2*nface+2],     xi_2[idx_xi+2*NY6*nface]), y_2[idx_y]  \
    ); \
} while(0)

#endif
