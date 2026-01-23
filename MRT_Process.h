#ifndef MRT_Process_FILE
#define MRT_Process_FILE
#include "MRT_Matrix.h"

//=============================================================================
// D2Q9 平衡態分佈函數 (Equilibrium Distribution Function)
// 輸入需要: rho_local, uy, uz (密度與速度)
// 輸出: F0_eq ~ F8_eq
//
// D2Q9 速度方向:
// F0: (0, 0)      w0 = 4/9
// F1: (+1, 0)     w1 = 1/9   (+Y)
// F2: (0, +1)     w2 = 1/9   (+Z)
// F3: (-1, 0)     w3 = 1/9   (-Y)
// F4: (0, -1)     w4 = 1/9   (-Z)
// F5: (+1, +1)    w5 = 1/36  (+Y, +Z)
// F6: (-1, +1)    w6 = 1/36  (-Y, +Z)
// F7: (-1, -1)    w7 = 1/36  (-Y, -Z)
// F8: (+1, -1)    w8 = 1/36  (+Y, -Z)
//=============================================================================
#define Equilibrium \
    double udot = uy*uy + uz*uz; \
    double F0_eq = (4.0/9.0)  * rho_local * (1.0 - 1.5*udot); \
    double F1_eq = (1.0/9.0)  * rho_local * (1.0 + 3.0*uy + 4.5*uy*uy - 1.5*udot); \
    double F2_eq = (1.0/9.0)  * rho_local * (1.0 + 3.0*uz + 4.5*uz*uz - 1.5*udot); \
    double F3_eq = (1.0/9.0)  * rho_local * (1.0 - 3.0*uy + 4.5*uy*uy - 1.5*udot); \
    double F4_eq = (1.0/9.0)  * rho_local * (1.0 - 3.0*uz + 4.5*uz*uz - 1.5*udot); \
    double F5_eq = (1.0/36.0) * rho_local * (1.0 + 3.0*(uy+uz) + 4.5*(uy+uz)*(uy+uz) - 1.5*udot); \
    double F6_eq = (1.0/36.0) * rho_local * (1.0 + 3.0*(-uy+uz) + 4.5*(-uy+uz)*(-uy+uz) - 1.5*udot); \
    double F7_eq = (1.0/36.0) * rho_local * (1.0 + 3.0*(-uy-uz) + 4.5*(-uy-uz)*(-uy-uz) - 1.5*udot); \
    double F8_eq = (1.0/36.0) * rho_local * (1.0 + 3.0*(uy-uz) + 4.5*(uy-uz)*(uy-uz) - 1.5*udot)

//=============================================================================
// 分佈函數向量矩 (Moment Vector)
// m = M * f
// 輸入需要: F0_in ~ F8_in, M[9][9]
// 輸出: m0 ~ m8
//=============================================================================
#define m_vector \
    m0 = F0_in*M[0][0] + F1_in*M[0][1] + F2_in*M[0][2] + F3_in*M[0][3] + F4_in*M[0][4] + F5_in*M[0][5] + F6_in*M[0][6] + F7_in*M[0][7] + F8_in*M[0][8]; \
    m1 = F0_in*M[1][0] + F1_in*M[1][1] + F2_in*M[1][2] + F3_in*M[1][3] + F4_in*M[1][4] + F5_in*M[1][5] + F6_in*M[1][6] + F7_in*M[1][7] + F8_in*M[1][8]; \
    m2 = F0_in*M[2][0] + F1_in*M[2][1] + F2_in*M[2][2] + F3_in*M[2][3] + F4_in*M[2][4] + F5_in*M[2][5] + F6_in*M[2][6] + F7_in*M[2][7] + F8_in*M[2][8]; \
    m3 = F0_in*M[3][0] + F1_in*M[3][1] + F2_in*M[3][2] + F3_in*M[3][3] + F4_in*M[3][4] + F5_in*M[3][5] + F6_in*M[3][6] + F7_in*M[3][7] + F8_in*M[3][8]; \
    m4 = F0_in*M[4][0] + F1_in*M[4][1] + F2_in*M[4][2] + F3_in*M[4][3] + F4_in*M[4][4] + F5_in*M[4][5] + F6_in*M[4][6] + F7_in*M[4][7] + F8_in*M[4][8]; \
    m5 = F0_in*M[5][0] + F1_in*M[5][1] + F2_in*M[5][2] + F3_in*M[5][3] + F4_in*M[5][4] + F5_in*M[5][5] + F6_in*M[5][6] + F7_in*M[5][7] + F8_in*M[5][8]; \
    m6 = F0_in*M[6][0] + F1_in*M[6][1] + F2_in*M[6][2] + F3_in*M[6][3] + F4_in*M[6][4] + F5_in*M[6][5] + F6_in*M[6][6] + F7_in*M[6][7] + F8_in*M[6][8]; \
    m7 = F0_in*M[7][0] + F1_in*M[7][1] + F2_in*M[7][2] + F3_in*M[7][3] + F4_in*M[7][4] + F5_in*M[7][5] + F6_in*M[7][6] + F7_in*M[7][7] + F8_in*M[7][8]; \
    m8 = F0_in*M[8][0] + F1_in*M[8][1] + F2_in*M[8][2] + F3_in*M[8][3] + F4_in*M[8][4] + F5_in*M[8][5] + F6_in*M[8][6] + F7_in*M[8][7] + F8_in*M[8][8]

//=============================================================================
// 平衡態分佈函數的矩 (Equilibrium Moment Vector)
// meq = M * feq
// 輸入需要: F0_eq ~ F8_eq (需先呼叫 Equilibrium), M[9][9]
// 輸出: meq0 ~ meq8
//=============================================================================
#define meq \
    meq0 = F0_eq*M[0][0] + F1_eq*M[0][1] + F2_eq*M[0][2] + F3_eq*M[0][3] + F4_eq*M[0][4] + F5_eq*M[0][5] + F6_eq*M[0][6] + F7_eq*M[0][7] + F8_eq*M[0][8] ; \
    meq1 = F0_eq*M[1][0] + F1_eq*M[1][1] + F2_eq*M[1][2] + F3_eq*M[1][3] + F4_eq*M[1][4] + F5_eq*M[1][5] + F6_eq*M[1][6] + F7_eq*M[1][7] + F8_eq*M[1][8] ; \
    meq2 = F0_eq*M[2][0] + F1_eq*M[2][1] + F2_eq*M[2][2] + F3_eq*M[2][3] + F4_eq*M[2][4] + F5_eq*M[2][5] + F6_eq*M[2][6] + F7_eq*M[2][7] + F8_eq*M[2][8] ; \
    meq3 = F0_eq*M[3][0] + F1_eq*M[3][1] + F2_eq*M[3][2] + F3_eq*M[3][3] + F4_eq*M[3][4] + F5_eq*M[3][5] + F6_eq*M[3][6] + F7_eq*M[3][7] + F8_eq*M[3][8] ; \
    meq4 = F0_eq*M[4][0] + F1_eq*M[4][1] + F2_eq*M[4][2] + F3_eq*M[4][3] + F4_eq*M[4][4] + F5_eq*M[4][5] + F6_eq*M[4][6] + F7_eq*M[4][7] + F8_eq*M[4][8] ; \
    meq5 = F0_eq*M[5][0] + F1_eq*M[5][1] + F2_eq*M[5][2] + F3_eq*M[5][3] + F4_eq*M[5][4] + F5_eq*M[5][5] + F6_eq*M[5][6] + F7_eq*M[5][7] + F8_eq*M[5][8] ; \
    meq6 = F0_eq*M[6][0] + F1_eq*M[6][1] + F2_eq*M[6][2] + F3_eq*M[6][3] + F4_eq*M[6][4] + F5_eq*M[6][5] + F6_eq*M[6][6] + F7_eq*M[6][7] + F8_eq*M[6][8] ; \
    meq7 = F0_eq*M[7][0] + F1_eq*M[7][1] + F2_eq*M[7][2] + F3_eq*M[7][3] + F4_eq*M[7][4] + F5_eq*M[7][5] + F6_eq*M[7][6] + F7_eq*M[7][7] + F8_eq*M[7][8] ; \
    meq8 = F0_eq*M[8][0] + F1_eq*M[8][1] + F2_eq*M[8][2] + F3_eq*M[8][3] + F4_eq*M[8][4] + F5_eq*M[8][5] + F6_eq*M[8][6] + F7_eq*M[8][7] + F8_eq*M[8][8] ; 

//=============================================================================
// MRT 碰撞 + Guo Forcing Scheme (正確版本)
// 
// Guo forcing 完整公式:
// F_i = (1 - 0.5*omega) * w_i * [ 3*(e_i - u)·F + 9*(e_i·u)(e_i·F) ]
// 
// 展開 (只有 Y 方向外力 Fy):
// F_i = (1 - 0.5*omega) * w_i * [ 3*(e_iy - ux)*Fy + 9*(e_iy*ux + e_iz*uy)*(e_iy*Fy) ]
//
// D2Q9 速度方向 (Y,Z):
//   F0: e=( 0, 0)  F1: e=(+1, 0)  F2: e=( 0,+1)  F3: e=(-1, 0)  F4: e=( 0,-1)
//   F5: e=(+1,+1)  F6: e=(-1,+1)  F7: e=(-1,-1)  F8: e=(+1,-1)
//=============================================================================
#define collision	\
    double omega_eff = omega_7; \
    double guo_factor = (1.0 - 0.5 * omega_eff); \
    double Fy = Force[0]; \
    double ux = v1_raw; \
    double uy = w1_raw; \
    /* F0: e=(0,0), e_iy=0, e_iz=0 */ \
    /* 3*(0-ux)*Fy + 9*(0)*0 = -3*ux*Fy */ \
    double guo0 = guo_factor * (4.0/9.0)  * (3.0 * ((0.0 - ux) * Fy)); \
    /* F1: e=(+1,0), e_iy=+1, e_iz=0 */ \
    /* 3*(1-ux)*Fy + 9*(1*ux+0*uy)*(1*Fy) = 3*(1-ux)*Fy + 9*ux*Fy */ \
    double guo1 = guo_factor * (1.0/9.0)  * (3.0 * ((1.0 - ux) * Fy) + 9.0 * (ux) * (Fy)); \
    /* F2: e=(0,+1), e_iy=0, e_iz=+1 */ \
    /* 3*(0-ux)*Fy + 9*(0*ux+1*uy)*(0*Fy) = -3*ux*Fy + 0 */ \
    double guo2 = guo_factor * (1.0/9.0)  * (3.0 * ((0.0 - ux) * Fy)); \
    /* F3: e=(-1,0), e_iy=-1, e_iz=0 */ \
    /* 3*(-1-ux)*Fy + 9*(-1*ux+0*uy)*(-1*Fy) = 3*(-1-ux)*Fy + 9*ux*Fy */ \
    double guo3 = guo_factor * (1.0/9.0)  * (3.0 * ((-1.0 - ux) * Fy) + 9.0 * (ux) * (Fy)); \
    /* F4: e=(0,-1), e_iy=0, e_iz=-1 */ \
    /* 3*(0-ux)*Fy + 9*(0*ux-1*uy)*(0*Fy) = -3*ux*Fy + 0 */ \
    double guo4 = guo_factor * (1.0/9.0)  * (3.0 * ((0.0 - ux) * Fy)); \
    /* F5: e=(+1,+1), e_iy=+1, e_iz=+1 */ \
    /* 3*(1-ux)*Fy + 9*(1*ux+1*uy)*(1*Fy) = 3*(1-ux)*Fy + 9*(ux+uy)*Fy */ \
    double guo5 = guo_factor * (1.0/36.0) * (3.0 * ((1.0 - ux) * Fy) + 9.0 * (ux + uy) * (Fy)); \
    /* F6: e=(-1,+1), e_iy=-1, e_iz=+1 */ \
    /* 3*(-1-ux)*Fy + 9*(-1*ux+1*uy)*(-1*Fy) = 3*(-1-ux)*Fy + 9*(ux-uy)*Fy */ \
    double guo6 = guo_factor * (1.0/36.0) * (3.0 * ((-1.0 - ux) * Fy) + 9.0 * (ux - uy) * (Fy)); \
    /* F7: e=(-1,-1), e_iy=-1, e_iz=-1 */ \
    /* 3*(-1-ux)*Fy + 9*(-1*ux-1*uy)*(-1*Fy) = 3*(-1-ux)*Fy + 9*(ux+uy)*Fy */ \
    double guo7 = guo_factor * (1.0/36.0) * (3.0 * ((-1.0 - ux) * Fy) + 9.0 * (ux + uy) * (Fy)); \
    /* F8: e=(+1,-1), e_iy=+1, e_iz=-1 */ \
    /* 3*(1-ux)*Fy + 9*(1*ux-1*uy)*(1*Fy) = 3*(1-ux)*Fy + 9*(ux-uy)*Fy */ \
    double guo8 = guo_factor * (1.0/36.0) * (3.0 * ((1.0 - ux) * Fy) + 9.0 * (ux - uy) * (Fy)); \
    F0_in  = F0_in - M_I[0][0]*s0*(m0-meq0)- M_I[0][1]*s1*(m1-meq1)- M_I[0][2]*s2*(m2-meq2)- M_I[0][3]*s3*(m3-meq3)- M_I[0][4]*s4*(m4-meq4)- M_I[0][5]*s5*(m5-meq5)- M_I[0][6]*s6*(m6-meq6)- M_I[0][7]*s7*(m7-meq7)- M_I[0][8]*s8*(m8-meq8) + guo0; \
    F1_in  = F1_in - M_I[1][0]*s0*(m0-meq0)- M_I[1][1]*s1*(m1-meq1)- M_I[1][2]*s2*(m2-meq2)- M_I[1][3]*s3*(m3-meq3)- M_I[1][4]*s4*(m4-meq4)- M_I[1][5]*s5*(m5-meq5)- M_I[1][6]*s6*(m6-meq6)- M_I[1][7]*s7*(m7-meq7)- M_I[1][8]*s8*(m8-meq8) + guo1; \
    F2_in  = F2_in - M_I[2][0]*s0*(m0-meq0)- M_I[2][1]*s1*(m1-meq1)- M_I[2][2]*s2*(m2-meq2)- M_I[2][3]*s3*(m3-meq3)- M_I[2][4]*s4*(m4-meq4)- M_I[2][5]*s5*(m5-meq5)- M_I[2][6]*s6*(m6-meq6)- M_I[2][7]*s7*(m7-meq7)- M_I[2][8]*s8*(m8-meq8) + guo2; \
    F3_in  = F3_in - M_I[3][0]*s0*(m0-meq0)- M_I[3][1]*s1*(m1-meq1)- M_I[3][2]*s2*(m2-meq2)- M_I[3][3]*s3*(m3-meq3)- M_I[3][4]*s4*(m4-meq4)- M_I[3][5]*s5*(m5-meq5)- M_I[3][6]*s6*(m6-meq6)- M_I[3][7]*s7*(m7-meq7)- M_I[3][8]*s8*(m8-meq8) + guo3; \
    F4_in  = F4_in - M_I[4][0]*s0*(m0-meq0)- M_I[4][1]*s1*(m1-meq1)- M_I[4][2]*s2*(m2-meq2)- M_I[4][3]*s3*(m3-meq3)- M_I[4][4]*s4*(m4-meq4)- M_I[4][5]*s5*(m5-meq5)- M_I[4][6]*s6*(m6-meq6)- M_I[4][7]*s7*(m7-meq7)- M_I[4][8]*s8*(m8-meq8) + guo4; \
    F5_in  = F5_in - M_I[5][0]*s0*(m0-meq0)- M_I[5][1]*s1*(m1-meq1)- M_I[5][2]*s2*(m2-meq2)- M_I[5][3]*s3*(m3-meq3)- M_I[5][4]*s4*(m4-meq4)- M_I[5][5]*s5*(m5-meq5)- M_I[5][6]*s6*(m6-meq6)- M_I[5][7]*s7*(m7-meq7)- M_I[5][8]*s8*(m8-meq8) + guo5; \
    F6_in  = F6_in - M_I[6][0]*s0*(m0-meq0)- M_I[6][1]*s1*(m1-meq1)- M_I[6][2]*s2*(m2-meq2)- M_I[6][3]*s3*(m3-meq3)- M_I[6][4]*s4*(m4-meq4)- M_I[6][5]*s5*(m5-meq5)- M_I[6][6]*s6*(m6-meq6)- M_I[6][7]*s7*(m7-meq7)- M_I[6][8]*s8*(m8-meq8) + guo6; \
    F7_in  = F7_in - M_I[7][0]*s0*(m0-meq0)- M_I[7][1]*s1*(m1-meq1)- M_I[7][2]*s2*(m2-meq2)- M_I[7][3]*s3*(m3-meq3)- M_I[7][4]*s4*(m4-meq4)- M_I[7][5]*s5*(m5-meq5)- M_I[7][6]*s6*(m6-meq6)- M_I[7][7]*s7*(m7-meq7)- M_I[7][8]*s8*(m8-meq8) + guo7; \
    F8_in  = F8_in - M_I[8][0]*s0*(m0-meq0)- M_I[8][1]*s1*(m1-meq1)- M_I[8][2]*s2*(m2-meq2)- M_I[8][3]*s3*(m3-meq3)- M_I[8][4]*s4*(m4-meq4)- M_I[8][5]*s5*(m5-meq5)- M_I[8][6]*s6*(m6-meq6)- M_I[8][7]*s7*(m7-meq7)- M_I[8][8]*s8*(m8-meq8) + guo8; 

    #endif