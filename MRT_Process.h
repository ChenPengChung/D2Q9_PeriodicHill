#ifndef MRT_Process_FILE
#define MRT_Process_FILE
#include "MRT_Matrix.h"
//此篇為分佈函數向量矩，作為MRT算子中的一種分佈函數
//做完Moment以後，分佈函數的矩作為一維陣列，是一向量
#define m_vector \
    m0 = F0_in*M[0][0] + F1_in*M[0][1] + F2_in*M[0][2] + F3_in*M[0][3] + F4_in*M[0][4] + F5_in*M[0][5] + F6_in*M[0][6] + F7_in*M[0][7] + F8_in*M[0][8]; \
    m1 = F0_in*M[1][0] + F1_in*M[1][1] + F2_in*M[1][2] + F3_in*M[1][3] + F4_in*M[1][4] + F5_in*M[1][5] + F6_in*M[1][6] + F7_in*M[1][7] + F8_in*M[1][8]; \
    m2 = F0_in*M[2][0] + F1_in*M[2][1] + F2_in*M[2][2] + F3_in*M[2][3] + F4_in*M[2][4] + F5_in*M[2][5] + F6_in*M[2][6] + F7_in*M[2][7] + F8_in*M[2][8]; \
    m3 = F0_in*M[3][0] + F1_in*M[3][1] + F2_in*M[3][2] + F3_in*M[3][3] + F4_in*M[3][4] + F5_in*M[3][5] + F6_in*M[3][6] + F7_in*M[3][7] + F8_in*M[3][8]; \
    m4 = F0_in*M[4][0] + F1_in*M[4][1] + F2_in*M[4][2] + F3_in*M[4][3] + F4_in*M[4][4] + F5_in*M[4][5] + F6_in*M[4][6] + F7_in*M[4][7] + F8_in*M[4][8]; \
    m5 = F0_in*M[5][0] + F1_in*M[5][1] + F2_in*M[5][2] + F3_in*M[5][3] + F4_in*M[5][4] + F5_in*M[5][5] + F6_in*M[5][6] + F7_in*M[5][7] + F8_in*M[5][8]; \
    m6 = F0_in*M[6][0] + F1_in*M[6][1] + F2_in*M[6][2] + F3_in*M[6][3] + F4_in*M[6][4] + F5_in*M[6][5] + F6_in*M[6][6] + F7_in*M[6][7] + F8_in*M[6][8]; \
    m7 = F0_in*M[7][0] + F1_in*M[7][1] + F2_in*M[7][2] + F3_in*M[7][3] + F4_in*M[7][4] + F5_in*M[7][5] + F6_in*M[7][6] + F7_in*M[7][7] + F8_in*M[7][8]; \
    m8 = F0_in*M[8][0] + F1_in*M[8][1] + F2_in*M[8][2] + F3_in*M[8][3] + F4_in*M[8][4] + F5_in*M[8][5] + F6_in*M[8][6] + F7_in*M[8][7] + F8_in*M[8][8];  

#define meq	\
    meq0 = F0_eq*M[0][0] + F1_eq*M[0][1] + F2_eq*M[0][2] + F3_eq*M[0][3] + F4_eq*M[0][4] + F5_eq*M[0][5] + F6_eq*M[0][6] + F7_eq*M[0][7] + F8_eq*M[0][8] ; \
    meq1 = F0_eq*M[1][0] + F1_eq*M[1][1] + F2_eq*M[1][2] + F3_eq*M[1][3] + F4_eq*M[1][4] + F5_eq*M[1][5] + F6_eq*M[1][6] + F7_eq*M[1][7] + F8_eq*M[1][8] ; \
    meq2 = F0_eq*M[2][0] + F1_eq*M[2][1] + F2_eq*M[2][2] + F3_eq*M[2][3] + F4_eq*M[2][4] + F5_eq*M[2][5] + F6_eq*M[2][6] + F7_eq*M[2][7] + F8_eq*M[2][8] ; \
    meq3 = F0_eq*M[3][0] + F1_eq*M[3][1] + F2_eq*M[3][2] + F3_eq*M[3][3] + F4_eq*M[3][4] + F5_eq*M[3][5] + F6_eq*M[3][6] + F7_eq*M[3][7] + F8_eq*M[3][8] ; \
    meq4 = F0_eq*M[4][0] + F1_eq*M[4][1] + F2_eq*M[4][2] + F3_eq*M[4][3] + F4_eq*M[4][4] + F5_eq*M[4][5] + F6_eq*M[4][6] + F7_eq*M[4][7] + F8_eq*M[4][8] ; \
    meq5 = F0_eq*M[5][0] + F1_eq*M[5][1] + F2_eq*M[5][2] + F3_eq*M[5][3] + F4_eq*M[5][4] + F5_eq*M[5][5] + F6_eq*M[5][6] + F7_eq*M[5][7] + F8_eq*M[5][8] ; \
    meq6 = F0_eq*M[6][0] + F1_eq*M[6][1] + F2_eq*M[6][2] + F3_eq*M[6][3] + F4_eq*M[6][4] + F5_eq*M[6][5] + F6_eq*M[6][6] + F7_eq*M[6][7] + F8_eq*M[6][8] ; \
    meq7 = F0_eq*M[7][0] + F1_eq*M[7][1] + F2_eq*M[7][2] + F3_eq*M[7][3] + F4_eq*M[7][4] + F5_eq*M[7][5] + F6_eq*M[7][6] + F7_eq*M[7][7] + F8_eq*M[7][8] ; \
    meq8 = F0_eq*M[8][0] + F1_eq*M[8][1] + F2_eq*M[8][2] + F3_eq*M[8][3] + F4_eq*M[8][4] + F5_eq*M[8][5] + F6_eq*M[8][6] + F7_eq*M[8][7] + F8_eq*M[8][8] ; 

#define collision	\
    F0_in  = F0_in - M_I[0][0]*s0*(m0-meq0)- M_I[0][1]*s1*(m1-meq1)- M_I[0][2]*s2*(m2-meq2)- M_I[0][3]*s3*(m3-meq3)- M_I[0][4]*s4*(m4-meq4)- M_I[0][5]*s5*(m5-meq5)- M_I[0][6]*s6*(m6-meq6)- M_I[0][7]*s7*(m7-meq7)- M_I[0][8]*s8*(m8-meq8) ; \
    F1_in  = F1_in - M_I[1][0]*s0*(m0-meq0)- M_I[1][1]*s1*(m1-meq1)- M_I[1][2]*s2*(m2-meq2)- M_I[1][3]*s3*(m3-meq3)- M_I[1][4]*s4*(m4-meq4)- M_I[1][5]*s5*(m5-meq5)- M_I[1][6]*s6*(m6-meq6)- M_I[1][7]*s7*(m7-meq7)- M_I[1][8]*s8*(m8-meq8) + 9 * dt * Force[0] * (1/9.0);\
    F2_in  = F2_in - M_I[2][0]*s0*(m0-meq0)- M_I[2][1]*s1*(m1-meq1)- M_I[2][2]*s2*(m2-meq2)- M_I[2][3]*s3*(m3-meq3)- M_I[2][4]*s4*(m4-meq4)- M_I[2][5]*s5*(m5-meq5)- M_I[2][6]*s6*(m6-meq6)- M_I[2][7]*s7*(m7-meq7)- M_I[2][8]*s8*(m8-meq8) ; \
    F3_in  = F3_in - M_I[3][0]*s0*(m0-meq0)- M_I[3][1]*s1*(m1-meq1)- M_I[3][2]*s2*(m2-meq2)- M_I[3][3]*s3*(m3-meq3)- M_I[3][4]*s4*(m4-meq4)- M_I[3][5]*s5*(m5-meq5)- M_I[3][6]*s6*(m6-meq6)- M_I[3][7]*s7*(m7-meq7)- M_I[3][8]*s8*(m8-meq8) + 9 * dt * Force[0] * (-1/9.0);\
    F4_in  = F4_in - M_I[4][0]*s0*(m0-meq0)- M_I[4][1]*s1*(m1-meq1)- M_I[4][2]*s2*(m2-meq2)- M_I[4][3]*s3*(m3-meq3)- M_I[4][4]*s4*(m4-meq4)- M_I[4][5]*s5*(m5-meq5)- M_I[4][6]*s6*(m6-meq6)- M_I[4][7]*s7*(m7-meq7)- M_I[4][8]*s8*(m8-meq8) ; \
    F5_in  = F5_in - M_I[5][0]*s0*(m0-meq0)- M_I[5][1]*s1*(m1-meq1)- M_I[5][2]*s2*(m2-meq2)- M_I[5][3]*s3*(m3-meq3)- M_I[5][4]*s4*(m4-meq4)- M_I[5][5]*s5*(m5-meq5)- M_I[5][6]*s6*(m6-meq6)- M_I[5][7]*s7*(m7-meq7)- M_I[5][8]*s8*(m8-meq8) + 9 * dt * Force[0] * (1/36.0);\
    F6_in  = F6_in - M_I[6][0]*s0*(m0-meq0)- M_I[6][1]*s1*(m1-meq1)- M_I[6][2]*s2*(m2-meq2)- M_I[6][3]*s3*(m3-meq3)- M_I[6][4]*s4*(m4-meq4)- M_I[6][5]*s5*(m5-meq5)- M_I[6][6]*s6*(m6-meq6)- M_I[6][7]*s7*(m7-meq7)- M_I[6][8]*s8*(m8-meq8) + 9 * dt * Force[0] * (-1/36.0);\
    F7_in  = F7_in - M_I[7][0]*s0*(m0-meq0)- M_I[7][1]*s1*(m1-meq1)- M_I[7][2]*s2*(m2-meq2)- M_I[7][3]*s3*(m3-meq3)- M_I[7][4]*s4*(m4-meq4)- M_I[7][5]*s5*(m5-meq5)- M_I[7][6]*s6*(m6-meq6)- M_I[7][7]*s7*(m7-meq7)- M_I[7][8]*s8*(m8-meq8) + 9 * dt * Force[0] * (-1/36.0);\
    F8_in  = F8_in - M_I[8][0]*s0*(m0-meq0)- M_I[8][1]*s1*(m1-meq1)- M_I[8][2]*s2*(m2-meq2)- M_I[8][3]*s3*(m3-meq3)- M_I[8][4]*s4*(m4-meq4)- M_I[8][5]*s5*(m5-meq5)- M_I[8][6]*s6*(m6-meq6)- M_I[8][7]*s7*(m7-meq7)- M_I[8][8]*s8*(m8-meq8) + 9 * dt * Force[0] * (1/36.0); 
#endif