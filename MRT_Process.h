
#ifndef MRT_Process_FILE
#define MRT_Process_FILE
//此篇為分佈函數向量矩，作為MRT算子中的一種分佈函數
//做完Moment以後，分佈函數的矩作為一維陣列，是一向量
#define m_vector \
    m0 = F0_in*M[0][0] + F1_in*M[0][1] + F2_in*M[0][2] + F3_in*M[0][3] + F4_in*M[0][4] + F5_in*M[0][5] + F6_in*M[0][6] + F7_in*M[0][7] + F8_in*M[0][8] + F9_in*M[0][9]\
    m1 = F0_in*M[1][0] + F1_in*M[1][1] + F2_in*M[1][2] + F3_in*M[1][3] + F4_in*M[1][4] + F5_in*M[1][5] + F6_in*M[1][6] + F7_in*M[1][7] + F8_in*M[1][8] + F9_in*M[1][9]\
    m2 = F0_in*M[2][0] + F1_in*M[2][1] + F2_in*M[2][2] + F3_in*M[2][3] + F4_in*M[2][4] + F5_in*M[2][5] + F6_in*M[2][6] + F7_in*M[2][7] + F8_in*M[2][8] + F9_in*M[2][9]\
    m3 = F0_in*M[3][0] + F1_in*M[3][1] + F2_in*M[3][2] + F3_in*M[3][3] + F4_in*M[3][4] + F5_in*M[3][5] + F6_in*M[3][6] + F7_in*M[3][7] + F8_in*M[3][8] + F9_in*M[3][9]\
    m4 = F0_in*M[4][0] + F1_in*M[4][1] + F2_in*M[4][2] + F3_in*M[4][3] + F4_in*M[4][4] + F5_in*M[4][5] + F6_in*M[4][6] + F7_in*M[4][7] + F8_in*M[4][8] + F9_in*M[4][9]\
    m5 = F0_in*M[5][0] + F1_in*M[5][1] + F2_in*M[5][2] + F3_in*M[5][3] + F4_in*M[5][4] + F5_in*M[5][5] + F6_in*M[5][6] + F7_in*M[5][7] + F8_in*M[5][8] + F9_in*M[5][9]\
    m6 = F0_in*M[6][0] + F1_in*M[6][1] + F2_in*M[6][2] + F3_in*M[6][3] + F4_in*M[6][4] + F5_in*M[6][5] + F6_in*M[6][6] + F7_in*M[6][7] + F8_in*M[6][8] + F9_in*M[6][9]\
    m7 = F0_in*M[7][0] + F1_in*M[7][1] + F2_in*M[7][2] + F3_in*M[7][3] + F4_in*M[7][4] + F5_in*M[7][5] + F6_in*M[7][6] + F7_in*M[7][7] + F8_in*M[7][8] + F9_in*M[7][9]\
    m8 = F0_in*M[8][0] + F1_in*M[8][1] + F2_in*M[8][2] + F3_in*M[8][3] + F4_in*M[8][4] + F5_in*M[8][5] + F6_in*M[8][6] + F7_in*M[8][7] + F8_in*M[8][8] + F9_in*M[8][9]\
#endif