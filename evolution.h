#ifndef EVOLUTION_FILE
#define EVOLUTION_FILE

#include "interpolationHillISLBM.h"
#include "MRT_Process.h"
#include "MRT_Matrix.h"
//1.物理空間計算點的平均密度場的時變量最小值
double dRhoglobal(double F1_in, double F2_in, double F3_in, double F4_in, double F5_in, double F6_in, double F7_in, double F8_in,
                  double f1_old, double f2_old, double f3_old, double f4_old, double f5_old, double f6_old, double f7_old, double f8_old){
    double rho_local;
    rho_local = F1_in + F2_in + F3_in + F4_in + F5_in + F6_in + F7_in + F8_in
                - (f1_old + f2_old + f3_old + f4_old + f5_old + f6_old + f7_old + f8_old);
    return rho_local;
}

void stream_collide(
    //f_old:上一個時間步所更新的物理空間計算點的碰撞後插值前一般態分佈函數
     double *f0_old, double *f1_old, double *f2_old, double *f3_old, double *f4_old, double *f5_old, double *f6_old, double *f7_old, double *f8_old,
    //f_new:本時間步所更新的物理空間計算點的碰撞後插值前一般態分佈函數
    double *f0_new, double *f1_new, double *f2_new, double *f3_new, double *f4_new, double *f5_new, double *f6_new, double *f7_new, double *f8_new, 
    //Ｙ方向預配置連乘權重一維連續記憶體
    double *Y0_0,  double *Y0_1, double *Y0_2,  double *Y0_3,  double *Y0_4,  double *Y0_5,  double *Y0_6,  //(處理F1,F5,F8等等y方向插值問題)
    double *Y2_0,  double *Y2_1,  double *Y2_2,  double *Y2_3,  double *Y2_4,  double *Y2_5,  double *Y2_6, //(處理F3,F6,F7等等y方向插值問題)
    //Z方向預配置連乘權重一維連續記憶體
    double* XiF1_0, double* XiF1_1, double* XiF1_2, double* XiF1_3, double* XiF1_4, double* XiF1_5, double* XiF1_6,
    double* XiF2_0, double* XiF2_1, double* XiF2_2, double* XiF2_3, double* XiF2_4, double* XiF2_5, double* XiF2_6,
    double* XiF3_0, double* XiF3_1, double* XiF3_2, double* XiF3_3, double* XiF3_4, double* XiF3_5, double* XiF3_6,
    double* XiF4_0, double* XiF4_1, double* XiF4_2, double* XiF4_3, double* XiF4_4, double* XiF4_5, double* XiF4_6,
    double* XiF5_0, double* XiF5_1, double* XiF5_2, double* XiF5_3, double* XiF5_4, double* XiF5_5, double* XiF5_6,
    double* XiF6_0, double* XiF6_1, double* XiF6_2, double* XiF6_3, double* XiF6_4, double* XiF6_5, double* XiF6_6,
    double* XiF7_0, double* XiF7_1, double* XiF7_2, double* XiF7_3, double* XiF7_4, double* XiF7_5, double* XiF7_6,
    double* XiF8_0, double* XiF8_1, double* XiF8_2, double* XiF8_3, double* XiF8_4, double* XiF8_5, double* XiF8_6,
    //BFL邊界條件(q<0.5)v下的y方向預配置連乘權重一維連續記憶體
    double* YBFLF3_0, double* YBFLF3_1, double* YBFLF3_2, double* YBFLF3_3, double* YBFLF3_4, double* YBFLF3_5, double* YBFLF3_6,
    double* YBFLF1_0, double* YBFLF1_1, double* YBFLF1_2, double* YBFLF1_3, double* YBFLF1_4, double* YBFLF1_5, double* YBFLF1_6,
    double* YBFLF7_0, double* YBFLF7_1, double* YBFLF7_2, double* YBFLF7_3, double* YBFLF7_4, double* YBFLF7_5, double* YBFLF7_6,
    double* YBFLF8_0, double* YBFLF8_1, double* YBFLF8_2, double* YBFLF8_3, double* YBFLF8_4, double* YBFLF8_5, double* YBFLF8_6,
    //BFL邊界條件(q<0.5)v下的z方向預配置連乘權重一維連續記憶體
    double* XiBFLF3_0, double* XiBFLF3_1, double* XiBFLF3_2, double* XiBFLF3_3, double* XiBFLF3_4, double* XiBFLF3_5, double* XiBFLF3_6,
    double* XiBFLF1_0, double* XiBFLF1_1, double* XiBFLF1_2, double* XiBFLF1_3, double* XiBFLF1_4, double* XiBFLF1_5, double* XiBFLF1_6,
    double* XiBFLF7_0, double* XiBFLF7_1, double* XiBFLF7_2, double* XiBFLF7_3, double* XiBFLF7_4, double* XiBFLF7_5, double* XiBFLF7_6,
    double* XiBFLF8_0, double* XiBFLF8_1, double* XiBFLF8_2, double* XiBFLF8_3, double* XiBFLF8_4, double* XiBFLF8_5, double* XiBFLF8_6,
    //宏觀參數
    double *v,           double *w,           double *rho_d,       double *Force,  double *rho_modify,
    //BFL邊界條件無因次化距離q
    double *Q3_d,        double*Q4_d,         double *Q15_d,       double*Q16_d)    
    {

}

#endif
