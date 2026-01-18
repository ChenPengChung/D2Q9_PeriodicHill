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
    //Ｙ方向預配置連乘權重
    double *Y0_0,  double *Y0_1, double *Y0_2,  double *Y0_3,  double *Y0_4,  double *Y0_5,  double *Y0_6,  
    double *Y2_0,  double *Y2_1,  double *Y2_2,  double *Y2_3,  double *Y2_4,  double *Y2_5,  double *Y2_6,
    ){

}

#endif
