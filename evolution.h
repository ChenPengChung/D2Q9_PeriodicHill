#ifndef EVOLUTION_FILE
#define EVOLUTION_FILE

#include "interpolationHillISLBM.h"
#include "initialization.h"
#include "MRT_Process.h"
#include "MRT_Matrix.h"
#include "variables.h"
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
    double *Q1_h,        double*Q3_h,         double *Q5_h,       double*Q6_h){ //本程式碼不分主機端與裝置端變數，統一已_h結尾表示物理空間計算點變數
    
    // MRT 矩陣與鬆弛參數 (巨集展開後會宣告 M[9][9], M_I[9][9], s0~s8)
    Matrix;
    Matrix_Inverse;
    Relaxation; 
    int nface = NZ6 ;
    //1.函數內直接開始執行雙重for迴圈
for(int j = 3 ; j < NY6-3 ; j++){
    for(int k = 3 ; k < NZ6-3 ; k++){
        int idx_xi = j *NZ6 + k ;
        int idx ; //出現在 interpolationHillISLBM.h 巨集定義內的中間變數idx
         
        //宣告物理空間計算點的碰撞前插值後一般態分佈函數
        double F0_in,  F1_in,  F2_in,  F3_in,  F4_in,  F5_in,  F6_in,  F7_in,  F8_in ;   
        //MRT Variables
        double m0, m1, m2, m3, m4, m5, m6, m7, m8;
        double meq0, meq1, meq2, meq3, meq4, meq5, meq6, meq7, meq8;

        //xi方向預配置連乘權重一維連續記憶體的內插成員起始編號
        //也是xi_h[NZ6]在Z方向有設立bufferlayer沒有使用的證據
        int cell_z = k-3;
        if( k <= 6 ) cell_z = 3;
        if( k >= NZ6-7 ) cell_z = NZ6-10;

        
        //1.Interpolation and Streaming 
        F0_Intrpl7(f0_old, j, k);
        F1_Intrpl7(f1_old,j,k,j-3,cell_z,j,idx_xi,Y0_0,Y0_1,Y0_2,Y0_3,Y0_4,Y0_5,Y0_6,XiF1_0,XiF1_1,XiF1_2,XiF1_3,XiF1_4,XiF1_5,XiF1_6);
        F3_Intrpl7(f3_old,j,k,j-3,cell_z,j,idx_xi,Y2_0,Y2_1,Y2_2,Y2_3,Y2_4,Y2_5,Y2_6,XiF3_0,XiF3_1,XiF3_2,XiF3_3,XiF3_4,XiF3_5,XiF3_6);
        F2_Intrpl7(f2_old, j,k, j-3, cell_z, j, idx_xi, XiF2_0, XiF2_1, XiF2_2, XiF2_3, XiF2_4, XiF2_5, XiF2_6);
        F4_Intrpl7(f4_old, j, k, j-3, cell_z, j, idx_xi, XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6);
        Y_XI_Intrpl7(f5_old, F5_in, j, k, j-3, cell_z, j, idx_xi, Y0_0,Y0_1,Y0_2,Y0_3,Y0_4,Y0_5,Y0_6, XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6);
        Y_XI_Intrpl7(f6_old, F6_in, j, k, j-3, cell_z, j, idx_xi, Y2_0,Y2_1,Y2_2,Y2_3,Y2_4,Y2_5,Y2_6, XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6);
        Y_XI_Intrpl7(f7_old, F7_in, j, k, j-3, cell_z, j, idx_xi, Y2_0,Y2_1,Y2_2,Y2_3,Y2_4,Y2_5,Y2_6, XiF7_0, XiF7_1, XiF7_2, XiF7_3, XiF7_4, XiF7_5, XiF7_6);
        Y_XI_Intrpl7(f8_old, F8_in, j, k, j-3, cell_z, j, idx_xi, Y0_0,Y0_1,Y0_2,Y0_3,Y0_4,Y0_5,Y0_6, XiF8_0, XiF8_1, XiF8_2, XiF8_3, XiF8_4, XiF8_5, XiF8_6);
        //2.Special case of Streaming Step : Boundry Treatment
        //(1.)Halhf-way Bounce-Back Boundary Condition
        if( k == 3 ){
        F2_in = f4_old[idx_xi];
        F5_in = f7_old[idx_xi];
        F6_in = f8_old[idx_xi];
        }
        if( k == NZ6-4 ){
        F4_in = f2_old[idx_xi];
        F7_in = f5_old[idx_xi];
        F8_in = f6_old[idx_xi];
        }
        //2.BFL curvlinear boundary treatment
        //透過BFLInitialization已經寫入q值到矩陣當中。\
        //左丘邊界，更新F1
        if(IsLeftHill_Boundary_yPlus(y_global[j], z_global[j*NZ6+k])){
            double q1 = Q1_h[idx_xi] ; 
            if(q1<0.5){
                //透過內插與Streaming更新F1
                Y_XI_Intrpl7(f3_old, F1_in, j, k, j-3, cell_z, j, idx_xi, 
                YBFLF3_0,YBFLF3_1,YBFLF3_2,YBFLF3_3,YBFLF3_4,YBFLF3_5,YBFLF3_6,
                XiBFLF3_0, XiBFLF3_1, XiBFLF3_2, XiBFLF3_3, XiBFLF3_4, XiBFLF3_5, XiBFLF3_6);
            }
            if(q1>0.5){//做線性內插
                F1_in = (1.0/(2.0*q1))*f3_old[idx_xi] + ((2.0*q1-1.0)/(2.0*q1))*f1_old[idx_xi];
            }
        }
        //右丘邊界，更新F3
        if(IsRightHill_Boundary_yMinus(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F3的邊界計算點
            double q3 = Q3_h[idx_xi] ;
            if(q3<0.5){
                //透過內插與Streaming更新F3
                Y_XI_Intrpl7(f1_old, F3_in, j, k, j-3, cell_z, j, idx_xi, 
                YBFLF1_0,YBFLF1_1,YBFLF1_2,YBFLF1_3,YBFLF1_4,YBFLF1_5,YBFLF1_6,
                XiBFLF1_0, XiBFLF1_1, XiBFLF1_2, XiBFLF1_3, XiBFLF1_4, XiBFLF1_5, XiBFLF1_6);
            }
            if(q3>0.5){
                F3_in = (1.0/(2.0*q3))*f1_old[idx_xi] + ((2.0*q3-1.0)/(2.0*q3))*f3_old[idx_xi];
            }
        }
        //左丘邊界，更新F5
        if(IsLeftHill_Boundary_Diagonal45(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F5的邊界計算點
            double q5 = Q5_h[idx_xi] ;
            if(q5<0.5){
                //透過內插與Streaming更新F5
                Y_XI_Intrpl7(f7_old, F5_in, j, k, j-3, cell_z, j, idx_xi, 
                YBFLF7_0,YBFLF7_1,YBFLF7_2,YBFLF7_3,YBFLF7_4,YBFLF7_5,YBFLF7_6,
                XiBFLF7_0, XiBFLF7_1, XiBFLF7_2, XiBFLF7_3, XiBFLF7_4, XiBFLF7_5, XiBFLF7_6);
            }
            if(q5>0.5){
                F5_in = (1.0/(2.0*q5))*f7_old[idx_xi] + ((2.0*q5-1.0)/(2.0*q5))*f5_old[idx_xi];
            }
        }
        //右丘邊界，更新F6
        if(IsRightHill_Boundary_Diagonal135(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F6的邊界計算點
            double q6 = Q6_h[idx_xi] ;
            if(q6<0.5){
                //透過內插與Streaming更新F6
                Y_XI_Intrpl7(f8_old, F6_in, j, k, j-3, cell_z, j, idx_xi, 
                YBFLF8_0,YBFLF8_1,YBFLF8_2,YBFLF8_3,YBFLF8_4,YBFLF8_5,YBFLF8_6,
                XiBFLF8_0, XiBFLF8_1, XiBFLF8_2, XiBFLF8_3, XiBFLF8_4, XiBFLF8_5, XiBFLF8_6);
            }
            if(q6>0.5){
                F6_in = (1.0/(2.0*q6))*f8_old[idx_xi] + ((2.0*q6-1.0)/(2.0*q6))*f6_old[idx_xi];
            }
        }
        //3.質量修正
        F0_in = F0_in + rho_modify[0];
        //4.計算equilibirium distribution function 
        double rho_s = F0_in  + F1_in  + F2_in  + F3_in  + F4_in  + F5_in  + F6_in  + F7_in  + F8_in; 
        double v1 = (F1_in+ F5_in+ F8_in -( F3_in+F6_in+F7_in)) / rho_s ;
	    double w1 = (F2_in+ F5_in+ F6_in- ( F4_in+F7_in+F8_in)) / rho_s ;
        double udot = v1*v1 + w1*w1;
        const double F0_eq  = (4./9.)  *rho_s*(1.0-1.5*udot);
        const double F1_eq  = (1./9.)  *rho_s*(1.0+3.0*v1 +4.5*v1*v1-1.5*udot);
        const double F2_eq  = (1./9.)  *rho_s*(1.0+3.0*w1 +4.5*w1*w1-1.5*udot);
        const double F3_eq  = (1./9.)  *rho_s*(1.0-3.0*v1 +4.5*v1*v1-1.5*udot);
        const double F4_eq  = (1./9.)  *rho_s*(1.0-3.0*w1 +4.5*w1*w1-1.5*udot);
        const double F5_eq  = (1./36.) *rho_s*(1.0+3.0*( v1 +w1) +4.5*( v1 +w1)*( v1 +w1)-1.5*udot);
        const double F6_eq  = (1./36.) *rho_s*(1.0+3.0*(-v1 +w1) +4.5*(-v1 +w1)*(-v1 +w1)-1.5*udot);
        const double F7_eq  = (1./36.) *rho_s*(1.0+3.0*(-v1 -w1) +4.5*(-v1 -w1)*(-v1 -w1)-1.5*udot);
        const double F8_eq  = (1./36.) *rho_s*(1.0+3.0*( v1 -w1) +4.5*( v1 -w1)*( v1 -w1)-1.5*udot);
        //5.MRT Collision Process
        m_vector;
        meq;
        collision;

        //6.在此時間步更新物理間計算點的碰撞後插值後一般態分佈函數
        f0_new[idx_xi] = F0_in;
        f1_new[idx_xi] = F1_in;
        f2_new[idx_xi] = F2_in;
        f3_new[idx_xi] = F3_in;
        f4_new[idx_xi] = F4_in;
        f5_new[idx_xi] = F5_in;
        f6_new[idx_xi] = F6_in;
        f7_new[idx_xi] = F7_in;
        f8_new[idx_xi] = F8_in;

        //7.更新宏觀量
        rho_d[idx_xi] = rho_s;
        v[idx_xi] = v1;
        w[idx_xi] = w1;
}}}
//y方向週期邊界條件 
void 
#endif
