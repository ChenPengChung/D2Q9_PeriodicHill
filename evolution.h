#ifndef EVOLUTION_FILE
#define EVOLUTION_FILE
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "InterpolationHillISLBM.h"
#include "initialization.h"
#include "MRT_Process.h"
#include "MRT_Matrix.h"
#include "variables.h"
using namespace std ; 
//==========================================
//1.物理空間計算點的平均密度場的時變量最小值
//==========================================
double dRhoglobal(double F1_in, double F2_in, double F3_in, double F4_in, double F5_in, double F6_in, double F7_in, double F8_in,
                  double f1_old, double f2_old, double f3_old, double f4_old, double f5_old, double f6_old, double f7_old, double f8_old){
    double rho_local;
    rho_local = F1_in + F2_in + F3_in + F4_in + F5_in + F6_in + F7_in + F8_in
                - (f1_old + f2_old + f3_old + f4_old + f5_old + f6_old + f7_old + f8_old);
    return rho_local;
}
//==========================================
//2.更新密度修正量(各點相同) 
//輸入: double* rho_d 舊更新密度場 ; double* rho_modify :密度修正量 (各點相同)
//輸出: rho_modify[0] :新的密度修正量 (各點相同)
//==========================================
void ComputeMassCorrection(double* rho_d, double* rho_modify) {
    // 1. 計算當前全域密度總和（只計算計算區域）
    double rho_sum = 0.0;
    for(int j = 3; j < NY6-3; j++) {
        for(int k = 3; k < NZ6-3; k++) {
            int idx = j * NZ6 + k;
            rho_sum += rho_d[idx];
        }
    }
    // 2. 計算目標總密度（初始密度 × 網格數）
    double rho_initial = 1.0;
    int num_cells = (NY6 - 6) * (NZ6 - 6);
    double rho_target = rho_initial * (double)num_cells;
    // 3. 計算"場平均密修正量" 加到每一個場點的密度
    rho_modify[0] = (rho_target - rho_sum) / (double)num_cells;
}

//==========================================
//3.計算當前密度場場平均值
//輸入: double* rho_d 舊更新密度場 ; int step :時間步
//輸出: rho_avg 當前密度場場平均直
//==========================================
double CheckMassConservation(double* rho_d, int step) {
    // 1. 計算全域密度總和
    double rho_sum = 0.0;
    for(int j = 3; j < NY6-3; j++) {
        for(int k = 3; k < NZ6-3; k++) {
            int idx = j * NZ6 + k;
            rho_sum += rho_d[idx];
        }
    }
    // 2. 計算平均密度
    int num_cells = (NY6 - 6) * (NZ6 - 6);
    double rho_avg = rho_sum / (double)num_cells;
    return rho_avg;
}

void stream_collide(
    //f_old:上一個時間步所更新的物理空間計算點的碰撞後插值前一般態分佈函數
     double *f0_old, double *f1_old, double *f2_old, double *f3_old, double *f4_old, double *f5_old, double *f6_old, double *f7_old, double *f8_old,
    //f_new:本時間步所更新的物理空間計算點的碰撞後插值前一般態分佈函數
    double *f0_new, double *f1_new, double *f2_new, double *f3_new, double *f4_new, double *f5_new, double *f6_new, double *f7_new, double *f8_new, 
    //Ｙ方向預配置連乘權重一維連續記憶體
    double *Y0_0,  double *Y0_1, double *Y0_2,  //double *Y0_3,  double *Y0_4,  double *Y0_5,  double *Y0_6,  //(處理F1,F5,F8等等y方向插值問題)
    double *Y2_0,  double *Y2_1,  double *Y2_2,  //double *Y2_3,  double *Y2_4,  double *Y2_5,  double *Y2_6, //(處理F3,F6,F7等等y方向插值問題)
    //Z方向預配置連乘權重一維連續記憶體
    double* XiF1_0, double* XiF1_1, double* XiF1_2, double* XiF1_3, double* XiF1_4, double* XiF1_5, double* XiF1_6,
    double* XiF2_0, double* XiF2_1, double* XiF2_2, double* XiF2_3, double* XiF2_4, double* XiF2_5, double* XiF2_6,
    double* XiF3_0, double* XiF3_1, double* XiF3_2, double* XiF3_3, double* XiF3_4, double* XiF3_5, double* XiF3_6,
    double* XiF4_0, double* XiF4_1, double* XiF4_2, double* XiF4_3, double* XiF4_4, double* XiF4_5, double* XiF4_6,
    double* XiF5_0, double* XiF5_1, double* XiF5_2, double* XiF5_3, double* XiF5_4, double* XiF5_5, double* XiF5_6,
    double* XiF6_0, double* XiF6_1, double* XiF6_2, double* XiF6_3, double* XiF6_4, double* XiF6_5, double* XiF6_6,
    double* XiF7_0, double* XiF7_1, double* XiF7_2, double* XiF7_3, double* XiF7_4, double* XiF7_5, double* XiF7_6,
    double* XiF8_0, double* XiF8_1, double* XiF8_2, double* XiF8_3, double* XiF8_4, double* XiF8_5, double* XiF8_6,
    //宏觀參數
    double *v,           double *w,           double *rho_d,    double *rho_modify
){
    
    // MRT 矩陣與鬆弛參數 (巨集展開後會宣告 M[9][9], M_I[9][9], s0~s8)
    Matrix;
    Matrix_Inverse;
    Relaxation; 
    int nface = NZ6 ;
    //計算密度修正量(各點相同)
    ComputeMassCorrection( rho_d, rho_modify) ; //第一步:先更新密度場修正量(各點相同)(更新方法:與密度為1的場比較)
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
        //1.Interpolation and Streaming
        //使用預計算的 stencil 起點（基於來源點位置計算，與 RelationXi 一致）
        //
        int cell_z = k-3;
        if( k <= 6 )    cell_z = 3;
        if( k >= NZ6-7) cell_z = NZ6-10;
        int cell_y = j-1 ; 
        if( j <= 4 )    cell_y = 3 ;
        if( j >= NY6-5) cell_y = NY6-6 ;        
        F0_in = f0_old[idx_xi];
        
        //========== 邊界判斷 ==========
        bool is_bottom = (k <= 6);   // 擴大範圍以保護 7-point stencil
        bool is_top = (k >= NZ6-7);   // 擴大範圍以保護 7-point stencil
        bool is_left = (j <= 6);      // 擴大範圍以保護 3-point Y stencil
        bool is_right = (j >= NY6-7); // 擴大範圍以保護 3-point Y stencil
        
        // 真正的壁面點（用於 moving wall 計算）
        bool is_wall_bottom = (k == 3);
        bool is_wall_top = (k == NZ6-4);
        bool is_wall_left = (j == 3);
        bool is_wall_right = (j == NY6-4);
        
        //========== 左右邊界 (Wet Node - Non-Equilibrium Extrapolation Method) ==========
        // f_i(wall) = f_i^eq(wall) + [f_i - f_i^eq]_neighbor
        // 壁面速度 = 0，密度從鄰近點外推
        if(is_wall_left || is_wall_right){
            // 1. 取得鄰近點 (內部一點) 的資料
            int j_neighbor = is_wall_left ? 4 : (NY6-5);  // j=4 或 j=NY6-5
            int idx_neighbor = j_neighbor * NZ6 + k;
            
            // 鄰近點的分佈函數
            double f0_nb = f0_old[idx_neighbor];
            double f1_nb = f1_old[idx_neighbor];
            double f2_nb = f2_old[idx_neighbor];
            double f3_nb = f3_old[idx_neighbor];
            double f4_nb = f4_old[idx_neighbor];
            double f5_nb = f5_old[idx_neighbor];
            double f6_nb = f6_old[idx_neighbor];
            double f7_nb = f7_old[idx_neighbor];
            double f8_nb = f8_old[idx_neighbor];
            
            // 鄰近點的密度和速度
            double rho_nb = f0_nb + f1_nb + f2_nb + f3_nb + f4_nb + f5_nb + f6_nb + f7_nb + f8_nb;
            double v_nb = (f1_nb + f5_nb + f8_nb - f3_nb - f6_nb - f7_nb) / rho_nb;
            double w_nb = (f2_nb + f5_nb + f6_nb - f4_nb - f7_nb - f8_nb) / rho_nb;
            
            // 2. 計算鄰近點的平衡態
            double v_sq_nb = v_nb * v_nb + w_nb * w_nb;
            double feq0_nb = (4./9.)  * rho_nb * (1.0 - 1.5 * v_sq_nb);
            double feq1_nb = (1./9.)  * rho_nb * (1.0 + 3.0*v_nb + 4.5*v_nb*v_nb - 1.5*v_sq_nb);
            double feq2_nb = (1./9.)  * rho_nb * (1.0 + 3.0*w_nb + 4.5*w_nb*w_nb - 1.5*v_sq_nb);
            double feq3_nb = (1./9.)  * rho_nb * (1.0 - 3.0*v_nb + 4.5*v_nb*v_nb - 1.5*v_sq_nb);
            double feq4_nb = (1./9.)  * rho_nb * (1.0 - 3.0*w_nb + 4.5*w_nb*w_nb - 1.5*v_sq_nb);
            double feq5_nb = (1./36.) * rho_nb * (1.0 + 3.0*(v_nb+w_nb) + 4.5*(v_nb+w_nb)*(v_nb+w_nb) - 1.5*v_sq_nb);
            double feq6_nb = (1./36.) * rho_nb * (1.0 + 3.0*(-v_nb+w_nb) + 4.5*(-v_nb+w_nb)*(-v_nb+w_nb) - 1.5*v_sq_nb);
            double feq7_nb = (1./36.) * rho_nb * (1.0 + 3.0*(-v_nb-w_nb) + 4.5*(-v_nb-w_nb)*(-v_nb-w_nb) - 1.5*v_sq_nb);
            double feq8_nb = (1./36.) * rho_nb * (1.0 + 3.0*(v_nb-w_nb) + 4.5*(v_nb-w_nb)*(v_nb-w_nb) - 1.5*v_sq_nb);
            
            // 3. 計算鄰近點的非平衡態部分
            double fneq0 = f0_nb - feq0_nb;
            double fneq1 = f1_nb - feq1_nb;
            double fneq2 = f2_nb - feq2_nb;
            double fneq3 = f3_nb - feq3_nb;
            double fneq4 = f4_nb - feq4_nb;
            double fneq5 = f5_nb - feq5_nb;
            double fneq6 = f6_nb - feq6_nb;
            double fneq7 = f7_nb - feq7_nb;
            double fneq8 = f8_nb - feq8_nb;
            
            // 4. 壁面的密度用鄰近點密度（或可外推），壁面速度 = 0
            double rho_wall = rho_nb + rho_modify[0];  // 質量修正
            double v_wall = 0.0;
            double w_wall = 0.0;
            
            // 5. 計算壁面的平衡態（速度=0）
            double feq0_wall = (4./9.)  * rho_wall;
            double feq1_wall = (1./9.)  * rho_wall;
            double feq2_wall = (1./9.)  * rho_wall;
            double feq3_wall = (1./9.)  * rho_wall;
            double feq4_wall = (1./9.)  * rho_wall;
            double feq5_wall = (1./36.) * rho_wall;
            double feq6_wall = (1./36.) * rho_wall;
            double feq7_wall = (1./36.) * rho_wall;
            double feq8_wall = (1./36.) * rho_wall;
            
            // 6. Non-Equilibrium Extrapolation: f_wall = feq_wall + fneq_neighbor
            F0_in = feq0_wall + fneq0;
            F1_in = feq1_wall + fneq1;
            F2_in = feq2_wall + fneq2;
            F3_in = feq3_wall + fneq3;
            F4_in = feq4_wall + fneq4;
            F5_in = feq5_wall + fneq5;
            F6_in = feq6_wall + fneq6;
            F7_in = feq7_wall + fneq7;
            F8_in = feq8_wall + fneq8;
            
            // 存儲結果
            f0_new[idx_xi] = F0_in;
            f1_new[idx_xi] = F1_in;
            f2_new[idx_xi] = F2_in;
            f3_new[idx_xi] = F3_in;
            f4_new[idx_xi] = F4_in;
            f5_new[idx_xi] = F5_in;
            f6_new[idx_xi] = F6_in;
            f7_new[idx_xi] = F7_in;
            f8_new[idx_xi] = F8_in;
            
            rho_d[idx_xi] = rho_wall;
            v[idx_xi] = v_wall;
            w[idx_xi] = w_wall;
            continue;  // 跳到下一個點
        }
        
        //========== 正常 Streaming + Bounce-Back ==========
        // 注意：左右壁面 (j==3, j==NY6-4) 已用 Equilibrium Scheme 處理並 continue
        // 這裡只處理內部點和上下壁面
        
        //--- F1 (+Y方向): 來源在 j-1 ---
        // 如果 j=4，來源 j-1=3 是左壁面（Equilibrium Scheme），需要特殊處理
        if(j == 4){
            // 左壁面旁邊的點：用 bounce-back（壁面速度=0）
            F1_in = f3_old[idx_xi];
        } else if(j <= 6 || is_top || is_bottom){
            // 接近左壁面或上下壁面區域：簡單 Y 方向移動
            F1_in = f1_old[(j-1) * NZ6 + k];
        } else {
            F1_Intrpl3(f1_old,j,k,cell_y,cell_z ,j,idx_xi,Y0_0,Y0_1,Y0_2,XiF1_0,XiF1_1,XiF1_2,XiF1_3,XiF1_4,XiF1_5,XiF1_6);
        }
        
        //--- F3 (-Y方向): 來源在 j+1 ---
        // 如果 j=NY6-5，來源 j+1=NY6-4 是右壁面（Equilibrium Scheme），需要特殊處理
        if(j == NY6-5){
            // 右壁面旁邊的點：用 bounce-back（壁面速度=0）
            F3_in = f1_old[idx_xi];
        } else if(j >= NY6-7 || is_top || is_bottom){
            // 接近右壁面或上下壁面區域：簡單 Y 方向移動
            F3_in = f3_old[(j+1) * NZ6 + k];
        } else {
            F3_Intrpl3(f3_old,j,k,cell_y,cell_z ,j,idx_xi,Y2_0,Y2_1,Y2_2,XiF3_0,XiF3_1,XiF3_2,XiF3_3,XiF3_4,XiF3_5,XiF3_6);
        }
        
        //--- F2 (+Z方向): 來源在 k-1 ---
        // F2 來自下方，接近下壁面需要 bounce-back
        if(k == 4){
            // k=4 的 F2 來源在 k=3（下壁面），用靜止 bounce-back
            F2_in = f4_old[idx_xi];
        } else if(k <= 6 || is_left || is_right){
            // 接近下壁面或左右壁面區域：簡單 streaming
            F2_in = f2_old[j * NZ6 + (k-1)];
        } else {
            F2_Intrpl7(f2_old, j,k,cell_y, cell_z , j, idx_xi, XiF2_0, XiF2_1, XiF2_2, XiF2_3, XiF2_4, XiF2_5, XiF2_6);
        }
        
        //--- F4 (-Z方向): 來源在 k+1 ---
        // F4 來自上方，接近上壁面需要 bounce-back
        if(k == NZ6-5 || is_wall_top){
            // k=NZ6-5 的 F4 來源在 k=NZ6-4（上壁面）
            // F4 方向是 (0,-1)，e4·u_wall = 0，所以 moving bounce-back = 靜止 bounce-back
            F4_in = f2_old[idx_xi];
        } else if(k >= NZ6-7 || is_left || is_right){
            // 接近上壁面或左右壁面區域：簡單 streaming
            F4_in = f4_old[j * NZ6 + (k+1)];
        } else {
            F4_Intrpl7(f4_old, j, k,cell_y,cell_z , j, idx_xi, XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6);
        }
        
        //--- F5 (+Y,+Z方向): 來源在 (j-1, k-1) ---
        // F5 來自 (j-1, k-1)
        if(j == 4 || k == 4){
            // j=4: 來源 j-1=3 是左壁面，用 bounce-back
            // k=4: 來源 k-1=3 是下壁面，用 bounce-back
            F5_in = f7_old[idx_xi];
        } else if(j <= 6 || k <= 6 || k >= NZ6-7){
            // 接近左壁面或上下壁面：簡單 streaming 避免插值取到壁面
            int src_j = (j > 3) ? j - 1 : j;  // 保護不要取到左壁面
            int src_k = k - 1;
            F5_in = f5_old[src_j * NZ6 + src_k];
        } else {
            Y_XI_Intrpl3(f5_old, F5_in, j, k, cell_y , cell_z, j, idx_xi, Y0_0,Y0_1,Y0_2, XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6);
        }
        
        //--- F6 (-Y,+Z方向): 來源在 (j+1, k-1) ---
        // F6 來自 (j+1, k-1)
        if(j == NY6-5 || k == 4){
            // j=NY6-5: 來源 j+1=NY6-4 是右壁面(wet node)，用靜止 bounce-back
            // k=4: 來源 k-1=3 是下壁面，用 bounce-back
            F6_in = f8_old[idx_xi];
        } else if(j >= NY6-7 || k <= 6 || k >= NZ6-7){
            // 接近右壁面或上下壁面：簡單 streaming
            int src_j = (j < NY6-4) ? j + 1 : j;  // 保護不要取到壁面
            int src_k = k - 1;
            F6_in = f6_old[src_j * NZ6 + src_k];
        } else {
            Y_XI_Intrpl3(f6_old, F6_in, j, k, cell_y , cell_z, j, idx_xi, Y2_0,Y2_1,Y2_2, XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6);
        }
        
        //--- F7 (-Y,-Z方向): 來源在 (j+1, k+1) ---
        // F7 來自 (j+1, k+1)
        if(k == NZ6-5){
            // k=NZ6-5: F7 來源在 k+1=NZ6-4（上壁面），用 moving bounce-back
            double rho_wall = f0_old[idx_xi] + f1_old[idx_xi] + f2_old[idx_xi] + f3_old[idx_xi] + f4_old[idx_xi] 
                            + f5_old[idx_xi] + f6_old[idx_xi] + f7_old[idx_xi] + f8_old[idx_xi];
            F7_in = f5_old[idx_xi] - (1.0/6.0) * rho_wall * Uref;
        } else if(j == NY6-5){
            // j=NY6-5: 來源 j+1=NY6-4 是右壁面(wet node)，用靜止 bounce-back
            F7_in = f5_old[idx_xi];
        } else if(is_wall_top){
            // 上壁面本身：moving bounce-back
            double rho_wall = f0_old[idx_xi] + f1_old[idx_xi] + f2_old[idx_xi] + f3_old[idx_xi] + f4_old[idx_xi] 
                            + f5_old[idx_xi] + f6_old[idx_xi] + f7_old[idx_xi] + f8_old[idx_xi];
            F7_in = f5_old[idx_xi] - (1.0/6.0) * rho_wall * Uref;
        } else if(j >= NY6-7 || k >= NZ6-7){
            // 接近右壁面或上壁面：簡單 streaming
            int src_j = (j < NY6-4) ? j + 1 : j;  // 保護不要取到壁面
            int src_k = (k < NZ6-4) ? k + 1 : k;
            F7_in = f7_old[src_j * NZ6 + src_k];
        } else {
            Y_XI_Intrpl3(f7_old, F7_in, j, k, cell_y , cell_z, j, idx_xi, Y2_0,Y2_1,Y2_2, XiF7_0, XiF7_1, XiF7_2, XiF7_3, XiF7_4, XiF7_5, XiF7_6);
        }
        
        //--- F8 (+Y,-Z方向): 來源在 (j-1, k+1) ---
        // F8 來自 (j-1, k+1)
        if(k == NZ6-5){
            // k=NZ6-5: F8 來源在 k+1=NZ6-4（上壁面），用 moving bounce-back
            double rho_wall = f0_old[idx_xi] + f1_old[idx_xi] + f2_old[idx_xi] + f3_old[idx_xi] + f4_old[idx_xi] 
                            + f5_old[idx_xi] + f6_old[idx_xi] + f7_old[idx_xi] + f8_old[idx_xi];
            F8_in = f6_old[idx_xi] + (1.0/6.0) * rho_wall * Uref;
        } else if(j == 4){
            // j=4: 來源 j-1=3 是左壁面(wet node)，用靜止 bounce-back
            F8_in = f6_old[idx_xi];
        } else if(is_wall_top){
            // 上壁面本身：moving bounce-back
            double rho_wall = f0_old[idx_xi] + f1_old[idx_xi] + f2_old[idx_xi] + f3_old[idx_xi] + f4_old[idx_xi] 
                            + f5_old[idx_xi] + f6_old[idx_xi] + f7_old[idx_xi] + f8_old[idx_xi];
            F8_in = f6_old[idx_xi] + (1.0/6.0) * rho_wall * Uref;
        } else if(j <= 6 || k >= NZ6-7){
            // 接近左壁面或上壁面：簡單 streaming
            int src_j = (j > 3) ? j - 1 : j;  // 保護不要取到左壁面
            int src_k = (k < NZ6-4) ? k + 1 : k;  // 保護不要取到上壁面
            F8_in = f8_old[src_j * NZ6 + src_k];
        } else {
            Y_XI_Intrpl3(f8_old, F8_in, j, k, cell_y , cell_z, j, idx_xi, Y0_0,Y0_1,Y0_2, XiF8_0, XiF8_1, XiF8_2, XiF8_3, XiF8_4, XiF8_5, XiF8_6);
        }
        
        
    
        //3.質量修正
        F0_in = F0_in + rho_modify[0];
        
        //4.計算密度和速度
        double rho_s = F0_in  + F1_in  + F2_in  + F3_in  + F4_in  + F5_in  + F6_in  + F7_in  + F8_in; 
        
        // 除錯：在計算 rho_s 後立即檢查
        if(rho_s < 0.5 || rho_s > 2.0 || std::isnan(rho_s)) {
            static int early_error = 0;
            if(early_error < 5) {
                std::cout << "EARLY_ABNORMAL at j=" << j << " k=" << k 
                          << " rho_s=" << rho_s
                          << " is_left=" << is_left << " is_right=" << is_right
                          << " is_top=" << is_top << " is_bottom=" << is_bottom
                          << "\n  F0=" << F0_in << " F1=" << F1_in << " F2=" << F2_in 
                          << " F3=" << F3_in << " F4=" << F4_in
                          << " F5=" << F5_in << " F6=" << F6_in 
                          << " F7=" << F7_in << " F8=" << F8_in << std::endl;
                early_error++;
            }
        }
        
        double v1 = (F1_in+ F5_in+ F8_in -( F3_in+F6_in+F7_in)) / rho_s ;
	    double w1 = (F2_in+ F5_in+ F6_in -( F4_in+F7_in+F8_in)) / rho_s ;
        
        //計算平衡態分佈函數
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
        
        //5.MRT Collision Process (內部點和上下壁面)
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
        
        //邊界速度強制為壁面速度
        if(is_bottom){
            v[idx_xi] = 0.0;
            w[idx_xi] = 0.0;
        }
        if(is_top){
            v[idx_xi] = Uref;
            w[idx_xi] = 0.0;
        }

        // 除錯：檢查是否有異常值（暫時放寬閾值以觀察趨勢）
        double rho_local = F0_in + F1_in + F2_in + F3_in + F4_in + F5_in + F6_in + F7_in + F8_in;
        if(rho_local > 5.0 || rho_local < 0.0 || std::isnan(rho_local)) {
            static int error_count = 0;
            if(error_count < 10) {
                std::cout << "ABNORMAL at j=" << j << " k=" << k 
                          << " rho=" << rho_local
                          << " F0=" << F0_in << " F1=" << F1_in << " F2=" << F2_in 
                          << " F3=" << F3_in << " F4=" << F4_in
                          << " F5=" << F5_in << " F6=" << F6_in 
                          << " F7=" << F7_in << " F8=" << F8_in << std::endl;
                error_count++;
            }
        }
}}}


#endif
