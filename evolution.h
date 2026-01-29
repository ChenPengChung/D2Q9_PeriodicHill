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
// Mach 數監控與速度限制函數
//==========================================

// 檢查並限制單點速度，返回是否進行了限制
inline bool LimitVelocity(double& uy, double& uz, double rho_local) {
    double u_mag = std::sqrt(uy*uy + uz*uz);
    if(u_mag > U_max) {
        // 按比例縮放速度到最大允許值
        double scale = U_max / u_mag;
        uy *= scale;
        uz *= scale;
        return true;  // 進行了限制
    }
    return false;
}

// 重建平衡態分佈函數（當速度被限制後需要重建）
inline void RebuildEquilibrium(
    double& F0, double& F1, double& F2, double& F3, double& F4,
    double& F5, double& F6, double& F7, double& F8,
    double rho, double uy, double uz)
{
    double udot = uy*uy + uz*uz;
    F0 = (4.0/9.0)  * rho * (1.0 - 1.5*udot);
    F1 = (1.0/9.0)  * rho * (1.0 + 3.0*uy + 4.5*uy*uy - 1.5*udot);
    F2 = (1.0/9.0)  * rho * (1.0 + 3.0*uz + 4.5*uz*uz - 1.5*udot);
    F3 = (1.0/9.0)  * rho * (1.0 - 3.0*uy + 4.5*uy*uy - 1.5*udot);
    F4 = (1.0/9.0)  * rho * (1.0 - 3.0*uz + 4.5*uz*uz - 1.5*udot);
    F5 = (1.0/36.0) * rho * (1.0 + 3.0*(uy+uz) + 4.5*(uy+uz)*(uy+uz) - 1.5*udot);
    F6 = (1.0/36.0) * rho * (1.0 + 3.0*(-uy+uz) + 4.5*(-uy+uz)*(-uy+uz) - 1.5*udot);
    F7 = (1.0/36.0) * rho * (1.0 + 3.0*(-uy-uz) + 4.5*(-uy-uz)*(-uy-uz) - 1.5*udot);
    F8 = (1.0/36.0) * rho * (1.0 + 3.0*(uy-uz) + 4.5*(uy-uz)*(uy-uz) - 1.5*udot);
}

// 全域 Mach 數統計結構
struct MachStats {
    double max_Ma;           // 最大 Mach 數
    double avg_Ma;           // 平均 Mach 數
    int count_exceed_warning; // 超過警告閾值的點數
    int count_exceed_max;     // 超過最大限制的點數
    int j_max, k_max;         // 最大 Mach 數位置
};

// 計算全場 Mach 數統計（每隔一定步數呼叫）
MachStats ComputeMachStats(double* v_field, double* w_field, double* rho_field) {
    MachStats stats = {0.0, 0.0, 0, 0, 0, 0};
    double Ma_sum = 0.0;
    int count = 0;

    for(int j = 3; j < NY6-3; j++) {
        for(int k = 3; k < NZ6-3; k++) {
            int idx = j * NZ6 + k;
            double u_mag = std::sqrt(v_field[idx]*v_field[idx] + w_field[idx]*w_field[idx]);
            double Ma_local = u_mag / cs;

            Ma_sum += Ma_local;
            count++;

            if(Ma_local > stats.max_Ma) {
                stats.max_Ma = Ma_local;
                stats.j_max = j;
                stats.k_max = k;
            }

            if(Ma_local > Ma_warning) stats.count_exceed_warning++;
            if(Ma_local > Ma_max) stats.count_exceed_max++;
        }
    }

    stats.avg_Ma = Ma_sum / (double)count;
    return stats;
}

// 輸出 Mach 數診斷資訊
void PrintMachDiagnostics(int timestep, const MachStats& stats) {
    std::cout << "[t=" << timestep << "] Mach stats: "
              << "max=" << stats.max_Ma << " (j=" << stats.j_max << ",k=" << stats.k_max << "), "
              << "avg=" << stats.avg_Ma;

    if(stats.count_exceed_warning > 0) {
        std::cout << ", WARNING: " << stats.count_exceed_warning << " points > Ma_warning";
    }
    if(stats.count_exceed_max > 0) {
        std::cout << ", CRITICAL: " << stats.count_exceed_max << " points > Ma_max (limited)";
    }
    std::cout << std::endl;
}

// 啟動時檢查理論 Mach 數
void CheckTheoreticalMach() {
    double Ma_theoretical = Uref / cs;
    std::cout << "=======================================" << std::endl;
    std::cout << "Mach Number Check at Startup" << std::endl;
    std::cout << "=======================================" << std::endl;
    std::cout << "Re = " << Re << std::endl;
    std::cout << "tau = " << tau << std::endl;
    std::cout << "niu = " << niu << std::endl;
    std::cout << "Uref = " << Uref << std::endl;
    std::cout << "cs = " << cs << std::endl;
    std::cout << "Ma_theoretical = Uref/cs = " << Ma_theoretical << std::endl;
    std::cout << "Ma_max (allowed) = " << Ma_max << std::endl;
    std::cout << "Ma_warning = " << Ma_warning << std::endl;

    if(Ma_theoretical > Ma_max) {
        std::cout << "!!! CRITICAL WARNING !!!" << std::endl;
        std::cout << "Theoretical Mach number (" << Ma_theoretical
                  << ") exceeds Ma_max (" << Ma_max << ")!" << std::endl;
        std::cout << "Simulation is likely to diverge!" << std::endl;
        std::cout << "Suggestions:" << std::endl;
        std::cout << "  1. Reduce CFL (current: " << CFL << ")" << std::endl;
        std::cout << "  2. Increase grid resolution (current: NY=" << NY << ", NZ=" << NZ << ")" << std::endl;
        std::cout << "  3. Reduce tau (current: " << tau << ") - but keep tau > 0.5" << std::endl;
    } else if(Ma_theoretical > Ma_warning) {
        std::cout << "WARNING: Theoretical Mach number is high but within limit." << std::endl;
        std::cout << "Consider reducing CFL or increasing resolution for better accuracy." << std::endl;
    } else {
        std::cout << "OK: Theoretical Mach number is within safe range." << std::endl;
    }
    std::cout << "=======================================" << std::endl;
}

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

//==========================================
//4.計算區域質量守恆 - 質量通量淨流入流出差
// 對每個格點計算: |流出質量通量 - 流入質量通量|
// 質量通量 = ρ*v (Y方向) + ρ*w (Z方向)
// 使用中心差分: div(ρu) ≈ [(ρv)_{j+1} - (ρv)_{j-1}]/2 + [(ρw)_{k+1} - (ρw)_{k-1}]/2
//==========================================
double ComputeMassFluxError(double* rho_d, double* v_field, double* w_field) {
    double total_error = 0.0;
    
    for(int j = 4; j < NY6-4; j++) {
        for(int k = 4; k < NZ6-4; k++) {
            int idx = j * NZ6 + k;
            int idx_jp1 = (j+1) * NZ6 + k;  // j+1
            int idx_jm1 = (j-1) * NZ6 + k;  // j-1
            int idx_kp1 = j * NZ6 + (k+1);  // k+1
            int idx_km1 = j * NZ6 + (k-1);  // k-1
            
            // Y方向質量通量差 (ρv)_{j+1} - (ρv)_{j-1}
            double flux_y = (rho_d[idx_jp1] * v_field[idx_jp1]) - (rho_d[idx_jm1] * v_field[idx_jm1]);
            
            // Z方向質量通量差 (ρw)_{k+1} - (ρw)_{k-1}
            double flux_z = (rho_d[idx_kp1] * w_field[idx_kp1]) - (rho_d[idx_km1] * w_field[idx_km1]);
            
            // 質量守恆誤差 = |div(ρu)| = |∂(ρv)/∂y + ∂(ρw)/∂z|
            double div_rho_u = std::fabs(flux_y + flux_z);
            total_error += div_rho_u;
        }
    }
    
    // 回傳總質量通量誤差（除以格點數得到平均值）
    int num_cells = (NY6 - 8) * (NZ6 - 8);
    return total_error / (double)num_cells;
}

//y方向週期邊界條件 
void periodicSW(
    double* f0_new, double* f1_new, double* f2_new, double* f3_new, double* f4_new, double* f5_new, double* f6_new, double* f7_new, double* f8_new
    ,double* v, double* w, double* rho_d
){//SW:Stream-Wise
    //目的：複製 Y 方向的 buffer Layer，實現週期性邊界
    //只複製 Z 方向的有效計算區域 (k=3 ~ NZ6-4)
    for(int k = 3 ; k < NZ6-3 ; k++){
        for(int i = 0 ; i <= 2 ; i++){
            //右邊左側buffer layer
            int idx_right = (i+NY6-7)*NZ6 + k ;
            int buffer_left = i*NZ6 + k ;
            f0_new[buffer_left] = f0_new[idx_right] ;
            f1_new[buffer_left] = f1_new[idx_right] ;
            f2_new[buffer_left] = f2_new[idx_right] ;
            f3_new[buffer_left] = f3_new[idx_right] ;
            f4_new[buffer_left] = f4_new[idx_right] ;
            f5_new[buffer_left] = f5_new[idx_right] ;
            f6_new[buffer_left] = f6_new[idx_right] ;
            f7_new[buffer_left] = f7_new[idx_right] ;
            f8_new[buffer_left] = f8_new[idx_right] ;
            v[buffer_left] = v[idx_right] ;
            w[buffer_left] = w[idx_right] ;
            rho_d[buffer_left] = rho_d[idx_right] ;
            //左邊右側buffer layer
            int idx_left = (i+4)*NZ6 + k ;
            int buffer_right = (NY6-3+i)*NZ6 + k ;
            f0_new[buffer_right] = f0_new[idx_left] ;
            f1_new[buffer_right] = f1_new[idx_left] ;
            f2_new[buffer_right] = f2_new[idx_left] ;
            f3_new[buffer_right] = f3_new[idx_left] ;
            f4_new[buffer_right] = f4_new[idx_left] ;
            f5_new[buffer_right] = f5_new[idx_left] ;
            f6_new[buffer_right] = f6_new[idx_left] ;
            f7_new[buffer_right] = f7_new[idx_left] ;
            f8_new[buffer_right] = f8_new[idx_left] ;
            v[buffer_right] = v[idx_left] ;
            w[buffer_right] = w[idx_left] ;
            rho_d[buffer_right] = rho_d[idx_left] ;
        }
    }
}


void stream_collide(
    //f_old:上一個時間步所更新的物理空間計算點的碰撞後插值前一般態分佈函數
     double *f0_old, double *f1_old, double *f2_old, double *f3_old, double *f4_old, double *f5_old, double *f6_old, double *f7_old, double *f8_old,
    //f_new:本時間步所更新的物理空間計算點的碰撞後插值前一般態分佈函數
    double *f0_new, double *f1_new, double *f2_new, double *f3_new, double *f4_new, double *f5_new, double *f6_new, double *f7_new, double *f8_new, 
    //Ｙ方向預配置連乘權重一維連續記憶體
    double *Y0_0,  double *Y0_1, double *Y0_2,  //double *Y0_3,  double *Y0_4,  double *Y0_5,  double *Y0_6,  //(處理F1,F5,F8等等y方向插值問題)
    double *Y2_0,  double *Y2_1, double *Y2_2,  //double *Y2_3,  double *Y2_4,  double *Y2_5,  double *Y2_6, //(處理F3,F6,F7等等y方向插值問題)
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
    double* YBFLF3_0, double* YBFLF3_1, double* YBFLF3_2, //double* YBFLF3_3, double* YBFLF3_4, double* YBFLF3_5, double* YBFLF3_6,
    double* YBFLF1_0, double* YBFLF1_1, double* YBFLF1_2, //double* YBFLF1_3, double* YBFLF1_4, double* YBFLF1_5, double* YBFLF1_6,
    double* YBFLF7_0, double* YBFLF7_1, double* YBFLF7_2, //double* YBFLF7_3, double* YBFLF7_4, double* YBFLF7_5, double* YBFLF7_6,
    double* YBFLF8_0, double* YBFLF8_1, double* YBFLF8_2, //double* YBFLF8_3, double* YBFLF8_4, double* YBFLF8_5, double* YBFLF8_6,
    //BFL邊界條件(q<0.5)v下的z方向預配置連乘權重一維連續記憶體
    double* XiBFLF3_0, double* XiBFLF3_1, double* XiBFLF3_2, double* XiBFLF3_3, double* XiBFLF3_4, double* XiBFLF3_5, double* XiBFLF3_6,
    double* XiBFLF1_0, double* XiBFLF1_1, double* XiBFLF1_2, double* XiBFLF1_3, double* XiBFLF1_4, double* XiBFLF1_5, double* XiBFLF1_6,
    double* XiBFLF7_0, double* XiBFLF7_1, double* XiBFLF7_2, double* XiBFLF7_3, double* XiBFLF7_4, double* XiBFLF7_5, double* XiBFLF7_6,
    double* XiBFLF8_0, double* XiBFLF8_1, double* XiBFLF8_2, double* XiBFLF8_3, double* XiBFLF8_4, double* XiBFLF8_5, double* XiBFLF8_6,
    //宏觀參數
    double *v,           double *w,           double *rho_d,       double *Force,  double *rho_modify,
    //BFL邊界條件無因次化距離q
    double *Q1_h,        double*Q3_h,         double *Q5_h,       double*Q6_h
){ //本程式碼不分主機端與裝置端變數，統一已_h結尾表示物理空間計算點變數
    
    // MRT 矩陣與鬆弛參數 (巨集展開後會宣告 M[9][9], M_I[9][9], s0~s8)
    Matrix;
    Matrix_Inverse;
    Relaxation; 
    int nface = NZ6 ;
    //計算密度修正量(各點相同)
    ComputeMassCorrection( rho_d, rho_modify) ; //第一步:先更新密度場修正量(各點相同)(更新方法:與密度為1的場比較)
    //施加週期性邊界條件
    periodicSW(f0_old, f1_old, f2_old, f3_old, f4_old, f5_old, f6_old, f7_old, f8_old, v, w, rho_d);//因為用舊值做操作 
    //1.函數內直接開始執行雙重for迴圈
for(int j = 3 ; j < NY6-3 ; j++){
    for(int k = 3 ; k < NZ6-3 ; k++){
        int idx_xi = j *NZ6 + k ;
        int idx ; //出現在 interpolationHillISLBM.h 巨集定義內的中間變數idx
        //宣告物理空間計算點的碰撞前插值後一般態分佈函數
        double F0_in = 0.0, F1_in = 0.0, F2_in = 0.0, F3_in = 0.0, F4_in = 0.0, F5_in = 0.0, F6_in = 0.0, F7_in = 0.0, F8_in = 0.0;   
        //MRT Variables
        double m0, m1, m2, m3, m4, m5, m6, m7, m8;
        double meq0, meq1, meq2, meq3, meq4, meq5, meq6, meq7, meq8;


        //1.Interpolation and Streaming
        //使用預計算的 stencil 起點（基於來源點位置計算，與 RelationXi 一致）
        //
        int cell_z = k-3;
        if( k <= 6 )    cell_z = 3;
        if( k >= NZ6-7) cell_z = NZ6-10;
        /*F0_Intrpl7(f0_old, j, k);
        F1_Intrpl7(f1_old,j,k,j-3,cell_z,j,idx_xi,Y0_0,Y0_1,Y0_2,Y0_3,Y0_4,Y0_5,Y0_6,XiF1_0,XiF1_1,XiF1_2,XiF1_3,XiF1_4,XiF1_5,XiF1_6);
        F3_Intrpl7(f3_old,j,k,j-3,cell_z,j,idx_xi,Y2_0,Y2_1,Y2_2,Y2_3,Y2_4,Y2_5,Y2_6,XiF3_0,XiF3_1,XiF3_2,XiF3_3,XiF3_4,XiF3_5,XiF3_6);
        F2_Intrpl7(f2_old, j,k, j-3, cell_z, j, idx_xi, XiF2_0, XiF2_1, XiF2_2, XiF2_3, XiF2_4, XiF2_5, XiF2_6);
        F4_Intrpl7(f4_old, j, k, j-3, cell_z, j, idx_xi, XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6);
        Y_XI_Intrpl7(f5_old, F5_in, j, k, j-3, cell_z, j, idx_xi, Y0_0,Y0_1,Y0_2,Y0_3,Y0_4,Y0_5,Y0_6, XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6);
        Y_XI_Intrpl7(f6_old, F6_in, j, k, j-3, cell_z, j, idx_xi, Y2_0,Y2_1,Y2_2,Y2_3,Y2_4,Y2_5,Y2_6, XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6);
        Y_XI_Intrpl7(f7_old, F7_in, j, k, j-3, cell_z, j, idx_xi, Y2_0,Y2_1,Y2_2,Y2_3,Y2_4,Y2_5,Y2_6, XiF7_0, XiF7_1, XiF7_2, XiF7_3, XiF7_4, XiF7_5, XiF7_6);
        Y_XI_Intrpl7(f8_old, F8_in, j, k, j-3, cell_z, j, idx_xi, Y0_0,Y0_1,Y0_2,Y0_3,Y0_4,Y0_5,Y0_6, XiF8_0, XiF8_1, XiF8_2, XiF8_3, XiF8_4, XiF8_5, XiF8_6);*/
        
        //降階版本
        F0_in = f0_old[idx_xi];
        
        // 對 Z 方向邊界使用條件式處理：
        // 7-point Xi stencil 需要 k_start 到 k_start+6 範圍內的有效資料
        // 當 cellZ_F* + 6 >= NZ6 - 3 時，會存取到上邊界 buffer 區
        // 當 cellZ_F* < 3 時，會存取到下邊界 buffer 區
        // 使用 half-way bounce-back 代替插值以避免外推問題
        
        // 下邊界 (k <= 5): F2,F5,F6 的來源位置可能落在下邊界 buffer 區
        // 上邊界 (k >= NZ6-7=149): F4,F7,F8 的來源位置可能落在上邊界 buffer 區
        // 注意：F1,F3 只有 Y 方向偏移，Z 方向不變，但 stencil 仍可能越界
        
        // Y 方向邊界檢查（週期性邊界需要用 streaming，不用插值）
        // 擴大 Y 邊界區域：包含 j <= 5 和 j >= NY6-6，以避免插值存取到邊界異常值
        bool y_boundary = (j <= 3) || (j >= NY6-4);
        bool z_lower = (k <= streaming_lower);  // 恢復到合理的插值範圍 (原為 11,保守用 8)
        bool z_upper = (k >= streaming_upper);  // k >= 121
        
        if( z_lower || z_upper || y_boundary ) {
            // 邊界附近：使用簡單的 streaming 替代插值

            // F1,F3: Y方向 streaming，需要週期性 wrap
            int jm1 = (j > 3) ? (j-1) : (NY6-5) ;  // 週期性 wrap
            int jp1 = (j < NY6-4) ? (j+1) : 4 ;     // 週期性 wrap

            // F1,F3: 所有邊界區域都需要處理（Y方向 streaming）
            F1_in = f1_old[jm1*NZ6 + k];
            F3_in = f3_old[jp1*NZ6 + k];

            // 根據 Z 位置處理 F2,F4,F5,F6,F7,F8
            if( z_lower ) {
                // 下邊界：F2,F5,F6 用 bounce-back
                F2_in = f4_old[idx_xi];
                F5_in = f7_old[idx_xi];
                F6_in = f8_old[idx_xi];

                // F4,F7,F8: 向 -Z 方向，從 z+Δ 取值
                F4_in = f4_old[j*NZ6 + k+1];
                F7_in = f7_old[jp1*NZ6 + k+1];
                F8_in = f8_old[jm1*NZ6 + k+1];

            } else if( z_upper ) {
                // 上邊界：F4,F7,F8 用 bounce-back
                F4_in = f2_old[idx_xi];
                F7_in = f5_old[idx_xi];
                F8_in = f6_old[idx_xi];

                // F2,F5,F6: 向 +Z 方向，從 z-Δ 取值
                F2_in = f2_old[j*NZ6 + k-1];
                F5_in = f5_old[jm1*NZ6 + k-1];
                F6_in = f6_old[jp1*NZ6 + k-1];

            } else {
                // Y邊界但非Z邊界：使用簡單 streaming
                //此為週期性邊界條件
                F2_in = f2_old[j*NZ6 + k-1];
                F4_in = f4_old[j*NZ6 + k+1];
                F5_in = f5_old[jm1*NZ6 + k-1];
                F6_in = f6_old[jp1*NZ6 + k-1];
                F7_in = f7_old[jp1*NZ6 + k+1];
                F8_in = f8_old[jm1*NZ6 + k+1];
            }
            
        } else {
            // 內部區域：正常插值
            F1_Intrpl3(f1_old,j,k,CellZ_F1,j,idx_xi,Y0_0,Y0_1,Y0_2,XiF1_0,XiF1_1,XiF1_2,XiF1_3,XiF1_4,XiF1_5,XiF1_6);
            F3_Intrpl3(f3_old,j,k,CellZ_F1,j,idx_xi,Y2_0,Y2_1,Y2_2,XiF3_0,XiF3_1,XiF3_2,XiF3_3,XiF3_4,XiF3_5,XiF3_6);
            F2_Intrpl7(f2_old, j, k, CellZ_F2, j, idx_xi, XiF2_0, XiF2_1, XiF2_2, XiF2_3, XiF2_4, XiF2_5, XiF2_6);
            F4_Intrpl7(f4_old, j, k, CellZ_F4, j, idx_xi, XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6);
            Y_XI_Intrpl3(f5_old, F5_in, j, k, CellZ_F2, j, idx_xi, Y0_0,Y0_1,Y0_2, XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6);
            Y_XI_Intrpl3(f6_old, F6_in, j, k, CellZ_F2, j, idx_xi, Y2_0,Y2_1,Y2_2, XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6);
            Y_XI_Intrpl3(f7_old, F7_in, j, k, CellZ_F4, j, idx_xi, Y2_0,Y2_1,Y2_2, XiF7_0, XiF7_1, XiF7_2, XiF7_3, XiF7_4, XiF7_5, XiF7_6);
            Y_XI_Intrpl3(f8_old, F8_in, j, k, CellZ_F4, j, idx_xi, Y0_0,Y0_1,Y0_2, XiF8_0, XiF8_1, XiF8_2, XiF8_3, XiF8_4, XiF8_5, XiF8_6);
        }

        //2.Special case of Streaming Step : Boundry Treatment
        //(1.)Halhf-way Bounce-Back Boundary Condition
        // 注意：對於 k==3 和 k==NZ6-4，已在上面的插值條件式中處理
        // 這裡保留對 k==3 的額外確認（最下邊界層）
        if( k == 3 ){
            F2_in = f4_old[idx_xi];
            F5_in = f7_old[idx_xi];
            F6_in = f8_old[idx_xi];
        }
        // 對於 k==NZ6-4 (k==152)，已在 k >= NZ6-6 條件中處理
        
        //2.BFL curvlinear boundary treatment
        //透過BFLInitialization已經寫入q值到矩陣當中
        //BFL 用的 stencil 起點（基於目標點位置，因為 BFL 權重是獨立計算的）
        int cell_z_bfl = k-3;
        if( k <= 3 ) cell_z_bfl = 3;
        if( k >= NZ6-7 ) cell_z_bfl = NZ6-10;

        //左丘邊界，更新F1
        if(IsLeftHill_Boundary_yPlus(y_global[j], z_global[j*NZ6+k])){
            double q1 = Q1_h[idx_xi] ;
            if(q1<0.5 && q1 >= 0 ){
                /*//透過內插與Streaming更新F1
                Y_XI_Intrpl7(f3_old, F1_in, j, k, j-3, cell_z_bfl, j, idx_xi,
                YBFLF3_0,YBFLF3_1,YBFLF3_2,YBFLF3_3,YBFLF3_4,YBFLF3_5,YBFLF3_6,
                XiBFLF3_0, XiBFLF3_1, XiBFLF3_2, XiBFLF3_3, XiBFLF3_4, XiBFLF3_5, XiBFLF3_6);*/
                //改用一階精度線性內插
                F1_in = (2*q1)*f3_old[idx_xi] + (1.0 - 2.0*q1)*f3_old[idx_xi+NZ6];
                //這邊有一個問題：書上做插值的點為Lattice上的計算點但是這邊取的是物理空間計算點的點 
            }
            if(q1>0.5){//做線性內插
                F1_in = (1.0/(2.0*q1))*f3_old[idx_xi] + ((2.0*q1-1.0)/(2.0*q1))*f1_old[idx_xi];
                //問題一:遷移後在鄰點上的值不存在
                //問題二：若插分的點要用鄰點新值，會有2個問題，
                //第一，該點不在物理間網格點上，只在Lattice上，
                //第二：就算你用的到Lattice上的值，本程式碼與沒有對他做紀錄，只有暫態變數記錄且該值也還沒更新，因為掃描方向從左至右 .....
            }
        }
        //右丘邊界，更新F3
        if(IsRightHill_Boundary_yMinus(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F3的邊界計算點
            double q3 = Q3_h[idx_xi] ;
            if(q3<0.5 && q3 >= 0.0){
                /*//透過內插與Streaming更新F3
                Y_XI_Intrpl7(f1_old, F3_in, j, k, j-3, cell_z_bfl, j, idx_xi,
                YBFLF1_0,YBFLF1_1,YBFLF1_2,YBFLF1_3,YBFLF1_4,YBFLF1_5,YBFLF1_6,
                XiBFLF1_0, XiBFLF1_1, XiBFLF1_2, XiBFLF1_3, XiBFLF1_4, XiBFLF1_5, XiBFLF1_6);*/
                 //改用一階精度線性內插
                F3_in = (2*q3)*f1_old[idx_xi] + (1.0 - 2.0*q3)*f1_old[idx_xi-NZ6];
            }
            if(q3>0.5 ){
                F3_in = (1.0/(2.0*q3))*f1_old[idx_xi] + ((2.0*q3-1.0)/(2.0*q3))*f3_old[idx_xi];
            }
        }
        //左丘邊界，更新F5
        if(IsLeftHill_Boundary_Diagonal45(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F5的邊界計算點
            double q5 = Q5_h[idx_xi] ;
            if(q5<0.5 && q5 >= 0.0){
                /*//透過內插與Streaming更新F5
                Y_XI_Intrpl7(f7_old, F5_in, j, k, j-3, cell_z_bfl, j, idx_xi,
                YBFLF7_0,YBFLF7_1,YBFLF7_2,YBFLF7_3,YBFLF7_4,YBFLF7_5,YBFLF7_6,
                XiBFLF7_0, XiBFLF7_1, XiBFLF7_2, XiBFLF7_3, XiBFLF7_4, XiBFLF7_5, XiBFLF7_6);*/
                //取右上左下線性內插
                F5_in = (2*q5)*f7_old[idx_xi] + (1.0 - 2.0*q5)*f7_old[idx_xi+NZ6+1];//往右邊再往上
            }
            if(q5>0.5){
                F5_in = (1.0/(2.0*q5))*f7_old[idx_xi] + ((2.0*q5-1.0)/(2.0*q5))*f5_old[idx_xi];
            }
        }
        //右丘邊界，更新F6
        if(IsRightHill_Boundary_Diagonal135(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F6的邊界計算點
            double q6 = Q6_h[idx_xi] ;
            if(q6<0.5 && q6 >= 0.0){
                /*//透過內插與Streaming更新F6
                Y_XI_Intrpl7(f8_old, F6_in, j, k, j-1, cell_z_bfl, j, idx_xi,
                YBFLF8_0,YBFLF8_1,YBFLF8_2,YBFLF8_3,YBFLF8_4,YBFLF8_5,YBFLF8_6,
                XiBFLF8_0, XiBFLF8_1, XiBFLF8_2, XiBFLF8_3, XiBFLF8_4, XiBFLF8_5, XiBFLF8_6);*/
                //取左上右下線性內插
                //F6_in = (2*q6)*f8_old[idx_xi] + (1.0 - 2.0*q6)*f8_old[idx_xi-NZ6+1];//往左邊再往上
                Y_XI_Intrpl3(f8_old, F6_in, j, k, CellZ_F8, j, idx_xi, 
                YBFLF8_0,YBFLF8_1,YBFLF8_2 , 
                XiBFLF8_0, XiBFLF8_1, XiBFLF8_2, XiBFLF8_3, XiBFLF8_4, XiBFLF8_5, XiBFLF8_6);
            }
            if(q6>0.5){
                F6_in = (1.0/(2.0*q6))*f8_old[idx_xi] + ((2.0*q6-1.0)/(2.0*q6))*f6_old[idx_xi];
            }
        }

        //3.質量修正
        F0_in = F0_in + rho_modify[0];
        //診斷代碼：在 t=150 附近檢查異常值
        double rho_local = F0_in + F1_in + F2_in + F3_in + F4_in + F5_in + F6_in + F7_in + F8_in;
        if(rho_local > 2.0 || rho_local < 0.5 || std::isnan(rho_local)) {
            static int error_count = 0;
            if(error_count < 20) {
                bool in_lower = (k <= 15);
                bool in_upper = (k >= NZ6-35);
                bool in_ybnd = (j <= 4) || (j >= NY6-5);
                std::cout << "ABNORMAL at j=" << j << " k=" << k 
                          << " region=" << (in_lower ? "[bottom boundary region]" : (in_upper ? "[upper boundary region]" : (in_ybnd ? "[left and right boundary region]" : "[internal region]")))
                          << " rho=" << rho_local
                          << " F=[" << F0_in << "," << F1_in << "," << F2_in << "," << F3_in << "," << F4_in
                          << "," << F5_in << "," << F6_in << "," << F7_in << "," << F8_in << "]" << std::endl;
                error_count++;
            }
        }
        //4.計算equilibirium distribution function
        double rho_s = F0_in  + F1_in  + F2_in  + F3_in  + F4_in  + F5_in  + F6_in  + F7_in  + F8_in;
        double v1 = (F1_in+ F5_in+ F8_in -( F3_in+F6_in+F7_in)) / rho_s ;
	    double w1 = (F2_in+ F5_in+ F6_in -( F4_in+F7_in+F8_in)) / rho_s ;

        //4.5 Mach 數限制：若速度超過上限，縮放至安全範圍並重建分佈函數
        if(LimitVelocity(v1, w1, rho_s)) {
            // 速度被限制，需要用限制後的速度重建平衡態分佈函數
            RebuildEquilibrium(F0_in, F1_in, F2_in, F3_in, F4_in,
                              F5_in, F6_in, F7_in, F8_in, rho_s, v1, w1);
        }

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

//=============================================================================
// 累積 Y 方向平均速度（每個時間步呼叫）
// 輸入：v_field - Y方向速度場
// 輸出：Ub_sum_ptr - 累積的速度總和
//=============================================================================
void AccumulateUbulk(double* Volume_rate_k , double* v , double* z_global) {
    //計算同一個 k 值下，體積流率的總和
    for(int k = 3 ; k < NZ6-3 ; k++){
        //每一個k值計算彼此獨立
        Volume_rate_k[k] = 0.0 ;
        for(int j = 3 ; j  <= NY6-4 ; j++){
            double dz = (z_global[j*NZ6 + k +1 ] - z_global[j*NZ6 + k-1]) / 2.0 ;
            Volume_rate_k[k] += dz * v[j*NZ6+k] ; //同一個k值下．沿著StreamWise 的方向做疊加
        }
    }
}

//=============================================================================
// 修正外力項（每 NDTFRC 步呼叫一次）
// 使用比例控制器讓實際流速趨近目標流速 Uref
//
// 控制公式：F_new = F_old + β × (U_target - U_actual) × U_ref / L
//   - β = max(0.001, 3/Re)：控制增益
//   - U_target = Uref：目標流速
//   - U_actual = Ub_avg：實際平均流速
//
// 輸入：
//   Force - 外力陣列
//   Ub_sum_ptr - 累積的速度總和
//   NDTFRC - 累積的時間步數
// 輸出：
//   更新後的 Force[0]
//   重置 *Ub_sum_ptr = 0
//=============================================================================
void ModifyForcingTerm() {
    //計算入口區的截面速平均速度
    double Ub_avg = 0.0;
    for( int k = 3; k < NZ6-3; k++ ){
        Ub_avg = Ub_avg + Ub_avg_h[k];
        Ub_avg_h[k] = 0.0;
    }
    Ub_avg = Ub_avg / (double)((LZ-1.0))/NDTFRC;


    // 1. 計算時間與空間平均速度
    int num_cells = (NY6 - 6) * (NZ6 - 6);  // 計算區域的網格數
    double Ub_avg = (*Ub_sum_ptr) / (double)(num_cells * NDTFRC);

    // 2. 計算控制增益（低雷諾數時較大，高雷諾數時較小）
    double beta = std::fmax(0.001, 3.0 / (double)Re);
    //透過平均速度與參考速度修正外力項
    // 3. 調整外力（比例控制器）
    Force[0] = Force[0] + beta * (Uref - Ub_avg) * Uref / LZ;
    
    //2026.1.29 外力輸出加上限幅
    double F_max = 0.01 ; 
    Force[0] = std::fmin(Force[0], F_max);
    Force[0] = std::fmax(Force[0], -F_max);  

    // 4. 輸出監控資訊
    double ratio = Ub_avg / Uref;
    double Ma_actual = Ub_avg / cs;
    std::cout << "Force Update: Ub_avg=" << Ub_avg << ", Uref=" << Uref
              << ", Ub/Uref=" << ratio
              << ", Ma_actual=" << Ma_actual
              << ", Force=" << Force[0] << std::endl;

    // 5. 重置累加器
    *Ub_sum_ptr = 0.0;
}

//U_bulk_average 針對所有平面的平均的加總 git
#endif
