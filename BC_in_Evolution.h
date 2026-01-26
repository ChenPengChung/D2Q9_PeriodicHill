// =============================================================================
// BC_in_Evolution.h - 備用邊界條件處理函數
// 包含兩種版本的插值/streaming 處理策略
// =============================================================================
#ifndef BC_IN_EVOLUTION_H
#define BC_IN_EVOLUTION_H

#include "globalVariables.h"
#include "model.h"
#include "InterpolationHillISLBM.h"
#include "initializationTool.h"
#include <iostream>
// =============================================================================
// 版本一：統一邊界區域處理
// 使用固定邊界範圍判斷，在邊界區用 streaming，內部用插值
// 優點：簡單直觀，避免插值越界
// =============================================================================
inline void BC1(
    int j, int k, int idx_xi, int nface,
    double* f0_old, double* f1_old, double* f2_old, double* f3_old, double* f4_old,
    double* f5_old, double* f6_old, double* f7_old, double* f8_old,
    double& F0_in, double& F1_in, double& F2_in, double& F3_in, double& F4_in,
    double& F5_in, double& F6_in, double& F7_in, double& F8_in,
    // Y 插值權重
    double* Y0_0, double* Y0_1, double* Y0_2,
    double* Y2_0, double* Y2_1, double* Y2_2,
    // Xi 插值權重 (F1~F8)
    double* XiF1_0, double* XiF1_1, double* XiF1_2, double* XiF1_3, double* XiF1_4, double* XiF1_5, double* XiF1_6,
    double* XiF2_0, double* XiF2_1, double* XiF2_2, double* XiF2_3, double* XiF2_4, double* XiF2_5, double* XiF2_6,
    double* XiF3_0, double* XiF3_1, double* XiF3_2, double* XiF3_3, double* XiF3_4, double* XiF3_5, double* XiF3_6,
    double* XiF4_0, double* XiF4_1, double* XiF4_2, double* XiF4_3, double* XiF4_4, double* XiF4_5, double* XiF4_6,
    double* XiF5_0, double* XiF5_1, double* XiF5_2, double* XiF5_3, double* XiF5_4, double* XiF5_5, double* XiF5_6,
    double* XiF6_0, double* XiF6_1, double* XiF6_2, double* XiF6_3, double* XiF6_4, double* XiF6_5, double* XiF6_6,
    double* XiF7_0, double* XiF7_1, double* XiF7_2, double* XiF7_3, double* XiF7_4, double* XiF7_5, double* XiF7_6,
    double* XiF8_0, double* XiF8_1, double* XiF8_2, double* XiF8_3, double* XiF8_4, double* XiF8_5, double* XiF8_6
) {
    // CellZ 陣列對應關係（根據 initialization.h）:
    // F1, F3: z 無偏移，使用 CellZ_F1
    // F2, F5, F6: z-Δ，使用 CellZ_F2
    // F4, F7, F8: z+Δ，使用 CellZ_F4
    
    F0_in = f0_old[idx_xi];
    
    // Y 方向邊界檢查（週期性邊界需要用 streaming，不用插值）
    bool y_boundary = (j <= 3) || (j >= NY6-4);
    // Z 方向邊界：下邊界和上邊界
    bool z_lower = (k <= 11);
    bool z_upper = (k >= NZ6-10);
    
    if( z_lower || z_upper || y_boundary ) {
        // ====================================================================
        // 邊界附近：使用簡單的 streaming 替代插值
        // ====================================================================
        
        // F1,F3: Y方向 streaming，需要週期性 wrap
        int jm1 = (j > 3) ? (j-1) : (NY6-4);  // 週期性 wrap
        int jp1 = (j < NY6-4) ? (j+1) : 3;     // 週期性 wrap
        F1_in = f1_old[jm1*NZ6 + k];
        F3_in = f3_old[jp1*NZ6 + k];
        
        // 根據 Z 位置處理 F2,F4,F5,F6,F7,F8
        if( z_lower ) {
            // Z下邊界：F2,F5,F6 使用 bounce-back (因為來源在下邊界外)
            F2_in = f4_old[idx_xi];
            F5_in = f7_old[idx_xi];
            F6_in = f8_old[idx_xi];
            
            // F4,F7,F8: 向 -Z 方向，從 z+Δ 取值
            F4_in = f4_old[j*NZ6 + k+1];
            F7_in = f7_old[jp1*NZ6 + k+1];
            F8_in = f8_old[jm1*NZ6 + k+1];
            
        } else if( z_upper ) {
            // Z上邊界：F4,F7,F8 使用 bounce-back (因為來源在上邊界外)
            F4_in = f2_old[idx_xi];
            F7_in = f5_old[idx_xi];
            F8_in = f6_old[idx_xi];
            
            // F2,F5,F6: 向 +Z 方向，從 z-Δ 取值
            F2_in = f2_old[j*NZ6 + k-1];
            F5_in = f5_old[jm1*NZ6 + k-1];
            F6_in = f6_old[jp1*NZ6 + k-1];
            
        } else {
            // Y邊界但非Z邊界：使用簡單 streaming
            F2_in = f2_old[j*NZ6 + k-1];
            F4_in = f4_old[j*NZ6 + k+1];
            F5_in = f5_old[jm1*NZ6 + k-1];
            F6_in = f6_old[jp1*NZ6 + k-1];
            F7_in = f7_old[jp1*NZ6 + k+1];
            F8_in = f8_old[jm1*NZ6 + k+1];
        }
        
    } else {
        // ====================================================================
        // 內部區域：正常插值
        // ====================================================================
        // 使用正確的 CellZ 對應關係
        F1_Intrpl3(f1_old, j, k, CellZ_F1, j, idx_xi, Y0_0, Y0_1, Y0_2, XiF1_0, XiF1_1, XiF1_2, XiF1_3, XiF1_4, XiF1_5, XiF1_6);
        F3_Intrpl3(f3_old, j, k, CellZ_F1, j, idx_xi, Y2_0, Y2_1, Y2_2, XiF3_0, XiF3_1, XiF3_2, XiF3_3, XiF3_4, XiF3_5, XiF3_6);
        F2_Intrpl7(f2_old, j, k, CellZ_F2, j, idx_xi, XiF2_0, XiF2_1, XiF2_2, XiF2_3, XiF2_4, XiF2_5, XiF2_6);
        F4_Intrpl7(f4_old, j, k, CellZ_F4, j, idx_xi, XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6);
        // F5, F6: z-Δ 方向，使用 CellZ_F2
        Y_XI_Intrpl3(f5_old, F5_in, j, k, CellZ_F2, j, idx_xi, Y0_0, Y0_1, Y0_2, XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6);
        Y_XI_Intrpl3(f6_old, F6_in, j, k, CellZ_F2, j, idx_xi, Y2_0, Y2_1, Y2_2, XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6);
        // F7, F8: z+Δ 方向，使用 CellZ_F4
        Y_XI_Intrpl3(f7_old, F7_in, j, k, CellZ_F4, j, idx_xi, Y2_0, Y2_1, Y2_2, XiF7_0, XiF7_1, XiF7_2, XiF7_3, XiF7_4, XiF7_5, XiF7_6);
        Y_XI_Intrpl3(f8_old, F8_in, j, k, CellZ_F4, j, idx_xi, Y0_0, Y0_1, Y0_2, XiF8_0, XiF8_1, XiF8_2, XiF8_3, XiF8_4, XiF8_5, XiF8_6);
    }
}

// =============================================================================
// 版本一的 BFL 曲面邊界處理（配合 BC_Version1 使用）
// =============================================================================
inline void BC1_BFL(
    int j, int k, int idx_xi,
    double* f1_old, double* f3_old, double* f5_old, double* f6_old, double* f7_old, double* f8_old,
    double& F1_in, double& F2_in, double& F3_in, double& F5_in, double& F6_in,
    double* Q1_h, double* Q3_h, double* Q5_h, double* Q6_h,
    double* y_global, double* z_global
) {
    // Half-way bounce-back 底部邊界
    if( k == 3 ){
        F2_in = f5_old[idx_xi];  // 使用已傳入的 f4_old 值 (caller 處理)
        F5_in = f7_old[idx_xi];
        F6_in = f8_old[idx_xi];
    }
    
    // BFL 曲面邊界處理
    // 左丘邊界，更新F1
    if(IsLeftHill_Boundary_yPlus(y_global[j], z_global[j*NZ6+k])){
        double q1 = Q1_h[idx_xi];
        if(q1 < 0.5 && q1 >= 0){
            F1_in = (2*q1)*f3_old[idx_xi] + (1.0 - 2.0*q1)*f3_old[idx_xi+NZ6];
        }
        if(q1 > 0.5){
            F1_in = (1.0/(2.0*q1))*f3_old[idx_xi] + ((2.0*q1-1.0)/(2.0*q1))*f1_old[idx_xi];
        }
    }
    
    // 右丘邊界，更新F3
    if(IsRightHill_Boundary_yMinus(y_global[j], z_global[j*NZ6+k])){
        double q3 = Q3_h[idx_xi];
        if(q3 < 0.5 && q3 >= 0.0){
            F3_in = (2*q3)*f1_old[idx_xi] + (1.0 - 2.0*q3)*f1_old[idx_xi-NZ6];
        }
        if(q3 > 0.5){
            F3_in = (1.0/(2.0*q3))*f1_old[idx_xi] + ((2.0*q3-1.0)/(2.0*q3))*f3_old[idx_xi];
        }
    }
    
    // 左丘邊界，更新F5
    if(IsLeftHill_Boundary_Diagonal45(y_global[j], z_global[j*NZ6+k])){
        double q5 = Q5_h[idx_xi];
        if(q5 < 0.5 && q5 >= 0.0){
            F5_in = (2*q5)*f7_old[idx_xi] + (1.0 - 2.0*q5)*f7_old[idx_xi+NZ6+1];
        }
        if(q5 > 0.5){
            F5_in = (1.0/(2.0*q5))*f7_old[idx_xi] + ((2.0*q5-1.0)/(2.0*q5))*f5_old[idx_xi];
        }
    }
    
    // 右丘邊界，更新F6
    if(IsRightHill_Boundary_Diagonal135(y_global[j], z_global[j*NZ6+k])){
        double q6 = Q6_h[idx_xi];
        if(q6 < 0.5 && q6 >= 0.0){
            F6_in = (2*q6)*f8_old[idx_xi] + (1.0 - 2.0*q6)*f8_old[idx_xi-NZ6+1];
        }
        if(q6 > 0.5){
            F6_in = (1.0/(2.0*q6))*f8_old[idx_xi] + ((2.0*q6-1.0)/(2.0*q6))*f6_old[idx_xi];
        }
    }
}

// =============================================================================
// 版本二：每個分佈函數獨立處理
// 優先順序：BFL 曲面 → 底/頂邊界 bounce-back → 內部插值
// 優點：曲面邊界優先，更精確處理山丘區域
// =============================================================================
inline void BC2(
    int j, int k, int idx_xi, int nface,
    double* f0_old, double* f1_old, double* f2_old, double* f3_old, double* f4_old,
    double* f5_old, double* f6_old, double* f7_old, double* f8_old,
    double& F0_in, double& F1_in, double& F2_in, double& F3_in, double& F4_in,
    double& F5_in, double& F6_in, double& F7_in, double& F8_in,
    double* Q1_h, double* Q3_h, double* Q5_h, double* Q6_h,
    double* y_global, double* z_global,
    // Y 插值權重
    double* Y0_0, double* Y0_1, double* Y0_2,
    double* Y2_0, double* Y2_1, double* Y2_2,
    // Xi 插值權重
    double* XiF1_0, double* XiF1_1, double* XiF1_2, double* XiF1_3, double* XiF1_4, double* XiF1_5, double* XiF1_6,
    double* XiF2_0, double* XiF2_1, double* XiF2_2, double* XiF2_3, double* XiF2_4, double* XiF2_5, double* XiF2_6,
    double* XiF3_0, double* XiF3_1, double* XiF3_2, double* XiF3_3, double* XiF3_4, double* XiF3_5, double* XiF3_6,
    double* XiF4_0, double* XiF4_1, double* XiF4_2, double* XiF4_3, double* XiF4_4, double* XiF4_5, double* XiF4_6,
    double* XiF5_0, double* XiF5_1, double* XiF5_2, double* XiF5_3, double* XiF5_4, double* XiF5_5, double* XiF5_6,
    double* XiF6_0, double* XiF6_1, double* XiF6_2, double* XiF6_3, double* XiF6_4, double* XiF6_5, double* XiF6_6,
    double* XiF7_0, double* XiF7_1, double* XiF7_2, double* XiF7_3, double* XiF7_4, double* XiF7_5, double* XiF7_6,
    double* XiF8_0, double* XiF8_1, double* XiF8_2, double* XiF8_3, double* XiF8_4, double* XiF8_5, double* XiF8_6
) {
    // CellZ 陣列對應關係:
    // F1, F3: z 無偏移，使用 CellZ_F1
    // F2, F5, F6: z-Δ，使用 CellZ_F2
    // F4, F7, F8: z+Δ，使用 CellZ_F4
    
    F0_in = f0_old[idx_xi];
    
    // ========================================================================
    // F2: +Z 方向
    // ========================================================================
    if(k == 3) {
        F2_in = f4_old[idx_xi];  // half-way bounce-back
    } else {
        F2_Intrpl7(f2_old, j, k, CellZ_F2, j, idx_xi, XiF2_0, XiF2_1, XiF2_2, XiF2_3, XiF2_4, XiF2_5, XiF2_6);
    }
    
    // ========================================================================
    // F4: -Z 方向
    // ========================================================================
    if(k == NZ6-4) {
        F4_in = f2_old[idx_xi];  // half-way bounce-back
    } else {
        F4_Intrpl7(f4_old, j, k, CellZ_F4, j, idx_xi, XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6);
    }
    
    // ========================================================================
    // F1: +Y 方向 (含左丘 BFL)
    // ========================================================================
    if(k == 3) {
        F1_in = f3_old[idx_xi];  // 底部邊界 bounce-back
    } else if(IsLeftHill_Boundary_yPlus(y_global[j], z_global[j*NZ6+k]) && Q1_h[idx_xi] > 0.0) {
        double q1 = Q1_h[idx_xi];
        if(q1 < 0.5 && q1 >= 0) {
            F1_in = (2*q1)*f3_old[idx_xi] + (1.0 - 2.0*q1)*f3_old[idx_xi+NZ6];
        } else if(q1 >= 0.5) {
            F1_in = (1.0/(2.0*q1))*f3_old[idx_xi] + ((2.0*q1-1.0)/(2.0*q1))*f1_old[idx_xi];
        }
    } else {
        F1_Intrpl3(f1_old, j, k, CellZ_F1, j, idx_xi, Y0_0, Y0_1, Y0_2, XiF1_0, XiF1_1, XiF1_2, XiF1_3, XiF1_4, XiF1_5, XiF1_6);
    }
    
    // ========================================================================
    // F3: -Y 方向 (含右丘 BFL)
    // ========================================================================
    if(k == 3) {
        F3_in = f1_old[idx_xi];  // 底部邊界 bounce-back
    } else if(IsRightHill_Boundary_yMinus(y_global[j], z_global[j*NZ6+k]) && Q3_h[idx_xi] > 0.0) {
        double q3 = Q3_h[idx_xi];
        if(q3 < 0.5 && q3 >= 0.0) {
            F3_in = (2*q3)*f1_old[idx_xi] + (1.0 - 2.0*q3)*f1_old[idx_xi-NZ6];
        } else if(q3 >= 0.5) {
            F3_in = (1.0/(2.0*q3))*f1_old[idx_xi] + ((2.0*q3-1.0)/(2.0*q3))*f3_old[idx_xi];
        }
    } else {
        F3_Intrpl3(f3_old, j, k, CellZ_F1, j, idx_xi, Y2_0, Y2_1, Y2_2, XiF3_0, XiF3_1, XiF3_2, XiF3_3, XiF3_4, XiF3_5, XiF3_6);
    }
    
    // ========================================================================
    // F5: +Y, +Z 方向 (含左丘對角 BFL)
    // ========================================================================
    if(IsLeftHill_Boundary_Diagonal45(y_global[j], z_global[j*NZ6+k]) && Q5_h[idx_xi] > 0.0) {
        double q5 = Q5_h[idx_xi];
        if(q5 < 0.5 && q5 >= 0.0) {
            F5_in = (2*q5)*f7_old[idx_xi] + (1.0 - 2.0*q5)*f7_old[idx_xi+NZ6+1];
        }
        if(q5 > 0.5) {
            F5_in = (1.0/(2.0*q5))*f7_old[idx_xi] + ((2.0*q5-1.0)/(2.0*q5))*f5_old[idx_xi];
        }
    } else if(k == 3) {
        F5_in = f7_old[idx_xi];  // half-way bounce-back
    } else {
        Y_XI_Intrpl3(f5_old, F5_in, j, k, CellZ_F2, j, idx_xi, Y0_0, Y0_1, Y0_2, XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6);
    }
    
    // ========================================================================
    // F6: -Y, +Z 方向 (含右丘對角 BFL)
    // ========================================================================
    if(IsRightHill_Boundary_Diagonal135(y_global[j], z_global[j*NZ6+k]) && Q6_h[idx_xi] > 0.0) {
        double q6 = Q6_h[idx_xi];
        if(q6 < 0.5 && q6 >= 0.0) {
            F6_in = (2*q6)*f8_old[idx_xi] + (1.0 - 2.0*q6)*f8_old[idx_xi-NZ6+1];
        }
        if(q6 > 0.5) {
            F6_in = (1.0/(2.0*q6))*f8_old[idx_xi] + ((2.0*q6-1.0)/(2.0*q6))*f6_old[idx_xi];
        }
    } else if(k == 3) {
        F6_in = f8_old[idx_xi];  // half-way bounce-back
    } else {
        Y_XI_Intrpl3(f6_old, F6_in, j, k, CellZ_F2, j, idx_xi, Y2_0, Y2_1, Y2_2, XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6);
    }
    
    // ========================================================================
    // F7: -Y, -Z 方向
    // ========================================================================
    if(k == NZ6-4) {
        F7_in = f5_old[idx_xi];  // half-way bounce-back
    } else {
        Y_XI_Intrpl3(f7_old, F7_in, j, k, CellZ_F4, j, idx_xi, Y2_0, Y2_1, Y2_2, XiF7_0, XiF7_1, XiF7_2, XiF7_3, XiF7_4, XiF7_5, XiF7_6);
    }
    
    // ========================================================================
    // F8: +Y, -Z 方向
    // ========================================================================
    if(k == NZ6-4) {
        F8_in = f6_old[idx_xi];  // half-way bounce-back
    } else {
        Y_XI_Intrpl3(f8_old, F8_in, j, k, CellZ_F4, j, idx_xi, Y0_0, Y0_1, Y0_2, XiF8_0, XiF8_1, XiF8_2, XiF8_3, XiF8_4, XiF8_5, XiF8_6);
    }
}

// =============================================================================
// 診斷函數：檢查異常值
// =============================================================================
inline void BC_DiagnoseAbnormal(
    int j, int k,
    double F0_in, double F1_in, double F2_in, double F3_in, double F4_in,
    double F5_in, double F6_in, double F7_in, double F8_in
) {
    double rho_local = F0_in + F1_in + F2_in + F3_in + F4_in + F5_in + F6_in + F7_in + F8_in;
    if(rho_local > 2.0 || rho_local < 0.5 || std::isnan(rho_local)) {
        static int error_count = 0;
        if(error_count < 20) {
            bool in_lower = (k <= 15);
            bool in_upper = (k >= NZ6-35);
            bool in_ybnd = (j <= 4) || (j >= NY6-5);
            std::cout << "ABNORMAL at j=" << j << " k=" << k 
                      << " region=" << (in_lower ? "LOWER" : (in_upper ? "UPPER" : (in_ybnd ? "YBND" : "INTERNAL")))
                      << " rho=" << rho_local
                      << " F=[" << F0_in << "," << F1_in << "," << F2_in << "," << F3_in << "," << F4_in
                      << "," << F5_in << "," << F6_in << "," << F7_in << "," << F8_in << "]" << std::endl;
            error_count++;
        }
    }
}

#endif // BC_IN_EVOLUTION_H
