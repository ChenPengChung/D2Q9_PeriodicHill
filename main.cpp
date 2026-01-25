//=============================================================================
// main.cpp - D2Q9 Periodic Hill LBM 模擬主程式
// 編譯：g++ -O3 -o hill main.cpp -lm
// 執行：./hill
//=============================================================================
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>
using namespace std;
//=============================================================================
// [區塊 1] 引入參數定義（不依賴全域變數的標頭檔）
//=============================================================================
#include "variables.h"
#include "model.h"
#include "paraviewOutput.h"
//=============================================================================
// [區塊 2] 全域變數宣告
// 注意：必須在 #include "initialization.h" 之前宣告
//=============================================================================
//-----------------------------------------------------------------------------
// 2.1 流場宏觀量
//-----------------------------------------------------------------------------
double rho[NY6 * NZ6];           // 密度場
double v[NY6 * NZ6];             // Y方向速度（主流場方向）
double w[NY6 * NZ6];             // Z方向速度（垂直壁面方向）
//-----------------------------------------------------------------------------
// 2.2 分佈函數 (D2Q9)
// 使用雙緩衝區：f_old 讀取, f_new 寫入，每步交換
//-----------------------------------------------------------------------------
double f[9][NY6 * NZ6];          // 初始化用（InitialUsingDftFunc 使用）
double f_old[9][NY6 * NZ6];      // 上一時間步的分佈函數
double f_new[9][NY6 * NZ6];      // 當前時間步的分佈函數
//-----------------------------------------------------------------------------
// 2.3 外力項與質量修正
//-----------------------------------------------------------------------------
double Force[2];                 // 外力項 [0]=Fy, [1]=Fz
double rho_modify[1] = {0.0};    // 質量修正項
//-----------------------------------------------------------------------------
// 2.4 座標陣列
//-----------------------------------------------------------------------------
double y_global[NY6];            // Y方向物理座標
double z_global[NY6 * NZ6];      // (含山丘)離散化全域Z座標
double xi_h[NZ6];                // 無因次化Z座標
double nonuni_a = 0.0;           // 非均勻網格參數（初始化後設定一次）
//-----------------------------------------------------------------------------
// 2.5 Y方向插值權重 (6階Lagrange = 7點)
// YPara0: 用於 F1, F5, F8 (從 y-minSize 位置插值)
// YPara2: 用於 F3, F6, F7 (從 y+minSize 位置插值)
//-----------------------------------------------------------------------------
double* YPara0_h[3];
double* YPara2_h[3];
//-----------------------------------------------------------------------------
// 2.6 Xi方向插值權重 (每個速度方向各一組)
// 維度: [7個權重][NY6*NZ6個計算點]
//-----------------------------------------------------------------------------
double* XiParaF1_h[7];
double* XiParaF2_h[7];
double* XiParaF3_h[7];
double* XiParaF4_h[7];
double* XiParaF5_h[7];
double* XiParaF6_h[7];
double* XiParaF7_h[7];
double* XiParaF8_h[7];
//-----------------------------------------------------------------------------
// 2.7 BFL邊界條件：Y方向插值權重
// 命名規則：YBFLParaF*_h 中的 F* 表示「用來插值 F* 分佈函數」
// - YBFLParaF3: 插值 F3，用於更新 F1（左丘邊界）
// - YBFLParaF1: 插值 F1，用於更新 F3（右丘邊界）
// - YBFLParaF7: 插值 F7，用於更新 F5（左丘對角線）
// - YBFLParaF8: 插值 F8，用於更新 F6（右丘對角線）
//-----------------------------------------------------------------------------
double* YBFLParaF3_h[3];
double* YBFLParaF1_h[3];
double* YBFLParaF7_h[3];
double* YBFLParaF8_h[3];
//-----------------------------------------------------------------------------
// 2.8 BFL邊界條件：Xi方向插值權重
//-----------------------------------------------------------------------------
double* XiBFLParaF3_h[7];
double* XiBFLParaF1_h[7];
double* XiBFLParaF7_h[7];
double* XiBFLParaF8_h[7];
//-----------------------------------------------------------------------------
// 2.9 BFL邊界條件：無因次化距離 q
// q = 計算點到壁面的距離 / minSize
// 只有邊界計算點會被賦值，其他位置為 0
//-----------------------------------------------------------------------------
double Q1_h[NY6 * NZ6];          // F1 方向的 q 值（左丘 +Y）
double Q3_h[NY6 * NZ6];          // F3 方向的 q 值（右丘 -Y）
double Q5_h[NY6 * NZ6];          // F5 方向的 q 值（左丘 +Y+Z）
double Q6_h[NY6 * NZ6];          // F6 方向的 q 值（右丘 -Y+Z）
//-----------------------------------------------------------------------------
// 2.10 Stencil 起點索引（基於來源點位置預先計算）
// 用於解決 RelationXi 與 evolution.h 中 cell_z 計算不一致的問題
//-----------------------------------------------------------------------------
int CellZ_F1[3*NY6 * NZ6];         // F1 方向的 Z stencil 起點
int CellZ_F2[3*NY6 * NZ6];         // F2 方向的 Z stencil 起點
int CellZ_F3[3*NY6 * NZ6];         // F3 方向的 Z stencil 起點
int CellZ_F4[3*NY6 * NZ6];         // F4 方向的 Z stencil 起點
int CellZ_F5[3*NY6 * NZ6];         // F5 方向的 Z stencil 起點
int CellZ_F6[3*NY6 * NZ6];         // F6 方向的 Z stencil 起點
int CellZ_F7[3*NY6 * NZ6];         // F7 方向的 Z stencil 起點
int CellZ_F8[3*NY6 * NZ6];         // F8 方向的 Z stencil 起點
//-----------------------------------------------------------------------------
// 2.11 外力修正用變數
//-----------------------------------------------------------------------------
double Ub_sum = 0.0;             // 累積的平均速度
int force_update_count = 0;      // 累積的時間步數
const int NDTFRC = 10000;        // 每多少步修正一次外力
//-----------------------------------------------------------------------------
// 2.11 輸出控制變數
//-----------------------------------------------------------------------------
const int outputInterval_VTK = 1000;     // VTK 檔案輸出間隔（步數）
const int outputInterval_Stats = 1000;   // 終端統計輸出間隔（步數）
//=============================================================================
// [區塊 3] 引入依賴全域變數的標頭檔
//=============================================================================
#include "memoryAllocator.h"
#include "initializationTool.h"
#include "initialization.h"
#include "InterpolationHillISLBM.h"
#include "MRT_Matrix.h"
#include "MRT_Process.h"
#include "evolution.h"
//=============================================================================
// [區塊 4] 輔助函數
//=============================================================================
//-----------------------------------------------------------------------------
// 4.1 初始化所有陣列為零
//-----------------------------------------------------------------------------
void initializeArrays() {
    // 流場
    memset(rho, 0, sizeof(rho));
    memset(v, 0, sizeof(v));
    memset(w, 0, sizeof(w));
    // 分佈函數
    memset(f, 0, sizeof(f));
    memset(f_old, 0, sizeof(f_old));
    memset(f_new, 0, sizeof(f_new));
    // 座標
    memset(y_global, 0, sizeof(y_global));
    memset(z_global, 0, sizeof(z_global));
    memset(xi_h, 0, sizeof(xi_h));
    // BFL q 值
    memset(Q1_h, 0, sizeof(Q1_h));
    memset(Q3_h, 0, sizeof(Q3_h));
    memset(Q5_h, 0, sizeof(Q5_h));
    memset(Q6_h, 0, sizeof(Q6_h));
    // 外力
    Force[0] = 0.0;
    Force[1] = 0.0;
}

//-----------------------------------------------------------------------------
// 4.2 計算並輸出流場統計量
//-----------------------------------------------------------------------------
void printStatistics(int step) {
    double rho_sum = 0.0, rho_min = 1e10, rho_max = -1e10;
    double v_sum = 0.0, v_max = -1e10;
    double w_sum = 0.0;
    int count = 0;

    for(int j = 3; j < NY6-3; j++) {
        for(int k = 3; k < NZ6-3; k++) {
            int idx = j * NZ6 + k;
            rho_sum += rho[idx];
            rho_min = fmin(rho_min, rho[idx]);
            rho_max = fmax(rho_max, rho[idx]);
            v_sum += v[idx];
            v_max = fmax(v_max, fabs(v[idx]));
            w_sum += w[idx];
            count++;
        }
    }
    double rho_avg = rho_sum / count;
    double v_avg = v_sum / count;
    double w_avg = w_sum / count;
    ostringstream oss;
    oss << "Step " << setw(6) << step
        << ": rho=[" << fixed << setprecision(4)
        << rho_min << ", " << rho_avg << ", " << rho_max << "], "
        << "v_avg=" << scientific << setprecision(6) << v_avg
        << ", v_max=" << scientific << setprecision(6) << v_max
        << ", w_avg=" << scientific << setprecision(6) << w_avg
        << '\n';
    cout << oss.str();
}
//-----------------------------------------------------------------------------
// 4.3 交換 f_old 與 f_new
//-----------------------------------------------------------------------------
void swapDistributions() {
    for(int dir = 0; dir < 9; dir++) {
        for(int idx = 0; idx < NY6 * NZ6; idx++) {
            double temp = f_old[dir][idx];
            f_old[dir][idx] = f_new[dir][idx];
            f_new[dir][idx] = temp;
        }
    }
}
//=============================================================================
// [區塊 5] 主程式
//=============================================================================
int main() {
    //-------------------------------------------------------------------------
    // 5.1 初始化階段
    //-------------------------------------------------------------------------
    cout << "Initializing..." << endl ;
    // 步驟 1：配置記憶體
    AllocateAllWeightArrays();
    //步驟 2 :清零所有陣列
    initializeArrays();
    //步驟 3 : 建立網格
    cout << "Generating mesh..." << endl ;
    GenerateMesh_Y();
    GenerateMesh_Z();
    //步驟 4 :預計算插值權重
    cout << "Computing interpolation weights...." << endl ;
    GetIntrplParameter_Y();
    GetIntrplParameter_Xi();
    //步驟 5 : BFL 邊界初始化
    cout << "Initializing BFL boundary...." << endl ;
    BFLInitialization(Q1_h, Q3_h, Q5_h, Q6_h);
    //步驟 6 : 初始化流場與分佈函數
    cout << "Initializing flow field...." << endl ;
    InitialUsingDftFunc();
    //步驟 7 : 複製初始分佈函數到 f_old
    for(int dir = 0; dir < 9; dir++) {
        for(int idx = 0; idx < NY6 * NZ6; idx++) {
            f_old[dir][idx] = f[dir][idx];
        }
    }//f_old為上一個時間步所更新的碰撞後插值後一般態分佈函數 
    cout << "Initialization complete.." << endl ;
    //-------------------------------------------------------------------------
    // 5.2 時間迴圈
    //-------------------------------------------------------------------------
    cout << "Starting time loop..." << endl ;
    cout << "loop = " << loop << endl ;
    cout.flush();
    for(int t = 0; t < loop; t++) {
        if(t % 10 == 0) { 
            double rho_check = CheckMassConservation(rho, t);
            if(std::isnan(rho_check) || std::isinf(rho_check)) {
                cout << "DIVERGED at t = " << t << ", rho = " << rho_check << endl;
                return 1;
            }
            if(t < 200) cout << "t = " << t << ", rho = " << rho_check << endl; 
            cout.flush(); 
        }
        // 5.2.1 Stream + Collide
        stream_collide(
            // f_old (9個)
            f_old[0], f_old[1], f_old[2], f_old[3], f_old[4],
            f_old[5], f_old[6], f_old[7], f_old[8],
            // f_new (9個)
            f_new[0], f_new[1], f_new[2], f_new[3], f_new[4],
            f_new[5], f_new[6], f_new[7], f_new[8],
            // Y方向權重 YPara0 (7個)
            YPara0_h[0], YPara0_h[1], YPara0_h[2], 
            //YPara0_h[3], YPara0_h[4], YPara0_h[5], YPara0_h[6],
            // Y方向權重 YPara2 (7個)
            YPara2_h[0], YPara2_h[1], YPara2_h[2], 
            //YPara2_h[3], YPara2_h[4], YPara2_h[5], YPara2_h[6],
            // Xi方向權重 F1 (7個)
            XiParaF1_h[0], XiParaF1_h[1], XiParaF1_h[2], XiParaF1_h[3],
            XiParaF1_h[4], XiParaF1_h[5], XiParaF1_h[6],
            // Xi方向權重 F2 (7個)
            XiParaF2_h[0], XiParaF2_h[1], XiParaF2_h[2], XiParaF2_h[3],
            XiParaF2_h[4], XiParaF2_h[5], XiParaF2_h[6],
            // Xi方向權重 F3 (7個)
            XiParaF3_h[0], XiParaF3_h[1], XiParaF3_h[2], XiParaF3_h[3],
            XiParaF3_h[4], XiParaF3_h[5], XiParaF3_h[6],
            // Xi方向權重 F4 (7個)
            XiParaF4_h[0], XiParaF4_h[1], XiParaF4_h[2], XiParaF4_h[3],
            XiParaF4_h[4], XiParaF4_h[5], XiParaF4_h[6],
            // Xi方向權重 F5 (7個)
            XiParaF5_h[0], XiParaF5_h[1], XiParaF5_h[2], XiParaF5_h[3],
            XiParaF5_h[4], XiParaF5_h[5], XiParaF5_h[6],
            // Xi方向權重 F6 (7個)
            XiParaF6_h[0], XiParaF6_h[1], XiParaF6_h[2], XiParaF6_h[3],
            XiParaF6_h[4], XiParaF6_h[5], XiParaF6_h[6],
            // Xi方向權重 F7 (7個)
            XiParaF7_h[0], XiParaF7_h[1], XiParaF7_h[2], XiParaF7_h[3],
            XiParaF7_h[4], XiParaF7_h[5], XiParaF7_h[6],
            // Xi方向權重 F8 (7個)
            XiParaF8_h[0], XiParaF8_h[1], XiParaF8_h[2], XiParaF8_h[3],
            XiParaF8_h[4], XiParaF8_h[5], XiParaF8_h[6],
            // BFL Y方向權重 F3 (7個)
            YBFLParaF3_h[0], YBFLParaF3_h[1], YBFLParaF3_h[2], //YBFLParaF3_h[3], YBFLParaF3_h[4], YBFLParaF3_h[5], YBFLParaF3_h[6],
            // BFL Y方向權重 F1 (7個)
            YBFLParaF1_h[0], YBFLParaF1_h[1], YBFLParaF1_h[2], //YBFLParaF1_h[3], YBFLParaF1_h[4], YBFLParaF1_h[5], YBFLParaF1_h[6],
            // BFL Y方向權重 F7 (7個)
            YBFLParaF7_h[0], YBFLParaF7_h[1], YBFLParaF7_h[2], //YBFLParaF7_h[3], YBFLParaF7_h[4], YBFLParaF7_h[5], YBFLParaF7_h[6],
            // BFL Y方向權重 F8 (7個)
            YBFLParaF8_h[0], YBFLParaF8_h[1], YBFLParaF8_h[2], //YBFLParaF8_h[3],  YBFLParaF8_h[4], YBFLParaF8_h[5], YBFLParaF8_h[6],
            // BFL Xi方向權重 F3 (7個)
            XiBFLParaF3_h[0], XiBFLParaF3_h[1], XiBFLParaF3_h[2], XiBFLParaF3_h[3],
            XiBFLParaF3_h[4], XiBFLParaF3_h[5], XiBFLParaF3_h[6],
            // BFL Xi方向權重 F1 (7個)
            XiBFLParaF1_h[0], XiBFLParaF1_h[1], XiBFLParaF1_h[2], XiBFLParaF1_h[3],
            XiBFLParaF1_h[4], XiBFLParaF1_h[5], XiBFLParaF1_h[6],
            // BFL Xi方向權重 F7 (7個)
            XiBFLParaF7_h[0], XiBFLParaF7_h[1], XiBFLParaF7_h[2], XiBFLParaF7_h[3],
            XiBFLParaF7_h[4], XiBFLParaF7_h[5], XiBFLParaF7_h[6],
            // BFL Xi方向權重 F8 (7個)
            XiBFLParaF8_h[0], XiBFLParaF8_h[1], XiBFLParaF8_h[2], XiBFLParaF8_h[3],
            XiBFLParaF8_h[4], XiBFLParaF8_h[5], XiBFLParaF8_h[6],
            // 宏觀參數
            v, w, rho, Force, rho_modify,
            // BFL q值
            Q1_h, Q3_h, Q5_h, Q6_h
        );

        // 5.2.2 週期性邊界條件
        periodicSW(
            f_new[0], f_new[1], f_new[2], f_new[3], f_new[4],
            f_new[5], f_new[6], f_new[7], f_new[8],
            v, w, rho
        );

        // 5.2.3 累積平均速度（用於外力修正）
        AccumulateUbulk(v, &Ub_sum);
        force_update_count++;

        // 5.2.4 每 NDTFRC 步修正外力
        if(force_update_count >= NDTFRC) {
            ModifyForcingTerm(Force, &Ub_sum, NDTFRC);
            force_update_count = 0;
        }

        // 5.2.5 交換 f_old 與 f_new
        swapDistributions();

        //---------------------------------------------------------------------
        // 5.2.6 輸出數據
        //---------------------------------------------------------------------
        // VTK 檔案輸出（較大間隔，用於 ParaView 視覺化）
        if(t % outputInterval_VTK == 0) {
            cout << "==================================================================" << endl ; 
            cout << "\n=== time step t = " << t << " ===" << endl;
            cout << "Output VTK file and Statics Data..." << endl;
            OutputVTK(t, y_global, z_global, rho, v, w);
        }
        // 終端統計輸出（較小間隔，用於監控模擬進度）
        if(t % outputInterval_Stats == 0) {
            //每 outputInterval_Stats 步輸出一次平均密度
            cout << "++      ++     ++      ++      ++      ++    ++      ++     ++      ++" << endl ; 
            cout << "Time=" << t << setw(6) <<
            " ; Average Density=" << CheckMassConservation(rho,t) << setw(6) << 
            " ; Density Correction=" << fabs(rho_modify[0] )<< endl ;
            //printStatistics(t);
        }
    }
    //-------------------------------------------------------------------------
    // 5.3 結束
    //-------------------------------------------------------------------------
    cout << "Simulation complete." << endl ;
    printStatistics(loop);
    return 0;
}
