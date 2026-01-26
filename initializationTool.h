#ifndef INITIALIZATIONTOOL_FILE
#define INITIALIZATIONTOOL_FILE

/**
 * @file initializationTool.h
 * @brief 初始化工具函數集 - 用於週期性山丘流場的幾何計算與插值
 * * 本檔案包含：
 * * - 山丘函數反函數計算
 * * - 非均勻網格參數計算
 * * - Lagrange 插值函數
 * * - BFL 邊界條件判斷函數
 */
#include <cmath>
#include <fstream>
#include "variables.h"
#include "model.h"
using namespace std;

// * 半邊山丘寬度 (無因次化)
#define HillHalfWidth (54.0/28.0) 

//1.
/**
 * @brief 左半丘反函數 (數值法)
 * * 使用二分法求解左半丘的反函數，給定高度 z 返回對應的 y 座標
 * @param z 目標高度值
 * @return 對應的 y 座標
 * ! 注意：使用二分法，精度為 1e-12
 * ? 左半丘是遞減函數：y 越大，HillFunction(y) 越小
 * @see HillFunction(), HillFunction_Inverse_Right()
 */
double HillFunction_Inverse_Left(double z) {
    double y_low = 0.0, y_high = HillHalfWidth;
    double y_mid;

    // ? 左半丘是遞減函數: y 越大，HillFunction(y) 越小
    while (y_high - y_low > 1e-12) {
        y_mid = (y_low + y_high) / 2.0;
        if (HillFunction(y_mid) > z) {
            y_low = y_mid;
        } else {
            y_high = y_mid;
        }
    }
    return y_mid;
}

//2.
/**
 * @brief 右半丘反函數 (數值法)
 * * 使用二分法求解右半丘的反函數，給定高度 z 返回對應的 y 座標
 * @param z 目標高度值
 * @return 對應的 y 座標
 * ! 注意：使用二分法，精度為 1e-12
 * ? 右半丘是遞增函數：y 越大，HillFunction(y) 越大
 * @see HillFunction(), HillFunction_Inverse_Left()
 */
double HillFunction_Inverse_Right(double z) {
    double y_low = LY - HillHalfWidth, y_high = LY;
    double y_mid;
    // ? 右半丘是遞增函數: y 越大，HillFunction(y) 越大
    while (y_high - y_low > 1e-12) {
        y_mid = (y_low + y_high) / 2.0;
        if (HillFunction(y_mid) < z) {
            y_low = y_mid;
        } else {
            y_high = y_mid;
        }
    }
    return y_mid;
}

//3.
/**
 * @brief 雙曲正切非均勻網格座標轉換巨集
 * * 核心網格轉換公式
 * @param L 計算點所包圍的總長度
 * @param MinSize 粒子晶格大小 (最小網格尺寸)
 * @param a 非均勻網格伸縮參數 (0 < a < 1)
 * @param j 計算點編號
 * @param N 計算點為節點之網格數量
 * ! 注意：a 必須在 (0, 1) 範圍內
 * @return 第 j 個計算點的物理座標
 */
#define tanhFunction( L , MinSize , a, j, N)\
(           \
    L/2.0 + MinSize/2.0 + ((L/2.0)/a)*tanh((-1.0+2.0*(double)(j)/(double)(N))/2.0*log((1.0+a)/(1.0-a)))\
)//實際應用 : tanhFunction(LT,mimnSize,a,j,N) - 0.5minSize

// 將 xi_h（已扣掉 minSize/2.0 的值）反映回連續索引 j，用於建立 xi_h 與 j 的對應
inline double Inverse_tanh_index(double xi_val, double L, double MinSize, double a, int N) {
    const double val = xi_val + MinSize / 2.0;   // 還原到 tanhFunction 的輸出值
    const double center = L / 2.0 + MinSize / 2.0;
    const double scale = (L / 2.0) / a;
    const double log_term = log((1.0 + a) / (1.0 - a));
    double u = (val - center) / scale;           // tanh 的輸出
    if (u > 0.999999999) u = 0.999999999;
    if (u < -0.999999999) u = -0.999999999;
    double t = atanh(u);
    double j_cont = (double)N * (1.0 + (2.0 / log_term) * t) / 2.0;
    if (j_cont < 0.0) j_cont = 0.0;
    if (j_cont > (double)N) j_cont = (double)N;
    return j_cont;
}

//4.
/**
 * @brief 計算非均勻網格伸縮參數 a
 * * 使用二分法求解伸縮參數，使最小網格間距等於 minSize
 * @return 伸縮參數 a (0 < a < 1)
 * ! minSize 定義為以 LZ-(y=0.0 下的山坡高度) 為總長度均勻切割下的網格大小 × 0.6
 * TODO: 考慮改用 Newton-Raphson 加速收斂
 * @see tanhFunction
 */
double GetNonuniParameter() {
    double total = LZ - HillFunction( 0.0 ) - minSize;
    double a_temp[2] = {0.1, 1.0};
    double a_mid;

    double x_temp[2], dx;
    do{
        a_mid = (a_temp[0]+a_temp[1]) / 2.0;
        // * dx = Z方向y = 0.0 非均勻網格下的最小間距 (從y=0.0取)
        // ? minSize的定義為以LZ-(y = 0.0下的山坡高度)為總長度均勻切割下的網格大小*0.6
        // ! 判斷標準：最小的網格必須要大於minSize
        x_temp[0] = tanhFunction(total, minSize, a_mid, 0, (NZ6-7));
        x_temp[1] = tanhFunction(total, minSize, a_mid, 1, (NZ6-7));
        dx = x_temp[1] - x_temp[0];
        if( dx - minSize >= 0.0 ){
            a_temp[0] = a_mid;
        } else {
            a_temp[1] = a_mid;
        }
    } while (fabs( dx - minSize) > 1e-14 );
    return a_mid;
}


//5.
/**
 * @brief 六階 Lagrange 插值基底函數
 * * 計算在位置 pos 處，以 x_i 為插值節點的 Lagrange 基底函數值
 * @param pos 目標插值位置
 * @param x_i 當前插值節點
 * @param x1 第 1 個節點座標
 * @param x2 第 2 個節點座標
 * @param x3 第 3 個節點座標
 * @param x4 第 4 個節點座標
 * @param x5 第 5 個節點座標
 * @param x6 第 6 個節點座標
 * ! 注意：確保所有節點座標互不相同，避免除零錯誤
 * @return Lagrange 基底函數值
 * @see GetParameter_6th()
 */
//線性外插
double Extrapolation(double pos , double x1 , double x2 ){
    double a = ( (pos - x1)/(x2 - x1) ) ; 
    return a; 
}

//夾層
double Lagrange_2nd(double pos , double x_i , double x1 , double x2 ){
    double Lagrange = (pos - x1)/(x_i - x1)*(pos - x2)/(x_i - x2);
    return Lagrange; 
}
//內層
double Lagrange_6th(double pos , double x_i , double x1 , double x2 , double x3 , double x4 , double x5 , double x6){
    double Lagrange = (pos - x1)/(x_i - x1)*(pos - x2)/(x_i - x2)*(pos - x3)/(x_i - x3)*(pos - x4)/(x_i - x4)*(pos - x5)/(x_i - x5)*(pos - x6)/(x_i - x6);
    return Lagrange; 
}


void RelationXi(double pos_z , int j , int k ,  int* cell_z , double a){//double* RelazationXi 為輸出七點座標
    //pos_z可能為 k的Z値 + minSize 或 k - minSize
    if (k<3){
        //防呆設計
        cout << "Error: k value need to be larger than 3 !" << endl;
        exit(1); //k  = 3 進行差值，則山坡區域將產生外插問題
    }else{
        //因為在lAGRANGE 插值中，需要以相同z座標(有因次含山丘版本)最為基準，由左邊二點與右邊一點取y方向之三點格式內插
        //令內插目標點的z座標 為 pos_z
        //該Y座標的z方向計算點所包圍的長度為L
        //目標: [pos_z]->[j]->[k]
        //寫入三個 相鄰 y座標的Z座標起始點
        //===============================//
        //寫入cell_z[idx_xi+0] , cell_z[idx_xi+1] , cell_z[idx_xi+2]作為起始點，共同點目標內插分的Z_global相同
        //儲存索引使用原始 j，計算使用周期映射後的 j_calc
        const int j_store = j;  // 保留原始 j 用於儲存
        int j_calc = j;         // 用於計算 y_global 的週期映射
        if(j_calc<4) j_calc = j_calc + NY6-7 ;
        if(j_calc>NY6-5) j_calc = j_calc - (NY6-7) ;
        double pos_z0 = pos_z -  0.5*minSize - HillFunction(y_global[j_calc-1]);
        double pos_z1 = pos_z -  0.5*minSize - HillFunction(y_global[j_calc]);
        double pos_z2 = pos_z -  0.5*minSize - HillFunction(y_global[j_calc+1]);
        double L0 = LZ - HillFunction(y_global[j_calc-1]) - minSize;
        double L1 = LZ - HillFunction(y_global[j_calc]) - minSize;
        double L2 = LZ - HillFunction(y_global[j_calc+1]) - minSize;
        //計算該高度在不同y值上的編號
        double index_z0 = Inverse_tanh_index( pos_z0 , L0 , minSize , a , (NZ6-7) )+3;
        double index_z1 = Inverse_tanh_index( pos_z1 , L1 , minSize , a , (NZ6-7) )+3;
        double index_z2 = Inverse_tanh_index( pos_z2 , L2 , minSize , a , (NZ6-7) )+3;//a 為 非均勻網格伸縮參數
        int k_0 = 0, k_1 = 0, k_2 = 0; 
        //找最接近物理空間計算點 編號  //Danger Zone 起點判斷 
        //剩下區域用七點內插，所以結構為---
        //j-1:外插|三點內插|七點內插|三點內插|外插----
        //j  :外插|三點內插|七點內插|三點內插|外插----
        //j+1:外插|三點內插|七點內插|三點內插|外插----
        if (index_z0 < 3) {
            k_0 = 3; /*兩點外插出去*/
        } else if (index_z0 > NZ6-4) {
            k_0 = NZ6-5;
        } else if (index_z0 < 6) {
            k_0 = (int)round(index_z0);
        } else if (index_z0 > NZ6-7) {
            k_0 = (int)ceil(index_z0) - 6; // ceil - 2 - 4 -> ceil -6
        } else {
            k_0 = (int)floor(index_z0) - 3; // 正確地以 index 減 3 作為起點
        }

        if (index_z1 < 3) {
            k_1 = 3;
        } else if (index_z1 > NZ6-4) {
            k_1 = NZ6-5;
        } else if (index_z1 < 6) {
            k_1 = (int)round(index_z1);
        } else if (index_z1 > NZ6-7) {
            k_1 = (int)ceil(index_z1) - 6;
        } else {
            k_1 = (int)floor(index_z1) - 3;
        }

        if (index_z2 < 3) {
            k_2 = 3;
        } else if (index_z2 > NZ6-4) {
            k_2 = NZ6-5;
        } else if (index_z2 < 6) {
            k_2 = (int)round(index_z2);
        } else if (index_z2 > NZ6-7) {
            k_2 = (int)ceil(index_z2) - 6;
        } else {
            k_2 = (int)floor(index_z2) - 3; // -2為真正的起點 -4為統一往上編號所需
        }
        //寫入每個idx_xi的三個(Y座標) 起始內插成員Z方向起始點編號
        //使用原始 j_store 作為儲存索引，確保不會因週期映射而錯位
        cell_z[NZ6*(j_store)+k + 0 * NY6*NZ6] = k_0;
        cell_z[NZ6*j_store+k   + 1 * NY6*NZ6] = k_1;
        cell_z[NZ6*(j_store)+k + 2 * NY6*NZ6] = k_2;
        //相鄰y列上的Z座標起始編號 寫入完成 !!
    }
}

//6.
/**
 * @brief 產生六階 Lagrange 插值預配置權重陣列
 * * 計算並儲存七個插值節點的 Lagrange 權重係數
 * @param[out] Para 二維權重陣列 [7][now]，儲存七個節點的權重
 * @param[in] Position 目標插值位置
 * @param[in] phy 物理座標陣列指標
 * @param[in] now 當前計算點編號
 * @param[in] start 起始節點索引
 * ! 確保 Para 陣列已正確分配記憶體
 * @see Lagrange_6th()
 */
//降階版本
void GetParameter_2nd(
    double *Para_h[3],      double Position,
    double *Pos,            int i,              int n  )
{
    Para_h[0][i] = Lagrange_2nd(Position, Pos[n],   Pos[n+1], Pos[n+2]);
    Para_h[1][i] = Lagrange_2nd(Position, Pos[n+1], Pos[n],   Pos[n+2]);
    Para_h[2][i] = Lagrange_2nd(Position, Pos[n+2], Pos[n],   Pos[n+1]);
}
void GetParameter_2ndlow(double** XiPara , double pos_z ,  double* RelationXi , int r , int index_xi){
    const int layer_stride = NY6 * NZ6;
    const int base = index_xi + r * layer_stride;
    XiPara[0][base] = Lagrange_2nd(pos_z, RelationXi[0],  RelationXi[1],  RelationXi[2] ); 
    XiPara[1][base] = Lagrange_2nd(pos_z, RelationXi[1],  RelationXi[0],  RelationXi[2] ); 
    XiPara[2][base] = Lagrange_2nd(pos_z, RelationXi[2],  RelationXi[0],  RelationXi[1] );     
}
void GetParameter_2ndhigh(double** XiPara , double pos_z ,  double* RelationXi , int r , int index_xi){
    const int layer_stride = NY6 * NZ6;
    const int base = index_xi + r * layer_stride;
    XiPara[0][base] = Lagrange_2nd(pos_z, RelationXi[4],  RelationXi[5],  RelationXi[6] ); 
    XiPara[1][base] = Lagrange_2nd(pos_z, RelationXi[5],  RelationXi[4],  RelationXi[6] ); 
    XiPara[2][base] = Lagrange_2nd(pos_z, RelationXi[6],  RelationXi[4],  RelationXi[5] );     
}
void GetParameter_6th(
    double *Para_h[7],      double Position,
    double *Pos,            int i,              int n  )
{
    Para_h[0][i] = Lagrange_6th(Position, Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+5], Pos[n+6]);
    Para_h[1][i] = Lagrange_6th(Position, Pos[n+1], Pos[n],   Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+5], Pos[n+6]);
    Para_h[2][i] = Lagrange_6th(Position, Pos[n+2], Pos[n],   Pos[n+1], Pos[n+3], Pos[n+4], Pos[n+5], Pos[n+6]);
    Para_h[3][i] = Lagrange_6th(Position, Pos[n+3], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+4], Pos[n+5], Pos[n+6]);
    Para_h[4][i] = Lagrange_6th(Position, Pos[n+4], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+5], Pos[n+6]);
    Para_h[5][i] = Lagrange_6th(Position, Pos[n+5], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+6]);
    Para_h[6][i] = Lagrange_6th(Position, Pos[n+6], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+5]);
    
}
void GetParameter_6th2(double** XiPara , double pos_z ,  double* RelationXi , int r , int j , int k  , int index_z0){
    const int layer_stride = NY6 * NZ6;
    const int base = j*NZ6 + k + r * layer_stride;
    if (index_z0 < 3){
    //線性外插 
    XiPara[0][base] =  1-Extrapolation(pos_z , RelationXi[0] , RelationXi[1] ) ;
    XiPara[1][base] =  Extrapolation(pos_z , RelationXi[0] , RelationXi[1] ) ;
    for(int i = 2 ; i <=6 ; i++) {XiPara[i][base] = 0.0 ;} 
    }
    else if (index_z0 > NZ6-4){
    for(int i = 0 ; i <=4 ; i++) {XiPara[i][base] = 0.0 ;}     
    XiPara[5][base] =  Extrapolation(pos_z , RelationXi[6] , RelationXi[5] ) ;
    XiPara[6][base] =  1-Extrapolation(pos_z , RelationXi[6] , RelationXi[5] ) ;
    }  
    else if (index_z0 < 6){
    for(int i = 3 ; i <=6 ; i++) {XiPara[i][base] = 0.0 ;} 
    XiPara[0][base] =  Lagrange_2nd(pos_z, RelationXi[0],  RelationXi[1],  RelationXi[2] ); 
    XiPara[1][base] =  Lagrange_2nd(pos_z, RelationXi[1],  RelationXi[0],  RelationXi[2] ); 
    XiPara[2][base] =  Lagrange_2nd(pos_z, RelationXi[2],  RelationXi[0],  RelationXi[1] );
    } 
    else if (index_z0 > NZ6-7){
    for(int i = 0 ; i <=3 ; i++) {XiPara[i][base] = 0.0 ;} 
    XiPara[4][base] =  Lagrange_2nd(pos_z, RelationXi[4],  RelationXi[5],  RelationXi[6] ); 
    XiPara[5][base] =  Lagrange_2nd(pos_z, RelationXi[5],  RelationXi[4],  RelationXi[6] );
    XiPara[6][base] =  Lagrange_2nd(pos_z, RelationXi[6],  RelationXi[4],  RelationXi[5] );
    }
    else {
    XiPara[0][base] = Lagrange_6th(pos_z, RelationXi[0],  RelationXi[1],  RelationXi[2] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[1][base] = Lagrange_6th(pos_z, RelationXi[1],  RelationXi[0],  RelationXi[2] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[2][base] = Lagrange_6th(pos_z, RelationXi[2],  RelationXi[0],  RelationXi[1] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[3][base] = Lagrange_6th(pos_z, RelationXi[3],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[4][base] = Lagrange_6th(pos_z, RelationXi[4],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[5], RelationXi[6]); 
    XiPara[5][base] = Lagrange_6th(pos_z, RelationXi[5],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[4], RelationXi[6]); 
    XiPara[6][base] = Lagrange_6th(pos_z, RelationXi[6],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[4], RelationXi[5]);    
}}//pos_xi為換算過後的無因次化Z座標 

//7.0度去向邊界計算點
/**
 * @brief 判斷是否為左山丘的 +y 方向物理空間邊界計算點
 * * 考察任意物理空間計算點是否為曲面邊界的邊界點
 * ? 判斷標準：以曲面邊界為本體，往 +y 方向考察半格移動距離 (minSize)
 * ? 若移動後會進入山丘內部，則為邊界計算點
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 * @return true 為 +y 方向邊界計算點，可套用 BFL 邊界條件
 * @return false 非邊界計算點
 * ! 考察範圍：左下角方形區域 pos_y: [0, HillHalfWidth], pos_z: [0, 1]
 * @see Left_q_yPlus(), IsRightHill_Boundary_yMinus()
 */
bool IsLeftHill_Boundary_yPlus(double y , double z){
    // * 第一步先初步篩選範圍：左下角方形區域
    if( y < 0.0 || y > HillHalfWidth || z < HillFunction(y) || z > 1.0 ){
        return false; // ! 第一步篩選未通過
    }else{// ? 第二步判斷，往左移動會撞倒代表現在高度在移動後座標點之山丘下方
        if( z <= HillFunction(y- minSize) ){
            return true;
        }else{
            return false;
        }
    }
}

//8.180度去向邊界計算點
/**
 * @brief 判斷是否為右山丘的 -y 方向物理空間邊界計算點 
 * * 判斷往 -y 方向移動一個 minSize 後是否會進入右山丘內部
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 * @return true 為 -y 方向邊界計算點，可套用 BFL 邊界條件
 * @return false 非邊界計算點
 * ! 考察範圍：右下角方形區域 pos_y: [LY-HillHalfWidth, LY], pos_z: [0, 1]
 * @see Right_q_yMinus(), IsLeftHill_Boundary_yPlus()
 */
bool IsRightHill_Boundary_yMinus(double y , double z){
    // * 第一步先初步篩選範圍：右下角方形區域
    if( y < (LY - HillHalfWidth) || y > LY || z < HillFunction(y) || z > 1.0 ){
        return false; // ! 第一步篩選未通過
    }else{// ? 第二步判斷，往右移動會撞倒代表現在高度在移動後座標點之山丘下方 
        if( z <= HillFunction(y + minSize) ){
            return true;
        }else{
            return false;
        }
    }
}

//9.45度去向邊界計算點
/**
 * @brief 判斷是否為左山丘的 45 度斜向物理空間邊界計算點
 * * 判斷往 (-Y, -Z) 方向 (即 45 度斜下方) 移動後是否會進入左山丘內部
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 * @return true 為 45 度斜向邊界計算點，可套用 BFL 邊界條件
 * @return false 非邊界計算點
 * ! 移動距離為 minSize (Y 和 Z 方向各移動 minSize)
 * @see Left_q_Diagonal45(), IsRightHill_Boundary_Diagonal135()
 */
bool IsLeftHill_Boundary_Diagonal45(double y , double z){
    // * 第一步先初步篩選範圍：左下角方形區域
    if( y < 0.0 || y > HillHalfWidth || z < HillFunction(y) || z > 1.0 ){
        return false; // ! 第一步篩選未通過
    }else{
        // * 移動後的位置
        double y_new = y - minSize;
        double z_new = z - minSize;
        // 判斷移動後是否進入山丘內部 (z_new <= 曲面高度)
        if( z_new <= HillFunction(y_new) ){
            return true;
        }else{
            return false;
        }
    }
}

//10.135度去向邊界計算點*
/**
 * @brief 判斷是否為右山丘的 135 度斜向物理空間邊界計算點 
 * * 判斷往 (+Y, -Z) 方向 (即 135 度斜下方) 移動後是否會進入右山丘內部
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 * @return true 為 135 斜向邊界計算點，可套用 BFL 邊界條件
 * @return false 非邊界計算點
 * ! 移動距離為 minSize (Y 和 Z 方向各移動 minSize)
 * @see Right_q_Diagonal135(), IsLeftHill_Boundary_Diagonal45()
 */
bool IsRightHill_Boundary_Diagonal135(double y , double z){
    // * 第一步先初步篩選範圍：右下角方形區域
    if( y < (LY - HillHalfWidth) || y > LY || z < HillFunction(y) || z > 1.0 ){
        return false; // ! 第一步篩選未通過
    }else{
        // * 移動後的位置
        double y_new = y + minSize;
        double z_new = z - minSize;
        // 判斷移動後是否進入山丘內部 (z_new <= 曲面高度)
        if( z_new <= HillFunction(y_new) ){
            return true;
        }else{
            return false;
        }
    }
}
//11.-45度去向邊界計算點
/**
 * @brief 判斷是否為左山丘的 -45 度斜向物理空間邊界計算點
 *
 * * 判斷往 (-Y, +Z) 方向 (即 -45 度斜上方) 移動後是否會進入左山丘內部
 * * 對應到格子速度方向 (+Y, -Z) 的「來源點」是否落在固體內（山丘斜率 > 1 時可能發生）
 *
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 *
 * @return true 為 -45 度斜向邊界計算點，可套用 BFL 邊界條件
 * @return false 非邊界計算點
 *
 * ! 移動距離為 minSize (Y 和 Z 方向各移動 minSize)
 * @see IsRightHill_Boundary_DiagonalMinus135()
 */
bool IsLeftHill_Boundary_DiagonalMinus45(double y , double z){
    // * 第一步先初步篩選範圍：左下角方形區域
    if( y < 0.0 || y > HillHalfWidth || z < HillFunction(y) || z > 1.0 ){
        return false; // ! 第一步篩選未通過
    }else{
        // * 移動後的位置
        const double y_new = y - minSize;
        const double z_new = z + minSize;
        // 判斷移動後是否進入山丘內部 (z_new <= 曲面高度)
        if( z_new <= HillFunction(y_new) ){
            return true;
        }else{
            return false;
        }
    }
}

//12.-135度去向邊界計算點
/**
 * @brief 判斷是否為右山丘的 -135 度斜向物理空間邊界計算點
 *
 * * 判斷往 (+Y, +Z) 方向 (即 -135 度斜上方) 移動後是否會進入右山丘內部
 * * 對應到格子速度方向 (-Y, -Z) 的「來源點」是否落在固體內（山丘斜率 > 1 時可能發生）
 *
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 *
 * @return true 為 -135 度斜向邊界計算點，可套用 BFL 邊界條件
 * @return false 非邊界計算點
 *
 * ! 移動距離為 minSize (Y 和 Z 方向各移動 minSize)
 * @see IsLeftHill_Boundary_DiagonalMinus45()
 */
bool IsRightHill_Boundary_DiagonalMinus135(double y , double z){
    // * 第一步先初步篩選範圍：右下角方形區域
    if( y < (LY - HillHalfWidth) || y > LY || z < HillFunction(y) || z > 1.0 ){
        return false; // ! 第一步篩選未通過
    }else{
        // * 移動後的位置
        const double y_new = y + minSize;
        const double z_new = z + minSize;
        // 判斷移動後是否進入山丘內部 (z_new <= 曲面高度)
        if( z_new <= HillFunction(y_new) ){
            return true;
        }else{
            return false;
        }
    }
}

//13.
/**
 * @brief 計算左山丘 +y 方向的 BFL 邊界條件 q 值
 * * 計算 +y 方向邊界計算點到曲面壁面的無因次距離 q
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 * @return q 值 (0 < q < 1)，若非邊界點則返回 -1.0
 * ! q = |y - y_wall| / minSize，用於 BFL 插值邊界條件
 * @see IsLeftHill_Boundary_yPlus(), HillFunction_Inverse_Left()
 */
double Left_q_yPlus(double y , double z){
    if(!IsLeftHill_Boundary_yPlus(y, z)){
        return -1.0; // ! 非邊界點返回-1
    }// * 返回負值作為錯誤之邊界條件實施判據，因為距離沒有負值
    // * 計算 q 值：計算點到曲面的 Y 方向距離 / 格子大小
    return fabs(y - HillFunction_Inverse_Left(z)) / minSize;
}

//14.
/**
 * @brief 計算右山丘 -y 方向的 BFL 邊界條件 q 值
 * * 計算 -y 方向邊界計算點到曲面壁面的無因次距離 q
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 * @return q 值 (0 < q < 1)，若非邊界點則返回 -1.0
 * ! q = |y - y_wall| / minSize，用於 BFL 插值邊界條件
 * @see IsRightHill_Boundary_yMinus(), HillFunction_Inverse_Right()
 */
double Right_q_yMinus(double y , double z){
    if(!IsRightHill_Boundary_yMinus(y, z)){
        return -1.0; // ! 非邊界點返回-1
    }
    // * 計算 q 值：計算點到曲面的 Y 方向距離 / 格子大小
    return fabs(y - HillFunction_Inverse_Right(z)) / minSize;
}

//15.
/**
 * @brief 計算左山丘 45 度斜向的 BFL 邊界條件 q 值
 * * 使用二分法求解 45 度射線與山丘曲面的交點，計算無因次距離 q
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 * @return q 值 (0 < q < 1)，若非邊界點則返回 -1.0
 * ! 45 度射線方程：z = y + (z0 - y0)
 * @see IsLeftHill_Boundary_Diagonal45()
 */
double Left_q_Diagonal45(double y , double z){
    if(!IsLeftHill_Boundary_Diagonal45(y , z)){
        return -1.0 ; // ! 非邊界點返回-1
    }
    double y_target ;
    // * 利用區間搜尋法尋找y_target ; 
    double y_temp[2] = {y-minSize ,y} ; 
    //y[0] = y_minSize ; y[1] = y ; 
    double y_middle = 0.0 ;
    double z_middle = 0.0 ; //在45度射線上相對應y_middle的z座標
    double Length_z = 0.0 ;
    do{
        y_middle = (y_temp[0] + y_temp[1]) / 2.0;
        z_middle = y_middle + (z-y) ; 
        Length_z = fabs(z_middle - HillFunction(y_middle));
        if(HillFunction(y_middle) > z_middle){
            //將y_middle往右移動
            y_temp[0] = y_middle;
        }else{
            //將y_middle往左移動
            y_temp[1] = y_middle;
            
        }
    }while(Length_z >= 1e-12) ;
    y_target = y_middle ;
    return fabs(y - y_target) / minSize ;
}

//16.
/**
 * @brief 計算右山丘 135 度斜向的 BFL 邊界條件 q 值 
 * * 使用二分法求解 135 度射線與山丘曲面的交點，計算無因次距離 q
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 * @return q 值 (0 < q < 1)，若非邊界點則返回 -1.0
 * ! 135 度射線方程：z = -y + (z0 + y0)
 * @see IsRightHill_Boundary_Diagonal135()
 */
double Right_q_Diagonal135(double y , double z){
    if(!IsRightHill_Boundary_Diagonal135(y , z)){
        return -1.0 ; // ! 非邊界點返回-1
    }
    double y_target ;
    // * 利用區間搜尋法尋找y_target ;
    double y_temp[2] = {y , y+minSize} ;
    //y[0] = y ; y[1] = y+minSize ;
    double y_middle = 0.0 ;
    double z_middle = 0.0 ; //在135度射線上相對應y_middle的z座標
    double Length_z = 0.0 ;
    do{
        y_middle = (y_temp[0] + y_temp[1]) / 2.0;
        z_middle = -y_middle + (z+y) ; // 135度射線：z = -y + (z0+y0)
        Length_z = fabs(z_middle - HillFunction(y_middle));
        if(HillFunction(y_middle) > z_middle){
            //將y_middle往左移動
            y_temp[1] = y_middle;
        }else{
            //將y_middle往右移動
            y_temp[0] = y_middle;
        }
    }while(Length_z >= 1e-12) ;
    y_target = y_middle ;
    return fabs(y - y_target) / minSize ;
}

//17.
/**
 * @brief 計算左山丘 -45 度斜向的 BFL 邊界條件 q 值
 * * 使用二分法求解 -45 度射線與山丘曲面的交點，計算無因次距離 q
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 * @return q 值 (0 < q < 1)，若非邊界點則返回 -1.0
 * ! -45 度射線方程：z = -y + (z0 + y0)
 * @see IsLeftHill_Boundary_DiagonalMinus45()
 */
double Left_q_DiagonalMinus45(double y , double z){
    if(!IsLeftHill_Boundary_DiagonalMinus45(y , z)){
        return -1.0 ; // ! 非邊界點返回-1
    }
    double y_target ;
    // * 利用區間搜尋法尋找y_target ;
    double y_temp[2] = {y-minSize ,y} ;
    //y[0] = y-minSize ; y[1] = y ;
    double y_middle = 0.0 ;
    double z_middle = 0.0 ; //在-45度射線上相對應y_middle的z座標
    double Length_z = 0.0 ;
    do{
        y_middle = (y_temp[0] + y_temp[1]) / 2.0;
        z_middle = -y_middle + (z+y) ; // -45度射線：z = -y + (z0+y0)
        Length_z = fabs(z_middle - HillFunction(y_middle));
        if(HillFunction(y_middle) > z_middle){
            //將y_middle往左移動（左丘：y越小，HillFunction越大）
            y_temp[1] = y_middle;
        }else{
            //將y_middle往右移動
            y_temp[0] = y_middle;
        }
    }while(Length_z >= 1e-12) ;
    y_target = y_middle ;
    return fabs(y - y_target) / minSize ;
}

//18.
/**
 * @brief 計算右山丘 -135 度斜向的 BFL 邊界條件 q 值
 * * 使用二分法求解 -135 度射線與山丘曲面的交點，計算無因次距離 q
 * @param y Y 座標 (流向)
 * @param z Z 座標 (法向)
 * @return q 值 (0 < q < 1)，若非邊界點則返回 -1.0
 * ! -135 度射線方程：z = y + (z0 - y0)
 * @see IsRightHill_Boundary_DiagonalMinus135()
 */
double Right_q_DiagonalMinus135(double y , double z){
    if(!IsRightHill_Boundary_DiagonalMinus135(y , z)){
        return -1.0 ; // ! 非邊界點返回-1
    }
    double y_target ;
    // * 利用區間搜尋法尋找y_target ;
    double y_temp[2] = {y , y+minSize} ;
    //y[0] = y ; y[1] = y+minSize ;
    double y_middle = 0.0 ;
    double z_middle = 0.0 ; //在-135度射線上相對應y_middle的z座標
    double Length_z = 0.0 ;
    do{
        y_middle = (y_temp[0] + y_temp[1]) / 2.0;
        z_middle = y_middle + (z-y) ; // -135度射線：z = y + (z0-y0)
        Length_z = fabs(z_middle - HillFunction(y_middle));
        if(HillFunction(y_middle) > z_middle){
            //將y_middle往右移動（右丘：y越大，HillFunction越大）
            y_temp[0] = y_middle;
        }else{
            //將y_middle往左移動
            y_temp[1] = y_middle;
        }
    }while(Length_z >= 1e-12) ;
    y_target = y_middle ;
    return fabs(y - y_target) / minSize ;
}


#endif

