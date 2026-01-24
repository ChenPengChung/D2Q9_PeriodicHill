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
#include "variables.h"
using namespace std;

//1.
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

//2.
/**
 * @brief 計算非均勻網格伸縮參數 a
 * * 使用二分法求解伸縮參數，使最小網格間距等於 minSize
 * @return 伸縮參數 a (0 < a < 1)
 * ! minSize 定義為以 LZ-(y=0.0 下的山坡高度) 為總長度均勻切割下的網格大小 × 0.6
 * TODO: 考慮改用 Newton-Raphson 加速收斂
 * @see tanhFunction
 */
double GetNonuniParameter() {
    double total = LZ - 0.0  - minSize;
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


//3.
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
//降階版本
double Lagrange_2nd(double pos , double x_i , double x1 , double x2 ){
    double Lagrange = (pos - x1)/(x_i - x1)*(pos - x2)/(x_i - x2);
    return Lagrange; 
}

double Lagrange_6th(double pos , double x_i , double x1 , double x2 , double x3 , double x4 , double x5 , double x6){
    double Lagrange = (pos - x1)/(x_i - x1)*(pos - x2)/(x_i - x2)*(pos - x3)/(x_i - x3)*(pos - x4)/(x_i - x4)*(pos - x5)/(x_i - x5)*(pos - x6)/(x_i - x6);
    return Lagrange; 
}

//給我一個編號，產生該Y值所對應的七個無因次化座標
void RelationXi(int k , double L , double MinSize , double a , int N , double* RelationXi ){
    int j = k-3 ; 
    if (j <3)  j =3 ;
    if (j > N-4) j = N-4 ;
    RelationXi[0] = tanhFunction( L , MinSize , a, j-3 , N) - MinSize/2.0;
    RelationXi[1] = tanhFunction( L , MinSize , a, j-2 , N) - MinSize/2.0;
    RelationXi[2] = tanhFunction( L , MinSize , a, j-1 , N) - MinSize/2.0;
    RelationXi[3] = tanhFunction( L , MinSize , a, j , N) - MinSize/2.0;
    RelationXi[4] = tanhFunction( L , MinSize , a, j+1 , N) - MinSize/2.0;
    RelationXi[5] = tanhFunction( L , MinSize , a, j+2 , N) - MinSize/2.0;
    RelationXi[6] = tanhFunction( L , MinSize , a, j+3 , N) - MinSize/2.0;
}

//4.
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
void GetParameter_6th2(double** XiPara , double pos_z ,  double* RelationXi , int r , int index_xi){
    const int layer_stride = NY6 * NZ6;
    const int base = index_xi + r * layer_stride;
    XiPara[0][base] = Lagrange_6th(pos_z, RelationXi[0],  RelationXi[1],  RelationXi[2] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[1][base] = Lagrange_6th(pos_z, RelationXi[1],  RelationXi[0],  RelationXi[2] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[2][base] = Lagrange_6th(pos_z, RelationXi[2],  RelationXi[0],  RelationXi[1] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[3][base] = Lagrange_6th(pos_z, RelationXi[3],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[4][base] = Lagrange_6th(pos_z, RelationXi[4],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[5], RelationXi[6]); 
    XiPara[5][base] = Lagrange_6th(pos_z, RelationXi[5],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[4], RelationXi[6]); 
    XiPara[6][base] = Lagrange_6th(pos_z, RelationXi[6],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[4], RelationXi[5]);    
}//pos_xi為換算過後的無因次化Z座標 




#endif

