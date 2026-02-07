#ifndef VARIABLES_FILE
#define VARIABLES_FILE
//包含山坡之物理模型(矩形體)
#define     pi     3.14159265358979323846264338327950
//二維模型之長寬設定
#define     LY     (9.0)
#define     LZ     (3.036)
//分配之格子數量，計算點在網格中心點
#define     NY      128
#define     NZ      64
//加Buffer之計算區域下的總體網格數量
//Stream-Wise方向不考慮GPU切割
#define     NY6    (NY+7)
#define     NZ6    (NZ+6)
//設定Lattice大小，定義為最小物理網格大小的0.6倍，作為一格粒子的移動距離
//CFL number 為 速度/格子速度 ，在此為比較 晶格大小/最小物理網格大小
//所以最小物理大小的CFL才是晶格大小，為minSize
#define     CFL                 0.2
#define     minSize             ((LZ-1.0)/(NZ6-6)*CFL)
//非均勻網格之判斷式
//1 : Yes,  0 : No
#define     Uniform_In_Ydir     1
#define     Uniform_In_Zdir     0
//定義無因次化長度上限
#define LXi (10.0)
/**********Secondary Parameter********************/
#define dt (minSize) 
#define cs (1.0/1.732050807568877)
//=============================================================================
// MRT 鬆弛參數設定
// 
// LBM 穩定性限制：馬赫數 Ma = U/cs < 0.3，否則可壓縮性誤差過大
// 聲速 cs = 1/sqrt(3) ≈ 0.577
// 因此 Uref < 0.3 * 0.577 ≈ 0.173
//
// 使用較低的鬆弛參數 omega_7 來增加黏滯係數，提高穩定性
//=============================================================================
//模擬迴圈上限值
#define     loop        10000000
#define      Re         300                            // 目標雷諾數
#define     tau         0.6833                            // 提高 tau 以增加穩定性（原 1.7 → 1.95）
#define     niu         ((tau-0.5)/3.0*dt)
#define     Uref        (Re*niu)
#define     L_char      1.0                             // 特徵長度 (山坡高度 h)
#define     omega_7     (1 / tau)                    // 剪切鬆弛參數 ≈ 0.513（更保守）
#define     omega_2     1.19                            // 能量鬆弛（提高至 0.5，原 0.001）

#endif
