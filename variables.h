#ifndef VARIABLES_FILE
#define VARIABLES_FILE
//包含山坡之物理模型(矩形)
#define     pi     3.14159265358979323846264338327950
//二維模型之長寬設定
#define     LY     (9.0)
#define     LZ     (3.036)
//分配之格子數量，計算點在網格中心點
#define     NY      256 
#define     NZ      128
//加Buffer之計算區域下的總體網格數量
//Stream-Wise方向不考慮GPU切割
#define     NY6    (NY+7)
#define     NZ6    (NZ+6)
//設定Lattice大小，定義為最小物理網格大小的0.6倍，作為一格粒子的移動距離
//CFL number 為 速度/格子速度 ，在此為比較 晶格大小/最小物理網格大小
//所以最小物理大小的CFL才是晶格大小，為minSize
#define     CFL                 0.2  // Re=1000 時用 0.85 更安全
#define     minSize             ((LZ-1.0)/(NZ6-6)*CFL)
//非均勻網格之判斷式
//1 : Yes,  0 : No
#define     Uniform_In_Ydir     1
#define     Uniform_In_Zdir     0
//定義無因次化長度上限
#define LXi (10.0)

//模擬迴圈上限值
#define     loop        10000000
#define      Re         100                            // 目標雷諾數
#define     tau         1.2                            // 提高 tau 以增加穩定性（原 1.7 → 1.95）
#define     niu         ((tau-0.5)/3.0*dt)
#define     Uref        (Re*niu)
#define     L_char      1.0                             // 特徵長度 (山坡高度 h)
#define     omega_7     (1.0 / tau)                    // 剪切鬆弛參數 ≈ 0.513（更保守）
#define     omega_2     0.1                            // 能量鬆弛（提高至 0.5，原 0.001）
#define     interpolation_lower  (20)                // 七點內插下界 : 6 為最保守
#define     interpolation_upper  (NZ6-14)            // 七點內插上界 : NZ6-7 為最保守
#define     streaming_lower  (10)
#define     streaming_upper  (NZ6-11)
//=== 次要參數 ===/
#define     dt          minSize
#define     cs          (1.0/1.732050807568877)              // 晶格聲速 = 1/sqrt(3)

//=== Mach 數限制參數 ===//
// LBM 不可壓縮性要求：Ma < 0.3 (嚴格)，Ma < 0.1 (保守)
#define     Ma_max      0.3                                  // 最大允許 Mach 數
#define     Ma_warning  0.15                                 // 警告閾值
#define     U_max       (Ma_max * cs)                        // 最大允許速度 ≈ 0.173
#define     U_warning   (Ma_warning * cs)                    // 警告速度 ≈ 0.087

// 計算當前參數下的理論 Mach 數
// Ma_theoretical = Uref / cs = Uref * sqrt(3)
// 注意：若 Ma_theoretical > Ma_max，模擬很可能不穩定

#endif



