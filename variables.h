#ifndef VARIABLES_FILE
#define VARIABLES_FILE
//包含山坡之物理模型(矩形)
#define     pi     3.14159265358979323846264338327950
//二維模型之長寬設定
#define     LY     (9.0)
#define     LZ     (3.036)
//分配之格子數量，計算點在網格中心點
#define     NY      512 
#define     NZ      256
//加Buffer之計算區域下的總體網格數量
//Stream-Wise方向不考慮GPU切割
#define     NY6    (NY+7)
#define     NZ6    (NZ+6)
//設定Lattice大小，定義為最小物理網格大小的0.6倍，作為一格粒子的移動距離
//CFL number 為 速度/格子速度 ，在此為比較 晶格大小/最小物理網格大小
//所以最小物理大小的CFL才是晶格大小，為minSize
#define     CFL                 0.8  // Re=1000 時用 0.85 更安全
#define     minSize             ((LZ-1.0)/(NZ6-6)*CFL)
//非均勻網格之判斷式
//1 : Yes,  0 : No
#define     Uniform_In_Ydir     1
#define     Uniform_In_Zdir     0
//定義無因次化長度上限
#define LXi (10.0)

//模擬迴圈上限值
#define     loop        10000000
#define      Re         700                            // 目標雷諾數
#define     tau         0.6833                            // 提高 tau 以增加穩定性（原 1.7 → 1.95）
#define     niu         ((tau-0.5)/3.0*dt)
#define     Uref        (Re*niu)
#define     L_char      1.0                             // 特徵長度 (山坡高度 h)
#define     omega_7     (1.0 / tau)                    // 剪切鬆弛參數 ≈ 0.513（更保守）
#define     omega_2     0.1                            // 能量鬆弛（提高至 0.5，原 0.001）
// 統一邊界定義：streaming 和 interpolation 必須一致
// streaming_lower/upper: evolution.h 中用於判斷是否用 streaming 代替插值
// interpolation_lower/upper: initialization.h 中用於判斷使用幾點插值
// 重要：streaming_lower 必須 >= interpolation_lower 以確保一致性
#define     interpolation_lower  (6)                // 七點內插下界 (index_z < 此值用三點插值)
#define     interpolation_upper  (NZ6-7)            // 七點內插上界 (index_z > 此值用三點插值)

//=== 動態 Streaming 邊界參數（分階段漸進式擴大解析層）===//
// 初始值（保守，更大的 streaming 區域）
#define     streaming_lower_init     (25)            // 初始下界 (k <= 50 用 streaming)
#define     streaming_upper_init     (NZ6-19)        // 初始上界 (k >= NZ6-51 用 streaming)

// === 第一階段：開放七點插值區 (streaming → interpolation_lower) ===
#define     streaming_lower_phase1   (10)  // 第一階段目標: 25
#define     streaming_upper_phase1   (NZ6-12)  // 第一階段目標: NZ6-26
#define     phase1_start_time        (0)             // 第一階段開始
#define     phase1_end_time          (100000)        // 第一階段結束

// === 第二階段：開放三點插值緩衝區 (interpolation_lower → target) ===
#define     streaming_lower_target   (6)            // 最終目標下界
#define     streaming_upper_target   (NZ6-7)         // 最終目標上界
#define     phase2_start_time        (10000)        // 第二階段開始
#define     phase2_end_time          (20000)        // 第二階段結束

// 全域變數宣告（在 main.cpp 中定義）
extern int streaming_lower;  // 動態下界，由 UpdateStreamingBounds() 更新
extern int streaming_upper;  // 動態上界，由 UpdateStreamingBounds() 更新
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



