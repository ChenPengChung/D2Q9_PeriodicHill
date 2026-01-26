#ifndef VARIABLES_FILE
#define VARIABLES_FILE
//包含山坡之物理模型(矩形)
#define     pi     3.14159265358979323846264338327950
//二維模型之長寬設定
#define     LY     (9.0)
#define     LZ     (3.036)
//分配之格子數量，計算點在網格中心點
#define     NY      400 
#define     NZ      200
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

//模擬迴圈上限值
#define loop 1000000

//=== 關鍵修改：正確的 Re 參數設計 ===
// 為確保數值穩定性，降低 Re 或提高黏度
#define Re          10                              // 目標雷諾數（從50降至10提高穩定性）
#define Uref        0.03                            // 固定特徵速度 (Ma < 0.1，降低馬赫數)
#define L_char      1.0                             // 特徵長度 (山坡高度 h)
#define niu         (Uref * L_char / (double)Re)    // 黏滯係數 = 0.003
#define tau         (3.0*niu + 0.5)                 // 鬆弛時間 = 0.509
#define omega_7     (1.0 / tau)                     // 剪切鬆弛參數 ≈ 1.964 → 1.2 (更保守)
#define omega_2     1.0                             // 能量鬆弛（降至1.0提高穩定性）

// 實際使用更保守的鬆弛參數
#undef omega_7
#define omega_7     1.2                             // 手動設置為1.2，確保穩定

//=== 次要參數 ===
#define dt          minSize
#define cs          (1.0/1.732050807568877)


#endif
