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
//雷諾數
#define Re 5 
//模擬迴圈上限值
#define loop 1000000
/**********Secondary Parameter********************/
#define dt (minSize) 
#define cs (1.0/1.732050807568877)
//Parameters of periodic hills Using MRT operator    
#define     omega_2     1.2
#define     omega_7     1.2
#define     niu         (1/3.)*(1/omega_7 - 0.5 )
//以山坡高度流場的特徵長度，計算特徵速度
#define     Uref        (Re*niu)
#endif
