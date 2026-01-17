# D2Q9 Periodic Hill 專案開發記錄

## 專案概述
使用 Lattice Boltzmann Method (LBM) 模擬 Periodic Hill 流場
- 模型：D2Q9 + MRT 碰撞算子
- 雷諾數：Re = 50
- 網格：NY=128, NZ=64 (含 buffer 為 NY6=135, NZ6=70)
- 非均勻網格：Z 方向使用 tanh 分佈

---

## 2025-01-14 進度檢查

### 已完成模組

| 檔案 | 狀態 | 說明 |
|------|------|------|
| `variables.h` | ✅ | 物理參數與網格配置 |
| `model.h` | ✅ | 山坡幾何函數 HillFunction() |
| `MRT_Matrix.h` | ✅ | MRT 矩陣定義（2025-01-15 已修復） |
| `MRT_Process.h` | ⚠️ | m_vector 巨集，邏輯待補充 |
| `initializationTool.h` | ✅ | 非均勻網格、Lagrange 插值 |
| `InterpolationHillISLBM.h` | ✅ | 7點插值巨集（經檢查無語法錯誤） |

### 已修復問題 (2025-01-15)

1. **MRT_Matrix.h** - 二維陣列各列缺少逗號分隔符
   - `Matrix` 巨集：各列 `{...}` 後補上逗號
   - `Matrix_Inverse` 巨集：各列 `{...}` 後補上逗號

2. **InterpolationHillISLBM.h** - 經檢查第50行語法正確，無需修改

### 待實作項目

- [ ] 主程式 main.cu
- [ ] MRT 碰撞運算完整邏輯
- [ ] Streaming 步驟
- [ ] 邊界條件處理 (週期性邊界、壁面)
- [ ] 初始化流場
- [ ] 結果輸出

---

## 開發筆記

- `y_h` 為 main.cu 中的全域變數，initialization.h 包含後可直接讀寫
- 包含順序需在全域定義之後

### 待定義全域變數 (main.cu)

以下變數在 `initialization.h` 中使用，需在 `main.cu` 中宣告：

```cpp
double rho[NY6 * NZ6];           // 密度場
double v[NY6 * NZ6];             // Y方向速度
double w[NY6 * NZ6];             // Z方向速度
double f[9][NY6 * NZ6];          // 分佈函數 (D2Q9)
double Force[2];                 // 外力項 (Fy, Fz)
double xi_h[NZ6];                // 無因次化Z座標
```

**注意**：`#include "initialization.h"` 必須放在這些全域變數宣告之後

---

## 2025-01-14 BFL 邊界條件判斷式

### 山丘幾何結構
- 左半丘：`y ∈ [0, 54/28]` ≈ `[0, 1.9286]`
- 右半丘：`y ∈ [LY - 54/28, LY]` ≈ `[7.07, 9.0]`
- 週期性邊界連接後形成完整山丘

### 新增判斷函數 (initializationTool.h)

```c
//5. 判斷是否為 +Y 方向邊界計算點 (針對左半丘)
#define HillHalfWidth (54.0/28.0)

#define IsBoundary_LeftHill_PlusY(y, z) \
(                                        \
    (y) >= 0.0 &&                        \
    (y) + minSize <= HillHalfWidth &&    \
    (z) > HillFunction(y) &&             \
    (z) <= HillFunction((y) + minSize)   \
)

//6. 判斷是否為 -Y 方向邊界計算點 (針對右半丘)
#define IsBoundary_RightHill_MinusY(y, z) \
(                                          \
    (y) <= LY &&                           \
    (y) - minSize >= (LY - HillHalfWidth) && \
    (z) > HillFunction(y) &&               \
    (z) <= HillFunction((y) - minSize)     \
)
```

### 邏輯說明
- 左半丘：計算點在流體區，往 +Y 方向移動 minSize 會碰到曲面
- 右半丘：計算點在流體區，往 -Y 方向移動 minSize 會碰到曲面

---

## 2026-01-16 完整程式碼審查

### 程式碼結構總覽

專案包含 7 個標頭檔：
1. `variables.h` - 物理參數定義
2. `model.h` - 山坡幾何模型
3. `MRT_Matrix.h` - MRT 轉換矩陣
4. `MRT_Process.h` - MRT 碰撞處理
5. `initializationTool.h` - 初始化工具函數
6. `InterpolationHillISLBM.h` - ISLBM 插值巨集
7. `initialization.h` - 包裝檔

### 各模組詳細分析

#### ✅ variables.h (正常)
- 定義物理參數：Re=50, LY=9.0, LZ=3.036
- 網格設定：NY=128, NZ=64
- Buffer 網格：NY6=135, NZ6=70
- CFL=0.6, minSize 計算正確
- tau=0.6833, omega_2 和 omega_7 定義正確
- Uref = Re*niu

#### ✅ model.h (正常)
- HillFunction() 實作完整
- 左半丘：6 段三次多項式 (y ∈ [0, 54/28])
- 右半丘：6 段三次多項式 (y ∈ [LY-54/28, LY])
- 週期性邊界處理正確

#### ⚠️ MRT_Matrix.h (需檢查)
- **問題 1**：第 10 行定義 `s1 = omega_1` 但 `omega_1` 未在 variables.h 中定義
  - 根據註解應該是 free parameter，建議改為 `s1 = 1.0`
- **問題 2**：第 9、11、14、16 行缺少分號
- Matrix M[9][9] 定義正確（已於 2025-01-15 修復）
- Matrix_Inverse M_I[9][9] 定義正確（已於 2025-01-15 修復）

#### ❌ MRT_Process.h (嚴重錯誤)
- **問題 1**：陣列索引越界
  - 各行使用 `M[i][9]` 但 D2Q9 只有 9 個方向 (索引 0-8)
  - 正確應為 `M[i][8]`
- **問題 2**：各行缺少分號
- **問題 3**：變數命名錯誤
  - 使用 F9_in 但 D2Q9 只有 F0-F8
  - 應改為 F8_in
- **問題 4**：m_vector 巨集邏輯不完整
  - 僅計算 moment，未實作碰撞鬆弛和逆轉換

#### ✅ initializationTool.h (功能完整，有小錯誤)
- **問題**：第 6 行 `HillHalfWidth` 定義錯誤
  - 現在：`#define HillHalfWidth 9.0*(54.0/28.0)` ≈ 17.36
  - 正確：`#define HillHalfWidth (54.0/28.0)` ≈ 1.93

**已實作功能**：
- ✅ HillFunction_Inverse_Left/Right (二分搜尋法)
- ✅ tanhFunction 巨集 (非均勻網格)
- ✅ GetNonuniParameter() (計算伸縮參數 a，使物理空間計算點的最小距離等於晶格大小 minSize)
  - **目的**：調整 tanh 非均勻網格的伸縮參數 a，確保最小網格間距 = minSize (晶格大小)
  - **方法**：二分法搜尋 a ∈ (0.1, 1.0)
  - **關鍵判斷**：
    ```c
    x_temp[0] = tanhFunction(total, minSize, a_mid, 0, (NZ6-7));
    x_temp[1] = tanhFunction(total, minSize, a_mid, 1, (NZ6-7));
    dx = x_temp[1] - x_temp[0];  // 第 0 與第 1 個計算點的間距
    if( dx - minSize >= 0.0 )    // dx >= minSize 則 a 可以更大
        a_temp[0] = a_mid;
    else
        a_temp[1] = a_mid;
    // 收斂條件：|dx - minSize| < 1e-14
    ```
  - **物理意義**：讓壁面附近(k=0)的網格間距精確等於晶格大小，確保 ISLBM 插值精度
- ✅ Lagrange_6th() (6 階插值)
- ✅ GetParameter_6th() (預計算插值權重)
- ✅ IsLeftHill_Boundary_yPlus() (左丘 +Y 邊界判斷)
- ✅ IsRightHill_Boundary_yMinus() (右丘 -Y 邊界判斷)
- ✅ IsLeftHill_Boundary_Diagonal45() (左丘 45° 邊界判斷)
- ✅ IsRightHill_Boundary_Diagonal135() (右丘 135° 邊界判斷)
- ✅ Left_q_yPlus() (左丘 +Y BFL 距離)
- ✅ Right_q_yMinus() (右丘 -Y BFL 距離)
- ✅ Left_q_Diagonal45() (左丘 45° BFL 距離)
- ✅ Right_q_Diagonal135() (右丘 135° BFL 距離)

#### ✅ InterpolationHillISLBM.h (正常)
- 定義 7 點插值巨集 Intrpl7()
- F0 到 F8 的 2D 插值巨集
- 支援非均勻網格的 Lagrange 插值
- 語法檢查無誤（2025-01-15 已確認）

#### ✅ initialization.h (正常)
- 簡單包裝檔，包含 initializationTool.h

### 需修復的錯誤清單

#### 高優先級 (阻斷性錯誤)
1. **MRT_Process.h 陣列越界** - 必須立即修復
2. **MRT_Process.h 缺少分號** - 編譯錯誤
3. **initializationTool.h HillHalfWidth 錯誤** - 影響邊界判斷

#### 中優先級
4. **MRT_Matrix.h omega_1 未定義** - 可能導致編譯錯誤
5. **MRT_Matrix.h 缺少分號** - 編譯錯誤

#### 待實作功能
- [ ] MRT 碰撞運算完整邏輯 (鬆弛 + 逆轉換)
- [ ] Streaming 步驟
- [ ] 邊界條件處理邏輯
- [ ] 主程式 main.cu
- [ ] 初始化流場設定
- [ ] 結果輸出功能

### D2Q9 模型速度方向

```
方向編號與速度向量：
F0: (0, 0)
F1: (1, 0)   F2: (0, 1)
F3: (-1, 0)  F4: (0, -1)
F5: (1, 1)   F6: (-1, 1)
F7: (-1, -1) F8: (1, -1)
```

---

## 2026-01-16 錯誤修復完成

### 已修復的高優先級錯誤

#### 1. initializationTool.h:20 - HillHalfWidth 定義錯誤 ✅
- **問題**：`#define HillHalfWidth 9.0*(54.0/28.0)` 值為 17.36
- **修復**：改為 `#define HillHalfWidth (54.0/28.0)` ≈ 1.9286
- **影響**：修正所有 BFL 邊界判斷函數的計算範圍

#### 2. MRT_Process.h - 陣列越界錯誤 ✅
- **問題**：所有行使用 `M[i][9]` 和 `F9_in`，但 D2Q9 模型只有 0-8 索引
- **修復**：將 `M[i][9]` 改為 `M[i][8]`，`F9_in` 改為 `F8_in`
- **影響**：避免記憶體越界存取

#### 3. MRT_Process.h - 缺少分號 ✅
- **問題**：m0-m8 各行定義後缺少分號
- **修復**：每行結尾加上分號 `;`
- **影響**：修正巨集語法錯誤

### 已修復的中優先級錯誤

#### 4. MRT_Matrix.h:10 - omega_1 未定義 ✅
- **問題**：`s1 = omega_1` 但 omega_1 未在 variables.h 定義
- **修復**：改為 `s1 = omega_2` (根據 variables.h:36 定義)
- **影響**：使用正確的鬆弛參數

#### 5. MRT_Matrix.h - 缺少分號 ✅
- **問題**：Relaxation 巨集各行缺少分號
- **修復**：所有行結尾補上分號 `;`
- **影響**：修正巨集語法錯誤

### 修復總結

所有阻斷性編譯錯誤已修復，程式碼目前狀態：
- ✅ 語法錯誤：全部修復
- ✅ 陣列越界：已修正
- ✅ 未定義變數：已修正
- ⚠️ 功能完整性：MRT 碰撞邏輯仍待實作

---

## 2026-01-16 外力項實作討論

### 問題分析：Half-Way Correction Force

目前程式的力離散化**缺少 Half-Way (Guo) 修正**，存在以下問題：

#### 1. 缺少 Half-Way 修正項
**現狀**：力項只用常數係數，不含 `(1 - s/2)` 因子，也不含速度 `u` 相關項

**影響**：
- 力對速度的影響只有一階精度
- τ 偏離 1 或流速較大時會引入偏差
- 無法正確捕捉流體對外力的響應

#### 2. 力項分配不完整
**現狀**：只對含 Y 分量的速度方向加入外力
- D2Q9 模型中，只有 F1(+Y), F3(-Y), F5(+Y,+Z), F6(-Y,+Z), F7(-Y,-Z), F8(+Y,-Z) 加入力
- F2(+Z), F4(-Z) 沒有加入力項

**原因**：程式將 `Force[0]` 當成單一軸向的體力（沿 Y 方向），其他分量為 0

### Half-Way (Guo) 力模型修正方案

#### 修正 1：力增量 (Force Distribution)

在碰撞前，對每個方向的力項使用 Guo 公式：

```c
Fi = wi * (1 - sα/2) * [ (ei - u)/cs^2 + (ei·u) ei / cs^4 ] · F
```

**參數說明**：
- `wi`：D2Q9 權重係數
  - w0 = 4/9 (靜止方向)
  - w1,2,3,4 = 1/9 (直向：±Y, ±Z)
  - w5,6,7,8 = 1/36 (斜向：對角線)
- `sα`：對應 moment 的鬆弛率
  - BGK: 使用 `(1 - 1/(2τ))`
  - MRT: 使用對應矩的 `sα` 值
- `ei`：第 i 個速度方向
- `u`：宏觀速度
- `F`：體積力向量
- `cs = 1/√3`：晶格聲速

#### 修正 2：速度修正 (Velocity Correction)

計算宏觀速度時使用半步修正：

```c
u = (Σ ei fi + 0.5 * F * dt) / rho
```

或分兩步：
```c
u_raw = (Σ ei fi) / rho
u = u_raw + 0.5 * F * dt / rho
```

這樣 momentum 更新與力增量一致，保證二階精度。

### D2Q9 模型的力項分配

#### 速度方向與權重

```
方向  速度向量 (ey, ez)  權重 wi    是否含Y分量
F0    (0, 0)             4/9       否
F1    (1, 0)             1/9       是 (+Y)
F2    (0, 1)             1/9       否 (純+Z)
F3    (-1, 0)            1/9       是 (-Y)
F4    (0, -1)            1/9       否 (純-Z)
F5    (1, 1)             1/36      是 (+Y+Z)
F6    (-1, 1)            1/36      是 (-Y+Z)
F7    (-1, -1)           1/36      是 (-Y-Z)
F8    (1, -1)            1/36      是 (+Y-Z)
```

#### 目前實作問題

**單軸力假設**：
- 假設只有 Fy ≠ 0，Fz = 0
- 因此 F2, F4 (純 Z 方向) 的 `ei · F = 0`，不加力項
- 這是簡化假設，適用於純 Y 方向驅動流

**通用外力實作**：
- 如需支援任意方向外力 (Fy, Fz)，需對**所有 9 個方向**按 Guo 公式分配
- 包含 `(1 - sα/2)` 因子和速度相關項 `u`

### 待實作項目

- [ ] 實作 Guo 力模型的力增量計算
- [ ] 加入 `(1 - sα/2)` 修正因子
- [ ] 實作速度相關項 `(ei·u) ei / cs^4`
- [ ] 修正宏觀速度計算（半步修正）
- [ ] 支援向量外力 `(Fy, Fz)` 而非單軸
- [ ] 對所有 9 個方向正確分配力項

### 物理意義

**為何需要 Half-Way 修正**：
1. **二階精度**：確保離散 Navier-Stokes 方程的時間精度達到 O(dt²)
2. **正確動量守恆**：力增量與速度更新在半時間步一致
3. **減少格子效應**：避免力項引入的數值偏差

**參考文獻**：
- Guo, Z., Zheng, C., & Shi, B. (2002). Discrete lattice effects on the forcing term in the lattice Boltzmann method. *Physical Review E*, 65(4), 046308.

---

## 2026-01-17 Xi 插值權重筆記（GetIntrplParameter_Xi / GetXiParameter / GetParameterXi）

### 為什麼需要同時輸入 `pos_y` 與 `pos_z`？

補充：`GetIntrplParameter_Xi()` 本身是 `void` 無參數；真正需要 `pos_z`/`pos_y` 的是它在迴圈中呼叫的 `GetXiParameter(...)`（CPU 版對應 `GetParameterXi(...)`），用來把 `(y,z)` 轉成 `pos_xi` 並產生插值權重。

Periodic Hill 的底壁高度為 `HillFunction(y)`（隨 `y` 改變），因此從物理座標 `(y, z)` 轉成無因次化的 `xi` 時，需要同時知道：
- 這個位置的底壁高度 `HillFunction(pos_y)`
- 這個位置的垂直座標 `pos_z`

關鍵在於：局部通道高度與離壁距離都跟 `y` 有關，所以同一個 `pos_z` 在不同的 `pos_y` 會得到不同的 `pos_xi`。

`pos_xi` 的核心關係式（概念上就是你筆記裡的 `pos_z - Hill - 0.5*minSize`）：

```cpp
L      = LZ - HillFunction(pos_y) - minSize;
pos_xi = LXi * (pos_z - (HillFunction(pos_y) + minSize/2.0)) / L;
```

也因此在使用上，像 `y+` 與 `y-` 會對應到不同的 `HillFunction(y±minSize)`，所以即使 `pos_z` 一樣，`pos_xi` 仍會不同（權重需要分開算）。

### 兩個指標（`XiPara` / `Pos_xi`）的意義

以 `GetParameterXi(...)` / `GetXiParameter(...)` 這類函式來說，常見有兩個「指標類」參數：
- `XiPara`（或 GPU 版的 `double *XiPara*_h[7]`）：用來**存插值權重**。6th order 會有 7 個權重，所以第一維是 7；第二維是「要存到哪個網格點」的索引。
- `Pos_xi`（例如 `xi_h`）：`xi` 方向的**座標陣列**，提供給 `GetParameter_6th(...)` 去定位 stencil 並計算 7 點權重。

### 為什麼「矩陣第二個參數」會是 `NYD6*NZ6`？

`z_h` 在程式裡是以 `(j,k)` 的 2D 分佈使用，但通常用一維連續記憶體存放並 flatten：
- `index = j*NZ6 + k`
- 總點數 = `NYD6 * NZ6`

而 `xi` 插值權重會隨 `y` 改變（因為 `HillFunction(y)` 會影響 `L` 與 `pos_xi`），所以每個 `(j,k)` 都需要各自一組 7 點權重；因此每個權重陣列的長度就會是 `NYD6*NZ6`。

### 關鍵程式碼（從 `periodic hill_6thIBLBM` 搜出）

檔案：`periodic hill_6thIBLBM/50hill_4GPU/initialization.h`

```cpp
void GetXiParameter(
    double *XiPara_h[7],    double pos_z,       double pos_y,
    double *Pos_xi,         int IdxToStore,     int k  )
{
    double L = LZ - HillFunction(pos_y) - minSize;
    double pos_xi = LXi * (pos_z - (HillFunction(pos_y)+minSize/2.0)) / L;

    if( k >= 3 && k <= 6 ){
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, 3 );
    } else if ( k >= NZ6-7 && k <= NZ6-4 ) {
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, NZ6-10 );
    } else {
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, k-3 );
    }
}

void GetIntrplParameter_Xi() {
    for( int j = 3; j < NYD6-3; j++ ){
    for( int k = 3; k < NZ6-3;  k++ ){
        GetXiParameter( XiParaF3_h,  z_h[j*NZ6+k],         y_h[j]-minSize, xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF4_h,  z_h[j*NZ6+k],         y_h[j]+minSize, xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF5_h,  z_h[j*NZ6+k]-minSize, y_h[j],         xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF6_h,  z_h[j*NZ6+k]+minSize, y_h[j],         xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF15_h, z_h[j*NZ6+k]-minSize, y_h[j]-minSize, xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF16_h, z_h[j*NZ6+k]-minSize, y_h[j]+minSize, xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF17_h, z_h[j*NZ6+k]+minSize, y_h[j]-minSize, xi_h, j*NZ6+k, k );
        GetXiParameter( XiParaF18_h, z_h[j*NZ6+k]+minSize, y_h[j]+minSize, xi_h, j*NZ6+k, k );
    }}
}
```

（對照概念）檔案：`initialization.h`

```cpp
void GetParameterXi(double** XiPara , double pos_y , double pos_z , double* Pos_xi , double now , double start){
    double L = LZ - HillFunction( pos_y ) - minSize;
    double pos_xi = (LXi / L) * (pos_z - HillFunction( pos_y ) - minSize/2.0);
    ...
}
```
