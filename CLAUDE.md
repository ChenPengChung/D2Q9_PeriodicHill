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
- ✅ IsLeftHill_Boundary_DiagonalMinus45() (左丘 -45° 邊界判斷)
- ✅ IsRightHill_Boundary_DiagonalMinus135() (右丘 -135° 邊界判斷)
- ✅ Left_q_yPlus() (左丘 +Y BFL 距離)
- ✅ Right_q_yMinus() (右丘 -Y BFL 距離)
- ✅ Left_q_Diagonal45() (左丘 45° BFL 距離)
- ✅ Right_q_Diagonal135() (右丘 135° BFL 距離)
- ✅ Left_q_DiagonalMinus45() (左丘 -45° BFL 距離)
- ✅ Right_q_DiagonalMinus135() (右丘 -135° BFL 距離)

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

### 你修改版 `GetParameterXi(...)` 是否等價？（避免顯式無因次化）

你提出的改法核心是：先算相對位置比例 `rate`（等同 `xi/LXi`），再把同一個比例映回某個參考列 `y_ref = y_global[now_y]` 的實體座標 `pos_z_mapped`，最後直接在 `Pos_z` 上算 6th 權重。

在目前專案的 `tanhFunction(...)` 網格定義下，這樣做在數學上是**等價**的（權重會一致，差異只會來自浮點誤差），因為：
- 由定義可得 `z(y,k) - Hill(y) - minSize/2` 與 `xi(k)` 對同一組 `k` 是**線性比例縮放**（`tanhFunction` 對長度參數 `L` 是線性的）。
- 6th Lagrange 權重在座標做仿射變換（線性縮放 + 平移）下不變；用 `xi_h` 當座標算權重，或改用對應 `y_ref` 的 `z` 座標陣列算權重，本質上等價。

你版本中這兩行的物理意義也可以這樣看：
- `rate = (pos_z - Hill(pos_y) - minSize/2) / L(pos_y)` 其實就是 `rate = pos_xi / LXi`（只是你不再顯式寫出 `pos_xi`）。
- `pos_z_mapped = rate * L(y_ref) + Hill(y_ref) + minSize/2` 是把同一個 `rate` 映射回參考列的 `z`。

實作上注意兩點：
- 你程式片段中 `double pos_z = ...;` 會和參數 `pos_z` 重名，會編譯失敗；請改名成 `pos_z_mapped`（或類似）。
- `Pos_z` 必須是 **`y_ref = y_global[now_y]` 那一列**的 `z` 座標陣列（同樣的 flatten 索引），而且 `now_y/now_z` 的索引要跟 `Pos_z` 的儲存方式一致，才會等價。

### 關於曲面下「2D/3D 預配置連乘權重」是否會引入額外誤差？

你提出的直覺是合理的：如果用「平面直角座標」的想法去看，`y` 不同時 `Hill(y)` 不同，確實會導致同一個物理 `z` 對應到不同的無因次座標；看起來像是「多個點共用同一組 `xi` 權重」會產生誤差。

但這份程式（`GetXiParameter` + `F3/F4..._Intrpl7` 的寫法）其實是在做**曲線座標 (y, xi) 的張量積插值**，而不是在 (y, z) 的直角網格上插值：
- 網格資料點是 `(y_j, xi_k)` 的 tensor-product；對應到物理座標才會變成 `(y_j, z(y_j,xi_k))` 的「曲面網格」。
- 因此在做 2D/3D 插值時，內層先用同一個 `xi_dep`（由目標點 `(pos_y,pos_z)` 算出）在每個 `y_j` 上求 `f(y_j, xi_dep)`，外層再沿 `y` 插到 `y_dep`，這是標準的 separable/interpolation-in-computational-space 作法。
- 換句話說，程式**不是**假設不同 `y` 時「同一個 `z`」要共用同一套權重；而是用 `xi` 把不同 `y` 的垂直方向做正規化後，在同一個 `xi` 上插值。

仍然會有誤差，但主要是「插值截斷誤差」（跟網格解析度/階數/函數光滑度有關），不是因為權重被錯誤共用。

如果你想做更「直角物理空間」的做法（固定物理 `z` 在不同 `y` 上取值），那就需要對每個 `y` stencil 點各自算 `xi(y_j, z_dep)`，讓每一個 `y_j` 都用**不同的 `xi` 權重**做內層插值；那會變成非張量積（或更昂貴的）插值流程，成本與實作複雜度都會明顯上升。

### `xi_h[0..2]`（以及尾端 buffer）未初始化是否會有問題？

在原作者的 `periodic hill_6thIBLBM/50hill_4GPU` 中，`xi_h` 的初始化刻意只做「非 buffer」區間：
- `periodic hill_6thIBLBM/50hill_4GPU/initialization.h:130`～`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:136`
  - `for( int k = bfr; k < NZ6-bfr; k++ )`，其中 `bfr = 3`
  - 所以只會算 `k = 3 .. NZ6-4`
  - `k = 0,1,2` 與 `k = NZ6-3, NZ6-2, NZ6-1` 都不會被賦值（屬於 buffer/ghost layers）

對應程式碼（`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:110` 內的 `GenerateMesh_Z()`）：

```cpp
int bfr = 3;
...
for( int k = bfr; k < NZ6-bfr; k++ ){
    xi_h[k] = tanhFunction( LXi, minSize, a, (k-3), (NZ6-7) ) - minSize/2.0;
}
```

**對求解本身通常不會有問題**，原因是後續插值權重的預配置與 kernel 計算也同樣跳過 buffer 範圍：
- `GetIntrplParameter_Xi()` 只對 `k = 3 .. NZ6-4` 預先計算權重  
  - `periodic hill_6thIBLBM/50hill_4GPU/initialization.h:231`～`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:241`
- CUDA kernel 直接排除 `k <= 2` 與 `k >= NZ6-3` 的 cell（不會用到 buffer 位置的 `idx_xi = j*NZ6 + k`）  
  - `periodic hill_6thIBLBM/50hill_4GPU/evolution.h:80`（`if( ... k <= 2 || k >= NZ6-3 ) return;`）

對應程式碼（預配置只算 interior）：

```cpp
for( int j = 3; j < NYD6-3; j++ ){
for( int k = 3; k < NZ6-3;  k++ ){
    GetXiParameter( XiParaF3_h, z_h[j*NZ6+k], y_h[j]-minSize, xi_h, j*NZ6+k, k );
    ...
}}
```

對應程式碼（kernel 直接跳過 buffer）：

```cpp
if( i <= 2 || i >= NX6-3 || k <= 2 || k >= NZ6-3 ) return;
```

而且就算是在靠近下壁的 `k = 3..6`，權重也會用固定 stencil 起點 `n=3`，因此完全不需要碰到 `xi_h[0..2]`：
- `periodic hill_6thIBLBM/50hill_4GPU/initialization.h:192`～`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:198`

對應程式碼（靠近邊界固定用 `n=3`，所以用到的是 `xi_h[3..9]`）：

```cpp
if( k >= 3 && k <= 6 ){
    GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, 3 );
}
```

**但有一個「容易誤會」的副作用**：程式會把整個 `xi_h[0..NZ6-1]` 都輸出/拷貝，包含未初始化的 buffer 值：
- 輸出 `meshXi.DAT` 會把 `xi_h[0..NZ6-1]` 全部寫出來  
  - `periodic hill_6thIBLBM/50hill_4GPU/initialization.h:170`～`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:175`
- `cudaMemcpy(xi_d, xi_h, NZ6*sizeof(double), ...)` 也會把未初始化的 buffer 一起拷到 GPU  
  - `periodic hill_6thIBLBM/50hill_4GPU/initialization.h:178`

對應程式碼（輸出/拷貝包含 buffer）：

```cpp
for( int k = 0; k < NZ6; k++ ){
    fprintf( meshXi, "%.15lf\n", xi_h[k] );
}
...
CHECK_CUDA( cudaMemcpy(xi_d, xi_h, NZ6*sizeof(double), cudaMemcpyHostToDevice) );
```

因此：
- **模擬計算面**：通常安全（因為 buffer 不會被用到）。
- **後處理/除錯面**：`meshXi.DAT` 的前 3 個與最後 3 個值可能是隨機/垃圾值，容易造成誤解。

建議（非必要，但能避免踩雷）：
- 在輸出 `meshXi.DAT` 時只輸出 `k=3..NZ6-4`；或
- 明確把 `xi_h[0..2]`、`xi_h[NZ6-3..NZ6-1]` 設成有意義的值（例如複製邊界、或設成 `NAN` 讓除錯一眼看出來）。

### 為什麼 `XiParameter*_d` 不會越界、而 buffer 區也不會被用到？

雖然 `XiParaF*_d[i]` 的配置大小是 `NYD6*NZ6`（每個點一組權重），但程式在「產生」與「使用」兩個階段都把 buffer 排除掉：
- **產生階段**：只填 `k = 3 .. NZ6-4`（`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:231`～`periodic hill_6thIBLBM/50hill_4GPU/initialization.h:241`）
- **使用階段**：kernel 直接 `return` 掉 `k <= 2` 與 `k >= NZ6-3`（`periodic hill_6thIBLBM/50hill_4GPU/evolution.h:80`）

所以像你提到的最靠近壁面的有效點（例如 `k=3,4,5`，對應 `idx_xi = NZ6*j + 3/4/5`）：
- **不會越界**：`idx_xi` 仍在 `0 .. NYD6*NZ6-1` 範圍內。
- **會被使用**：因為 `k=3,4,5` 都不在 `k<=2` 的排除範圍內。
- **不依賴 `xi_h[0..2]`**：權重 stencil 會從 `n=3` 開始取（避免碰到 buffer 的 `xi_h[0..2]`）。

---

## 2026-01-18 stream_collide CUDA Kernel 呼叫機制

### stream_collide 不是 for 迴圈，而是 CUDA kernel

`stream_collide` 和 `stream_collide_Buffer` 是 **CUDA `__global__` kernel**，不使用傳統 CPU 的 for 迴圈。它們的「迴圈」是透過 **CUDA grid/block 配置**來實現平行計算。

### Grid/Block 配置

檔案：`periodic hill_6thIBLBM/50hill_4GPU/evolution.h:785`～`793`

```cpp
// stream_collide 用這組（主要計算區域）
dim3 griddim(  NX6/NT+1, NYD6, NZ6);
dim3 blockdim( NT, 1, 1);

// stream_collide_Buffer 用這組（Y方向邊界 buffer 區）
dim3 griddimBuf(NX6/NT+1, 1, NZ6);
dim3 blockdimBuf(NT, 4, 1);
```

### Kernel 呼叫位置

檔案：`periodic hill_6thIBLBM/50hill_4GPU/evolution.h`

```cpp
// evolution.h:795 - 處理 j=3 的 buffer 區
stream_collide_Buffer<<<griddimBuf, blockdimBuf, 0, stream1>>>(
    f_old[0], f_old[1], ..., f_old[18],
    f_new[0], f_new[1], ..., f_new[18],
    XPara0_d[0..6], XPara2_d[0..6],
    YPara0_d[0..6], YPara2_d[0..6],
    XiParaF3_d[0..6], XiParaF4_d[0..6], ...
    u, v, w, rho_d, Force_d, 3, ...  // 注意最後參數 3 表示 j=3
);

// evolution.h:823 - 處理 j=NYD6-7 的 buffer 區
stream_collide_Buffer<<<griddimBuf, blockdimBuf, 0, stream1>>>(
    ...
    u, v, w, rho_d, Force_d, NYD6-7, ...  // 注意最後參數 NYD6-7
);

// evolution.h:859 - 處理主要計算區域（interior）
stream_collide<<<griddim, blockdim, 0, stream0>>>(
    f_old[0], f_old[1], ..., f_old[18],
    f_new[0], f_new[1], ..., f_new[18],
    XPara0_d[0..6], XPara2_d[0..6],
    YPara0_d[0..6], YPara2_d[0..6],
    XiParaF3_d[0..6], XiParaF4_d[0..6], ...
    u, v, w, rho_d, Force_d, ...
);
```

### Kernel 內部的索引計算

在 kernel 內部，會用 `blockIdx`/`threadIdx` 計算出對應的 `(i, j, k)` 索引：

檔案：`periodic hill_6thIBLBM/50hill_4GPU/evolution.h:80` 附近

```cpp
int i = blockIdx.x * blockDim.x + threadIdx.x;  // X 方向
int j = blockIdx.y;                              // Y 方向
int k = blockIdx.z;                              // Z 方向

// 跳過 buffer 區域
if( i <= 2 || i >= NX6-3 || k <= 2 || k >= NZ6-3 ) return;
```

### CUDA kernel 等效於 CPU 三層 for 迴圈

上述 CUDA 配置等同於 CPU 上的三層 for 迴圈：

```cpp
// CPU 等效寫法（僅供理解，實際是 GPU 平行執行）
for(int i = 3; i < NX6-3; i++)      // X 方向（由 griddim.x * blockdim.x 展開）
for(int j = 0; j < NYD6; j++)        // Y 方向（由 griddim.y 展開）
for(int k = 3; k < NZ6-3; k++)       // Z 方向（由 griddim.z 展開）
{
    // stream + collide 計算邏輯
}
```

### 為什麼需要分開處理 Buffer 區？

`stream_collide_Buffer` 專門處理 Y 方向邊界（`j=3` 和 `j=NYD6-7`），原因：
1. 這些位置需要特殊的週期性邊界條件處理
2. 分開處理可以使用不同的 stream（`stream1`），與主計算（`stream0`）平行執行
3. `blockdimBuf` 的 `y=4` 表示一次處理 4 列（`j=3,4,5,6` 或 `j=NYD6-7, NYD6-6, NYD6-5, NYD6-4`）

### 執行順序

1. `stream_collide_Buffer`（j=3）在 `stream1` 上啟動
2. `stream_collide_Buffer`（j=NYD6-7）在 `stream1` 上排隊
3. `AccumulateUbulk` 在 `stream1` 上計算平均速度
4. `stream_collide`（主計算區）在 `stream0` 上啟動

由於使用不同的 CUDA stream，Buffer 處理與主計算可以部分重疊執行。

---

## 2026-01-18 BFL 邊界條件初始化筆記

### BFL 邊界處理的核心概念

**反彈方向（離開壁面的方向）才需要做邊界處理**，不是入射方向。

當粒子撞到壁面後會反彈：
- **入射方向**：粒子**進入**壁面的方向
- **反彈方向**：粒子**離開**壁面的方向 → **這個才需要 BFL 處理**

### 邊界判斷函數的設計邏輯

邊界判斷函數檢測的是：「從這個計算點往**入射方向**移動，會不會撞到壁面？」
- 如果會撞到 → 這個點是**反彈方向**的邊界計算點
- 該**反彈方向**的分佈函數需要用 BFL 更新

### D2Q9 方向與邊界處理對應關係

| 邊界判斷函數 | 入射方向（撞牆） | 反彈方向（離開） | 需要更新 | 用誰的值來插值 |
|-------------|-----------------|-----------------|---------|---------------|
| `IsLeftHill_Boundary_yPlus` | F3 (-Y) 撞左丘 | **F1 (+Y)** | F1 | F3 |
| `IsRightHill_Boundary_yMinus` | F1 (+Y) 撞右丘 | **F3 (-Y)** | F3 | F1 |
| `IsLeftHill_Boundary_Diagonal45` | F7 (-Y,-Z) 撞左丘 | **F5 (+Y,+Z)** | F5 | F7 |
| `IsRightHill_Boundary_Diagonal135` | F8 (+Y,-Z) 撞右丘 | **F6 (-Y,+Z)** | F6 | F8 |

### BFLInitialization 函數結構

檔案：`initialization.h:129`

```cpp
void BFLInitialization() {
    for(int j = 3; j < NY6-3; j++){
        for(int k = 3; k < NZ6-3; k++){

            // F1 的邊界處理：F3 入射撞左丘 → F1 反彈離開
            if(IsLeftHill_Boundary_yPlus(y_global[j], z_global[j*NZ6+k])){
                double q1 = Left_q_yPlus(y_global[j], z_global[j*NZ6+k]);
                double delta1 = minSize * (1.0 - 2.0*q1);
                // BFL 反彈點在 +Y 方向: y + delta, z 不變
                GetParameter_6th(YBFLParaF3_h, y_global[j]+delta1, ...);  // 用 F3 插值
                GetXiParameter(XiBFLParaF3_h, z_global[j*NZ6+k], y_global[j]+delta1, ...);
                Q1_h[j*NZ6+k] = q1;
            }

            // F3 的邊界處理：F1 入射撞右丘 → F3 反彈離開
            if(IsRightHill_Boundary_yMinus(y_global[j], z_global[j*NZ6+k])){
                double q3 = Right_q_yMinus(y_global[j], z_global[j*NZ6+k]);
                double delta3 = minSize * (1.0 - 2.0*q3);
                // BFL 反彈點在 -Y 方向: y - delta, z 不變
                GetParameter_6th(YBFLParaF1_h, y_global[j]-delta3, ...);  // 用 F1 插值
                GetXiParameter(XiBFLParaF1_h, z_global[j*NZ6+k], y_global[j]-delta3, ...);
                Q3_h[j*NZ6+k] = q3;
            }

            // F5 的邊界處理：F7 入射撞左丘 → F5 反彈離開
            if(IsLeftHill_Boundary_Diagonal45(y_global[j], z_global[j*NZ6+k])){
                double q5 = Left_q_Diagonal45(y_global[j], z_global[j*NZ6+k]);
                double delta5 = minSize * (1.0 - 2.0*q5);
                // BFL 反彈點在 (+Y,+Z) 方向: y + delta, z + delta
                GetParameter_6th(YBFLParaF7_h, y_global[j]+delta5, ...);  // 用 F7 插值
                GetXiParameter(XiBFLParaF7_h, z_global[j*NZ6+k]+delta5, y_global[j]+delta5, ...);
                Q5_h[j*NZ6+k] = q5;
            }

            // F6 的邊界處理：F8 入射撞右丘 → F6 反彈離開
            if(IsRightHill_Boundary_Diagonal135(y_global[j], z_global[j*NZ6+k])){
                double q6 = Right_q_Diagonal135(y_global[j], z_global[j*NZ6+k]);
                double delta6 = minSize * (1.0 - 2.0*q6);
                // BFL 反彈點在 (-Y,+Z) 方向: y - delta, z + delta
                GetParameter_6th(YBFLParaF8_h, y_global[j]-delta6, ...);  // 用 F8 插值
                GetXiParameter(XiBFLParaF8_h, z_global[j*NZ6+k]+delta6, y_global[j]-delta6, ...);
                Q6_h[j*NZ6+k] = q6;
            }
        }
    }
}
```

### delta 的計算

**delta** 是 BFL 反彈點相對於計算點的偏移量：

```cpp
delta = minSize * (1.0 - 2.0*q)
```

其中 `q` 是計算點到壁面的無因次距離（`0 < q < 1`）。

- 當 `q = 0`（計算點在壁面上）：`delta = minSize`
- 當 `q = 0.5`（計算點在中間）：`delta = 0`
- 當 `q = 1`（計算點離壁面一格）：`delta = -minSize`

對於斜向（45°/135°），若 `q` 是用這份程式的定義（`q = |Δy|/minSize`，等同沿 link 的無因次距離），則 `delta` 仍用 `delta = minSize*(1-2q)`；斜向的「實際位移長度」才會是 `sqrt(2.0)*delta`（但你在座標上做的是 `y±delta, z±delta`）。

### 權重陣列命名規則

`YBFLParaF3_h` 和 `XiBFLParaF3_h` 的命名意義：
- **F3** 代表這個權重陣列是用來**插值 F3 分佈函數**的
- 插值結果會用來**更新 F1**（F3 的反方向）

同理：
- `YBFLParaF1_h` → 插值 F1，更新 F3
- `YBFLParaF7_h` → 插值 F7，更新 F5
- `YBFLParaF8_h` → 插值 F8，更新 F6

---

## 重要程式碼索引

### initializationTool.h - 初始化工具函數

| 編號 | 函數名稱 | 行號 | 功能說明 |
|------|---------|------|----------|
| 1 | `HillFunction_Inverse_Left()` | 29-43 | 左半丘反函數（二分法） |
| 2 | `HillFunction_Inverse_Right()` | 55-68 | 右半丘反函數（二分法） |
| 3 | `tanhFunction` | 82-85 | 雙曲正切非均勻網格座標轉換巨集 |
| 4 | `GetNonuniParameter()` | 96-117 | 計算非均勻網格伸縮參數 a |
| 5 | `Lagrange_6th()` | 135-138 | 六階 Lagrange 插值基底函數 |
| 6 | `GetParameter_6th()` | 151-159 | 產生六階 Lagrange 插值預配置權重 |
| 7 | `IsLeftHill_Boundary_yPlus()` | 174-185 | 判斷 F1 (+Y) 邊界計算點 |
| 8 | `IsRightHill_Boundary_yMinus()` | 198-209 | 判斷 F3 (-Y) 邊界計算點 |
| 9 | `IsLeftHill_Boundary_Diagonal45()` | 222-237 | 判斷 F5 (+Y,+Z) 邊界計算點 |
| 10 | `IsRightHill_Boundary_Diagonal135()` | 250-265 | 判斷 F6 (-Y,+Z) 邊界計算點 |
| 11 | `IsLeftHill_Boundary_DiagonalMinus45()` | 282-297 | 判斷 F8 邊界（斜率>1時用） |
| 12 | `IsRightHill_Boundary_DiagonalMinus135()` | 315-330 | 判斷 F7 邊界（斜率>1時用） |
| 13 | `Left_q_yPlus()` | 342-348 | 計算左丘 +Y 方向 q 值 |
| 14 | `Right_q_yMinus()` | 360-366 | 計算右丘 -Y 方向 q 值 |
| 15 | `Left_q_Diagonal45()` | 378-404 | 計算左丘 45° 方向 q 值 |
| 16 | `Right_q_Diagonal135()` | 416-441 | 計算右丘 135° 方向 q 值 |
| 17 | `Left_q_DiagonalMinus45()` | 453-478 | 計算左丘 -45° 方向 q 值 |
| 18 | `Right_q_DiagonalMinus135()` | 490-515 | 計算右丘 -135° 方向 q 值 |

### initialization.h - 初始化主函數

| 函數名稱 | 行號 | 功能說明 |
|---------|------|----------|
| `InitialUsingDftFunc()` | 7-31 | 初始化分佈函數與外力項 |
| `GenerateMesh_Y()` | 34-47 | 建立 Y 方向均勻網格 |
| `GenerateMesh_Z()` | 49-74 | 建立 Z 方向非均勻網格（含山丘） |
| `GetXiParameter()` | 75-90 | 計算 Xi 方向插值權重 |
| `GetIntrplParameter_Y()` | 92-97 | 預配置 Y 方向插值權重 |
| `GetIntrplParameter_Xi()` | 100-127 | 預配置 Xi 方向插值權重（所有 F1-F8） |
| `BFLInitialization()` | 129-183 | BFL 邊界條件初始化（權重與 q 值） |

### 關鍵變數對照表

| 變數名稱 | 用途 |
|---------|------|
| `y_global[NY6]` | Y 方向物理座標陣列 |
| `z_global[NY6*NZ6]` | Z 方向物理座標陣列（2D flatten） |
| `xi_h[NZ6]` | 無因次化 Z 座標（不含山丘影響） |
| `XiParaF*_h[7]` | F* 方向的 Xi 插值權重（7點） |
| `YBFLParaF*_h[7]` | BFL 用的 Y 方向插值權重 |
| `XiBFLParaF*_h[7]` | BFL 用的 Xi 方向插值權重 |
| `Q*_h[]` | 各方向的 BFL q 值（計算點到壁面距離）|

---

## 2026-01-18 未初始化清單（主程式撰寫時要統整）

以下是「在目前程式流程中，不一定會被寫入」的資料區域；若後續程式會讀到它們，需在配置/初始化階段做 `memset`、填 0、或填 `NAN`/sentinel。

### 網格/座標

- `xi_h` 的 buffer 區：`initialization.h:58` 只寫 `k=3..NZ6-4` → `xi_h[0..2]` 與 `xi_h[NZ6-3..NZ6-1]` 未初始化。
- `z_global` 的 buffer 區：`initialization.h:67` 只寫 `k=3..NZ6-4`，另外只補 `k=2` 與 `k=NZ6-3`（`initialization.h:71-72`）→ `k=0,1,NZ6-2,NZ6-1` 未初始化（每個 `j` 都一樣）。
- `Force[1]`：`InitialUsingDftFunc()` 只賦值 `Force[0]`（`initialization.h:30`）→ 若你把 `Force` 當成 `(Fy,Fz)`，則 `Force[1]` 需要明確初始化。

### 插值權重（非 BFL）

- `YPara0_h[*][i]`、`YPara2_h[*][i]`：`GetIntrplParameter_Y()` 只算 `i=3..NY6-4`（`initialization.h:92-97`）→ buffer `i<=2` 與 `i>=NY6-3` 未初始化。
- `XiParaF1_h..XiParaF8_h`：`GetIntrplParameter_Xi()` 只算 `j=3..NYD6-4`、`k=3..NZ6-4`（`initialization.h:100-127`）→ 任何在這個範圍外的 `(j,k)` 權重未初始化。

### BFL 權重與 q 值

- `Q1_h/Q3_h/Q5_h/Q6_h`：只在對應 `Is*Boundary*()` 成立時才寫入（`initialization.h:140-180`）→ 非邊界點位置未初始化。
- `YBFLParaF3_h/YBFLParaF1_h/YBFLParaF7_h/YBFLParaF8_h` 與 `XiBFLParaF3_h/XiBFLParaF1_h/XiBFLParaF7_h/XiBFLParaF8_h`：同上，只在邊界點填權重 → 非邊界點位置未初始化。

### 目前在 repo 內「找不到宣告」的符號（只存在於筆記，尚未在可編譯的程式檔出現）

`initialization.h` 直接使用但在程式碼裡尚未搜尋到宣告（你說「假設已宣告」的那批）：
- 流場/分佈：`rho[]`, `v[]`, `w[]`, `f[9][]`, `Force[]`
- 網格座標：`y_global[]`, `z_global[]`, `xi_h[]`, 以及你在權重生成時用的 `y_h[]`, `z_h[]`
- 一般插值權重：`YPara0_h`, `YPara2_h`, `XiParaF1_h..XiParaF8_h`
- BFL：`Q1_h`, `Q3_h`, `Q5_h`, `Q6_h`, `YBFLParaF3_h/YBFLParaF1_h/YBFLParaF7_h/YBFLParaF8_h`, `XiBFLParaF3_h/XiBFLParaF1_h/XiBFLParaF7_h/XiBFLParaF8_h`

### D2Q9 速度方向定義

```
     F6(-1,+1)  F2(0,+1)  F5(+1,+1)
              \    |    /
               \   |   /
     F3(-1,0) ←― F0 ―→ F1(+1,0)
               /   |   \
              /    |    \
     F7(-1,-1)  F4(0,-1)  F8(+1,-1)
```

| 方向 | 速度向量 (ey, ez) | 權重 |
|------|------------------|------|
| F0 | (0, 0) | 4/9 |
| F1 | (+1, 0) | 1/9 |
| F2 | (0, +1) | 1/9 |
| F3 | (-1, 0) | 1/9 |
| F4 | (0, -1) | 1/9 |
| F5 | (+1, +1) | 1/36 |
| F6 | (-1, +1) | 1/36 |
| F7 | (-1, -1) | 1/36 |
| F8 | (+1, -1) | 1/36 |
