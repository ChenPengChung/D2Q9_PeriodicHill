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

---

## 2025-01-20 插值與非均勻網格快取調整（AI 協助）

### 新增/變更
- 新增全域 `nonuni_a`：`globalVariables.h` 宣告，`main.cpp` 定義，`GenerateMesh_Z` 計算一次後共用，避免多次呼叫 `GetNonuniParameter()`。
- `GenerateMesh_Z` 內部：使用 `nonuni_a` 生成 `xi_h`、`z_global`，邏輯不變，僅移除重複求 a。
- 新增 `BuildXiWeights`（`initialization.h`）：對同一 `(j,k)` 一次性計算七列的 `RelationXi`/`GetParameter_6th2`，寫入對應 `XiPara*_h[?][index_xi + row*NZ6]`；內含 y 週期處理，使用預先的 `nonuni_a`，避免 8 倍重算。
- `GetIntrplParameter_Xi`：八個速度方向改呼叫 `BuildXiWeights`，原 `GetXiParameter` 保留但不再使用。
- 修正 `GetParameter_6th2` 第五列使用 `pos_z` 的錯誤，統一用同一組 `pos_z2`。

#### 主要程式碼片段

`globalVariables.h` / `main.cpp`：新增共用的非均勻參數 `nonuni_a`。
```cpp
// globalVariables.h
extern double nonuni_a;

// main.cpp
double nonuni_a = 0.0;
```

`initialization.h`：在產網格時計算一次 a，之後都用同一數值。
```cpp
void GenerateMesh_Z() {
    ...
    nonuni_a = GetNonuniParameter();      // 只算一次
    const double a = nonuni_a;
    xi_h[k] = tanhFunction(LXi, minSize, a, (k-3), (NZ6-7)) - minSize/2.0;
    ...
}
```



`initialization.h`：`BuildXiWeights` 將同一 `(j,k)` 的七列權重一次算完，供 8 個方向共用。
```cpp
inline void BuildXiWeights(double* XiPara_h[7], double pos_z, double pos_y,
                           int index_xi, int j, int k) {
    int jj = j;
    if (jj < 3) jj += NY6 - 7;
    if (jj > NY6 - 4) jj -= (NY6 - 7);
    double L = LZ - HillFunction(pos_y) - minSize;
    double pos_xi = pos_z - (HillFunction(pos_y) + minSize/2.0);
    double j_cont = Inverse_tanh_index(pos_xi, L, minSize, nonuni_a, (NZ6-7));

    auto fillRow = [&](int rowOffset, int storeRow) {
        int rowIdx = jj + rowOffset;
        if (rowIdx < 0) rowIdx += NY6;
        if (rowIdx >= NY6) rowIdx -= NY6;
        double H = HillFunction(y_global[rowIdx]);
        double Lrow = LZ - H - minSize;
        double pos_z_row = tanhFunction(Lrow, minSize, nonuni_a, j_cont, (NZ6-7)) - minSize/2.0;
        double rel[7];
        RelationXi(j_cont, Lrow, minSize, nonuni_a, (NZ6-7), rel);
        GetParameter_6th2(XiPara_h, pos_z_row, rel, storeRow, index_xi);
    };
    fillRow(-3,0); fillRow(-2,1); fillRow(-1,2); fillRow(0,3);
    fillRow(1,4);  fillRow(2,5);  fillRow(3,6);
}
```

`GetIntrplParameter_Xi`：八個方向共用 `BuildXiWeights`，避免重算。
```cpp
int idx = j * NZ6 + k;
BuildXiWeights(XiParaF1_h, z_global[idx],         y_global[j]-minSize, idx, j, k);
BuildXiWeights(XiParaF2_h, z_global[idx]-minSize, y_global[j],         idx, j, k);
BuildXiWeights(XiParaF3_h, z_global[idx],         y_global[j]+minSize, idx, j, k);
BuildXiWeights(XiParaF4_h, z_global[idx]+minSize, y_global[j],         idx, j, k);
BuildXiWeights(XiParaF5_h, z_global[idx]-minSize, y_global[j]-minSize, idx, j, k);
BuildXiWeights(XiParaF6_h, z_global[idx]-minSize, y_global[j]+minSize, idx, j, k);
BuildXiWeights(XiParaF7_h, z_global[idx]+minSize, y_global[j]+minSize, idx, j, k);
BuildXiWeights(XiParaF8_h, z_global[idx]+minSize, y_global[j]-minSize, idx, j, k);
```

`GetParameter_6th2` 呼叫修正：`(j+1)` 那列使用同一組 `pos_z2`。
```cpp
RelationXi(j_cont, LT, minSize, a, (NZ6-7), RelationXi_4);
GetParameter_6th2(XiPara_h, pos_z2, RelationXi_4, 4, index_xi);
```

### 未改動的邏輯
- 網格生成公式、HillFunction、Y 向插值、BFL 初始化、MRT/演化主流程均未改。
- 既有的 `GetXiParameter` 仍在檔案中，可刪可留，不影響現行流程。
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

---

## 2026-01-18 MRT 矩陣巨集使用方式

### 巨集定義位置

檔案：`MRT_Matrix.h`

### 三個主要巨集

| 巨集名稱 | 展開後宣告 | 用途 |
|---------|-----------|------|
| `Matrix` | `double M[9][9] = {...}` | 分佈函數 → 矩空間變換矩陣 |
| `Matrix_Inverse` | `double M_I[9][9] = {...}` | 矩空間 → 分佈函數逆變換矩陣 |
| `Relaxation` | `double s0=0; ... double s8=omega_7;` | 鬆弛參數（對角矩陣元素） |

### evolution.h 中的使用方式

檔案：`evolution.h:64-68`

```cpp
// MRT 矩陣與鬆弛參數 (巨集展開後會宣告 M[9][9], M_I[9][9], s0~s8)
Matrix;
Matrix_Inverse;
Relaxation;
```

**注意**：
- 巨集已包含變數宣告，**不需要**預先宣告 `s0~s8`
- 展開後即可直接使用 `M[i][j]`、`M_I[i][j]`、`s0~s8`

### 錯誤示範（會導致重複宣告）

```cpp
// ❌ 錯誤：s0~s8 會被重複宣告
double s0, s1, s2, s3, s4, s5, s6, s7, s8;  // 手動宣告
Relaxation;  // 巨集內又宣告一次 → 編譯錯誤
```

### 正確示範

```cpp
// ✅ 正確：直接使用巨集
Matrix;           // → double M[9][9] = {...};
Matrix_Inverse;   // → double M_I[9][9] = {...};
Relaxation;       // → double s0=0; double s1=omega_2; ... double s8=omega_7;

// 現在可以直接使用 M, M_I, s0~s8
```

### 常見錯誤：巨集名稱拼錯

| 錯誤寫法 | 正確寫法 |
|---------|---------|
| `Inverse_Matrix;` | `Matrix_Inverse;` |
| `Relaxtion;` | `Relaxation;` |

### 鬆弛參數物理意義

```
s0 = 0        // 密度守恆（不鬆弛）
s1 = omega_2  // 能量鬆弛
s2 = 1.0      // free parameter
s3 = 0        // Y動量守恆（不鬆弛）
s4 = 1.0      // free parameter
s5 = 0        // Z動量守恆（不鬆弛）
s6 = 1.0      // free parameter
s7 = omega_7  // 剪切應力鬆弛（與黏滯係數相關）
s8 = omega_7  // 剪切應力鬆弛（與黏滯係數相關）
```

其中 `omega_2` 和 `omega_7` 定義在 `variables.h:36-37`：

```cpp
#define omega_2  1/(niu/9.0 + 0.5)
#define omega_7  1/(niu/3.0 + 0.5)
```

### MRT 碰撞流程（搭配 MRT_Process.h）

```cpp
// 1. 宣告矩陣與參數
Matrix;
Matrix_Inverse;
Relaxation;

// 2. 計算矩 m = M * f （使用 m_vector 巨集）
m_vector;  // → m0, m1, ..., m8

// 3. 計算平衡態矩 meq = M * feq （使用 meq 巨集）
meq;  // → meq0, meq1, ..., meq8

// 4. 碰撞：f = f - M^(-1) * S * (m - meq) （使用 collision 巨集）
collision;  // → F0_in, F1_in, ..., F8_in 更新
```

### CUDA Constant Memory 版本（可選）

如果要在 CUDA kernel 中使用，可以啟用 constant memory 版本以提升效能：

```cpp
// 在 main.cu 開頭
#define USE_CUDA_CONSTANT
#include "MRT_Matrix.h"

// 初始化
initRelaxationToGPU();  // 將鬆弛參數傳到 GPU

// 在 kernel 中直接使用
__global__ void kernel() {
    // 使用 d_M[i][j], d_M_I[i][j], d_S[i]
    // 從 constant memory 讀取，不需要每次重新建立陣列
}
```

---

## 2026-01-19 evolution.h Equilibrium 檢查與修復

### Equilibrium 分佈函數驗證結果

檔案：`evolution.h:153-165`

#### ✅ 正確的部分

**1. 密度計算 (Line 153)**
```cpp
double rho_s = F0_in + F1_in + F2_in + F3_in + F4_in + F5_in + F6_in + F7_in + F8_in;
```
符合 D2Q9 密度定義：`ρ = Σ fi`

**2. 速度計算 (Lines 154-155)**
```cpp
double v1 = (F1_in + F5_in + F8_in - (F3_in + F6_in + F7_in)) / rho_s;  // Y方向速度
double w1 = (F2_in + F5_in + F6_in - (F4_in + F7_in + F8_in)) / rho_s;  // Z方向速度
```
符合 D2Q9 動量定義：
- `ρuy = F1 + F5 + F8 - F3 - F6 - F7`（含 +Y 分量的方向 - 含 -Y 分量的方向）
- `ρuz = F2 + F5 + F6 - F4 - F7 - F8`（含 +Z 分量的方向 - 含 -Z 分量的方向）

**3. 平衡態分佈函數公式 (Lines 157-165)**

所有 9 個方向的權重與公式都正確，符合標準 D2Q9 平衡態：
```
f_eq[i] = w[i] * ρ * (1 + 3(e·u) + 4.5(e·u)² - 1.5u²)
```

| 方向 | 權重 | 速度向量 | 公式驗證 |
|-----|------|---------|---------|
| F0 | 4/9 | (0,0) | ✅ `1 - 1.5*udot` |
| F1 | 1/9 | (+1,0) | ✅ `1 + 3v + 4.5v² - 1.5*udot` |
| F2 | 1/9 | (0,+1) | ✅ `1 + 3w + 4.5w² - 1.5*udot` |
| F3 | 1/9 | (-1,0) | ✅ `1 - 3v + 4.5v² - 1.5*udot` |
| F4 | 1/9 | (0,-1) | ✅ `1 - 3w + 4.5w² - 1.5*udot` |
| F5 | 1/36 | (+1,+1) | ✅ `1 + 3(v+w) + 4.5(v+w)² - 1.5*udot` |
| F6 | 1/36 | (-1,+1) | ✅ `1 + 3(-v+w) + 4.5(-v+w)² - 1.5*udot` |
| F7 | 1/36 | (-1,-1) | ✅ `1 + 3(-v-w) + 4.5(-v-w)² - 1.5*udot` |
| F8 | 1/36 | (+1,-1) | ✅ `1 + 3(v-w) + 4.5(v-w)² - 1.5*udot` |

### 已修復的錯誤

#### 1. Line 8: 中文冒號語法錯誤 ✅
- **問題**：`#include "variables.h"：`（中文全形冒號）
- **修復**：刪除冒號 → `#include "variables.h"`

#### 2. Line 56: 迴圈範圍錯誤 ✅
- **問題**：`for(int j = 3 ; j <= NZ6-4 ; j++)`（j 應該是 Y 方向，用 NY6）
- **修復**：`for(int j = 3 ; j < NY6-3 ; j++)`

#### 3. Line 57: 迴圈範圍修正 ✅
- **問題**：`for(int k = 3 ; k <= NZ6-4 ; k++)`
- **修復**：`for(int k = 3 ; k < NZ6-3 ; k++)`

#### 4. Lines 108, 121, 134, 147: BFL q>0.5 運算子優先順序錯誤 ✅

**問題**：整數除法導致錯誤結果
```cpp
// 錯誤：1/2*q1 = (1/2)*q1 = 0*q1 = 0（整數除法）
F1_in = (1/2*q1)*f3_old[idx_xi] + ((2*q1-1)/2*q1)*f1_old[idx_xi];
```

**修復**：使用浮點數除法
```cpp
// 正確：1.0/(2.0*q1) 得到正確的分數
F1_in = (1.0/(2.0*q1))*f3_old[idx_xi] + ((2.0*q1-1.0)/(2.0*q1))*f1_old[idx_xi];
```

**BFL 線性插值公式說明**：
當 `q > 0.5` 時，使用線性插值：
```
F_反彈 = (1/(2q)) * F_入射(壁面點) + ((2q-1)/(2q)) * F_反彈(計算點)
```

### 修復後的程式碼狀態

| 行號 | 修復內容 |
|-----|---------|
| 8 | 刪除中文冒號 |
| 56 | `j < NY6-3` |
| 57 | `k < NZ6-3` |
| 108 | `1.0/(2.0*q1)` 和 `(2.0*q1-1.0)/(2.0*q1)` |
| 121 | `1.0/(2.0*q3)` 和 `(2.0*q3-1.0)/(2.0*q3)` |
| 134 | `1.0/(2.0*q5)` 和 `(2.0*q5-1.0)/(2.0*q5)` |
| 147 | `1.0/(2.0*q6)` 和 `(2.0*q6-1.0)/(2.0*q6)` |

### Equilibrium 與 MRT_Process.h 的關係

`MRT_Process.h` 中的 `Equilibrium` 巨集與 `evolution.h` 中手寫的平衡態計算功能相同：

**MRT_Process.h 巨集版本**：
```cpp
#define Equilibrium \
    double udot = uy*uy + uz*uz; \
    double F0_eq = (4.0/9.0)  * rho_local * (1.0 - 1.5*udot); \
    ...
```
- 輸入：`rho_local`, `uy`, `uz`
- 輸出：`F0_eq ~ F8_eq`

**evolution.h 手寫版本**（Lines 153-165）：
```cpp
double rho_s = F0_in + ... + F8_in;
double v1 = (...) / rho_s;
double w1 = (...) / rho_s;
double udot = v1*v1 + w1*w1;
const double F0_eq = (4./9.) * rho_s * (1.0-1.5*udot);
...
```
- 直接從分佈函數計算密度與速度
- 然後計算平衡態

兩種方式數學上等價，`evolution.h` 的版本更適合在 stream-collide 流程中使用，因為它直接從 post-streaming 的分佈函數開始計算。

---

## 2026-01-19 外力項動態修正機制

### 功能目的

在週期性驅動流中，初始外力是根據 Poiseuille 流解析解估算的，但由於：
1. Periodic Hill 的幾何形狀造成額外壓力損失
2. 非均勻網格的數值誤差
3. MRT 碰撞的數值黏滯

實際流速可能偏離目標值 `Uref`。**動態外力調整**使用比例控制器讓實際流速趨近目標流速。

### 控制器公式

```
F_new = F_old + β × (U_target - U_actual) × U_ref / L
```

| 參數 | 意義 | 數值 |
|------|------|------|
| `β` | 控制增益 | `max(0.001, 3/Re)` |
| `U_target` | 目標流速 | `Uref` |
| `U_actual` | 實際平均流速 | `Ub_avg` |
| `L` | 特徵長度 | `LZ` |

**控制邏輯**：
- 當 `U_actual < U_target` → 誤差為正 → 增加外力 → 加速流體
- 當 `U_actual > U_target` → 誤差為負 → 減少外力 → 減速流體

### 新增函數

#### 1. `AccumulateUbulk()` (evolution.h:235-242)

```cpp
void AccumulateUbulk(double* v_field, double* Ub_sum_ptr) {
    for(int j = 3; j < NY6-3; j++) {
        for(int k = 3; k < NZ6-3; k++) {
            int idx = j * NZ6 + k;
            *Ub_sum_ptr += v_field[idx];
        }
    }
}
```

- **功能**：累積 Y 方向平均速度
- **呼叫時機**：每個時間步
- **輸入**：`v_field` - Y 方向速度場
- **輸出**：`Ub_sum_ptr` - 累積的速度總和（透過指標更新）

#### 2. `ModifyForcingTerm()` (evolution.h:261-278)

```cpp
void ModifyForcingTerm(double* Force, double* Ub_sum_ptr, int NDTFRC) {
    // 1. 計算時間與空間平均速度
    int num_cells = (NY6 - 6) * (NZ6 - 6);
    double Ub_avg = (*Ub_sum_ptr) / (double)(num_cells * NDTFRC);

    // 2. 計算控制增益
    double beta = fmax(0.001, 3.0 / (double)Re);

    // 3. 調整外力
    Force[0] = Force[0] + beta * (Uref - Ub_avg) * Uref / LZ;

    // 4. 輸出監控資訊
    printf("Force Update: Ub_avg = %.6f, Uref = %.6f, Force = %.5e\n",
           Ub_avg, Uref, Force[0]);

    // 5. 重置累加器
    *Ub_sum_ptr = 0.0;
}
```

- **功能**：使用比例控制器調整外力
- **呼叫時機**：每 `NDTFRC` 步（建議 10000 步）
- **輸入**：
  - `Force` - 外力陣列
  - `Ub_sum_ptr` - 累積的速度總和
  - `NDTFRC` - 累積的時間步數
- **輸出**：
  - 更新後的 `Force[0]`
  - 重置 `*Ub_sum_ptr = 0`

### 主程式使用範例

```cpp
// main.cpp
#include "evolution.h"

int main() {
    // ... 初始化 ...

    double Ub_sum = 0.0;           // 累積的平均速度
    int force_update_count = 0;    // 累積的時間步數
    const int NDTFRC = 10000;      // 每多少步修正一次

    for(int t = 0; t < loop; t++) {
        // 1. Stream + Collide
        stream_collide(...);

        // 2. 週期性邊界
        periodicSW(...);

        // 3. 累積平均速度
        AccumulateUbulk(v, &Ub_sum);
        force_update_count++;

        // 4. 每 NDTFRC 步修正外力
        if(force_update_count >= NDTFRC) {
            ModifyForcingTerm(Force, &Ub_sum, NDTFRC);
            force_update_count = 0;
        }

        // 5. 交換 f_old 與 f_new 指標
        std::swap(f0_old, f0_new);
        std::swap(f1_old, f1_new);
        // ... 其他方向 ...
    }
    return 0;
}
```

### 參數建議

| 參數 | 建議值 | 說明 |
|------|--------|------|
| `NDTFRC` | 10000 | 太小會震盪，太大會收斂慢 |
| `β` 下限 | 0.001 | 避免高 Re 時控制增益過小 |
| `β` 上限 | 3/Re | 低 Re 時較大，高 Re 時較小 |

### 與原始 GPU 版本的差異

| 項目 | 原始 GPU 版本 | 目前 CPU 版本 |
|------|--------------|--------------|
| 速度累積 | GPU kernel + cudaMemcpy | CPU 直接迴圈 |
| 多 GPU 同步 | MPI_Reduce + MPI_Bcast | 不需要 |
| 記憶體 | Ub_avg_d[NX6*NZ6] | Ub_sum (純量) |

---

## 2026-01-19 週期性邊界條件 periodicSW()

### 功能說明

`periodicSW()` 實現 Y 方向（Stream-Wise，主流場方向）的週期性邊界條件。

### 網格配置

```
Y 方向索引 (j):
  0   1   2  |  3   4   5  ...  NY6-6  NY6-5  NY6-4  |  NY6-3  NY6-2  NY6-1
  ←buffer→  |  ←────────── 計算區域 ──────────────→  |  ←───buffer───→
```

### 複製邏輯

#### 左側 Buffer 填充（從右邊複製）
```cpp
int idx_right = (i+NY6-6)*NZ6 + k;   // 來源：j = NY6-6, NY6-5, NY6-4
int buffer_left = i*NZ6 + k;         // 目標：j = 0, 1, 2
```

#### 右側 Buffer 填充（從左邊複製）
```cpp
int idx_left = (i+3)*NZ6 + k;        // 來源：j = 3, 4, 5
int buffer_right = (i+NY6-3)*NZ6+k;  // 目標：j = NY6-3, NY6-2, NY6-1
```

### 複製內容

- 分佈函數：`f0_new ~ f8_new`
- 宏觀量：`v`, `w`, `rho_d`

### 已修正問題

將 `k` 的範圍從 `0 ~ NZ6` 改為 `3 ~ NZ6-3`，只複製 Z 方向的有效計算區域。

---

## 2026-01-20 流場條紋問題診斷

### 問題描述

模擬結果出現**斜向條紋圖案**：
- 條紋方向：約 45° 斜向
- 條紋週期：約 7-8 格（與 7 點 stencil 相關）
- 狀態：數值穩定，沒有發散，但明顯是系統性誤差

### 流場圖像特徵分析

```
觀察到的現象：
1. 條紋沿 45° 方向 → 暗示 Y 和 Z 方向有某種耦合問題
2. 條紋週期約 7-8 格 → 可能與 6 階插值的 7 點 stencil 有關
3. 數值穩定 → 不是發散問題，是系統性的權重/索引錯誤
```

### 已排除的問題

經過詳細檢查，以下項目**確認無誤**：

#### ✅ 平衡態分佈函數計算 (evolution.h:202-214)
- 密度計算正確：`rho_s = Σ Fi_in`
- 速度計算正確：`v1 = (F1+F5+F8-F3-F6-F7)/rho_s`
- 平衡態公式正確：符合標準 D2Q9 公式

#### ✅ MRT 碰撞巨集 (MRT_Process.h)
- m_vector 巨集：正確計算 `m = M * f`
- meq 巨集：正確計算 `meq = M * feq`
- collision 巨集：正確計算 `f = f - M_I * S * (m - meq)`

#### ✅ 週期邊界條件 (evolution.h:237-276)
- 左 buffer `(0,1,2)` ← 從 `(NY6-6, NY6-5, NY6-4)` 複製 ✅
- 右 buffer `(NY6-3, NY6-2, NY6-1)` ← 從 `(3, 4, 5)` 複製 ✅

#### ✅ 初始化 (initialization.h:10-34)
- 全場初始化為 `rho=1, v=0, w=0` ✅
- 分佈函數用平衡態初始化 ✅

#### ✅ Lagrange 插值權重計算 (initializationTool.h:154-162)
- `GetParameter_6th()` 計算正確的 Lagrange 基底函數

### 可能的問題來源

#### 🔴 疑點 1：`GetXiParameter` 的 stencil 起點判斷使用目標點的 k

**檔案**：`initialization.h:79-93`

```cpp
void GetXiParameter(
    double** XiPara_h,    double pos_z,       double pos_y,
    double *Pos_xi,         int IdxToStore,     int k  )
{
    ...
    if( k >= 3 && k <= 6 ){
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, 3 );
    } else if ( k >= NZ6-7 && k <= NZ6-4 ) {
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, NZ6-10 );
    } else {
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, k-3 );
    }
}
```

**問題**：傳入的 `k` 是**目標點**的索引，但對於 streaming 來說，應該用**來源點**的位置來決定 stencil 起點。

例如 F2 (+Z)：
- 目標點 `(j, k)`
- 來源點 `(j, k-1)`（Z 方向往下一格）
- 程式用 `k`（目標點）判斷 stencil 起點
- 應該用來源點的位置來判斷

**潛在影響**：當來源點和目標點跨越 stencil 邊界時（如 `k=7` 目標但來源點應該用 `start=3`），會使用錯誤的 stencil 起點。

#### 🔴 疑點 2：`F1_Intrpl7` / `F3_Intrpl7` 的插值順序與權重可能不匹配

**檔案**：`interpolationHillISLBM.h:14-36`

`F1_Intrpl7` 和 `F3_Intrpl7` 的插值順序：
```
外層: Y 方向插值 (y_0..y_6)
內層: Xi 方向插值 (xi_0..xi_6)
```

`Y_XI_Intrpl7` (用於 F5-F8) 的插值順序：
```
外層: Xi 方向插值 (xi_0..xi_6)
內層: Y 方向插值 (y_0..y_6)
```

**兩者順序相反**！雖然理論上 2D 張量積插值結果與順序無關，但如果權重計算方式與使用方式不一致，可能產生問題。

#### 🔴 疑點 3：CLAUDE.md 參考版本與現有程式的 Xi 權重計算不一致

**CLAUDE.md:438-441 的參考版本**：
```cpp
GetXiParameter( XiParaF3_h,  z_h[j*NZ6+k],         y_h[j]-minSize, xi_h, j*NZ6+k, k );
GetXiParameter( XiParaF4_h,  z_h[j*NZ6+k],         y_h[j]+minSize, xi_h, j*NZ6+k, k );
GetXiParameter( XiParaF5_h,  z_h[j*NZ6+k]-minSize, y_h[j],         xi_h, j*NZ6+k, k );
```

**現有程式 (initialization.h:115-129)**：
```cpp
// F1: pos_y = y_global[j]-minSize
GetXiParameter( XiParaF1_h,  z_global[j*NZ6+k],         y_global[j]-minSize, xi_h, j*NZ6+k, k );
// F5: pos_y = y_global[j]-minSize  ← 與參考版本不同！參考版本是 y_h[j]
GetXiParameter( XiParaF5_h,  z_global[j*NZ6+k]-minSize, y_global[j]-minSize, xi_h, j*NZ6+k, k );
```

**差異**：
- 參考版本 F5 用 `y_h[j]`（不偏移）
- 現有程式 F5 用 `y_global[j]-minSize`（有偏移）

**物理意義**：F5 (+Y,+Z) 的來源點是 `(y-Δ, z-Δ)`，所以 `pos_y = y-Δ` 是合理的。但這與參考版本不同，需要確認哪個是正確的。

### 建議的診斷步驟

#### 步驟 1：輸出權重總和驗證

在某個特定 `(j, k)` 位置，輸出所有 7 個權重並檢查總和是否為 1：

```cpp
// 加在 GetIntrplParameter_Xi() 結尾
if(j == 64 && k == 35) {
    double sum = 0.0;
    for(int i = 0; i < 7; i++) {
        sum += XiParaF1_h[i][j*NZ6+k];
        printf("XiParaF1_h[%d] = %e\n", i, XiParaF1_h[i][j*NZ6+k]);
    }
    printf("Sum of weights = %f (should be 1.0)\n", sum);
}
```

#### 步驟 2：輸出 pos_xi 範圍驗證

確認 `pos_xi` 在合理範圍內（0 到 LXi）：

```cpp
// 加在 GetXiParameter() 內
if(pos_xi < 0 || pos_xi > LXi) {
    printf("WARNING: pos_xi = %f out of range [0, %f] at k=%d\n", pos_xi, LXi, k);
}
```

#### 步驟 3：簡化測試 - 改用均勻網格

暫時把 Z 方向改成均勻網格，看條紋是否消失：

```cpp
// 在 GenerateMesh_Z() 中，暫時改用：
for( int k = 3; k < NZ6-3; k++ ){
    z_global[j*NZ6+k] = HillFunction(y_global[j]) + (LZ - HillFunction(y_global[j])) * (k-3) / (NZ6-7);
}
```

#### 步驟 4：單方向測試

暫時只執行 F0（靜止）和 F2/F4（純 Z 方向），關閉 Y 方向和斜向方向：

```cpp
// 在 stream_collide 中暫時註解掉：
// F1_Intrpl7(...);
// F3_Intrpl7(...);
// Y_XI_Intrpl7(f5_old, ...);
// Y_XI_Intrpl7(f6_old, ...);
// Y_XI_Intrpl7(f7_old, ...);
// Y_XI_Intrpl7(f8_old, ...);

// 改為直接複製：
F1_in = f1_old[idx_xi];
F3_in = f3_old[idx_xi];
F5_in = f5_old[idx_xi];
F6_in = f6_old[idx_xi];
F7_in = f7_old[idx_xi];
F8_in = f8_old[idx_xi];
```

如果條紋消失，問題在 Y 方向或斜向插值；如果仍有條紋，問題在 Z 方向插值。

### 待修復項目清單

| 優先級 | 項目 | 狀態 |
|-------|------|------|
| 高 | 確認 `GetXiParameter` 的 stencil 起點判斷邏輯 | ⏳ 待驗證 |
| 高 | 比對 CLAUDE.md 參考版本與現有程式的 Xi 權重計算 | ⏳ 待確認 |
| 中 | 驗證 `F1_Intrpl7` / `F3_Intrpl7` 的插值順序正確性 | ⏳ 待驗證 |
| 低 | 檢查 `cell_z` 計算與 `GetXiParameter` 的一致性 | ⏳ 待確認 |

### 相關檔案索引

| 檔案 | 關鍵行號 | 內容 |
|------|---------|------|
| `interpolationHillISLBM.h` | 14-24 | `F1_Intrpl7` 巨集 |
| `interpolationHillISLBM.h` | 26-36 | `F3_Intrpl7` 巨集 |
| `interpolationHillISLBM.h` | 46-56 | `Y_XI_Intrpl7` 巨集 |
| `initialization.h` | 79-93 | `GetXiParameter` 函數 |
| `initialization.h` | 104-130 | `GetIntrplParameter_Xi` 函數 |
| `evolution.h` | 123-131 | 插值巨集呼叫 |
| `initializationTool.h` | 154-162 | `GetParameter_6th` 函數 |

### 補充：D2Q9 Streaming 方向與來源點對應

| 速度方向 | 速度向量 | Streaming 來源點 | 備註 |
|---------|---------|-----------------|------|
| F0 | (0, 0) | 原地 | 不需插值 |
| F1 | (+1, 0) | (y-Δ, z) | Y 方向偏移 |
| F2 | (0, +1) | (y, z-Δ) | Z 方向偏移 |
| F3 | (-1, 0) | (y+Δ, z) | Y 方向偏移 |
| F4 | (0, -1) | (y, z+Δ) | Z 方向偏移 |
| F5 | (+1, +1) | (y-Δ, z-Δ) | 斜向偏移 |
| F6 | (-1, +1) | (y+Δ, z-Δ) | 斜向偏移 |
| F7 | (-1, -1) | (y+Δ, z+Δ) | 斜向偏移 |
| F8 | (+1, -1) | (y-Δ, z+Δ) | 斜向偏移 |

## 2026-01-20 程式碼改進 

```cpp
//給我一個編號，產生該Y值所對應的七個無因次化座標
void RelationXi(double nonintindex, double L , double MinSize , double a , int N , double* RelationXi){
    int i = (int)floor(nonintindex);
    RelationXi[0] = tanhFunction( L , MinSize , a, i-3 , N) - MinSize/2.0;
    RelationXi[1] = tanhFunction( L , MinSize , a, i-2 , N) - MinSize/2.0;
    RelationXi[2] = tanhFunction( L , MinSize , a, i-1 , N) - MinSize/2.0;
    RelationXi[3] = tanhFunction( L , MinSize , a, i , N) - MinSize/2.0;
    RelationXi[4] = tanhFunction( L , MinSize , a, i+1 , N) - MinSize/2.0;
    RelationXi[5] = tanhFunction( L , MinSize , a, i+2 , N) - MinSize/2.0;
    RelationXi[6] = tanhFunction( L , MinSize , a, i+3 , N) - MinSize/2.0;
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
    XiPara[0][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[0],  RelationXi[1],  RelationXi[2] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[1][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[1],  RelationXi[0],  RelationXi[2] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[2][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[2],  RelationXi[0],  RelationXi[1] , RelationXi[3], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[3][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[3],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[4], RelationXi[5], RelationXi[6]); 
    XiPara[4][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[4],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[5], RelationXi[6]); 
    XiPara[5][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[5],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[4], RelationXi[6]); 
    XiPara[6][index_xi+r*NZ6] = Lagrange_6th(pos_z, RelationXi[6],  RelationXi[0],  RelationXi[1] , RelationXi[2], RelationXi[3], RelationXi[4], RelationXi[5]);    
}//pos_xi為換算過後的無因次化Z座標 




void GetXiParameter(double* XiPara_h[7], double pos_z, double pos_y, int index_xi, int j, int k ) 
{   
    //越界防呆
    if(j<3) j = j + NY6-7 ; 
    if(j>NY6-4) j = j - (NY6-7) ;
    //Z方向係數編號系統(第二格)  = (j+relation)*NZ6 + k  ; relation:1~7
    double L = LZ - HillFunction(pos_y) - minSize;
    //每一個y位置而言每一個計算間都不一樣 
    double pos_xi = LXi * (pos_z - (HillFunction(pos_y)+minSize/2.0)) / L;
    double a = GetNonuniParameter();
     //求解目標點的編號(映射回均勻系統)
    double j_cont = Inverse_tanh_index( pos_xi , L , minSize , GetNonuniParameter() , (NZ6-7) );
    //Z方向成員起始編號
    int cell_z1 = k-3;
    if( k <= 6 ) cell_z1 = 3;
    if( k >= NZ6-7 ) cell_z1 = NZ6-10; 
    //給我一個double 給出相對應七個內插的無因次化座標
    double RelationXi_0[7] ; //(j-3)
    double LT = LZ - HillFunction(y_global[j-3]) - minSize;
    double pos_z =  tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_0);//寫入第一套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_0 , 0 , index_xi); //XiPara_h[0][index_xi + 0*NZ6] ~ XiPara_h[6][index_xi + 0*NZ6]
    
    double RelationXi_1[7] ; //(j-2)
    LT = LZ - HillFunction(y_global[j-2]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_1);//寫入第二套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_1 , 1 , index_xi); //XiPara_h[0][index_xi + 1*NZ6] ~ XiPara_h[6][index_xi + 1*NZ6]
    
    double RelationXi_2[7] ; //(j-1)
    LT = LZ - HillFunction(y_global[j-1]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_2);//寫入第三套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_2 , 2 , index_xi); //XiPara_h[0][index_xi + 2*NZ6] ~ XiPara_h[6][index_xi + 2*NZ6]
    
    double RelationXi_3[7] ; //(j)
    LT = LZ - HillFunction(y_global[j]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_3);//寫入第四套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_3 , 3 , index_xi); //XiPara_h[0][index_xi + 3*NZ6] ~ XiPara_h[6][index_xi + 3*NZ6]

    double RelationXi_4[7] ; //(j+1)
    LT = LZ - HillFunction(y_global[j+1]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_4);//寫入第五套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_4 , 4 , index_xi); //XiPara_h[0][index_xi + 4*NZ6] ~ XiPara_h[6][index_xi + 4*NZ6]

    double RelationXi_5[7] ; //(j+2)
    LT = LZ - HillFunction(y_global[j+2]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_5);//寫入第六套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_5 , 5 , index_xi); //XiPara_h[0][index_xi + 5*NZ6] ~ XiPara_h[6][index_xi + 5*NZ6]
   
    double RelationXi_6[7] ; //(j+3)
    LT = LZ - HillFunction(y_global[j+3]) - minSize;
    pos_z = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_6);//寫入第七套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z , RelationXi_6 , 6 , index_xi); //XiPara_h[0][index_xi + 6*NZ6] ~ XiPara_h[6][index_xi + 6*NZ6]
}
```
---

## 2026-01-23 Xi 權重覆寫問題（今日討論）

### 問題來源（關鍵程式段）

1) `initializationTool.h` / `GetParameter_6th2(...)`

```cpp
XiPara[2][index_xi + r*NZ6] = Lagrange_6th(...);
```

2) `initialization.h` / `GetXiParameter(...)`

```cpp
GetParameter_6th2(XiPara_h, pos_z2, RelationXi_0, 0, index_xi);
GetParameter_6th2(XiPara_h, pos_z2, RelationXi_1, 1, index_xi);
...
GetParameter_6th2(XiPara_h, pos_z2, RelationXi_6, 6, index_xi);
```

3) `initialization.h` / `GetIntrplParameter_Xi()` 與 `BFLInitialization()`

```cpp
GetXiParameter(XiParaF1_h, ...);
...
GetXiParameter(XiParaF8_h, ...);
```

### 覆寫是什麼意思
- 不同空間點寫到同一個權重位置，後寫覆蓋先寫（不是越界）。

### 為什麼會覆寫（例子）
因為 `index_xi = j*NZ6 + k`，實際寫入位置是：

```
write_idx = index_xi + r*NZ6 = (j + r)*NZ6 + k
```

所以會發生：

```
(j=3, k=3, r=1) -> write_idx = (4,3)
(j=4, k=3, r=0) -> write_idx = (4,3)
```

因此 `(j=3, r=1)` 的權重會被 `(j=4, r=0)` 覆蓋。

### 結論
- 目前的 `index_xi + r*NZ6` 只是把 r 映射成 j 位移，等同把 7 層權重疊在同一張 2D 格子上，會系統性覆寫。


## 2026-01-23 為避免記憶體爆炸，Y方向只用三點格式降階
1.先調整 RelaxationXi 陣列 

---

## 2026-01-24 邊界處理修正與數值穩定性工作

### 問題診斷過程

#### 1. 初始問題：程式發散 (t=10 開始出現 nan)
- **症狀**：`rho = -54151.9` (t=10)，迅速發散
- **除錯輸出**：在 k=152 (上邊界 NZ6-4) 發現 F4 為負值 (-1.38293)
- **根本原因**：F4 (-Z方向) 的來源位置是 z+Δ，當 k=152 時會超出計算域

#### 2. cellZ_F* 計算問題
- **問題**：`GetXiParameter` 中 `cell_z1 = k-3` 使用的是目標點 k，而非來源點
- **修正**：改為基於 `j_cont`（映射後的連續座標）計算
```cpp
// 舊版
int cell_z1 = k-3;

// 新版  
int k_source = static_cast<int>(std::round(j_cont)) + 3;
int cell_z1 = k_source - 3;
if(cell_z1 < 0) cell_z1 = 0;
if(cell_z1 > NZ6-7) cell_z1 = NZ6-7;
```

#### 3. 邊界區域的插值外推問題
- **問題**：7-point stencil 在邊界附近會存取 buffer 區域，導致 Lagrange 外推產生極端權重
- **影響範圍**：
  - 下邊界 (k ≤ 5)：F2, F5, F6 的來源可能越界
  - 上邊界 (k ≥ NZ6-6=150)：F4, F7, F8 的來源可能越界
  - F1, F3 的 Xi stencil 也會受影響

### 最終解決方案：邊界區域使用簡化 streaming

```cpp
// evolution.h 中的邊界處理邏輯
if( k <= 5 ) {
    // 下邊界附近：使用簡單 streaming 和 bounce-back
    F1_in = f1_old[(j-1)*NZ6 + k];
    F3_in = f3_old[(j+1)*NZ6 + k];
    F2_in = f4_old[idx_xi];  // bounce-back
    F5_in = f7_old[idx_xi];  // bounce-back
    F6_in = f8_old[idx_xi];  // bounce-back
    F4_Intrpl7(...);  // 正常插值（遠離下邊界）
    F7_in = f7_old[(j+1)*NZ6 + k+1];
    F8_in = f8_old[(j-1)*NZ6 + k+1];
    
} else if( k >= NZ6-6 ) {
    // 上邊界附近：使用簡單 streaming 和 bounce-back
    F1_in = f1_old[(j-1)*NZ6 + k];
    F3_in = f3_old[(j+1)*NZ6 + k];
    F2_in = f2_old[j*NZ6 + k-1];      // 簡單 streaming
    F5_in = f5_old[(j-1)*NZ6 + k-1];  // 簡單 streaming
    F6_in = f6_old[(j+1)*NZ6 + k-1];  // 簡單 streaming
    F4_in = f2_old[idx_xi];  // bounce-back
    F7_in = f5_old[idx_xi];  // bounce-back
    F8_in = f6_old[idx_xi];  // bounce-back
    
} else {
    // 內部區域：正常 ISLBM 插值
    F1_Intrpl3(...); F3_Intrpl3(...);
    F2_Intrpl7(...); F4_Intrpl7(...);
    Y_XI_Intrpl3(...);  // F5~F8
}
```

### 修改的檔案清單

| 檔案 | 修改內容 |
|------|----------|
| `initialization.h` | `GetXiParameter` 中 cell_z1 計算改用 j_cont |
| `evolution.h` | 新增邊界區域條件判斷，使用簡化 streaming 替代插值 |

### 測試結果
- ✅ 程式穩定運行至 t=1000+
- ✅ 密度收斂到 ~0.988
- ⚠️ 出現 **Checkerboard instability**（棋盤格振盪）

### Checkerboard 問題分析

**觀察**：
- Y 方向七階插值 → 條紋較粗
- Y 方向三階插值 → 條紋較密

**可能原因**：
1. 奇偶解耦 (odd-even decoupling)
2. 非均勻網格的數值色散
3. Lagrange 插值的振盪特性

**待嘗試解決方案**：
- [ ] Z 方向也降為三階插值
- [ ] 調整 MRT 鬆弛參數
- [ ] 加入數值濾波器

### 下次工作方向

1. **Z 方向降階測試**
   - 新增 `F2_Intrpl3` 和 `F4_Intrpl3` 巨集
   - 修改 `evolution.h` 使用三階版本
   
2. **評估降階效果**
   - 若 checkerboard 消失且結果合理 → 繼續使用
   - 若精度不足 → 考慮其他抑制方法

### 目前配置摘要

| 項目 | 設定值 |
|------|--------|
| Y 方向插值 | 3 階 Lagrange |
| Z 方向插值 | 7 階 Lagrange |
| 下邊界處理 (k≤5) | 簡化 streaming + bounce-back |
| 上邊界處理 (k≥150) | 簡化 streaming + bounce-back |
| 內部區域 | ISLBM 插值 |
| Re | 1 |
| CFL | 0.2 |
| 網格 | NY=200, NZ=150 (含 buffer: NY6=207, NZ6=156) |

---

## 2026-01-28 重大突破：動態 Streaming 邊界與 Re=500 穩定運行

### 今日核心成果

| 項目 | 狀態 | 說明 |
|------|------|------|
| **Re=500 穩定運行** | ✅ | 成功跑到 t=39000+ 無發散 |
| **動態 Streaming 邊界** | ✅ | 漸進式擴大解析層設計完成 |
| **質量通量誤差** | ✅ | 改用 div(ρu) 計算局部質量守恆 |

---

### 1. 重大發現：Streaming Layer 與 Interpolation 邊界的獨立性

#### 關鍵洞察

之前錯誤認為 `streaming_lower` 必須 >= `interpolation_lower`，但實際上**兩者是獨立的**！

| 參數 | 作用 | 計算時機 | 可否動態調整 |
|------|------|----------|--------------|
| `interpolation_lower/upper` | 決定 XiPara[] 預存的插值權重類型（3點或7點）| **初始化時**預計算 | 否（需重算權重）|
| `streaming_lower/upper` | 決定是否跳過插值改用 streaming | **每個時間步**即時判斷 | ✅ 可以動態調整 |

#### 邏輯說明

```cpp
// interpolation_lower = 25 意味著：
// k < 25：XiPara 存「三點插值」權重
// k >= 25：XiPara 存「七點插值」權重

// streaming_lower 意味著：
// k <= streaming_lower：跳過插值，用 streaming（最保守）
// k > streaming_lower：讀取 XiPara[] 做插值
```

#### 重要修正

**之前的錯誤認知**：`streaming_lower >= interpolation_lower`（必須）

**正確理解**：`streaming_lower` 可以 < `interpolation_lower`！
- 這樣三點插值區就能作為**穩定緩衝區**發揮作用
- 漸進開放順序：streaming → 三點插值 → 七點插值

---

### 2. 動態 Streaming 邊界設計

#### 設計原理

```
時間進展 →

t=0:      |===streaming (k≤50)===|---七點插值---|===streaming===|
                                                 
t=50000:  |==streaming (k≤30)==|--三點--|--七點--|--三點--|==streaming==|

t=100000: |=s(k≤10)=|--三點--|-------七點插值 (完整)-------|--三點--|=s=|
                     ↑                                      ↑
                緩衝區開放！                            緩衝區開放！
```

#### 最終解析層結構 (t >= 100000)

```
k=0    k=10       k=25                    k=236      k=255  k=262
|======|==========|========================|==========|======|
stream   三點插值        七點插值          三點插值    stream
 (10層)  (15層緩衝)       (主解析區)         (19層緩衝)  (7層)
```

#### 參數設定（variables.h）

```cpp
//=== 動態 Streaming 邊界參數（漸進式擴大解析層）===//
// 固定的插值邊界（決定權重類型，初始化時計算）
#define     interpolation_lower  (25)                // 七點內插下界
#define     interpolation_upper  (NZ6-26)            // 七點內插上界

// 初始值（保守，更大的 streaming 區域）
#define     streaming_lower_init     (50)            // 初始下界 (k <= 50 用 streaming)
#define     streaming_upper_init     (NZ6-51)        // 初始上界 (k >= NZ6-51 用 streaming)

// 目標值（比 interpolation 更激進，開放三點插值緩衝區）
#define     streaming_lower_target   (10)            // 目標下界
#define     streaming_upper_target   (NZ6-7)         // 目標上界

// 過渡時間設定
#define     ramp_start_time          (0)             // 開始漸進的時間步
#define     ramp_end_time            (100000)        // 完成漸進的時間步

// 全域變數宣告（在 main.cpp 中定義）
extern int streaming_lower;  // 動態下界
extern int streaming_upper;  // 動態上界
```

#### 更新函數實作（main.cpp）

```cpp
//-----------------------------------------------------------------------------
// 2.12 動態 Streaming 邊界（漸進式擴大解析層）
//-----------------------------------------------------------------------------
int streaming_lower = streaming_lower_init;  // 動態下界，初始為保守值
int streaming_upper = streaming_upper_init;  // 動態上界，初始為保守值

// 使用 tanh 平滑過渡更新 streaming 邊界
void UpdateStreamingBounds(int t) {
    if (t >= ramp_end_time) {
        // 過渡完成，使用目標值
        streaming_lower = streaming_lower_target;
        streaming_upper = streaming_upper_target;
    } else if (t <= ramp_start_time) {
        // 尚未開始，使用初始值
        streaming_lower = streaming_lower_init;
        streaming_upper = streaming_upper_init;
    } else {
        // 過渡期：使用 tanh 平滑過渡
        double progress = (double)(t - ramp_start_time) / (ramp_end_time - ramp_start_time);
        // tanh 平滑：將 [0,1] 映射到 [0,1]，但中間過渡更平滑
        double smooth_ratio = 0.5 * (1.0 + tanh(6.0 * (progress - 0.5)));
        
        // 計算當前邊界值
        streaming_lower = streaming_lower_init - 
            (int)(smooth_ratio * (streaming_lower_init - streaming_lower_target));
        streaming_upper = streaming_upper_init + 
            (int)(smooth_ratio * (streaming_upper_target - streaming_upper_init));
    }
}
```

#### 時間迴圈中的呼叫（main.cpp）

```cpp
for(t = 0; t < loop; t++) {
    // 更新動態 streaming 邊界（漸進式擴大解析層）
    UpdateStreamingBounds(t);
    
    // ... 其餘時間迴圈邏輯 ...
}
```

#### 輸出監控（main.cpp）

```cpp
// 輸出當前 streaming 邊界（顯示解析層擴大進度）
cout << "[Streaming Bounds] lower=" << streaming_lower 
     << " upper=" << streaming_upper 
     << " (target: " << streaming_lower_target << "/" << streaming_upper_target << ")" << endl;
```

---

### 3. 質量通量誤差計算改進（evolution.h）

#### 舊版：密度和方法

```cpp
// 舊版 - 簡單密度和
double ComputeMaxLocalMassError(double* rho_d) {
    double max_error = 0.0;
    for(int j = 4; j < NY6-4; j++) {
        for(int k = streaming_lower+1; k < streaming_upper-1; k++) {
            double local_mass = 0.0;
            for(int dj = -1; dj <= 1; dj++) {
                for(int dk = -1; dk <= 1; dk++) {
                    local_mass += rho_d[(j + dj) * NZ6 + (k + dk)];
                }
            }
            double error = std::fabs(local_mass - 9.0) / 9.0;
            if(error > max_error) max_error = error;
        }
    }
    return max_error;
}
```

#### 新版：質量通量散度方法

```cpp
//==========================================
//4.計算區域質量守恆 - 質量通量淨流入流出差
// 對每個格點計算: |流出質量通量 - 流入質量通量|
// 質量通量 = ρ*v (Y方向) + ρ*w (Z方向)
// 使用中心差分: div(ρu) ≈ [(ρv)_{j+1} - (ρv)_{j-1}]/2 + [(ρw)_{k+1} - (ρw)_{k-1}]/2
//==========================================
double ComputeMassFluxError(double* rho_d, double* v_field, double* w_field) {
    double total_error = 0.0;
    
    for(int j = 4; j < NY6-4; j++) {
        for(int k = 4; k < NZ6-4; k++) {
            int idx = j * NZ6 + k;
            int idx_jp1 = (j+1) * NZ6 + k;  // j+1
            int idx_jm1 = (j-1) * NZ6 + k;  // j-1
            int idx_kp1 = j * NZ6 + (k+1);  // k+1
            int idx_km1 = j * NZ6 + (k-1);  // k-1
            
            // Y方向質量通量差 (ρv)_{j+1} - (ρv)_{j-1}
            double flux_y = (rho_d[idx_jp1] * v_field[idx_jp1]) - (rho_d[idx_jm1] * v_field[idx_jm1]);
            
            // Z方向質量通量差 (ρw)_{k+1} - (ρw)_{k-1}
            double flux_z = (rho_d[idx_kp1] * w_field[idx_kp1]) - (rho_d[idx_km1] * w_field[idx_km1]);
            
            // 質量守恆誤差 = |div(ρu)| = |∂(ρv)/∂y + ∂(ρw)/∂z|
            double div_rho_u = std::fabs(flux_y + flux_z);
            total_error += div_rho_u;
        }
    }
    
    // 回傳總質量通量誤差（除以格點數得到平均值）
    int num_cells = (NY6 - 8) * (NZ6 - 8);
    return total_error / (double)num_cells;
}
```

---

### 4. Re=500 穩定運行記錄

#### 模擬參數

| 參數 | 值 | 說明 |
|------|-----|------|
| Re | 500 | 雷諾數 |
| tau | 0.6833 | 鬆弛時間 |
| CFL | 0.8 | CFL 數 |
| NY | 512 | Y 方向網格數 |
| NZ | 256 | Z 方向網格數 |
| NY6 | 519 | 含 buffer 的 Y 網格數 |
| NZ6 | 262 | 含 buffer 的 Z 網格數 |
| Ma_theoretical | 0.336666 | 理論 Mach 數（超過 Ma_max=0.3）|

#### 運行輸出（t=39000 時）

```
Time=39000 ; Average Density=1 ; Density Correction=... ; Mass Flux Err=...
[Streaming Bounds] lower=25 upper=236 (target: 10/255)
[t=39000] Mach stats: max=0.3 (j=...,k=...), avg=0.17...
```

#### 流場分析結果

```
============================================================
Global Velocity Statistics (t=34000~39000)
============================================================
    Time     Uy_max     Uy_min    Uy_mean     Uz_max     Uz_min
------------------------------------------------------------
   34000    0.17320   -0.01115    0.10711    0.05660   -0.04676
   35000    0.17320   -0.01402    0.10701    0.05144   -0.04878
   36000    0.17320   -0.01578    0.10697    0.04710   -0.05055
   37000    0.17320   -0.01486    0.10728    0.04853   -0.05171
   38000    0.17320   -0.01069    0.10776    0.05344   -0.05275
   39000    0.17320   -0.01112    0.10809    0.05523   -0.05398

============================================================
Recirculation Zone Analysis
============================================================
k= 12: Backflow region Y=[0.65, 7.77], 406 points
k= 25: Backflow region Y=[0.70, 3.87], 181 points
k= 42: Backflow region Y=[1.21, 2.14], 54 points
k= 64: No backflow detected
```

#### 關鍵觀察

1. **Uy_max = 0.17320** 被 Ma_max 限制住
2. **回流區明確存在**：靠近底部有大範圍回流（符合 Periodic Hill 物理特徵）
3. **密度穩定在 ~1.0**，質量守恆良好

---

### 5. 從 Edit6 初始版本到今天的演進歷程

#### Edit6 初始版本特點
- 自適應性內插權重使用全域變數賦值
- `streaming_lower/upper` 為靜態 `#define` 常數
- 固定邊界，無法動態調整

#### 今日改進

| 改動項目 | 舊版 | 新版 |
|----------|------|------|
| streaming 邊界 | `#define streaming_lower (25)` 靜態 | `extern int streaming_lower` 動態全域變數 |
| 邊界更新 | 無 | `UpdateStreamingBounds(t)` 隨時間 tanh 平滑過渡 |
| 質量誤差計算 | 密度和方法 | 質量通量散度 div(ρu) |
| 解析層策略 | 固定 | 漸進式擴大（streaming → 三點 → 七點）|

---

### 6. 修改的檔案完整清單

#### variables.h 修改

```cpp
// === 原本 (靜態) ===
#define     streaming_lower      (25)
#define     streaming_upper      (NZ6-26)

// === 改為 (動態) ===
#define     interpolation_lower  (25)
#define     interpolation_upper  (NZ6-26)

//=== 動態 Streaming 邊界參數（漸進式擴大解析層）===//
#define     streaming_lower_init     (50)
#define     streaming_upper_init     (NZ6-51)
#define     streaming_lower_target   (10)
#define     streaming_upper_target   (NZ6-7)
#define     ramp_start_time          (0)
#define     ramp_end_time            (100000)

extern int streaming_lower;
extern int streaming_upper;
```

#### main.cpp 修改

```cpp
// === 新增全域變數與更新函數 ===
int streaming_lower = streaming_lower_init;
int streaming_upper = streaming_upper_init;

void UpdateStreamingBounds(int t) {
    if (t >= ramp_end_time) {
        streaming_lower = streaming_lower_target;
        streaming_upper = streaming_upper_target;
    } else if (t <= ramp_start_time) {
        streaming_lower = streaming_lower_init;
        streaming_upper = streaming_upper_init;
    } else {
        double progress = (double)(t - ramp_start_time) / (ramp_end_time - ramp_start_time);
        double smooth_ratio = 0.5 * (1.0 + tanh(6.0 * (progress - 0.5)));
        streaming_lower = streaming_lower_init - 
            (int)(smooth_ratio * (streaming_lower_init - streaming_lower_target));
        streaming_upper = streaming_upper_init + 
            (int)(smooth_ratio * (streaming_upper_target - streaming_upper_init));
    }
}

// === 時間迴圈中呼叫 ===
for(t = 0; t < loop; t++) {
    UpdateStreamingBounds(t);  // 新增
    // ...
}

// === 輸出監控 ===
cout << "[Streaming Bounds] lower=" << streaming_lower 
     << " upper=" << streaming_upper 
     << " (target: " << streaming_lower_target << "/" << streaming_upper_target << ")" << endl;
```

#### evolution.h 修改

```cpp
// === 改為質量通量誤差 ===
double ComputeMassFluxError(double* rho_d, double* v_field, double* w_field) {
    double total_error = 0.0;
    for(int j = 4; j < NY6-4; j++) {
        for(int k = 4; k < NZ6-4; k++) {
            int idx = j * NZ6 + k;
            int idx_jp1 = (j+1) * NZ6 + k;
            int idx_jm1 = (j-1) * NZ6 + k;
            int idx_kp1 = j * NZ6 + (k+1);
            int idx_km1 = j * NZ6 + (k-1);
            double flux_y = (rho_d[idx_jp1] * v_field[idx_jp1]) - (rho_d[idx_jm1] * v_field[idx_jm1]);
            double flux_z = (rho_d[idx_kp1] * w_field[idx_kp1]) - (rho_d[idx_km1] * w_field[idx_km1]);
            double div_rho_u = std::fabs(flux_y + flux_z);
            total_error += div_rho_u;
        }
    }
    int num_cells = (NY6 - 8) * (NZ6 - 8);
    return total_error / (double)num_cells;
}

// === main.cpp 呼叫改為 ===
" ; Mass Flux Err=" << ComputeMassFluxError(rho, v, w) << endl;
```

---

### 7. 下一步工作方向

1. **測試漸進式邊界效果**
   - 觀察 streaming bounds 隨時間變化的輸出
   - 確認在過渡期間模擬穩定

2. **挑戰更高雷諾數**
   - 若 Re=500 穩定後，嘗試 Re=700, Re=1000
   - 可能需要調整 ramp_end_time 或初始邊界值

3. **評估三點插值緩衝區效果**
   - 比較有/無緩衝區的穩定性差異
   - 分析流場精度影響

---

## 2026-01-28 修正：分階段過渡解決 t=53000 崩潰問題

### 問題描述

模擬在 t≈53000 時連續崩潰兩次。分析發現：

```
t=53000 時：
- progress = 53000/100000 = 0.53
- smooth_ratio ≈ 0.5 * (1 + tanh(6*(0.53-0.5))) ≈ 0.59
- streaming_lower = 50 - 0.59*(50-10) ≈ 26
```

此時 `streaming_lower ≈ 26` 正好接近 `interpolation_lower = 25`，這是一個**臨界過渡點**！

### 問題根因

原設計使用單階段過渡：
- streaming 50 → 10（跨越了三點插值區和七點插值區的邊界）
- 在 streaming_lower 接近 interpolation_lower 時，系統從 streaming 直接跳到使用七點插值
- 這個突變可能導致數值不穩定

### 解決方案：分階段過渡

#### 設計原理

```
=== 第一階段 (t=0 ~ 100000) ===
streaming 50 → 25（開放七點插值區）
   
t=0:      |====streaming (k≤50)====|---七點插值---|====streaming====|
t=100000: |==streaming (k≤25)==|-------七點插值-------|==streaming==|
                               ↑                       ↑
                          interpolation_lower    interpolation_upper

=== 第二階段 (t=100000 ~ 200000) ===
streaming 25 → 10（開放三點插值緩衝區）

t=100000: |==streaming (k≤25)==|-------七點插值-------|==streaming==|
t=200000: |=s(k≤10)=|--三點--|-------七點插值-------|--三點--|=s=|
              ↑    ↑                                 ↑      ↑
           target  interpolation_lower        interpolation_upper
```

#### 參數修改（variables.h）

```cpp
//=== 動態 Streaming 邊界參數（分階段漸進式擴大解析層）===//
// 初始值（保守，更大的 streaming 區域）
#define     streaming_lower_init     (50)            // 初始下界 (k <= 50 用 streaming)
#define     streaming_upper_init     (NZ6-51)        // 初始上界 (k >= NZ6-51 用 streaming)

// === 第一階段：開放七點插值區 (streaming → interpolation_lower) ===
#define     streaming_lower_phase1   (interpolation_lower)  // 第一階段目標: 25
#define     streaming_upper_phase1   (interpolation_upper)  // 第一階段目標: NZ6-26
#define     phase1_start_time        (0)             // 第一階段開始
#define     phase1_end_time          (100000)        // 第一階段結束

// === 第二階段：開放三點插值緩衝區 (interpolation_lower → target) ===
#define     streaming_lower_target   (10)            // 最終目標下界
#define     streaming_upper_target   (NZ6-7)         // 最終目標上界
#define     phase2_start_time        (100000)        // 第二階段開始
#define     phase2_end_time          (200000)        // 第二階段結束
```

#### 更新函數實作（main.cpp）

```cpp
// 使用 tanh 平滑過渡更新 streaming 邊界（分階段版本）
// 第一階段 (t=0~100000)：streaming 50→25（開放七點插值區）
// 第二階段 (t=100000~200000)：streaming 25→10（開放三點插值緩衝區）
void UpdateStreamingBounds(int t) {
    if (t >= phase2_end_time) {
        // 第二階段完成，使用最終目標值
        streaming_lower = streaming_lower_target;
        streaming_upper = streaming_upper_target;
    } else if (t >= phase2_start_time) {
        // 第二階段：從 interpolation 邊界過渡到最終目標（開放三點插值區）
        double progress = (double)(t - phase2_start_time) / (phase2_end_time - phase2_start_time);
        double smooth_ratio = 0.5 * (1.0 + tanh(6.0 * (progress - 0.5)));
        
        streaming_lower = streaming_lower_phase1 - 
            (int)(smooth_ratio * (streaming_lower_phase1 - streaming_lower_target));
        streaming_upper = streaming_upper_phase1 + 
            (int)(smooth_ratio * (streaming_upper_target - streaming_upper_phase1));
    } else if (t >= phase1_start_time) {
        // 第一階段：從初始值過渡到 interpolation 邊界（開放七點插值區）
        double progress = (double)(t - phase1_start_time) / (phase1_end_time - phase1_start_time);
        double smooth_ratio = 0.5 * (1.0 + tanh(6.0 * (progress - 0.5)));
        
        streaming_lower = streaming_lower_init - 
            (int)(smooth_ratio * (streaming_lower_init - streaming_lower_phase1));
        streaming_upper = streaming_upper_init + 
            (int)(smooth_ratio * (streaming_upper_phase1 - streaming_upper_init));
    } else {
        // 尚未開始，使用初始值
        streaming_lower = streaming_lower_init;
        streaming_upper = streaming_upper_init;
    }
}
```

### 設計優點

1. **避免跨越臨界點**：第一階段結束在 interpolation_lower，第二階段才繼續向下
2. **充分穩定時間**：每階段各 100000 步，讓系統充分適應
3. **保留 tanh 平滑過渡**：避免階梯跳變帶來的數值衝擊
4. **漸進式開放**：先開放七點插值區（精度高），穩定後再開放三點緩衝區

### 預期時間線

| 時間 | streaming_lower | streaming_upper | 開放區域 |
|------|-----------------|-----------------|----------|
| t=0 | 50 | 211 | streaming only |
| t=50000 | 30~35 | 225~230 | 七點插值區部分開放 |
| t=100000 | 25 | 236 | 七點插值區完全開放 |
| t=150000 | 16~20 | 245~250 | 三點緩衝區部分開放 |
| t=200000 | 10 | 255 | 完全開放（最終狀態） |

---

### Git Commit 記錄

```
commit e5f3c1b
Author: ChenPengChung
Date:   2026-01-28

feat: 重大發現 - 動態 Streaming 邊界漸進式擴大解析層

## 核心發現：Streaming Layer 與 Interpolation 邊界的獨立性
- interpolation_lower/upper: 決定 XiPara[] 預存什麼類型的插值權重（初始化時）
- streaming_lower/upper: 決定是否跳過插值改用 streaming（每步即時判斷，可動態調整）
- 重要修正：streaming_lower 可以 < interpolation_lower，讓三點插值區作為穩定緩衝

## 漸進式設計
- 使用 tanh 平滑過渡函數避免階梯跳變
- 初始值: streaming_lower=50, streaming_upper=NZ6-51
- 目標值: streaming_lower=10, streaming_upper=NZ6-7
- 過渡時間: t=0 ~ t=100000

## 其他改進
- 質量通量誤差計算改用 div(ρu) 公式

3 files changed, 81 insertions(+), 18 deletions(-)
```




目前的問題是，當我CFL下降，本來邊界不墜震盪，都會出問題，相較於ＣＦＬ＝0.8 另外，當時間步推進時，Ma會提高最後平均超過0.3時，整個結果開始出現震盪，因此開始結果不可信，有什麼方法可以抑制馬赫數換衝，難道目前限速氣是唯一解方，但是限速器會不正確地影想外力修正，也並非好方法，木前只做了家溺往格256->512 但是根本馬赫數飆高導致結果不可信的問題仍然無法解決ThiukHard (不要改code先討論)


```