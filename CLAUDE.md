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
