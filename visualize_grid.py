#!/usr/bin/env python3
"""
visualize_grid.py - 輸出初始化下的二維非均勻網格結構 VTK 檔案
用於 ParaView 可視化週期性山丘流場的非均勻網格

生成檔案：
  output/grid_structure.vtk       - STRUCTURED_GRID 格式（含網格品質指標）
  output/grid_quadcells.vtk       - UNSTRUCTURED_GRID 格式（每個 cell 獨立，便於檢視）

使用方法：
  python visualize_grid.py

ParaView 可視化步驟：
  1. 開啟 ParaView → File → Open → 選擇 output/grid_structure.vtk
  2. 點擊 Apply
  3. 在 Representation 下拉選單選 "Surface With Edges" 或 "Wireframe" 即可看到網格線
  4. 可用 Color By → CellSize / AspectRatio 檢查網格品質
"""

import numpy as np
import os
import math

# =============================================================================
# 參數定義（對應 variables.h）
# =============================================================================
LY = 9.0
LZ = 3.036
NY = 32
NZ = 16
NY6 = NY + 7    # 39
NZ6 = NZ + 6    # 22
CFL = 0.65
minSize = (LZ - 1.0) / (NZ6 - 6) * CFL
print(f"[設定] NY={NY}, NZ={NZ}, NY6={NY6}, NZ6={NZ6}, minSize={minSize:.8f}")
BUFFER = 3

# =============================================================================
# 山丘函數（對應 model.h 的 HillFunction）
# =============================================================================
def HillFunction(Y):
    """週期性山丘幾何函數，輸入 Y 座標，回傳山丘高度"""
    if Y < 0.0:
        Yb = Y + LY
    elif Y > LY:
        Yb = Y - LY
    else:
        Yb = Y

    model = 0.0
    # 左半丘 (6 段多項式)
    if Yb <= (54.0/28.0) * (9.0/54.0):
        model = (1.0/28.0) * min(28.0,
            28.0 + 0.006775070969851*(Yb*28)*(Yb*28)
            - 0.0021245277758000*(Yb*28)*(Yb*28)*(Yb*28))
    elif Yb <= (54.0/28.0) * (14.0/54.0):
        model = (1.0/28.0) * (25.07355893131
            + 0.9754803562315*(Yb*28)
            - 0.1016116352781*(Yb*28)*(Yb*28)
            + 0.001889794677828*(Yb*28)*(Yb*28)*(Yb*28))
    elif Yb <= (54.0/28.0) * (20.0/54.0):
        model = (1.0/28.0) * (25.79601052357
            + 0.8206693007457*(Yb*28)
            - 0.09055370274339*(Yb*28)*(Yb*28)
            + 0.001626510569859*(Yb*28)*(Yb*28)*(Yb*28))
    elif Yb <= (54.0/28.0) * (30.0/54.0):
        model = (1.0/28.0) * (40.46435022819
            - 1.379581654948*(Yb*28)
            + 0.019458845041284*(Yb*28)*(Yb*28)
            - 0.0002070318932190*(Yb*28)*(Yb*28)*(Yb*28))
    elif Yb <= (54.0/28.0) * (40.0/54.0):
        model = (1.0/28.0) * (17.92461334664
            + 0.8743920332081*(Yb*28)
            - 0.05567361123058*(Yb*28)*(Yb*28)
            + 0.0006277731764683*(Yb*28)*(Yb*28)*(Yb*28))
    elif Yb <= (54.0/28.0) * (54.0/54.0):
        model = (1.0/28.0) * max(0.0,
            56.39011190988
            - 2.010520359035*(Yb*28)
            + 0.01644919857549*(Yb*28)*(Yb*28)
            + 0.00002674976141766*(Yb*28)*(Yb*28)*(Yb*28))

    # 右半丘 (對稱)
    Yr = LY - Yb
    if Yr >= 0 and Yr <= (54.0/28.0):
        if Yr <= (54.0/28.0) * (9.0/54.0):
            model = (1.0/28.0) * min(28.0,
                28.0 + 0.006775070969851*(Yr*28)*(Yr*28)
                - 0.0021245277758000*(Yr*28)*(Yr*28)*(Yr*28))
        elif Yr <= (54.0/28.0) * (14.0/54.0):
            model = (1.0/28.0) * (25.07355893131
                + 0.9754803562315*(Yr*28)
                - 0.1016116352781*(Yr*28)*(Yr*28)
                + 0.001889794677828*(Yr*28)*(Yr*28)*(Yr*28))
        elif Yr <= (54.0/28.0) * (20.0/54.0):
            model = (1.0/28.0) * (25.79601052357
                + 0.8206693007457*(Yr*28)
                - 0.09055370274339*(Yr*28)*(Yr*28)
                + 0.001626510569859*(Yr*28)*(Yr*28)*(Yr*28))
        elif Yr <= (54.0/28.0) * (30.0/54.0):
            model = (1.0/28.0) * (40.46435022819
                - 1.379581654948*(Yr*28)
                + 0.019458845041284*(Yr*28)*(Yr*28)
                - 0.0002070318932190*(Yr*28)*(Yr*28)*(Yr*28))
        elif Yr <= (54.0/28.0) * (40.0/54.0):
            model = (1.0/28.0) * (17.92461334664
                + 0.8743920332081*(Yr*28)
                - 0.05567361123058*(Yr*28)*(Yr*28)
                + 0.0006277731764683*(Yr*28)*(Yr*28)*(Yr*28))
        elif Yr <= (54.0/28.0) * (54.0/54.0):
            model = (1.0/28.0) * max(0.0,
                56.39011190988
                - 2.010520359035*(Yr*28)
                + 0.01644919857549*(Yr*28)*(Yr*28)
                + 0.00002674976141766*(Yr*28)*(Yr*28)*(Yr*28))

    return model


# =============================================================================
# 非均勻網格函數（對應 initializationTool.h）
# =============================================================================
def tanhFunction(L, MinSize, a, j, N):
    """雙曲正切非均勻網格座標轉換"""
    return (L / 2.0 + MinSize / 2.0 +
            (L / 2.0 / a) * math.tanh(
                (-1.0 + 2.0 * j / N) / 2.0 * math.log((1.0 + a) / (1.0 - a))))


def GetNonuniParameter():
    """二分法求解非均勻網格伸縮參數 a"""
    total = LZ - HillFunction(0.0) - minSize
    a_low, a_high = 0.1, 1.0

    for _ in range(200):
        a_mid = (a_low + a_high) / 2.0
        x0 = tanhFunction(total, minSize, a_mid, 0, NZ6 - 7)
        x1 = tanhFunction(total, minSize, a_mid, 1, NZ6 - 7)
        dx = x1 - x0
        if dx - minSize >= 0.0:
            a_low = a_mid
        else:
            a_high = a_mid
        if abs(dx - minSize) < 1e-14:
            break
    return a_mid


# =============================================================================
# 網格生成（對應 initialization.h 的 GenerateMesh_Y / GenerateMesh_Z）
# =============================================================================
def generate_mesh():
    """生成完整的非均勻網格座標陣列"""
    # Y 方向均勻網格
    y_global = np.zeros(NY6)
    dy = LY / (NY6 - 2 * BUFFER - 1)
    for i in range(NY6):
        y_global[i] = dy * (i - BUFFER)

    # Z 方向非均勻網格
    nonuni_a = GetNonuniParameter()
    print(f"非均勻網格參數 a = {nonuni_a:.10f}")
    print(f"最小網格大小 minSize = {minSize:.10f}")

    z_global = np.zeros((NY6, NZ6))

    for j in range(NY6):
        hill_h = HillFunction(y_global[j])
        total = LZ - hill_h - minSize
        for k in range(BUFFER, NZ6 - BUFFER):
            z_global[j, k] = tanhFunction(total, minSize, nonuni_a,
                                          k - BUFFER, NZ6 - 7) + hill_h
        z_global[j, 2] = hill_h
        z_global[j, NZ6 - 3] = LZ

    return y_global, z_global, nonuni_a


# =============================================================================
# VTK 輸出 - STRUCTURED_GRID 格式（使用 POINT_DATA，與 flow VTK 相同格式）
# =============================================================================
def output_structured_grid_vtk(y_global, z_global, filename):
    """
    輸出 STRUCTURED_GRID 格式 VTK（使用 POINT_DATA）
    格式與 flow_000000.vtk 完全一致，確保 ParaView 可正確開啟
    ParaView: Representation → 'Surface With Edges' 即可看到網格線
    """
    ny_out = NY6 - 6  # 去掉 buffer
    nz_out = NZ6 - 6
    npoints = ny_out * nz_out

    dy = y_global[BUFFER + 1] - y_global[BUFFER]

    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Non-uniform Grid Structure for Periodic Hill\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write(f"DIMENSIONS {ny_out} {nz_out} 1\n")
        f.write(f"POINTS {npoints} double\n")

        # 點座標（與 flow VTK 相同迴圈順序: k外迴圈, j內迴圈）
        for k in range(BUFFER, NZ6 - BUFFER):
            for j in range(BUFFER, NY6 - BUFFER):
                f.write(f"{y_global[j]} {z_global[j, k]} 0.0\n")

        # --- POINT_DATA（與 flow VTK 相同格式）---
        f.write(f"\nPOINT_DATA {npoints}\n")

        # 1. 山丘高度（標記固/流體）
        f.write("\nSCALARS HillHeight double 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(BUFFER, NZ6 - BUFFER):
            for j in range(BUFFER, NY6 - BUFFER):
                f.write(f"{HillFunction(y_global[j])}\n")

        # 2. 節點處 Z 座標值（方便確認非均勻分佈）
        f.write("\nSCALARS Z_coordinate double 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(BUFFER, NZ6 - BUFFER):
            for j in range(BUFFER, NY6 - BUFFER):
                f.write(f"{z_global[j, k]}\n")

        # 3. 局部 dz（Z方向網格間距，差分近似到節點）
        f.write("\nSCALARS LocalDZ double 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(BUFFER, NZ6 - BUFFER):
            for j in range(BUFFER, NY6 - BUFFER):
                if k < NZ6 - BUFFER - 1:
                    dz = z_global[j, k + 1] - z_global[j, k]
                else:
                    dz = z_global[j, k] - z_global[j, k - 1]
                f.write(f"{dz}\n")

        # 4. 長寬比 Aspect Ratio（差分近似到節點）
        f.write("\nSCALARS AspectRatio double 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(BUFFER, NZ6 - BUFFER):
            for j in range(BUFFER, NY6 - BUFFER):
                if k < NZ6 - BUFFER - 1:
                    dz = z_global[j, k + 1] - z_global[j, k]
                else:
                    dz = z_global[j, k] - z_global[j, k - 1]
                if dz > 1e-15:
                    ar = max(dy / dz, dz / dy)
                else:
                    ar = 1e10
                f.write(f"{ar}\n")

        # 5. IsSolid（節點是否在山丘內部：1=固體，0=流體）
        f.write("\nSCALARS IsSolid int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(BUFFER, NZ6 - BUFFER):
            for j in range(BUFFER, NY6 - BUFFER):
                hill_h = HillFunction(y_global[j])
                is_solid = 1 if z_global[j, k] < hill_h else 0
                f.write(f"{is_solid}\n")

    print(f"已輸出 STRUCTURED_GRID VTK: {filename}")


# =============================================================================
# VTK 輸出 - UNSTRUCTURED_GRID 格式（每個 cell 獨立 quad）
# =============================================================================
def output_unstructured_grid_vtk(y_global, z_global, filename):
    """
    輸出 UNSTRUCTURED_GRID 格式 VTK
    每個 cell 為獨立四邊形，更清楚展示非均勻網格結構
    """
    ny_out = NY6 - 6
    nz_out = NZ6 - 6
    npoints = ny_out * nz_out
    ncells = (ny_out - 1) * (nz_out - 1)

    dy = y_global[BUFFER + 1] - y_global[BUFFER]

    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Non-uniform Grid Unstructured Quads\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {npoints} double\n")

        # 寫入所有節點（與 structured grid 相同順序: k外, j內）
        for k in range(BUFFER, NZ6 - BUFFER):
            for j in range(BUFFER, NY6 - BUFFER):
                f.write(f"{y_global[j]} {z_global[j, k]} 0.0\n")

        # 寫入 cell 連接 (Quad = VTK type 9)
        f.write(f"\nCELLS {ncells} {ncells * 5}\n")
        for k in range(nz_out - 1):
            for j in range(ny_out - 1):
                p0 = k * ny_out + j           # (j,   k)
                p1 = k * ny_out + (j + 1)     # (j+1, k)
                p2 = (k + 1) * ny_out + (j + 1)  # (j+1, k+1)
                p3 = (k + 1) * ny_out + j     # (j,   k+1)
                f.write(f"4 {p0} {p1} {p2} {p3}\n")

        f.write(f"\nCELL_TYPES {ncells}\n")
        for _ in range(ncells):
            f.write("9\n")  # VTK_QUAD = 9

        # --- POINT_DATA ---
        f.write(f"\nPOINT_DATA {npoints}\n")

        # IsSolid
        f.write("SCALARS IsSolid int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(BUFFER, NZ6 - BUFFER):
            for j in range(BUFFER, NY6 - BUFFER):
                hill_h = HillFunction(y_global[j])
                is_solid = 1 if z_global[j, k] < hill_h else 0
                f.write(f"{is_solid}\n")

        # Z 座標
        f.write("\nSCALARS Z_coordinate double 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(BUFFER, NZ6 - BUFFER):
            for j in range(BUFFER, NY6 - BUFFER):
                f.write(f"{z_global[j, k]}\n")

        # 局部 dz
        f.write("\nSCALARS LocalDZ double 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(BUFFER, NZ6 - BUFFER):
            for j in range(BUFFER, NY6 - BUFFER):
                if k < NZ6 - BUFFER - 1:
                    dz = z_global[j, k + 1] - z_global[j, k]
                else:
                    dz = z_global[j, k] - z_global[j, k - 1]
                f.write(f"{dz}\n")

        # 長寬比
        f.write("\nSCALARS AspectRatio double 1\n")
        f.write("LOOKUP_TABLE default\n")
        for k in range(BUFFER, NZ6 - BUFFER):
            for j in range(BUFFER, NY6 - BUFFER):
                if k < NZ6 - BUFFER - 1:
                    dz = z_global[j, k + 1] - z_global[j, k]
                else:
                    dz = z_global[j, k] - z_global[j, k - 1]
                if dz > 1e-15:
                    ar = max(dy / dz, dz / dy)
                else:
                    ar = 1e10
                f.write(f"{ar}\n")

    print(f"已輸出 UNSTRUCTURED_GRID VTK: {filename}")


# =============================================================================
# 網格資訊統計
# =============================================================================
def print_grid_statistics(y_global, z_global):
    """輸出網格統計資訊"""
    dy = y_global[BUFFER + 1] - y_global[BUFFER]
    print("\n" + "=" * 60)
    print("非均勻網格統計資訊")
    print("=" * 60)
    print(f"計算域: Y ∈ [0, {LY}],  Z ∈ [0, {LZ}]")
    print(f"網格數: NY={NY} (含buffer: {NY6}), NZ={NZ} (含buffer: {NZ6})")
    print(f"Buffer 層數: {BUFFER}")
    print(f"有效計算網格: {NY6-6} × {NZ6-6} = {(NY6-6)*(NZ6-6)} 個節點")
    print(f"\nY 方向 (均勻): dy = {dy:.8f}")
    print(f"minSize (lattice size) = {minSize:.8f}")
    print(f"CFL = {CFL}")

    # Z 方向統計
    dz_all = []
    for j in range(BUFFER, NY6 - BUFFER):
        for k in range(BUFFER, NZ6 - BUFFER - 1):
            dz = z_global[j, k + 1] - z_global[j, k]
            dz_all.append(dz)
    dz_arr = np.array(dz_all)
    print(f"\nZ 方向 (非均勻):")
    print(f"  最小 dz = {dz_arr.min():.8f}")
    print(f"  最大 dz = {dz_arr.max():.8f}")
    print(f"  平均 dz = {dz_arr.mean():.8f}")
    print(f"  拉伸比  = {dz_arr.max() / dz_arr.min():.2f}")

    # 長寬比統計
    ar_all = []
    for dz in dz_all:
        if dz > 1e-15:
            ar_all.append(max(dy / dz, dz / dy))
    ar_arr = np.array(ar_all)
    print(f"\n最長寬比 (Aspect Ratio):")
    print(f"  最小 AR = {ar_arr.min():.2f}")
    print(f"  最大 AR = {ar_arr.max():.2f}")
    print(f"  平均 AR = {ar_arr.mean():.2f}")
    print("=" * 60)


# =============================================================================
# 主程式
# =============================================================================
if __name__ == "__main__":
    os.makedirs("output", exist_ok=True)

    print("正在生成非均勻網格...")
    y_global, z_global, nonuni_a = generate_mesh()

    print_grid_statistics(y_global, z_global)

    print("\n正在輸出 VTK 檔案...")
    output_structured_grid_vtk(y_global, z_global, "output/grid_structure.vtk")
    output_unstructured_grid_vtk(y_global, z_global, "output/grid_quadcells.vtk")

    print("\n完成！請使用 ParaView 開啟 VTK 檔案：")
    print("  1. File → Open → output/grid_structure.vtk")
    print("  2. Apply")
    print("  3. Representation → 'Surface With Edges' 或 'Wireframe'")
    print("  4. Color By → CellSize_Z / AspectRatio 可觀察網格品質")
    print("\n或開啟 grid_quadcells.vtk（UNSTRUCTURED_GRID 格式）：")
    print("  - Color By → IsSolid 可分辨流體/固體區域")
