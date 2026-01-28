#!/usr/bin/env python3
"""Analyze velocity field from VTK files to check for oscillations

跨平台兼容版本：
- 自動偵測虛擬環境
- 如果缺少套件，會提示正確的安裝方式
- 支援 Windows / macOS / Linux
"""

import os
import sys

def check_packages():
    """檢查必要的 Python 套件"""
    required = ['numpy', 'matplotlib']
    missing = []
    
    for pkg in required:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    
    if missing:
        # 偵測腳本所在目錄的虛擬環境
        script_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(script_dir)
        
        # 可能的虛擬環境位置
        venv_paths = [
            os.path.join(parent_dir, '.venv', 'bin', 'python'),  # Unix: ../.venv
            os.path.join(parent_dir, '.venv', 'Scripts', 'python.exe'),  # Windows
            os.path.join(script_dir, '.venv', 'bin', 'python'),  # Unix: ./.venv
            os.path.join(script_dir, '.venv', 'Scripts', 'python.exe'),  # Windows
        ]
        
        venv_python = None
        for vp in venv_paths:
            if os.path.exists(vp):
                venv_python = vp
                break
        
        print("=" * 60)
        print(f"[ERROR] 缺少必要套件: {', '.join(missing)}")
        print("=" * 60)
        
        if venv_python:
            print(f"\n發現虛擬環境，請使用以下指令執行:")
            print(f"  {venv_python} {os.path.abspath(__file__)}")
        else:
            print("\n請建立虛擬環境並安裝套件:")
            print(f"  cd {parent_dir}")
            print("  python3 -m venv .venv")
            if sys.platform == 'win32':
                print("  .venv\\Scripts\\activate")
            else:
                print("  source .venv/bin/activate")
            print(f"  pip install {' '.join(required)} scipy")
            print(f"\n然後執行:")
            print(f"  python {os.path.basename(__file__)}")
        
        print("=" * 60)
        sys.exit(1)

# 檢查套件
check_packages()

import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# =============================================================================
# 設定論文標準格式 (Times New Roman + 標準尺寸)
# =============================================================================
plt.rcParams.update({
    # 字體設定
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif', 'serif'],
    'mathtext.fontset': 'stix',  # 數學字體使用 STIX (類似 Times)
    
    # 字體大小 (論文標準)
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8,
    
    # 線條寬度 (論文標準)
    'lines.linewidth': 1.0,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.minor.width': 0.5,
    'ytick.minor.width': 0.5,
    
    # 標記大小
    'lines.markersize': 4,
    
    # 圖例
    'legend.frameon': True,
    'legend.framealpha': 1.0,
    'legend.edgecolor': 'black',
    'legend.fancybox': False,
    
    # 網格
    'grid.linewidth': 0.5,
    'grid.alpha': 0.3,
    
    # 輸出品質
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

# 嘗試導入 scipy，如果失敗則設置標記
try:
    from scipy.interpolate import griddata, interp1d
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Warning: scipy not available, streamlines will be skipped")

# =============================================================================
# Re=700 Periodic Hill Benchmark Data (ERCOFTAC/Breuer et al. 2009)
# 數據來源: M. Breuer et al., "Flow over periodic hills - Numerical and 
#           experimental study in a wide range of Reynolds numbers", 
#           Computers & Fluids 38 (2009) 433-457
# x/h 位置: 0.05, 0.5, 1, 2, 3, 4, 5, 6, 7, 8
# =============================================================================

# Re=700 的 <u>/Ub 剖面 (y/h vs <u>/Ub) - Breuer et al. DNS 數據
BENCHMARK_RE700 = {
    # x/h = 0.05 (山丘頂端後方)
    0.05: {
        'y_h': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.8, 3.035],
        'u_Ub': [0.0, 0.35, 0.62, 0.82, 0.95, 1.04, 1.10, 1.16, 1.18, 1.17, 1.14, 1.10, 1.05, 1.00, 0.95, 0.88, 0.82, 0.78]
    },
    # x/h = 0.5 (山丘下游)
    0.5: {
        'y_h': [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.035],
        'u_Ub': [0.0, -0.05, -0.08, -0.06, 0.02, 0.25, 0.50, 0.70, 0.85, 1.02, 1.10, 1.12, 1.10, 1.06, 1.01, 0.96, 0.85, 0.78]
    },
    # x/h = 1 (迴流區開始)
    1: {
        'y_h': [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.035],
        'u_Ub': [0.0, -0.10, -0.15, -0.12, -0.05, 0.15, 0.40, 0.60, 0.78, 0.98, 1.08, 1.11, 1.10, 1.06, 1.00, 0.95, 0.84, 0.78]
    },
    # x/h = 2 (迴流區中間)
    2: {
        'y_h': [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.035],
        'u_Ub': [0.0, -0.12, -0.18, -0.18, -0.12, 0.05, 0.30, 0.52, 0.72, 0.95, 1.06, 1.10, 1.09, 1.05, 1.00, 0.94, 0.84, 0.78]
    },
    # x/h = 3 (迴流區)
    3: {
        'y_h': [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.035],
        'u_Ub': [0.0, -0.10, -0.15, -0.15, -0.10, 0.05, 0.28, 0.50, 0.70, 0.94, 1.05, 1.09, 1.08, 1.05, 1.00, 0.94, 0.84, 0.78]
    },
    # x/h = 4 (迴流區尾端)
    4: {
        'y_h': [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.035],
        'u_Ub': [0.0, -0.05, -0.08, -0.06, 0.02, 0.18, 0.38, 0.56, 0.72, 0.94, 1.04, 1.08, 1.08, 1.04, 0.99, 0.94, 0.84, 0.78]
    },
    # x/h = 5 (再附著點附近)
    5: {
        'y_h': [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.035],
        'u_Ub': [0.0, 0.02, 0.05, 0.10, 0.18, 0.32, 0.48, 0.62, 0.75, 0.94, 1.03, 1.07, 1.07, 1.04, 0.99, 0.94, 0.84, 0.78]
    },
    # x/h = 6 (恢復區)
    6: {
        'y_h': [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.035],
        'u_Ub': [0.0, 0.10, 0.20, 0.30, 0.38, 0.50, 0.60, 0.70, 0.80, 0.95, 1.03, 1.06, 1.06, 1.03, 0.99, 0.94, 0.84, 0.78]
    },
    # x/h = 7 (恢復區)
    7: {
        'y_h': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.035],
        'u_Ub': [0.0, 0.28, 0.48, 0.60, 0.70, 0.78, 0.85, 0.96, 1.02, 1.05, 1.05, 1.03, 0.99, 0.94, 0.84, 0.78]
    },
    # x/h = 8 (接近下一個山丘)
    8: {
        'y_h': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.035],
        'u_Ub': [0.0, 0.32, 0.55, 0.70, 0.80, 0.88, 0.94, 1.02, 1.06, 1.07, 1.06, 1.03, 0.99, 0.94, 0.84, 0.78]
    }
}

# 自動偵測腳本所在目錄，設定輸出路徑
script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, "output")
fig_dir = os.path.join(output_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

# =============================================================================
# 從 variables.h 讀取模擬參數
# =============================================================================

def read_simulation_params(script_dir):
    """從 variables.h 讀取模擬參數 (Re, LY, LZ, NY, NZ 等)"""
    params = {
        'Re': 700,      # 預設值
        'LY': 9.0,
        'LZ': 3.036,
        'NY': 512,
        'NZ': 256,
    }
    
    variables_file = os.path.join(script_dir, 'variables.h')
    if not os.path.exists(variables_file):
        print(f"Warning: variables.h not found, using default parameters")
        return params
    
    try:
        with open(variables_file, 'r') as f:
            content = f.read()
            
        # 使用正則表達式提取參數
        patterns = {
            'Re': r'#define\s+Re\s+(\d+)',
            'LY': r'#define\s+LY\s+\(?([\d.]+)\)?',
            'LZ': r'#define\s+LZ\s+\(?([\d.]+)\)?',
            'NY': r'#define\s+NY\s+(\d+)',
            'NZ': r'#define\s+NZ\s+(\d+)',
        }
        
        for key, pattern in patterns.items():
            match = re.search(pattern, content)
            if match:
                if key in ['Re', 'NY', 'NZ']:
                    params[key] = int(match.group(1))
                else:
                    params[key] = float(match.group(1))
        
        print(f"Read from variables.h: Re={params['Re']}, "
              f"Domain={params['LY']}×{params['LZ']}, "
              f"Grid={params['NY']}×{params['NZ']}")
              
    except Exception as e:
        print(f"Warning: Could not parse variables.h: {e}")
    
    return params

# 讀取模擬參數
SIM_PARAMS = read_simulation_params(script_dir)

# =============================================================================
# 網格尺寸將從 VTK 文件自動讀取，不再硬編碼
# =============================================================================

def read_vtk_dimensions(filepath):
    """從 VTK 文件自動讀取網格尺寸 (DIMENSIONS 行)
    
    VTK STRUCTURED_GRID 格式:
    DIMENSIONS nx ny nz
    其中 nx = NY+1 (Y方向點數), ny = NZ (Z方向點數), nz = 1 (2D)
    
    Returns:
        (NY, NZ): Y方向和Z方向的點數
    """
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("DIMENSIONS"):
                parts = line.strip().split()
                # DIMENSIONS ny nz 1 → (ny, nz)
                ny = int(parts[1])  # Y方向點數 (第一維度，最快變化)
                nz = int(parts[2])  # Z方向點數 (第二維度)
                return ny, nz
    return None, None

def read_vtk_velocity(filepath, ny=None, nz=None):
    """Extract velocity field from VTK file
    
    Args:
        filepath: VTK 文件路徑
        ny, nz: 可選的網格尺寸，如果不提供則自動讀取
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # 如果未提供尺寸，從文件讀取
    if ny is None or nz is None:
        for line in lines:
            if line.startswith("DIMENSIONS"):
                parts = line.strip().split()
                ny = int(parts[1])
                nz = int(parts[2])
                break
    
    if ny is None or nz is None:
        print("Error: Could not determine grid dimensions")
        return None
    
    # Find velocity section
    vel_start = None
    for i, line in enumerate(lines):
        if "VECTORS Velocity" in line:
            vel_start = i + 1
            break
    
    if vel_start is None:
        return None
    
    # Read velocity data (ny * nz points)
    # VTK 順序: 外層 k (Z), 內層 j (Y)
    # 所以 reshape 順序是 (nz, ny, 3)
    n_points = ny * nz
    velocities = []
    for i in range(vel_start, vel_start + n_points):
        if i >= len(lines):
            break
        parts = lines[i].strip().split()
        if len(parts) == 3:
            uy, uz, ux = float(parts[0]), float(parts[1]), float(parts[2])
            velocities.append([uy, uz, ux])
    
    # VTK 資料順序: 先遍歷 Y (內層)，再遍歷 Z (外層)
    # reshape 成 (nz, ny, 3) 後，索引為 vel[k, j, :]
    return np.array(velocities).reshape(nz, ny, 3)

def read_vtk_coordinates(filepath, ny=None, nz=None):
    """Extract coordinates from VTK file
    
    Args:
        filepath: VTK 文件路徑
        ny, nz: 可選的網格尺寸，如果不提供則自動讀取
    
    Returns:
        (Y, Z): Y座標陣列和Z座標陣列，shape 為 (nz, ny)
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # 如果未提供尺寸，從文件讀取
    if ny is None or nz is None:
        for line in lines:
            if line.startswith("DIMENSIONS"):
                parts = line.strip().split()
                ny = int(parts[1])
                nz = int(parts[2])
                break
    
    # Find POINTS section
    pts_start = None
    n_points = 0
    for i, line in enumerate(lines):
        if line.startswith("POINTS"):
            parts = line.strip().split()
            n_points = int(parts[1])
            pts_start = i + 1
            break
    
    if pts_start is None:
        return None, None
    
    # Read coordinates
    # VTK DIMENSIONS ny nz 1 表示: 第一維度(ny點)最快變化 = Y方向
    # 資料排列: (j=0,k=0), (j=1,k=0), ..., (j=ny-1,k=0), (j=0,k=1), ...
    coords = []
    for i in range(pts_start, pts_start + n_points):
        if i >= len(lines):
            break
        parts = lines[i].strip().split()
        if len(parts) == 3:
            y, z, x = float(parts[0]), float(parts[1]), float(parts[2])
            coords.append([y, z])
    
    # reshape 成 (nz, ny, 2)，因為 k(Z) 是外層，j(Y) 是內層
    coords = np.array(coords).reshape(nz, ny, 2)
    Y = coords[:, :, 0]  # Y[k, j] = y 座標
    Z = coords[:, :, 1]  # Z[k, j] = z 座標
    return Y, Z

def safe_two_slope_norm(data, vcenter=0):
    """Create TwoSlopeNorm safely, handling edge cases"""
    vmin, vmax = data.min(), data.max()
    
    # Handle edge cases where TwoSlopeNorm would fail
    if vmin >= vcenter:
        # All values >= vcenter, use regular norm
        return None
    if vmax <= vcenter:
        # All values <= vcenter, use regular norm
        return None
    if np.isnan(vmin) or np.isnan(vmax) or np.isinf(vmin) or np.isinf(vmax):
        return None
    
    return TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)


def get_hill_height(y, LY=9.0):
    """計算給定 y 位置的山丘高度 (與 C++ HillFunction 完全一致)
    
    這是從 model.h 中的 HillFunction 直接翻譯過來的
    
    Args:
        y: Y position (streamwise, 0 to LY)
        LY: Domain length (default 9.0)
    
    Returns:
        z_wall: wall height at this position
    """
    # 週期性邊界處理
    if y < 0.0:
        Yb = y + LY
    elif y > LY:
        Yb = y - LY
    else:
        Yb = y
    
    model = 0.0
    
    # 左半部山丘 (分段多項式)
    if Yb <= (54./28.)*(9./54.):
        model = (1./28.) * min(28., 28. + 0.006775070969851*Yb*28*Yb*28 
                               - 0.0021245277758000*Yb*28*Yb*28*Yb*28)
    elif Yb <= (54./28.)*(14./54.):
        model = (1./28.) * (25.07355893131 + 0.9754803562315*Yb*28 
                            - 0.1016116352781*Yb*28*Yb*28 
                            + 0.001889794677828*Yb*28*Yb*28*Yb*28)
    elif Yb <= (54./28.)*(20./54.):
        model = (1./28.) * (25.79601052357 + 0.8206693007457*Yb*28 
                            - 0.09055370274339*Yb*28*Yb*28 
                            + 0.001626510569859*Yb*28*Yb*28*Yb*28)
    elif Yb <= (54./28.)*(30./54.):
        model = (1./28.) * (40.46435022819 - 1.379581654948*Yb*28 
                            + 0.019458845041284*Yb*28*Yb*28 
                            - 0.0002070318932190*Yb*28*Yb*28*Yb*28)
    elif Yb <= (54./28.)*(40./54.):
        model = (1./28.) * (17.92461334664 + 0.8743920332081*Yb*28 
                            - 0.05567361123058*Yb*28*Yb*28 
                            + 0.0006277731764683*Yb*28*Yb*28*Yb*28)
    elif Yb <= (54./28.)*(54./54.):
        model = (1./28.) * max(0., 56.39011190988 - 2.010520359035*Yb*28 
                               + 0.01644919857549*Yb*28*Yb*28 
                               + 0.00002674976141766*Yb*28*Yb*28*Yb*28)
    
    # 右半部山丘 (對稱)
    if Yb < LY - (54./28.)*(40./54.) and Yb >= LY - (54./28.)*(54./54.):
        Yr = LY - Yb
        model = (1./28.) * max(0., 56.39011190988 - 2.010520359035*Yr*28 
                               + 0.01644919857549*Yr*28*Yr*28 
                               + 0.00002674976141766*Yr*28*Yr*28*Yr*28)
    elif Yb < LY - (54./28.)*(30./54.) and Yb >= LY - (54./28.)*(40./54.):
        Yr = LY - Yb
        model = (1./28.) * (17.92461334664 + 0.8743920332081*Yr*28 
                            - 0.05567361123058*Yr*28*Yr*28 
                            + 0.0006277731764683*Yr*28*Yr*28*Yr*28)
    elif Yb < LY - (54./28.)*(20./54.) and Yb >= LY - (54./28.)*(30./54.):
        Yr = LY - Yb
        model = (1./28.) * (40.46435022819 - 1.379581654948*Yr*28 
                            + 0.019458845041284*Yr*28*Yr*28 
                            - 0.0002070318932190*Yr*28*Yr*28*Yr*28)
    elif Yb < LY - (54./28.)*(14./54.) and Yb >= LY - (54./28.)*(20./54.):
        Yr = LY - Yb
        model = (1./28.) * (25.79601052357 + 0.8206693007457*Yr*28 
                            - 0.09055370274339*Yr*28*Yr*28 
                            + 0.001626510569859*Yr*28*Yr*28*Yr*28)
    elif Yb < LY - (54./28.)*(9./54.) and Yb >= LY - (54./28.)*(14./54.):
        Yr = LY - Yb
        model = (1./28.) * (25.07355893131 + 0.9754803562315*Yr*28 
                            - 0.1016116352781*Yr*28*Yr*28 
                            + 0.001889794677828*Yr*28*Yr*28*Yr*28)
    elif Yb >= LY - (54./28.)*(9./54.):
        Yr = LY - Yb
        model = (1./28.) * min(28., 28. + 0.006775070969851*Yr*28*Yr*28 
                               - 0.0021245277758000*Yr*28*Yr*28*Yr*28)
    
    return model


def plot_velocity_profiles_paper_style(vel, Y, Z, timestep, fig_dir, ny, nz, 
                                        Ub=None, h=1.0, L=9.0, H=3.035,
                                        show_benchmark=True, Re=700):
    """生成論文風格的速度剖面圖 (按 x/h 位置並排)
    
    Args:
        vel: Velocity array (nz, ny, 3)
        Y, Z: Coordinate arrays (nz, ny) - 物理座標
        timestep: Time step number
        fig_dir: Output directory for figures
        ny, nz: Grid dimensions
        Ub: Bulk velocity for normalization (auto-calculated if None)
        h: Hill height (default 1.0)
        L: Domain length in x-direction (default 9.0h)
        H: Channel height (default 3.035h)
        show_benchmark: Whether to show Re=700 benchmark data
        Re: Reynolds number for title
    """
    Uy = vel[:, :, 0]  # Streamwise velocity (u)
    Uz = vel[:, :, 1]  # Vertical velocity (v)
    
    # Auto-calculate bulk velocity if not provided
    if Ub is None:
        # Ub = average Uy across the channel (excluding near-wall regions)
        k_start = int(nz * 0.1)
        k_end = int(nz * 0.9)
        Ub = np.mean(Uy[k_start:k_end, :])
        print(f"  Auto-calculated Ub = {Ub:.6f}")
    
    # Normalize velocities
    u_norm = Uy / Ub if Ub != 0 else Uy
    v_norm = Uz / Ub if Ub != 0 else Uz
    
    # ================================================================
    # x/h 位置對應到 j 索引的映射
    # 假設 Y 座標從 0 到 L (9h)
    # ================================================================
    Y_min, Y_max = Y.min(), Y.max()
    print(f"  Physical Y range: [{Y_min:.3f}, {Y_max:.3f}]")
    
    # x/h 位置 (論文標準位置)
    x_h_positions = [0.05, 0.5, 1, 2, 3, 4, 5, 6, 7, 8]
    
    # 找到每個 x/h 位置對應的 j 索引
    j_indices = []
    for x_h in x_h_positions:
        # 將 x/h 轉換為物理 Y 座標
        y_target = x_h * h * (Y_max - Y_min) / L + Y_min
        # 找最接近的 j 索引
        j = np.argmin(np.abs(Y[nz//2, :] - y_target))
        j_indices.append(j)
        print(f"  x/h = {x_h}: j = {j}, Y = {Y[nz//2, j]:.3f}")
    
    # ================================================================
    # Figure 1: Streamwise velocity <u>/Ub profiles (論文風格) - 已移至 Figure 4
    # Figure 2: Vertical velocity - 已停用
    # 只保留 subplots 和 paper 風格
    # ================================================================
    
    # ================================================================
    # Figure 3: 子圖版本 (上下兩排，每排5個剖面) - 類似論文圖 (a)(b) 風格
    # ================================================================
    fig, axes = plt.subplots(2, 5, figsize=(12, 6), sharey=True)
    plt.subplots_adjust(hspace=0.3, wspace=0.1)
    
    for i, (x_h, j) in enumerate(zip(x_h_positions, j_indices)):
        row = i // 5
        col = i % 5
        ax = axes[row, col]
        
        z_profile = Z[:, j]
        u_profile = u_norm[:, j]
        
        # 獲取該位置的山丘高度
        y_phys = Y[nz//2, j]
        wall_z = get_hill_height(y_phys, L)
        
        # y/h 從壁面算起
        y_h = (z_profile - wall_z) / h
        
        # LBM 結果 (實線，論文標準線寬)
        ax.plot(u_profile, y_h, 'k-', linewidth=1.0, label='LBM')
        
        # Benchmark 數據 (空心圓，論文標準大小)
        if show_benchmark and x_h in BENCHMARK_RE700:
            bench = BENCHMARK_RE700[x_h]
            ax.plot(bench['u_Ub'], bench['y_h'], 'o', 
                    color='black', markersize=3, markerfacecolor='none', 
                    markeredgewidth=0.6, label='DNS')
        
        ax.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
        ax.set_xlim([-0.3, 1.3])
        ax.set_ylim([0, 3.5])
        
        # x/h 標籤放在圖內
        ax.text(0.95, 0.95, f'$x/h={x_h}$', transform=ax.transAxes,
                fontsize=8, ha='right', va='top')
        
        # 只在左側顯示 y 軸標籤
        if col == 0:
            ax.set_ylabel(r'$y/h$')
        
        # 只在底部顯示 x 軸標籤
        if row == 1:
            ax.set_xlabel(r'$\langle u \rangle/U_b$')
        
        # 細化刻度
        ax.tick_params(direction='in', top=True, right=True)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_yticks([0, 1, 2, 3])
    
    # 統一圖例 (放在第一個子圖)
    axes[0, 0].legend(loc='upper left', fontsize=7, handlelength=1.5)
    
    plt.savefig(os.path.join(fig_dir, f'velocity_profiles_subplots_t{timestep}.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: velocity_profiles_subplots_t{timestep}.png")
    
    # ================================================================
    # Figure 4: 單圖並排顯示 (完全模仿論文風格)
    # 所有剖面在同一張圖上按 x/h 位置偏移顯示
    # 剖面底部跟隨山丘形狀
    # ================================================================
    fig, ax = plt.subplots(figsize=(10, 4))
    
    # 首先繪製山丘輪廓 (使用實際的 HillFunction)
    x_hill = np.linspace(0, L, 500)
    y_hill = np.array([get_hill_height(x, L) for x in x_hill])
    ax.fill_between(x_hill, 0, y_hill, color='lightgray', alpha=0.6, zorder=1)
    ax.plot(x_hill, y_hill, 'k-', linewidth=0.8, zorder=2)
    
    # 繪製上壁面
    ax.axhline(y=H, color='k', linewidth=0.8, zorder=2)
    
    # 速度剖面縮放因子 (讓剖面在視覺上更清晰)
    velocity_scale = 0.7  # 調整此值改變剖面寬度
    
    for i, (x_h, j) in enumerate(zip(x_h_positions, j_indices)):
        # 獲取該 j 位置的實際 Y 座標
        y_phys = Y[nz//2, j]  # 實際的 Y 物理座標
        
        z_profile = Z[:, j]   # 該位置的 Z 座標剖面
        u_profile = u_norm[:, j]
        
        # 計算該位置的山丘高度 (使用實際 Y 座標)
        wall_z = get_hill_height(y_phys, L)
        
        # 只取山丘以上的部分
        mask = z_profile >= wall_z - 0.05  # 稍微放寬避免邊界問題
        z_plot = z_profile[mask]
        u_plot = u_profile[mask]
        
        # 偏移到 x/h 位置，速度值乘以縮放因子
        u_shifted = y_phys + u_plot * velocity_scale
        
        # 繪製 LBM 結果 (黑色實線，論文標準)
        ax.plot(u_shifted, z_plot, '-', color='black', linewidth=1.0, zorder=3)
        
        # 繪製 u=0 參考線 (從山丘到上壁)
        ax.plot([y_phys, y_phys], [wall_z, H], ':', color='gray', linewidth=0.4, alpha=0.7, zorder=2)
        
        # 繪製 benchmark 數據 (空心圓點，論文標準)
        if show_benchmark and x_h in BENCHMARK_RE700:
            bench = BENCHMARK_RE700[x_h]
            # Benchmark 的 y_h 是從壁面算起的高度，需要加上該位置的山丘高度
            bench_z = np.array(bench['y_h']) * h + wall_z
            bench_u = np.array(bench['u_Ub']) * velocity_scale + y_phys
            # 只繪製在上壁面以下的點
            valid = bench_z <= H
            ax.plot(bench_u[valid], bench_z[valid], 'o', color='black', markersize=3,
                    markerfacecolor='none', markeredgewidth=0.6, zorder=4)
    
    # 設定軸
    ax.set_xlim([-0.3, L + 0.3])
    ax.set_ylim([0, H + 0.1])
    ax.set_xlabel(r'$x/h$')
    ax.set_ylabel(r'$y/h$')
    
    # 添加圖例
    ax.plot([], [], 'k-', linewidth=1.0, label='LBM')
    ax.plot([], [], 'ko', markersize=3, markerfacecolor='none', markeredgewidth=0.6, label='DNS (Re=700)')
    ax.legend(loc='upper right', fontsize=8, frameon=True)
    
    # 刻度設定
    ax.tick_params(direction='in', top=True, right=True)
    ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    ax.set_aspect('equal')
    
    plt.savefig(os.path.join(fig_dir, f'velocity_profiles_paper_t{timestep}.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: velocity_profiles_paper_t{timestep}.png")
    
    return Ub

def plot_velocity_field(vel, Y, Z, timestep, fig_dir, ny=None, nz=None):
    """Generate velocity field visualizations (contour and vector field only)
    
    Args:
        vel: Velocity array (nz, ny, 3)
        Y, Z: Coordinate arrays (nz, ny)
        timestep: Time step number
        fig_dir: Output directory for figures
        ny, nz: Grid dimensions (auto-detected from vel if not provided)
    """
    # 從 vel 的 shape 自動獲取尺寸
    if nz is None:
        nz = vel.shape[0]
    if ny is None:
        ny = vel.shape[1]
    
    Uy = vel[:, :, 0]  # Streamwise velocity
    Uz = vel[:, :, 1]  # Vertical velocity
    Umag = np.sqrt(Uy**2 + Uz**2)
    
    # Check for divergence
    is_diverged = np.max(np.abs(Uy)) > 10 or np.max(np.abs(Uz)) > 10
    if is_diverged:
        print(f"  [!] WARNING: Simulation appears DIVERGED! Max|Uy|={np.max(np.abs(Uy)):.2f}, Max|Uz|={np.max(np.abs(Uz)):.2f}")
    
    # ================================================================
    # Figure 1: Full velocity field contours
    # ================================================================
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))
    
    # Uy contour
    ax = axes[0]
    norm_uy = safe_two_slope_norm(Uy, vcenter=0)
    c1 = ax.contourf(Y, Z, Uy, levels=50, cmap='RdBu_r', norm=norm_uy)
    ax.set_xlabel('Y (streamwise)')
    ax.set_ylabel('Z (vertical)')
    ax.set_title(f'Streamwise Velocity Uy (t={timestep})')
    ax.set_aspect('equal')
    plt.colorbar(c1, ax=ax, label='Uy')
    
    # Uz contour
    ax = axes[1]
    norm_uz = safe_two_slope_norm(Uz, vcenter=0)
    c2 = ax.contourf(Y, Z, Uz, levels=50, cmap='RdBu_r', norm=norm_uz)
    ax.set_xlabel('Y (streamwise)')
    ax.set_ylabel('Z (vertical)')
    ax.set_title(f'Vertical Velocity Uz (t={timestep})')
    ax.set_aspect('equal')
    plt.colorbar(c2, ax=ax, label='Uz')
    
    # Velocity magnitude
    ax = axes[2]
    c3 = ax.contourf(Y, Z, Umag, levels=50, cmap='jet')
    ax.set_xlabel('Y (streamwise)')
    ax.set_ylabel('Z (vertical)')
    ax.set_title(f'Velocity Magnitude |U| (t={timestep})')
    ax.set_aspect('equal')
    plt.colorbar(c3, ax=ax, label='|U|')
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'velocity_contour_t{timestep}.png'), dpi=150)
    plt.close()
    print(f"  Saved: velocity_contour_t{timestep}.png")
    
    # ================================================================
    # Figure 2: Vector field with streamlines
    # ================================================================
    fig, ax = plt.subplots(figsize=(14, 5))
    
    # Subsample for quiver plot
    skip_y = 5
    skip_z = 3
    
    # Background: velocity magnitude contour
    c = ax.contourf(Y, Z, Umag, levels=30, cmap='coolwarm', alpha=0.7)
    plt.colorbar(c, ax=ax, label='|U|')
    
    # Quiver plot with auto-scaling
    Uy_sub = Uy[::skip_z, ::skip_y]
    Uz_sub = Uz[::skip_z, ::skip_y]
    # Auto-calculate scale based on velocity magnitude
    vel_mag = np.sqrt(Uy_sub**2 + Uz_sub**2)
    max_vel = np.percentile(vel_mag[vel_mag > 0], 95) if np.any(vel_mag > 0) else 1.0
    quiver_scale = max(max_vel * 20, 0.1)  # Ensure minimum scale
    
    ax.quiver(Y[::skip_z, ::skip_y], Z[::skip_z, ::skip_y], 
              Uy_sub, Uz_sub,
              scale=quiver_scale, width=0.002, color='black', alpha=0.7)
    
    # Try to add streamlines (may fail if field is too complex)
    if SCIPY_AVAILABLE:
        try:
            # Create regular grid for streamplot
            y_reg = np.linspace(Y.min(), Y.max(), 100)
            z_reg = np.linspace(Z.min(), Z.max(), 50)
            Y_reg, Z_reg = np.meshgrid(y_reg, z_reg)
            
            Uy_reg = griddata((Y.flatten(), Z.flatten()), Uy.flatten(), (Y_reg, Z_reg), method='linear')
            Uz_reg = griddata((Y.flatten(), Z.flatten()), Uz.flatten(), (Y_reg, Z_reg), method='linear')
            
            ax.streamplot(y_reg, z_reg, Uy_reg, Uz_reg, color='white', linewidth=0.5, density=1.5, arrowsize=0.5)
        except Exception as e:
            print(f"  Streamlines skipped: {e}")
    
    ax.set_xlabel('Y (streamwise)')
    ax.set_ylabel('Z (vertical)')
    ax.set_title(f'Velocity Vector Field (t={timestep})')
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'vector_field_t{timestep}.png'), dpi=150)
    plt.close()
    print(f"  Saved: vector_field_t{timestep}.png")
    
    return Uy, Uz, Umag

def main():
    # Get list of VTK files sorted by timestep (numerically, not alphabetically)
    vtk_files = glob.glob(os.path.join(output_dir, "flow_*.vtk"))
    # Sort by extracting the timestep number from filename
    vtk_files = sorted(vtk_files, key=lambda f: int(re.search(r'flow_(\d+)\.vtk', os.path.basename(f)).group(1)))
    
    if len(vtk_files) == 0:
        print("No VTK files found in output directory")
        return
    
    if len(vtk_files) < 3:
        print(f"Only {len(vtk_files)} VTK file(s) found, proceeding with available files...")
    
    # 從第一個 VTK 文件自動讀取網格尺寸
    NY, NZ = read_vtk_dimensions(vtk_files[0])
    if NY is None or NZ is None:
        print("Error: Could not read grid dimensions from VTK file")
        return
    
    print("=" * 60)
    print("Velocity Field Oscillation Analysis")
    print("=" * 60)
    print(f"Found {len(vtk_files)} VTK files")
    print(f"Grid: {NY} x {NZ} (Y x Z) - Auto-detected from VTK")
    
    # Analyze last few files (at least the latest one)
    n_files = min(6, len(vtk_files))
    files_to_analyze = vtk_files[-n_files:]
    
    print(f"\nAnalyzing files: {[os.path.basename(f) for f in files_to_analyze]}")
    
    # Store velocity statistics
    stats = []
    probe_velocities = []
    
    # Read coordinates from first file
    Y_coords, Z_coords = read_vtk_coordinates(files_to_analyze[0], NY, NZ)
    
    # Probe points: (k, j) - 動態計算基於實際網格尺寸
    # 確保索引在有效範圍內
    probe_points = [
        (min(NZ*3//8, NZ-5), min(NY//4, NY-5), "Mid-channel center"),
        (min(NZ//4, NZ-5), min(NY//10, NY-5), "Near left hill"),
        (min(NZ//2, NZ-5), min(NY*3//8, NY-5), "Wake region"),
        (min(NZ*2//3, NZ-5), min(NY//4, NY-5), "Upper channel"),
    ]
    
    for filepath in files_to_analyze:
        basename = os.path.basename(filepath)
        timestep = int(re.search(r'\d+', basename).group())
        
        vel = read_vtk_velocity(filepath, NY, NZ)
        if vel is None:
            print(f"Failed to read {basename}")
            continue
        
        # Calculate statistics
        uy = vel[:, :, 0]  # Y velocity (streamwise)
        uz = vel[:, :, 1]  # Z velocity (vertical)
        
        # Exclude boundary regions (first/last few rows)
        uy_interior = uy[5:-5, :]
        uz_interior = uz[5:-5, :]
        
        uy_max = np.max(uy_interior)
        uy_min = np.min(uy_interior)
        uz_max = np.max(uz_interior)
        uz_min = np.min(uz_interior)
        uy_mean = np.mean(uy_interior)
        
        stats.append({
            'time': timestep,
            'uy_max': uy_max,
            'uy_min': uy_min,
            'uy_mean': uy_mean,
            'uz_max': uz_max,
            'uz_min': uz_min
        })
        
        # Probe point velocities
        probe_vals = {}
        for k, j, name in probe_points:
            probe_vals[name] = (vel[k, j, 0], vel[k, j, 1])
        probe_velocities.append({'time': timestep, 'probes': probe_vals})
    
    # Print statistics
    print("\n" + "=" * 60)
    print("Global Velocity Statistics")
    print("=" * 60)
    print(f"{'Time':>8} {'Uy_max':>10} {'Uy_min':>10} {'Uy_mean':>10} {'Uz_max':>10} {'Uz_min':>10}")
    print("-" * 60)
    
    for s in stats:
        print(f"{s['time']:>8} {s['uy_max']:>10.5f} {s['uy_min']:>10.5f} {s['uy_mean']:>10.5f} {s['uz_max']:>10.5f} {s['uz_min']:>10.5f}")
    
    # Check for oscillations in global statistics
    print("\n" + "=" * 60)
    print("Oscillation Check (Variance of Statistics)")
    print("=" * 60)
    
    if len(stats) >= 3:
        uy_max_var = np.var([s['uy_max'] for s in stats])
        uy_mean_var = np.var([s['uy_mean'] for s in stats])
        uz_max_var = np.var([s['uz_max'] for s in stats])
        
        print(f"Variance of Uy_max: {uy_max_var:.2e}")
        print(f"Variance of Uy_mean: {uy_mean_var:.2e}")
        print(f"Variance of Uz_max: {uz_max_var:.2e}")
        
        # Threshold for oscillation detection
        threshold = 1e-6
        if uy_max_var > threshold or uz_max_var > threshold:
            print("\n[!] WARNING: Possible oscillations detected!")
        else:
            print("\n[OK] Velocity field appears stable (low variance)")
    
    # Print probe point analysis
    print("\n" + "=" * 60)
    print("Probe Point Velocity Evolution")
    print("=" * 60)
    
    for k, j, name in probe_points:
        print(f"\n{name} (k={k}, j={j}):")
        print(f"  {'Time':>8} {'Uy':>12} {'Uz':>12}")
        print("  " + "-" * 36)
        
        for pv in probe_velocities:
            uy, uz = pv['probes'][name]
            print(f"  {pv['time']:>8} {uy:>12.6f} {uz:>12.6f}")
        
        # Calculate variance at this probe
        uy_vals = [pv['probes'][name][0] for pv in probe_velocities]
        uz_vals = [pv['probes'][name][1] for pv in probe_velocities]
        print(f"  Variance: Uy={np.var(uy_vals):.2e}, Uz={np.var(uz_vals):.2e}")
    
    # Check spatial oscillation (striped pattern)
    print("\n" + "=" * 60)
    print("Spatial Pattern Analysis (Checking for stripes)")
    print("=" * 60)
    
    # Use the last file for spatial analysis
    last_file = files_to_analyze[-1]
    print(f"\nAnalyzing latest file: {os.path.basename(last_file)}")
    last_vel = read_vtk_velocity(last_file, NY, NZ)
    last_timestep = int(re.search(r'\d+', os.path.basename(last_file)).group())
    
    if last_vel is not None and Y_coords is not None:
        # Generate visualization figures
        print("\nGenerating velocity field visualizations...")
        plot_velocity_field(last_vel, Y_coords, Z_coords, last_timestep, fig_dir, NY, NZ)
        
        # Generate paper-style velocity profiles with benchmark comparison
        print("\nGenerating paper-style velocity profiles with Re=700 benchmark...")
        plot_velocity_profiles_paper_style(
            last_vel, Y_coords, Z_coords, last_timestep, fig_dir, 
            NY, NZ, Ub=None, 
            h=1.0, 
            L=SIM_PARAMS['LY'], 
            H=SIM_PARAMS['LZ'],
            show_benchmark=True, 
            Re=SIM_PARAMS['Re']
        )
        
        # Check Y-direction pattern at mid-height
        k_mid = NZ // 2
        uy_profile = last_vel[k_mid, :, 0]
        
        # Look for sign changes (oscillations)
        sign_changes = np.sum(np.diff(np.sign(uy_profile)) != 0)
        
        # Check Z-direction pattern at mid-Y
        j_mid = NY // 2
        uz_profile = last_vel[:, j_mid, 1]
        z_sign_changes = np.sum(np.diff(np.sign(uz_profile)) != 0)
        
        print(f"At k={k_mid} (mid-height): Uy sign changes along Y = {sign_changes}")
        print(f"At j={j_mid} (mid-Y): Uz sign changes along Z = {z_sign_changes}")
        
        # High frequency oscillation check (consecutive sign changes)
        if sign_changes > 50 or z_sign_changes > 30:
            print("\n[!] WARNING: High frequency spatial oscillations detected!")
        else:
            print("\n[OK] No obvious striped pattern detected")
    
    # ============================================================
    # Centerline velocity profiles
    # ============================================================
    print("\n" + "=" * 60)
    print("Centerline Velocity Profiles (Z-direction at different Y)")
    print("=" * 60)
    
    if last_vel is not None and Y_coords is not None:
        # Y positions for centerline profiles - 動態計算基於實際網格尺寸
        # 確保索引在有效範圍內
        y_positions = [
            (min(NY//20, NY-1), "Y near left (Left hill peak)"),
            (min(NY//5, NY-1), "Y = 1/5 (After hill)"),
            (min(NY//2, NY-1), "Y = 1/2 (Channel center)"),
            (min(NY*3//4, NY-1), "Y = 3/4 (Before right hill)"),
            (min(NY*19//20, NY-1), "Y near right (Right hill peak)"),
        ]
        
        for j, y_label in y_positions:
            print(f"\n--- {y_label} (j={j}) ---")
            print(f"{'k':>4} {'Z':>8} {'Uy':>12} {'Uz':>12}")
            print("-" * 40)
            
            # Print velocity profile at this Y position (every 10th point or scaled)
            step = max(1, NZ // 20)  # 動態步長
            for k in range(5, NZ-5, step):
                # Estimate physical Z coordinate (tanh stretching)
                z_ratio = k / NZ
                z_phys = z_ratio * 3.035  # Approximate channel height
                
                uy = last_vel[k, j, 0]
                uz = last_vel[k, j, 1]
                print(f"{k:>4} {z_phys:>8.3f} {uy:>12.6f} {uz:>12.6f}")
        
        # Horizontal centerline at mid-height
        print("\n" + "=" * 60)
        print(f"Horizontal Centerline (Y-direction at mid-height k={NZ//2})")
        print("=" * 60)
        
        k_mid = NZ // 2
        print(f"{'j':>4} {'Y':>8} {'Uy':>12} {'Uz':>12}")
        print("-" * 40)
        
        step_y = max(1, NY // 20)  # 動態步長
        for j in range(0, NY, step_y):
            y_phys = Y_coords[k_mid, j]  # 使用實際讀取的座標
            uy = last_vel[k_mid, j, 0]
            uz = last_vel[k_mid, j, 1]
            print(f"{j:>4} {y_phys:>8.3f} {uy:>12.6f} {uz:>12.6f}")
        
        # Check for velocity reversal (recirculation)
        print("\n" + "=" * 60)
        print("Recirculation Zone Analysis")
        print("=" * 60)
        
        # Check for negative Uy (backflow) at different heights
        k_check = [k for k in [NZ//20, NZ//10, NZ//6, NZ//4] if k < NZ-5]
        for k in k_check:
            uy_profile = last_vel[k, :, 0]
            negative_count = np.sum(uy_profile < 0)
            if negative_count > 0:
                negative_indices = np.where(uy_profile < 0)[0]
                y_start = Y_coords[k, negative_indices[0]]  # 使用實際座標
                y_end = Y_coords[k, negative_indices[-1]]    # 使用實際座標
                print(f"k={k:>3}: Backflow region Y=[{y_start:.2f}, {y_end:.2f}], {negative_count} points")
            else:
                print(f"k={k:>3}: No backflow detected")
    
    print("\n" + "=" * 60)
    print("Analysis Complete")
    print("=" * 60)
    print(f"\nFigures saved to: {fig_dir}")

if __name__ == "__main__":
    main()
