#!/usr/bin/env python3
"""Analyze velocity field from VTK files to check for oscillations"""

import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# 嘗試導入 scipy，如果失敗則設置標記
try:
    from scipy.interpolate import griddata
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Warning: scipy not available, streamlines will be skipped")

output_dir = r"c:\Users\88697.CHENPENGCHUNG12\D2Q9_PeriodicHill\output"
fig_dir = os.path.join(output_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

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

def plot_velocity_field(vel, Y, Z, timestep, fig_dir, ny=None, nz=None):
    """Generate comprehensive velocity field visualizations
    
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
    
    # ================================================================
    # Figure 2: Centerline velocity profiles
    # ================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Vertical centerline at different Y positions - 動態計算索引
    y_indices = [max(0, ny//20), ny//5, ny//2, ny*3//4, min(ny-1, ny*19//20)]
    y_labels = [f'j={j}' for j in y_indices]
    colors = plt.cm.viridis(np.linspace(0, 1, len(y_indices)))
    
    # Uy profile along Z
    ax = axes[0, 0]
    for j, label, color in zip(y_indices, y_labels, colors):
        z_profile = Z[:, j]
        uy_profile = Uy[:, j]
        ax.plot(uy_profile, z_profile, label=label, color=color, linewidth=1.5)
    ax.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Uy (streamwise velocity)')
    ax.set_ylabel('Z (height)')
    ax.set_title('Vertical Profiles of Uy')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Uz profile along Z
    ax = axes[0, 1]
    for j, label, color in zip(y_indices, y_labels, colors):
        z_profile = Z[:, j]
        uz_profile = Uz[:, j]
        ax.plot(uz_profile, z_profile, label=label, color=color, linewidth=1.5)
    ax.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Uz (vertical velocity)')
    ax.set_ylabel('Z (height)')
    ax.set_title('Vertical Profiles of Uz')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Horizontal centerline at different Z positions
    # Ensure indices are within bounds
    k_max = nz - 10
    k_indices = [k for k in [10, 30, 50, 75, 100, 130] if k < k_max]
    if not k_indices:  # 如果網格太小
        k_indices = [nz // 4, nz // 2, 3 * nz // 4]
    k_labels = [f'k={k}' for k in k_indices]
    colors = plt.cm.plasma(np.linspace(0, 1, len(k_indices)))
    
    # Uy profile along Y
    ax = axes[1, 0]
    for k, label, color in zip(k_indices, k_labels, colors):
        y_profile = Y[k, :]
        uy_profile = Uy[k, :]
        ax.plot(y_profile, uy_profile, label=label, color=color, linewidth=1.5)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Y (streamwise position)')
    ax.set_ylabel('Uy (streamwise velocity)')
    ax.set_title('Horizontal Profiles of Uy')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Uz profile along Y
    ax = axes[1, 1]
    for k, label, color in zip(k_indices, k_labels, colors):
        y_profile = Y[k, :]
        uz_profile = Uz[k, :]
        ax.plot(y_profile, uz_profile, label=label, color=color, linewidth=1.5)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Y (streamwise position)')
    ax.set_ylabel('Uz (vertical velocity)')
    ax.set_title('Horizontal Profiles of Uz')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'Centerline Velocity Profiles (t={timestep})', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'centerline_profiles_t{timestep}.png'), dpi=150)
    plt.close()
    
    # ================================================================
    # Figure 3: Vector field with streamlines
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
    
    # ================================================================
    # Figure 4: Oscillation check - consecutive point differences
    # ================================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Y-direction differences at mid-height (使用動態 nz)
    k_mid = nz // 2
    ax = axes[0, 0]
    uy_diff = np.diff(Uy[k_mid, :])
    ax.plot(Y[k_mid, :-1], uy_diff, 'b-', linewidth=1)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Y')
    ax.set_ylabel('ΔUy')
    ax.set_title(f'Uy Differences Along Y (k={k_mid})')
    ax.grid(True, alpha=0.3)
    
    ax = axes[0, 1]
    uz_diff = np.diff(Uz[k_mid, :])
    ax.plot(Y[k_mid, :-1], uz_diff, 'r-', linewidth=1)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Y')
    ax.set_ylabel('ΔUz')
    ax.set_title(f'Uz Differences Along Y (k={k_mid})')
    ax.grid(True, alpha=0.3)
    
    # Z-direction differences at mid-Y (使用動態 ny)
    j_mid = ny // 2
    ax = axes[1, 0]
    uy_diff_z = np.diff(Uy[:, j_mid])
    ax.plot(Z[:-1, j_mid], uy_diff_z, 'b-', linewidth=1)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Z')
    ax.set_ylabel('ΔUy')
    ax.set_title(f'Uy Differences Along Z (j={j_mid})')
    ax.grid(True, alpha=0.3)
    
    ax = axes[1, 1]
    uz_diff_z = np.diff(Uz[:, j_mid])
    ax.plot(Z[:-1, j_mid], uz_diff_z, 'r-', linewidth=1)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Z')
    ax.set_ylabel('ΔUz')
    ax.set_title(f'Uz Differences Along Z (j={j_mid})')
    ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'Spatial Oscillation Check (t={timestep})', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'oscillation_check_t{timestep}.png'), dpi=150)
    plt.close()
    
    print(f"  Saved figures to {fig_dir}/")
    
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
