#!/usr/bin/env python3
"""Analyze velocity field from VTK files to check for oscillations"""

import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

output_dir = r"c:\Users\88697.CHENPENGCHUNG12\D2Q9_PeriodicHill\output"
fig_dir = os.path.join(output_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

# Grid dimensions from VTK
NY = 201  # Y points
NZ = 150  # Z points

# Physical dimensions
Ly = 9.0   # Channel length in Y
Lz = 3.035 # Channel height in Z

def read_vtk_velocity(filepath):
    """Extract velocity field from VTK file"""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find velocity section
    vel_start = None
    for i, line in enumerate(lines):
        if "VECTORS Velocity" in line:
            vel_start = i + 1
            break
    
    if vel_start is None:
        return None
    
    # Read velocity data (NY * NZ points)
    n_points = NY * NZ
    velocities = []
    for i in range(vel_start, vel_start + n_points):
        parts = lines[i].strip().split()
        if len(parts) == 3:
            uy, uz, ux = float(parts[0]), float(parts[1]), float(parts[2])
            velocities.append([uy, uz, ux])
    
    return np.array(velocities).reshape(NZ, NY, 3)

def read_vtk_coordinates(filepath):
    """Extract coordinates from VTK file"""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
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
    coords = []
    for i in range(pts_start, pts_start + n_points):
        parts = lines[i].strip().split()
        if len(parts) == 3:
            y, z, x = float(parts[0]), float(parts[1]), float(parts[2])
            coords.append([y, z])
    
    coords = np.array(coords).reshape(NZ, NY, 2)
    Y = coords[:, :, 0]
    Z = coords[:, :, 1]
    return Y, Z

def plot_velocity_field(vel, Y, Z, timestep, fig_dir):
    """Generate comprehensive velocity field visualizations"""
    
    Uy = vel[:, :, 0]  # Streamwise velocity
    Uz = vel[:, :, 1]  # Vertical velocity
    Umag = np.sqrt(Uy**2 + Uz**2)
    
    # ================================================================
    # Figure 1: Full velocity field contours
    # ================================================================
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))
    
    # Uy contour
    ax = axes[0]
    norm_uy = TwoSlopeNorm(vmin=Uy.min(), vcenter=0, vmax=Uy.max())
    c1 = ax.contourf(Y, Z, Uy, levels=50, cmap='RdBu_r', norm=norm_uy)
    ax.set_xlabel('Y (streamwise)')
    ax.set_ylabel('Z (vertical)')
    ax.set_title(f'Streamwise Velocity Uy (t={timestep})')
    ax.set_aspect('equal')
    plt.colorbar(c1, ax=ax, label='Uy')
    
    # Uz contour
    ax = axes[1]
    norm_uz = TwoSlopeNorm(vmin=Uz.min(), vcenter=0, vmax=Uz.max())
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
    
    # Vertical centerline at different Y positions
    y_indices = [10, 50, 100, 150, 190]
    y_labels = ['Y=0.45 (Hill)', 'Y=2.25', 'Y=4.50 (Center)', 'Y=6.75', 'Y=8.55 (Hill)']
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
    k_indices = [10, 30, 50, 75, 100, 130]
    k_labels = ['k=10 (Near wall)', 'k=30', 'k=50', 'k=75 (Mid)', 'k=100', 'k=130 (Upper)']
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
    
    # Quiver plot
    ax.quiver(Y[::skip_z, ::skip_y], Z[::skip_z, ::skip_y], 
              Uy[::skip_z, ::skip_y], Uz[::skip_z, ::skip_y],
              scale=2, width=0.002, color='black', alpha=0.7)
    
    # Try to add streamlines (may fail if field is too complex)
    try:
        # Create regular grid for streamplot
        y_reg = np.linspace(Y.min(), Y.max(), 100)
        z_reg = np.linspace(Z.min(), Z.max(), 50)
        Y_reg, Z_reg = np.meshgrid(y_reg, z_reg)
        
        from scipy.interpolate import griddata
        Uy_reg = griddata((Y.flatten(), Z.flatten()), Uy.flatten(), (Y_reg, Z_reg), method='linear')
        Uz_reg = griddata((Y.flatten(), Z.flatten()), Uz.flatten(), (Y_reg, Z_reg), method='linear')
        
        ax.streamplot(y_reg, z_reg, Uy_reg, Uz_reg, color='white', linewidth=0.5, density=1.5, arrowsize=0.5)
    except:
        pass  # Skip streamlines if scipy not available or interpolation fails
    
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
    
    # Y-direction differences at mid-height
    k_mid = 75
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
    
    # Z-direction differences at mid-Y
    j_mid = 100
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
    # Get list of VTK files sorted by timestep
    vtk_files = sorted(glob.glob(os.path.join(output_dir, "flow_*.vtk")))
    
    if len(vtk_files) < 3:
        print("Not enough VTK files for analysis")
        return
    
    print("=" * 60)
    print("Velocity Field Oscillation Analysis")
    print("=" * 60)
    print(f"Found {len(vtk_files)} VTK files")
    print(f"Grid: {NY} x {NZ} (Y x Z)")
    
    # Analyze last few files
    files_to_analyze = vtk_files[-6:]
    
    print(f"\nAnalyzing files: {[os.path.basename(f) for f in files_to_analyze]}")
    
    # Store velocity statistics
    stats = []
    probe_velocities = []
    
    # Read coordinates from first file
    Y_coords, Z_coords = read_vtk_coordinates(files_to_analyze[0])
    
    # Probe points: (k, j) - mid-channel, near hill, in wake
    probe_points = [
        (75, 100, "Mid-channel center"),
        (50, 30, "Near left hill"),
        (100, 150, "Wake region"),
        (130, 100, "Upper channel"),
    ]
    
    for filepath in files_to_analyze:
        basename = os.path.basename(filepath)
        timestep = int(re.search(r'\d+', basename).group())
        
        vel = read_vtk_velocity(filepath)
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
            print("\n⚠️  WARNING: Possible oscillations detected!")
        else:
            print("\n✓ Velocity field appears stable (low variance)")
    
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
    last_vel = read_vtk_velocity(files_to_analyze[-1])
    last_timestep = int(re.search(r'\d+', os.path.basename(files_to_analyze[-1])).group())
    
    if last_vel is not None and Y_coords is not None:
        # Generate visualization figures
        print("\nGenerating velocity field visualizations...")
        plot_velocity_field(last_vel, Y_coords, Z_coords, last_timestep, fig_dir)
        
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
            print("\n⚠️  WARNING: High frequency spatial oscillations detected!")
        else:
            print("\n✓ No obvious striped pattern detected")
    
    # ============================================================
    # Centerline velocity profiles
    # ============================================================
    print("\n" + "=" * 60)
    print("Centerline Velocity Profiles (Z-direction at different Y)")
    print("=" * 60)
    
    if last_vel is not None and Y_coords is not None:
        # Y positions for centerline profiles (in lattice units)
        # Periodic hill: Y from 0 to ~9, hill peaks at Y~0.5 and Y~8.5
        y_positions = [
            (10, "Y=0.45 (Left hill peak)"),
            (50, "Y=2.25 (After hill)"),
            (100, "Y=4.50 (Channel center)"),
            (150, "Y=6.75 (Before right hill)"),
            (190, "Y=8.55 (Right hill peak)"),
        ]
        
        for j, y_label in y_positions:
            print(f"\n--- {y_label} (j={j}) ---")
            print(f"{'k':>4} {'Z':>8} {'Uy':>12} {'Uz':>12}")
            print("-" * 40)
            
            # Print velocity profile at this Y position (every 10th point)
            for k in range(5, NZ-5, 10):
                # Estimate physical Z coordinate (tanh stretching)
                z_ratio = k / NZ
                z_phys = z_ratio * 3.035  # Approximate channel height
                
                uy = last_vel[k, j, 0]
                uz = last_vel[k, j, 1]
                print(f"{k:>4} {z_phys:>8.3f} {uy:>12.6f} {uz:>12.6f}")
        
        # Horizontal centerline at mid-height
        print("\n" + "=" * 60)
        print("Horizontal Centerline (Y-direction at mid-height k=75)")
        print("=" * 60)
        
        k_mid = 75
        print(f"{'j':>4} {'Y':>8} {'Uy':>12} {'Uz':>12}")
        print("-" * 40)
        
        for j in range(0, NY, 10):
            y_phys = j * 0.045  # Grid spacing in Y
            uy = last_vel[k_mid, j, 0]
            uz = last_vel[k_mid, j, 1]
            print(f"{j:>4} {y_phys:>8.3f} {uy:>12.6f} {uz:>12.6f}")
        
        # Check for velocity reversal (recirculation)
        print("\n" + "=" * 60)
        print("Recirculation Zone Analysis")
        print("=" * 60)
        
        # Check for negative Uy (backflow) at different heights
        for k in [10, 20, 30, 50]:
            uy_profile = last_vel[k, :, 0]
            negative_count = np.sum(uy_profile < 0)
            if negative_count > 0:
                negative_indices = np.where(uy_profile < 0)[0]
                y_start = negative_indices[0] * 0.045
                y_end = negative_indices[-1] * 0.045
                print(f"k={k:>3}: Backflow region Y=[{y_start:.2f}, {y_end:.2f}], {negative_count} points")
            else:
                print(f"k={k:>3}: No backflow detected")
    
    print("\n" + "=" * 60)
    print("Analysis Complete")
    print("=" * 60)

if __name__ == "__main__":
    main()
