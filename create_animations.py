#!/usr/bin/env python3
"""
批次處理VTK檔案，生成向量場和等值線圖片並製作動畫
Batch process VTK files to generate vector field and contour images and create animations
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image
import re

# 嘗試載入 scipy（流線圖需要將非均勻網格插值到均勻網格）
try:
    from scipy.interpolate import griddata
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Warning: scipy not available, streamlines will use raw grid (may be less smooth)")

# 設置全局字體為Times New Roman加粗
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['figure.titleweight'] = 'bold'

def check_packages():
    """檢查必要的套件"""
    required = {
        'numpy': 'numpy',
        'matplotlib': 'matplotlib',
    }
    
    missing = []
    for module, package in required.items():
        try:
            __import__(module)
        except ImportError:
            missing.append(package)
    
    if missing:
        print("=" * 60)
        print(f"[ERROR] 缺少必要套件: {', '.join(missing)}")
        print("=" * 60)
        print("\n請使用以下指令安裝:")
        print(f"  pip install {' '.join(missing)}")
        print("=" * 60)
        sys.exit(1)

check_packages()


def read_vtk_structured_grid(filename):
    """讀取 VTK Structured Grid 格式的檔案"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 解析維度
    dims = None
    points = []
    velocity = []
    density = []
    velocity_magnitude = []
    
    reading_points = False
    reading_velocity = False
    reading_density = False
    reading_vel_mag = False
    
    for i, line in enumerate(lines):
        if 'DIMENSIONS' in line:
            dims = list(map(int, line.split()[1:4]))
        elif 'POINTS' in line and 'double' in line:
            reading_points = True
            reading_velocity = False
            reading_density = False
            reading_vel_mag = False
            continue
        elif 'POINT_DATA' in line:
            reading_points = False
            continue
        elif 'VECTORS Velocity' in line:
            reading_velocity = True
            reading_density = False
            reading_vel_mag = False
            reading_points = False
            continue
        elif 'SCALARS Density' in line:
            reading_velocity = False
            reading_density = True
            reading_vel_mag = False
            reading_points = False
            continue
        elif 'SCALARS VelocityMagnitude' in line:
            reading_velocity = False
            reading_density = False
            reading_vel_mag = True
            reading_points = False
            continue
        elif 'LOOKUP_TABLE' in line:
            continue
        elif reading_points:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    points.append([float(parts[0]), float(parts[1]), float(parts[2])])
                except ValueError:
                    pass
        elif reading_velocity:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    velocity.append([float(parts[0]), float(parts[1]), float(parts[2])])
                except ValueError:
                    reading_velocity = False
        elif reading_density:
            parts = line.strip().split()
            for val in parts:
                try:
                    density.append(float(val))
                except ValueError:
                    pass
        elif reading_vel_mag:
            parts = line.strip().split()
            for val in parts:
                try:
                    velocity_magnitude.append(float(val))
                except ValueError:
                    pass
    
    points = np.array(points)
    velocity = np.array(velocity)
    density = np.array(density)
    velocity_magnitude = np.array(velocity_magnitude)
    
    if dims is None:
        raise ValueError("無法從VTK檔案中讀取DIMENSIONS")
    
    nx, ny, nz = dims
    
    # Reshape數據 (對於2D模擬，nz通常是1)
    if len(points) == nx * ny * nz:
        points = points.reshape((ny, nx, 3))
    
    if len(velocity) == nx * ny * nz:
        velocity = velocity.reshape((ny, nx, 3))
    
    if len(density) == nx * ny * nz:
        density = density.reshape((ny, nx))
    
    if len(velocity_magnitude) == nx * ny * nz:
        velocity_magnitude = velocity_magnitude.reshape((ny, nx))
    
    return {
        'dims': dims,
        'points': points,
        'velocity': velocity,
        'density': density,
        'velocity_magnitude': velocity_magnitude
    }


def plot_2d_data(data):
    """提取2D數據用於繪圖"""
    velocity = data['velocity']
    density = data['density']
    velocity_magnitude = data['velocity_magnitude']
    points = data['points']
    dims = data['dims']
    nx, ny, nz = dims
    
    # 對於2D模擬，數據已經是2D的
    u = velocity[:, :, 0]  # x方向速度
    v = velocity[:, :, 1]  # y方向速度
    rho = density[:, :]     # 密度
    
    # 如果有預先計算的速度大小，使用它；否則計算
    if velocity_magnitude.size > 0:
        speed = velocity_magnitude[:, :]
    else:
        speed = np.sqrt(u**2 + v**2)
    
    # 計算馬赫數 (Ma = u/cs, cs = 1/sqrt(3) 在格子單位中)
    cs = 1.0 / np.sqrt(3.0)  # 格子聲速
    mach = speed / cs
    
    # 從points提取x, y座標
    x_coords = points[:, :, 0]
    y_coords = points[:, :, 1]
    
    return x_coords, y_coords, u, v, speed, rho, mach


def create_vector_field_image(X, Y, u, v, speed, mach, filename, title="Velocity Field"):
    """創建向量場圖片"""
    fig, ax = plt.subplots(figsize=(14, 7))
    
    # 繪製速度大小的contour作為背景
    contourf = ax.contourf(X, Y, speed, levels=30, cmap='jet', alpha=0.7)
    cbar = plt.colorbar(contourf, ax=ax, label='Speed')
    cbar.ax.tick_params(labelsize=12)
    
    # 繪製向量場（降低密度以便觀察）- 使用黑色箭頭
    skip_x = max(1, X.shape[1] // 40)  # 調整箭頭密度
    skip_y = max(1, X.shape[0] // 20)
    
    quiver = ax.quiver(X[::skip_y, ::skip_x], Y[::skip_y, ::skip_x], 
                       u[::skip_y, ::skip_x], v[::skip_y, ::skip_x],
                       scale=None, scale_units='xy',
                       color='black', alpha=0.8, width=0.002)
    
    # 計算並顯示馬赫數統計信息
    ma_max = np.max(mach)
    ma_avg = np.mean(mach)
    ma_limit = 0.3  # 不可壓縮流的馬赫數限制
    
    # 在圖片底部添加統計信息（放大字體，使用加粗）- 增加與主體的距離
    info_text = f'Ma max: {ma_max:.4f}   |   Ma avg: {ma_avg:.4f}   |   Ma limit: {ma_limit:.2f}'
    ax.text(0.5, -0.18, info_text, transform=ax.transAxes,
            fontsize=14, weight='bold', verticalalignment='top', horizontalalignment='center',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9, pad=0.8))
    
    ax.set_xlabel('X', fontsize=13, weight='bold')
    ax.set_ylabel('Y', fontsize=13, weight='bold')
    ax.set_title(title, fontsize=15, weight='bold', pad=10)
    ax.set_aspect('equal')
    ax.tick_params(labelsize=11)
    
    # 調整子圖邊距 - 增加主體大小，減少留白
    plt.subplots_adjust(left=0.07, right=0.96, top=0.95, bottom=0.15)
    plt.savefig(filename, dpi=150, bbox_inches='tight', pad_inches=0.25)
    plt.close()
    print(f"已保存向量場圖片: {filename} (Ma_max={ma_max:.4f})")


def create_contour_image(X, Y, field, mach, filename, title="Contour Plot", levels=30, cmap='jet', field_name='Value'):
    """創建等值線圖片"""
    fig, ax = plt.subplots(figsize=(14, 7))
    
    # 繪製filled contour
    contourf = ax.contourf(X, Y, field, levels=levels, cmap=cmap)
    cbar = plt.colorbar(contourf, ax=ax, label=field_name)
    cbar.ax.tick_params(labelsize=12)
    
    # 添加contour線
    contour = ax.contour(X, Y, field, levels=levels, colors='black', 
                         linewidths=0.5, alpha=0.3)
    
    # 計算並顯示馬赫數統計信息
    ma_max = np.max(mach)
    ma_avg = np.mean(mach)
    ma_limit = 0.3  # 不可壓縮流的馬赫數限制
    
    # 在圖片底部添加統計信息（放大字體，使用加粗）- 增加與主體的距離
    info_text = f'Ma max: {ma_max:.4f}   |   Ma avg: {ma_avg:.4f}   |   Ma limit: {ma_limit:.2f}'
    ax.text(0.5, -0.18, info_text, transform=ax.transAxes,
            fontsize=14, weight='bold', verticalalignment='top', horizontalalignment='center',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9, pad=0.8))
    
    ax.set_xlabel('X', fontsize=13, weight='bold')
    ax.set_ylabel('Y', fontsize=13, weight='bold')
    ax.set_title(title, fontsize=15, weight='bold', pad=10)
    ax.set_aspect('equal')
    ax.tick_params(labelsize=11)
    
    # 調整子圖邊距 - 增加主體大小，減少留白
    plt.subplots_adjust(left=0.07, right=0.96, top=0.95, bottom=0.15)
    plt.savefig(filename, dpi=150, bbox_inches='tight', pad_inches=0.25)
    plt.close()
    print(f"已保存等值線圖片: {filename}")


def create_streamline_image(X, Y, u, v, speed, mach, filename, title="Streamlines"):
    """創建流線圖（streamline plot）"""
    fig, ax = plt.subplots(figsize=(14, 7))

    # 背景：速度大小等值線
    contourf = ax.contourf(X, Y, speed, levels=30, cmap='coolwarm', alpha=0.8)
    cbar = plt.colorbar(contourf, ax=ax, label='Speed')
    cbar.ax.tick_params(labelsize=12)

    # 流線圖需要均勻網格；若有 scipy 則插值到均勻網格，否則直接用原始資料
    if SCIPY_AVAILABLE:
        # 建立均勻網格
        n_y = min(200, X.shape[1])
        n_z = min(100, X.shape[0])
        x_reg = np.linspace(X.min(), X.max(), n_y)
        y_reg = np.linspace(Y.min(), Y.max(), n_z)
        X_reg, Y_reg = np.meshgrid(x_reg, y_reg)

        u_reg = griddata((X.flatten(), Y.flatten()), u.flatten(),
                         (X_reg, Y_reg), method='linear')
        v_reg = griddata((X.flatten(), Y.flatten()), v.flatten(),
                         (X_reg, Y_reg), method='linear')

        # 處理 NaN（邊界外的點填 0）
        u_reg = np.nan_to_num(u_reg, nan=0.0)
        v_reg = np.nan_to_num(v_reg, nan=0.0)

        ax.streamplot(x_reg, y_reg, u_reg, v_reg,
                      color='black', linewidth=0.8, density=2.0,
                      arrowsize=1.0, arrowstyle='->')
    else:
        # 無 scipy 時直接使用原始網格（streamplot 要求嚴格遞增的 1D 座標）
        try:
            x_1d = X[0, :]
            y_1d = Y[:, 0]
            ax.streamplot(x_1d, y_1d, u, v,
                          color='black', linewidth=0.8, density=2.0,
                          arrowsize=1.0, arrowstyle='->')
        except Exception as e:
            print(f"  流線繪製失敗 (需要 scipy 處理非均勻網格): {e}")

    # 馬赫數統計
    ma_max = np.max(mach)
    ma_avg = np.mean(mach)
    ma_limit = 0.3

    info_text = f'Ma max: {ma_max:.4f}   |   Ma avg: {ma_avg:.4f}   |   Ma limit: {ma_limit:.2f}'
    ax.text(0.5, -0.18, info_text, transform=ax.transAxes,
            fontsize=14, weight='bold', verticalalignment='top', horizontalalignment='center',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9, pad=0.8))

    ax.set_xlabel('X', fontsize=13, weight='bold')
    ax.set_ylabel('Y', fontsize=13, weight='bold')
    ax.set_title(title, fontsize=15, weight='bold', pad=10)
    ax.set_aspect('equal')
    ax.tick_params(labelsize=11)

    plt.subplots_adjust(left=0.07, right=0.96, top=0.95, bottom=0.15)
    plt.savefig(filename, dpi=150, bbox_inches='tight', pad_inches=0.25)
    plt.close()
    print(f"已保存流線圖: {filename} (Ma_max={ma_max:.4f})")


def batch_process_vtk_files(vtk_dir='output', output_dir='output_images'):
    """批次處理所有VTK檔案"""
    # 創建輸出目錄
    vector_dir = os.path.join(output_dir, 'vectors')
    contour_speed_dir = os.path.join(output_dir, 'contour_speed')
    contour_density_dir = os.path.join(output_dir, 'contour_density')
    streamline_dir = os.path.join(output_dir, 'streamlines')

    os.makedirs(vector_dir, exist_ok=True)
    os.makedirs(contour_speed_dir, exist_ok=True)
    os.makedirs(contour_density_dir, exist_ok=True)
    os.makedirs(streamline_dir, exist_ok=True)
    
    # 獲取所有VTK檔案並排序
    vtk_files = sorted(glob.glob(os.path.join(vtk_dir, 'flow_*.vtk')))
    
    if not vtk_files:
        print(f"在 {vtk_dir} 目錄中沒有找到VTK檔案")
        return []
    
    print(f"找到 {len(vtk_files)} 個VTK檔案")
    
    processed_files = []
    
    for i, vtk_file in enumerate(vtk_files):
        basename = os.path.splitext(os.path.basename(vtk_file))[0]
        print(f"\n處理 [{i+1}/{len(vtk_files)}]: {basename}")
        
        try:
            # 讀取VTK數據
            data = read_vtk_structured_grid(vtk_file)
            
            # 提取2D數據
            X, Y, u, v, speed, rho, mach = plot_2d_data(data)
            
            # 生成向量場圖片
            vector_filename = os.path.join(vector_dir, f'{basename}_vector.png')
            create_vector_field_image(X, Y, u, v, speed, mach, vector_filename, 
                                     title=f'Velocity Field - {basename}')
            
            # 生成速度等值線圖片
            contour_speed_filename = os.path.join(contour_speed_dir, f'{basename}_speed.png')
            create_contour_image(X, Y, speed, mach, contour_speed_filename,
                                title=f'Speed Contour - {basename}', cmap='jet', field_name='Speed')
            
            # 生成密度等值線圖片
            contour_density_filename = os.path.join(contour_density_dir, f'{basename}_density.png')
            create_contour_image(X, Y, rho, mach, contour_density_filename,
                                title=f'Density Contour - {basename}', cmap='viridis', field_name='Density')

            # 生成流線圖
            streamline_filename = os.path.join(streamline_dir, f'{basename}_streamline.png')
            create_streamline_image(X, Y, u, v, speed, mach, streamline_filename,
                                    title=f'Streamlines - {basename}')

            processed_files.append({
                'vtk': vtk_file,
                'vector': vector_filename,
                'speed': contour_speed_filename,
                'density': contour_density_filename,
                'streamline': streamline_filename
            })
            
        except Exception as e:
            print(f"處理 {vtk_file} 時發生錯誤: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    return processed_files


def create_animation_from_images(image_files, output_file, fps=10, title_prefix=""):
    """從圖片序列創建動畫 - 直接使用圖片不添加額外留白"""
    if not image_files:
        print("沒有圖片可以製作動畫")
        return
    
    print(f"\n正在創建動畫: {output_file}")
    
    # 使用PIL直接讀取圖片並創建GIF
    images = []
    for img_file in image_files:
        img = Image.open(img_file)
        images.append(img)
    
    # 保存為GIF動畫
    images[0].save(
        output_file,
        save_all=True,
        append_images=images[1:],
        duration=int(1000/fps),  # 每幀持續時間（毫秒）
        loop=0  # 無限循環
    )
    
    print(f"動畫已保存: {output_file}")


def main():
    """主程序"""
    print("=" * 60)
    print("VTK檔案批次處理 - 生成向量場和等值線動畫")
    print("=" * 60)
    
    # 設置參數
    vtk_dir = 'output'
    output_dir = 'output_images'
    fps = 10     # 動畫幀率
    
    # 批次處理VTK檔案
    print(f"\n開始批次處理VTK檔案...")
    processed_files = batch_process_vtk_files(vtk_dir, output_dir)
    
    if not processed_files:
        print("沒有成功處理任何檔案")
        return
    
    print(f"\n成功處理 {len(processed_files)} 個檔案")
    
    # 創建動畫
    print("\n" + "=" * 60)
    print("創建動畫...")
    print("=" * 60)
    
    # 向量場動畫
    vector_images = [f['vector'] for f in processed_files]
    vector_animation = os.path.join(output_dir, 'animation_vector_field.gif')
    create_animation_from_images(vector_images, vector_animation, fps=fps, 
                                 title_prefix="Velocity Field")
    
    # 速度等值線動畫
    speed_images = [f['speed'] for f in processed_files]
    speed_animation = os.path.join(output_dir, 'animation_speed_contour.gif')
    create_animation_from_images(speed_images, speed_animation, fps=fps,
                                title_prefix="Speed Contour")
    
    # 密度等值線動畫
    density_images = [f['density'] for f in processed_files]
    density_animation = os.path.join(output_dir, 'animation_density_contour.gif')
    create_animation_from_images(density_images, density_animation, fps=fps,
                                title_prefix="Density Contour")

    # 流線圖動畫
    streamline_images = [f['streamline'] for f in processed_files]
    streamline_animation = os.path.join(output_dir, 'animation_streamlines.gif')
    create_animation_from_images(streamline_images, streamline_animation, fps=fps,
                                title_prefix="Streamlines")

    print("\n" + "=" * 60)
    print("完成！")
    print("=" * 60)
    print(f"\n生成的圖片位於: {output_dir}")
    print(f"  - 向量場圖片: {os.path.join(output_dir, 'vectors')}")
    print(f"  - 速度等值線: {os.path.join(output_dir, 'contour_speed')}")
    print(f"  - 密度等值線: {os.path.join(output_dir, 'contour_density')}")
    print(f"  - 流線圖:     {os.path.join(output_dir, 'streamlines')}")
    print(f"\n生成的動畫:")
    print(f"  - {vector_animation}")
    print(f"  - {speed_animation}")
    print(f"  - {density_animation}")
    print(f"  - {streamline_animation}")
    print("=" * 60)


if __name__ == '__main__':
    main()
