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
import re

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


def read_vtk_structured_points(filename):
    """讀取 VTK Structured Points 格式的檔案"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 解析維度
    dims = None
    spacing = None
    origin = None
    velocity = []
    density = []
    
    reading_velocity = False
    reading_density = False
    
    for i, line in enumerate(lines):
        if 'DIMENSIONS' in line:
            dims = list(map(int, line.split()[1:4]))
        elif 'SPACING' in line:
            spacing = list(map(float, line.split()[1:4]))
        elif 'ORIGIN' in line:
            origin = list(map(float, line.split()[1:4]))
        elif 'VECTORS Velocity' in line:
            reading_velocity = True
            reading_density = False
            continue
        elif 'SCALARS Density' in line or 'SCALARS' in line:
            reading_velocity = False
            reading_density = True
            continue
        elif 'LOOKUP_TABLE' in line:
            continue
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
    
    velocity = np.array(velocity)
    density = np.array(density)
    
    if dims is None:
        raise ValueError("無法從VTK檔案中讀取DIMENSIONS")
    
    nx, ny, nz = dims
    
    # Reshape數據
    if len(velocity) == nx * ny * nz:
        velocity = velocity.reshape((nz, ny, nx, 3))
    
    if len(density) == nx * ny * nz:
        density = density.reshape((nz, ny, nx))
    
    return {
        'dims': dims,
        'spacing': spacing,
        'origin': origin,
        'velocity': velocity,
        'density': density
    }


def plot_2d_slice(data, z_slice=0, output_dir='output_images'):
    """繪製2D切面的向量場和等值線圖"""
    velocity = data['velocity']
    density = data['density']
    dims = data['dims']
    nx, ny, nz = dims
    
    if z_slice >= nz:
        z_slice = nz // 2
    
    # 提取2D切面
    u = velocity[z_slice, :, :, 0]  # x方向速度
    v = velocity[z_slice, :, :, 1]  # y方向速度
    rho = density[z_slice, :, :]     # 密度
    
    # 計算速度大小
    speed = np.sqrt(u**2 + v**2)
    
    # 創建座標網格
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    
    return X, Y, u, v, speed, rho


def create_vector_field_image(X, Y, u, v, speed, filename, title="Velocity Field"):
    """創建向量場圖片"""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # 繪製速度大小的contour作為背景
    contourf = ax.contourf(X, Y, speed, levels=20, cmap='jet', alpha=0.6)
    plt.colorbar(contourf, ax=ax, label='Speed')
    
    # 繪製向量場（降低密度以便觀察）
    skip = max(1, X.shape[1] // 30)  # 調整箭頭密度
    quiver = ax.quiver(X[::skip, ::skip], Y[::skip, ::skip], 
                       u[::skip, ::skip], v[::skip, ::skip],
                       scale=None, scale_units='xy', angles='xy',
                       color='white', alpha=0.8, width=0.003)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(title)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"已保存向量場圖片: {filename}")


def create_contour_image(X, Y, field, filename, title="Contour Plot", levels=20, cmap='jet'):
    """創建等值線圖片"""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # 繪製filled contour
    contourf = ax.contourf(X, Y, field, levels=levels, cmap=cmap)
    plt.colorbar(contourf, ax=ax, label='Value')
    
    # 添加contour線
    contour = ax.contour(X, Y, field, levels=levels, colors='black', 
                         linewidths=0.5, alpha=0.3)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(title)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"已保存等值線圖片: {filename}")


def batch_process_vtk_files(vtk_dir='output', output_dir='output_images', z_slice=0):
    """批次處理所有VTK檔案"""
    # 創建輸出目錄
    vector_dir = os.path.join(output_dir, 'vectors')
    contour_speed_dir = os.path.join(output_dir, 'contour_speed')
    contour_density_dir = os.path.join(output_dir, 'contour_density')
    
    os.makedirs(vector_dir, exist_ok=True)
    os.makedirs(contour_speed_dir, exist_ok=True)
    os.makedirs(contour_density_dir, exist_ok=True)
    
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
            data = read_vtk_structured_points(vtk_file)
            
            # 提取2D切面
            X, Y, u, v, speed, rho = plot_2d_slice(data, z_slice=z_slice)
            
            # 生成向量場圖片
            vector_filename = os.path.join(vector_dir, f'{basename}_vector.png')
            create_vector_field_image(X, Y, u, v, speed, vector_filename, 
                                     title=f'Velocity Field - {basename}')
            
            # 生成速度等值線圖片
            contour_speed_filename = os.path.join(contour_speed_dir, f'{basename}_speed.png')
            create_contour_image(X, Y, speed, contour_speed_filename,
                                title=f'Speed Contour - {basename}', cmap='jet')
            
            # 生成密度等值線圖片
            contour_density_filename = os.path.join(contour_density_dir, f'{basename}_density.png')
            create_contour_image(X, Y, rho, contour_density_filename,
                                title=f'Density Contour - {basename}', cmap='viridis')
            
            processed_files.append({
                'vtk': vtk_file,
                'vector': vector_filename,
                'speed': contour_speed_filename,
                'density': contour_density_filename
            })
            
        except Exception as e:
            print(f"處理 {vtk_file} 時發生錯誤: {e}")
            continue
    
    return processed_files


def create_animation_from_images(image_files, output_file, fps=10, title_prefix=""):
    """從圖片序列創建動畫"""
    if not image_files:
        print("沒有圖片可以製作動畫")
        return
    
    print(f"\n正在創建動畫: {output_file}")
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis('off')
    
    # 讀取第一張圖片
    first_img = plt.imread(image_files[0])
    im = ax.imshow(first_img)
    
    title_text = ax.text(0.5, 1.05, '', transform=ax.transAxes, 
                         ha='center', fontsize=12)
    
    def update(frame):
        img = plt.imread(image_files[frame])
        im.set_array(img)
        title_text.set_text(f'{title_prefix} - Frame {frame+1}/{len(image_files)}')
        return [im, title_text]
    
    anim = FuncAnimation(fig, update, frames=len(image_files), 
                        interval=1000//fps, blit=True)
    
    # 保存為GIF
    writer = PillowWriter(fps=fps)
    anim.save(output_file, writer=writer, dpi=150)
    plt.close()
    
    print(f"動畫已保存: {output_file}")


def main():
    """主程序"""
    print("=" * 60)
    print("VTK檔案批次處理 - 生成向量場和等值線動畫")
    print("=" * 60)
    
    # 設置參數
    vtk_dir = 'output'
    output_dir = 'output_images'
    z_slice = 0  # 使用z=0的切面（對於2D模擬，通常只有一個切面）
    fps = 10     # 動畫幀率
    
    # 批次處理VTK檔案
    print(f"\n開始批次處理VTK檔案...")
    processed_files = batch_process_vtk_files(vtk_dir, output_dir, z_slice)
    
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
    
    print("\n" + "=" * 60)
    print("完成！")
    print("=" * 60)
    print(f"\n生成的圖片位於: {output_dir}")
    print(f"  - 向量場圖片: {os.path.join(output_dir, 'vectors')}")
    print(f"  - 速度等值線: {os.path.join(output_dir, 'contour_speed')}")
    print(f"  - 密度等值線: {os.path.join(output_dir, 'contour_density')}")
    print(f"\n生成的動畫:")
    print(f"  - {vector_animation}")
    print(f"  - {speed_animation}")
    print(f"  - {density_animation}")
    print("=" * 60)


if __name__ == '__main__':
    main()
