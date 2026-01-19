#ifndef PARAVIEWOUTPUT_H
#define PARAVIEWOUTPUT_H

/**
 * @file paraviewOutput.h
 * @brief ParaView 可視化輸出函數
 *
 * 支援格式：
 * - VTK Legacy Structured Grid (.vtk)
 * - CSV + Python 腳本
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include "variables.h"
#include "model.h"

//=============================================================================
// VTK Structured Grid 格式輸出（推薦）
//=============================================================================

/**
 * @brief 輸出 VTK Structured Grid 格式（ParaView 可直接讀取）
 *
 * @param timestep 時間步
 * @param y_global Y 方向座標陣列
 * @param z_global Z 方向座標陣列（2D flatten）
 * @param rho 密度場
 * @param v Y 方向速度
 * @param w Z 方向速度
 * @param outputDir 輸出目錄（預設為 "output"）
 *
 * 使用方式：
 * 在主迴圈中每隔 N 步呼叫一次
 *
 * 範例：
 * if(timestep % 1000 == 0) {
 *     OutputVTK(timestep, y_global, z_global, rho, v, w);
 * }
 */
inline void OutputVTK(
    int timestep,
    const double* y_global,
    const double* z_global,
    const double* rho,
    const double* v,
    const double* w,
    const std::string& outputDir = "output"
) {
    // 1. 建立輸出目錄（如果不存在）
    #ifdef _WIN32
        system(("mkdir " + outputDir + " 2>nul").c_str());
    #else
        system(("mkdir -p " + outputDir).c_str());
    #endif

    // 2. 建立檔案名稱
    std::ostringstream filename;
    filename << outputDir << "/flow_" << std::setfill('0') << std::setw(6) << timestep << ".vtk";

    // 3. 開啟檔案
    std::ofstream file(filename.str());
    if(!file.is_open()) {
        std::cerr << "錯誤：無法開啟檔案 " << filename.str() << std::endl;
        return;
    }

    std::cout << "輸出 VTK 檔案: " << filename.str() << std::endl;

    // 4. 寫入 VTK 標頭
    file << "# vtk DataFile Version 3.0\n";
    file << "D2Q9 Periodic Hill Flow at timestep " << timestep << "\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";

    // 5. 寫入網格尺寸（只輸出內部計算區域，不含 buffer）
    const int ny_out = NY6 - 6;  // 去掉 buffer
    const int nz_out = NZ6 - 6;
    file << "DIMENSIONS " << ny_out << " " << nz_out << " 1\n";
    file << "POINTS " << (ny_out * nz_out) << " double\n";

    // 6. 寫入網格點座標（2D 流場，z 座標設為 0）
    for(int k = 3; k < NZ6 - 3; ++k) {
        for(int j = 3; j < NY6 - 3; ++j) {
            int idx = j * NZ6 + k;
            file << std::scientific << std::setprecision(10)
                 << y_global[j] << " "
                 << z_global[idx] << " "
                 << "0.0\n";
        }
    }

    // 7. 寫入點數據
    file << "\nPOINT_DATA " << (ny_out * nz_out) << "\n";

    //-------------------------------------------------------------------------
    // 7.1 密度場（標量）
    //-------------------------------------------------------------------------
    file << "\nSCALARS Density double 1\n";
    file << "LOOKUP_TABLE default\n";
    for(int k = 3; k < NZ6 - 3; ++k) {
        for(int j = 3; j < NY6 - 3; ++j) {
            int idx = j * NZ6 + k;
            file << std::scientific << std::setprecision(10) << rho[idx] << "\n";
        }
    }

    //-------------------------------------------------------------------------
    // 7.2 速度場（向量）
    //-------------------------------------------------------------------------
    file << "\nVECTORS Velocity double\n";
    for(int k = 3; k < NZ6 - 3; ++k) {
        for(int j = 3; j < NY6 - 3; ++j) {
            int idx = j * NZ6 + k;
            file << std::scientific << std::setprecision(10)
                 << v[idx] << " "    // Y 方向速度
                 << w[idx] << " "    // Z 方向速度
                 << "0.0\n";         // 2D 流場，第三個分量為 0
        }
    }

    //-------------------------------------------------------------------------
    // 7.3 速度大小（標量）
    //-------------------------------------------------------------------------
    file << "\nSCALARS VelocityMagnitude double 1\n";
    file << "LOOKUP_TABLE default\n";
    for(int k = 3; k < NZ6 - 3; ++k) {
        for(int j = 3; j < NY6 - 3; ++j) {
            int idx = j * NZ6 + k;
            double vmag = std::sqrt(v[idx]*v[idx] + w[idx]*w[idx]);
            file << std::scientific << std::setprecision(10) << vmag << "\n";
        }
    }

    //-------------------------------------------------------------------------
    // 7.4 渦度（標量，2D 情況下是 ∂w/∂y - ∂v/∂z）
    //-------------------------------------------------------------------------
    file << "\nSCALARS Vorticity double 1\n";
    file << "LOOKUP_TABLE default\n";
    for(int k = 3; k < NZ6 - 3; ++k) {
        for(int j = 3; j < NY6 - 3; ++j) {
            // 中央差分計算渦度
            double dwdy = 0.0, dvdz = 0.0;

            if(j > 3 && j < NY6 - 4) {
                int idx_yp = (j+1) * NZ6 + k;
                int idx_ym = (j-1) * NZ6 + k;
                double dy = y_global[j+1] - y_global[j-1];
                dwdy = (w[idx_yp] - w[idx_ym]) / dy;
            }

            if(k > 3 && k < NZ6 - 4) {
                int idx_zp = j * NZ6 + (k+1);
                int idx_zm = j * NZ6 + (k-1);
                double dz = z_global[idx_zp] - z_global[idx_zm];
                dvdz = (v[idx_zp] - v[idx_zm]) / dz;
            }

            double vorticity = dwdy - dvdz;
            file << std::scientific << std::setprecision(10) << vorticity << "\n";
        }
    }

    file.close();
    std::cout << "✓ VTK 輸出完成：" << filename.str() << std::endl;
}

//=============================================================================
// 輸出統計資訊（可選）
//=============================================================================

/**
 * @brief 輸出流場統計資訊
 *
 * @param timestep 時間步
 * @param rho 密度場
 * @param v Y 方向速度
 * @param w Z 方向速度
 * @param logFile 日誌檔案名稱
 */
inline void OutputStatistics(
    int timestep,
    const double* rho,
    const double* v,
    const double* w,
    const std::string& logFile = "output/statistics.log"
) {
    // 計算統計量
    double rho_min = 1e10, rho_max = -1e10, rho_avg = 0.0;
    double v_min = 1e10, v_max = -1e10, v_avg = 0.0;
    double w_min = 1e10, w_max = -1e10, w_avg = 0.0;
    double vmag_max = 0.0;

    int count = 0;
    for(int j = 3; j < NY6 - 3; ++j) {
        for(int k = 3; k < NZ6 - 3; ++k) {
            int idx = j * NZ6 + k;

            rho_min = std::min(rho_min, rho[idx]);
            rho_max = std::max(rho_max, rho[idx]);
            rho_avg += rho[idx];

            v_min = std::min(v_min, v[idx]);
            v_max = std::max(v_max, v[idx]);
            v_avg += v[idx];

            w_min = std::min(w_min, w[idx]);
            w_max = std::max(w_max, w[idx]);
            w_avg += w[idx];

            double vmag = std::sqrt(v[idx]*v[idx] + w[idx]*w[idx]);
            vmag_max = std::max(vmag_max, vmag);

            count++;
        }
    }

    rho_avg /= count;
    v_avg /= count;
    w_avg /= count;

    // 寫入日誌檔案
    std::ofstream file;
    if(timestep == 0) {
        file.open(logFile);
        file << "# Timestep, Rho_avg, Rho_min, Rho_max, V_avg, V_max, W_avg, W_max, Vmag_max\n";
    } else {
        file.open(logFile, std::ios::app);
    }

    if(file.is_open()) {
        file << timestep << " "
             << std::scientific << std::setprecision(6)
             << rho_avg << " " << rho_min << " " << rho_max << " "
             << v_avg << " " << v_max << " "
             << w_avg << " " << w_max << " "
             << vmag_max << "\n";
        file.close();
    }

    // 印出到終端
    std::cout << "統計 [t=" << timestep << "]: "
              << "ρ=" << rho_avg << " "
              << "V=" << v_avg << " "
              << "Vmag_max=" << vmag_max << std::endl;
}

//=============================================================================
// 輸出山丘幾何（一次性輸出，用於 ParaView 顯示）
//=============================================================================

/**
 * @brief 輸出山丘幾何（VTK PolyData 格式）
 *
 * 用於在 ParaView 中顯示山丘邊界
 */
inline void OutputHillGeometry(const std::string& filename = "output/hill.vtk") {
    std::ofstream file(filename);
    if(!file.is_open()) {
        std::cerr << "錯誤：無法開啟檔案 " << filename << std::endl;
        return;
    }

    // 生成山丘表面點
    const int npoints = 200;
    file << "# vtk DataFile Version 3.0\n";
    file << "Hill Geometry\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    file << "POINTS " << (npoints * 2) << " double\n";

    // 左半丘
    for(int i = 0; i < npoints; ++i) {
        double y = (54.0/28.0) * i / (npoints - 1);
        double z = HillFunction(y);
        file << y << " " << z << " 0.0\n";
    }

    // 右半丘
    for(int i = 0; i < npoints; ++i) {
        double y = LY - (54.0/28.0) + (54.0/28.0) * i / (npoints - 1);
        double z = HillFunction(y);
        file << y << " " << z << " 0.0\n";
    }

    // 線段連接
    file << "\nLINES 2 " << (npoints * 2 + 2) << "\n";
    file << npoints << " ";
    for(int i = 0; i < npoints; ++i) {
        file << i << " ";
    }
    file << "\n";

    file << npoints << " ";
    for(int i = 0; i < npoints; ++i) {
        file << (npoints + i) << " ";
    }
    file << "\n";

    file.close();
    std::cout << "✓ 山丘幾何輸出完成：" << filename << std::endl;
}

//=============================================================================
// 整合輸出函數（建議使用）
//=============================================================================

/**
 * @brief 整合輸出函數（推薦在主迴圈中使用）
 *
 * @param timestep 時間步
 * @param y_global Y 方向座標
 * @param z_global Z 方向座標
 * @param rho 密度場
 * @param v Y 方向速度
 * @param w Z 方向速度
 * @param outputVTK 是否輸出 VTK（預設 true）
 * @param outputStats 是否輸出統計（預設 true）
 *
 * 使用範例：
 * if(timestep % 1000 == 0) {
 *     OutputFlowField(timestep, y_global, z_global, rho, v, w);
 * }
 */
inline void OutputFlowField(
    int timestep,
    const double* y_global,
    const double* z_global,
    const double* rho,
    const double* v,
    const double* w,
    bool outputVTK = true,
    bool outputStats = true
) {
    if(outputVTK) {
        OutputVTK(timestep, y_global, z_global, rho, v, w);
    }

    if(outputStats) {
        OutputStatistics(timestep, rho, v, w);
    }
}

//=============================================================================
// ParaView 批次處理腳本生成（Python）
//=============================================================================

/**
 * @brief 生成 ParaView Python 腳本（自動載入並渲染）
 */
inline void GenerateParaViewScript(const std::string& filename = "output/visualize.py") {
    std::ofstream file(filename);
    if(!file.is_open()) {
        std::cerr << "錯誤：無法開啟檔案 " << filename << std::endl;
        return;
    }

    file << "# ParaView Python Script\n";
    file << "# 使用方式: pvpython visualize.py\n";
    file << "# 或在 ParaView 中: Tools > Python Shell > Run Script\n\n";

    file << "from paraview.simple import *\n\n";

    file << "# 載入 VTK 檔案序列\n";
    file << "reader = LegacyVTKReader(FileNames=['flow_*.vtk'])\n";
    file << "reader.UpdatePipeline()\n\n";

    file << "# 建立視圖\n";
    file << "renderView = CreateView('RenderView')\n";
    file << "renderView.ViewSize = [1920, 1080]\n";
    file << "renderView.Background = [1, 1, 1]  # 白色背景\n\n";

    file << "# 顯示數據\n";
    file << "display = Show(reader, renderView)\n\n";

    file << "# 設定顏色映射（速度大小）\n";
    file << "ColorBy(display, ('POINTS', 'VelocityMagnitude'))\n";
    file << "display.RescaleTransferFunctionToDataRange()\n\n";

    file << "# 設定色條\n";
    file << "colorBar = GetScalarBar(GetColorTransferFunction('VelocityMagnitude'), renderView)\n";
    file << "colorBar.Title = 'Velocity Magnitude'\n";
    file << "colorBar.ComponentTitle = ''\n";
    file << "display.SetScalarBarVisibility(renderView, True)\n\n";

    file << "# 載入山丘幾何\n";
    file << "hill = LegacyVTKReader(FileNames=['hill.vtk'])\n";
    file << "hillDisplay = Show(hill, renderView)\n";
    file << "hillDisplay.LineWidth = 3.0\n";
    file << "hillDisplay.AmbientColor = [0, 0, 0]\n\n";

    file << "# 調整視角\n";
    file << "renderView.ResetCamera()\n\n";

    file << "# 渲染\n";
    file << "Render()\n\n";

    file << "# 保存截圖\n";
    file << "SaveScreenshot('flow_visualization.png', renderView, ImageResolution=[1920, 1080])\n\n";

    file << "print('視覺化完成！')\n";

    file.close();
    std::cout << "✓ ParaView 腳本生成完成：" << filename << std::endl;
}

#endif // PARAVIEWOUTPUT_H
