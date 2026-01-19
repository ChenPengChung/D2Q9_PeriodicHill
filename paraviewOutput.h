#ifndef PARAVIEWOUTPUT_H
#define PARAVIEWOUTPUT_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include "variables.h"

//=============================================================================
// VTK 輸出（密度 + 速度場）
//=============================================================================

void OutputVTK(
    int timestep,
    const double* y_global,
    const double* z_global,
    const double* rho,
    const double* v,
    const double* w,
    const std::string& outputDir = "output"
) {
    // 建立輸出目錄
    #ifdef _WIN32
        system(("mkdir " + outputDir + " 2>nul").c_str());
    #else
        system(("mkdir -p " + outputDir).c_str());
    #endif

    // 建立檔案名稱
    std::ostringstream filename;
    filename << outputDir << "/flow_" << std::setfill('0') << std::setw(6) << timestep << ".vtk";

    std::ofstream file(filename.str());
    if(!file.is_open()) {
        std::cerr << "錯誤：無法開啟檔案 " << filename.str() << std::endl;
        return;
    }

    // VTK 標頭
    file << "# vtk DataFile Version 3.0\n";
    file << "LBM Flow at timestep " << timestep << "\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";

    // 網格尺寸（去掉 buffer）
    const int ny_out = NY6 - 6;
    const int nz_out = NZ6 - 6;
    file << "DIMENSIONS " << ny_out << " " << nz_out << " 1\n";
    file << "POINTS " << (ny_out * nz_out) << " double\n";

    // 網格點座標
    for(int k = 3; k < NZ6 - 3; ++k) {
        for(int j = 3; j < NY6 - 3; ++j) {
            int idx = j * NZ6 + k;
            file << y_global[j] << " " << z_global[idx] << " 0.0\n";
        }
    }

    // 點數據
    file << "\nPOINT_DATA " << (ny_out * nz_out) << "\n";

    // 密度場
    file << "\nSCALARS Density double 1\n";
    file << "LOOKUP_TABLE default\n";
    for(int k = 3; k < NZ6 - 3; ++k) {
        for(int j = 3; j < NY6 - 3; ++j) {
            file << rho[j * NZ6 + k] << "\n";
        }
    }

    // 速度場（向量）
    file << "\nVECTORS Velocity double\n";
    for(int k = 3; k < NZ6 - 3; ++k) {
        for(int j = 3; j < NY6 - 3; ++j) {
            int idx = j * NZ6 + k;
            file << v[idx] << " " << w[idx] << " 0.0\n";
        }
    }

    // 速度大小
    file << "\nSCALARS VelocityMagnitude double 1\n";
    file << "LOOKUP_TABLE default\n";
    for(int k = 3; k < NZ6 - 3; ++k) {
        for(int j = 3; j < NY6 - 3; ++j) {
            int idx = j * NZ6 + k;
            double vmag = std::sqrt(v[idx]*v[idx] + w[idx]*w[idx]);
            file << vmag << "\n";
        }
    }

    file.close();
    std::cout << "✓ VTK 輸出: " << filename.str() << std::endl;
}

//=============================================================================
// 統計資訊輸出
//=============================================================================

void OutputStatistics(
    int timestep,
    const double* rho,
    const double* v,
    const double* w,
    const std::string& logFile = "output/statistics.log"
) {
    double rho_avg = 0.0, v_avg = 0.0, w_avg = 0.0;
    double v_max = 0.0, w_max = 0.0, vmag_max = 0.0;
    int count = 0;

    // 計算統計量
    for(int j = 3; j < NY6 - 3; ++j) {
        for(int k = 3; k < NZ6 - 3; ++k) {
            int idx = j * NZ6 + k;
            
            rho_avg += rho[idx];
            v_avg += v[idx];
            w_avg += w[idx];
            
            v_max = std::max(v_max, std::abs(v[idx]));
            w_max = std::max(w_max, std::abs(w[idx]));
            
            double vmag = std::sqrt(v[idx]*v[idx] + w[idx]*w[idx]);
            vmag_max = std::max(vmag_max, vmag);
            
            count++;
        }
    }

    rho_avg /= count;
    v_avg /= count;
    w_avg /= count;

    // 寫入日誌
    std::ofstream file;
    if(timestep == 0) {
        file.open(logFile);
        file << "# Timestep Rho_avg V_avg W_avg V_max W_max Vmag_max\n";
    } else {
        file.open(logFile, std::ios::app);
    }

    if(file.is_open()) {
        file << timestep << " "
             << rho_avg << " " << v_avg << " " << w_avg << " "
             << v_max << " " << w_max << " " << vmag_max << "\n";
        file.close();
    }

    // 終端輸出
    std::cout << "[t=" << timestep << "] "
              << "ρ=" << rho_avg << " "
              << "V_max=" << vmag_max << std::endl;
}

#endif // PARAVIEWOUTPUT_H