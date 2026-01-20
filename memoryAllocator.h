#ifndef MEMORYALLOCATOR_H
#define MEMORYALLOCATOR_H
#include <vector>
#include <iostream>
#include "variables.h"
using namespace std ; 
//=============================================================================
// 全域 buffer 容器（使用 std::vector 自動管理記憶體）
//=============================================================================

// Y 方向權重 buffer
vector<double> YPara0_buf;
vector<double> YPara2_buf;

// Xi 方向權重 buffer (F1~F8)
vector<double> XiParaF1_buf;
std::vector<double> XiParaF2_buf;
std::vector<double> XiParaF3_buf;
std::vector<double> XiParaF4_buf;
std::vector<double> XiParaF5_buf;
std::vector<double> XiParaF6_buf;
std::vector<double> XiParaF7_buf;
std::vector<double> XiParaF8_buf;

// BFL Y 方向權重 buffer
std::vector<double> YBFLParaF1_buf;
std::vector<double> YBFLParaF3_buf;
std::vector<double> YBFLParaF7_buf;
std::vector<double> YBFLParaF8_buf;

// BFL Xi 方向權重 buffer
std::vector<double> XiBFLParaF1_buf;
std::vector<double> XiBFLParaF3_buf;
std::vector<double> XiBFLParaF7_buf;
std::vector<double> XiBFLParaF8_buf;

//=============================================================================
// 記憶體配置函數
//=============================================================================

/**
 * @brief 為單個權重陣列配置連續記憶體並設定指標
 *
 * @param pointerArray 指標陣列 (例如 YPara0_h)
 * @param buffer vector buffer
 * @param numWeights 權重數量（通常是 7）
 * @param numPoints 每個權重對應的點數
 */
inline void AllocateWeightArray(
    double* pointerArray[7],
    std::vector<double>& buffer,
    int numWeights,
    int numPoints
) {
    // 1. 配置連續記憶體
    buffer.resize(numWeights * numPoints, 0.0);

    // 2. 設定指標，讓 pointerArray[p] 指向對應的記憶體段
    for(int p = 0; p < numWeights; ++p) {
        pointerArray[p] = buffer.data() + p * numPoints;
    }
}

/**
 * @brief 配置所有權重陣列的記憶體（主函數）
 *
 * 呼叫順序：
 * 1. 在 main() 開始時呼叫此函數
 * 2. 之後就可以使用 YPara0_h[p][i] 的方式訪問
 */
inline void AllocateAllWeightArrays() {
    const int lenY = NY6;
    const int lenXi = (NY6+7) * NZ6; //Z 方向曲面插值需要 

    std::cout << "=========================================" << std::endl;
    std::cout << "  Allocate memory for XiPara[7][NY6*NZ6] \ YPara[7][(NY6+7) * NZ6]" << std::endl;
    std::cout << "=========================================" << std::endl;

    // 外部宣告的指標陣列（假設已在 main.cpp 中宣告）
    extern double* YPara0_h[7];
    extern double* YPara2_h[7];
    extern double* XiParaF1_h[7];
    extern double* XiParaF2_h[7];
    extern double* XiParaF3_h[7];
    extern double* XiParaF4_h[7];
    extern double* XiParaF5_h[7];
    extern double* XiParaF6_h[7];
    extern double* XiParaF7_h[7];
    extern double* XiParaF8_h[7];
    extern double* YBFLParaF1_h[7];
    extern double* YBFLParaF3_h[7];
    extern double* YBFLParaF7_h[7];
    extern double* YBFLParaF8_h[7];
    extern double* XiBFLParaF1_h[7];
    extern double* XiBFLParaF3_h[7];
    extern double* XiBFLParaF7_h[7];
    extern double* XiBFLParaF8_h[7];

    //-------------------------------------------------------------------------
    // Y 方向插值權重
    //-------------------------------------------------------------------------
    std::cout << "[1/4] Y-direction interpolation weights..." << std::endl;
    AllocateWeightArray(YPara0_h, YPara0_buf, 7, lenY);
    AllocateWeightArray(YPara2_h, YPara2_buf, 7, lenY);
    std::cout << " Configuration: 2 arrays × 7 weights × " << lenY << " points " << std::endl;
    std::cout << " Size: " << (2 * 7 * lenY * sizeof(double)) / 1024.0 << " KB" << std::endl;

    //-------------------------------------------------------------------------
    // Xi 方向插值權重 (F1~F8)
    //-------------------------------------------------------------------------
    std::cout << "[2/4] Xi-direction interpolation weight (F1~F8)..." << std::endl;
    AllocateWeightArray(XiParaF1_h, XiParaF1_buf, 7, lenXi);
    AllocateWeightArray(XiParaF2_h, XiParaF2_buf, 7, lenXi);
    AllocateWeightArray(XiParaF3_h, XiParaF3_buf, 7, lenXi);
    AllocateWeightArray(XiParaF4_h, XiParaF4_buf, 7, lenXi);
    AllocateWeightArray(XiParaF5_h, XiParaF5_buf, 7, lenXi);
    AllocateWeightArray(XiParaF6_h, XiParaF6_buf, 7, lenXi);
    AllocateWeightArray(XiParaF7_h, XiParaF7_buf, 7, lenXi);
    AllocateWeightArray(XiParaF8_h, XiParaF8_buf, 7, lenXi);
    std::cout << " Configuration: 8 arrays × 7 weights × " << lenXi << " points" << std::endl;
    std::cout << " Size: " << (8 * 7 * lenXi * sizeof(double)) / (1024.0 * 1024.0) << " MB" << std::endl;

    //-------------------------------------------------------------------------
    // BFL Y 方向權重
    //-------------------------------------------------------------------------
    std::cout << "[3/4] BFL Y-direction interpolation weights..." << std::endl;
    AllocateWeightArray(YBFLParaF1_h, YBFLParaF1_buf, 7, lenY);
    AllocateWeightArray(YBFLParaF3_h, YBFLParaF3_buf, 7, lenY);
    AllocateWeightArray(YBFLParaF7_h, YBFLParaF7_buf, 7, lenY);
    AllocateWeightArray(YBFLParaF8_h, YBFLParaF8_buf, 7, lenY);
    std::cout << " Configuration: 4 arrays × 7 weights × " << lenY << " points" << std::endl;
    std::cout << " Size: " << (4 * 7 * lenY * sizeof(double)) / 1024.0 << " KB" << std::endl;

    //-------------------------------------------------------------------------
    // BFL Xi 方向權重
    //-------------------------------------------------------------------------
    std::cout << "[4/4] Xi-direction interpolation weight..." << std::endl;
    AllocateWeightArray(XiBFLParaF1_h, XiBFLParaF1_buf, 7, lenXi);
    AllocateWeightArray(XiBFLParaF3_h, XiBFLParaF3_buf, 7, lenXi);
    AllocateWeightArray(XiBFLParaF7_h, XiBFLParaF7_buf, 7, lenXi);
    AllocateWeightArray(XiBFLParaF8_h, XiBFLParaF8_buf, 7, lenXi);
    std::cout << " Configuration: 4 arrays × 7 weights × " << lenXi << " points" << std::endl;
    std::cout << " Size: " << (4 * 7 * lenXi * sizeof(double)) / (1024.0 * 1024.0) << " MB" << std::endl;

    //-------------------------------------------------------------------------
    // 統計總計
    //-------------------------------------------------------------------------
    double totalMB = (
        2 * 7 * lenY +      // Y 方向
        8 * 7 * lenXi +     // Xi 方向
        4 * 7 * lenY +      // BFL Y
        4 * 7 * lenXi       // BFL Xi
    ) * sizeof(double) / (1024.0 * 1024.0);

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "Total memory usage: " << totalMB << " MB" << std::endl;
    std::cout << "=======================================" << std::endl;
    std::cout << "Memory allocation complete!" << std::endl;
    std::cout << "==========================================" << std::endl;
}

//=============================================================================
// 記憶體釋放（vector 自動釋放，但提供清空功能）
//=============================================================================

/**
 * @brief 清空所有 buffer（釋放記憶體）
 *
 * 注意：vector 會在程式結束時自動釋放，此函數是可選的
 */
inline void FreeAllWeightArrays() {
    std::cout << "Release weighted array memory..." << std::endl;

    // 清空所有 buffer（vector 會自動釋放記憶體）
    YPara0_buf.clear();
    YPara0_buf.shrink_to_fit();
    YPara2_buf.clear();
    YPara2_buf.shrink_to_fit();

    XiParaF1_buf.clear();
    XiParaF1_buf.shrink_to_fit();
    XiParaF2_buf.clear();
    XiParaF2_buf.shrink_to_fit();
    XiParaF3_buf.clear();
    XiParaF3_buf.shrink_to_fit();
    XiParaF4_buf.clear();
    XiParaF4_buf.shrink_to_fit();
    XiParaF5_buf.clear();
    XiParaF5_buf.shrink_to_fit();
    XiParaF6_buf.clear();
    XiParaF6_buf.shrink_to_fit();
    XiParaF7_buf.clear();
    XiParaF7_buf.shrink_to_fit();
    XiParaF8_buf.clear();
    XiParaF8_buf.shrink_to_fit();

    YBFLParaF1_buf.clear();
    YBFLParaF1_buf.shrink_to_fit();
    YBFLParaF3_buf.clear();
    YBFLParaF3_buf.shrink_to_fit();
    YBFLParaF7_buf.clear();
    YBFLParaF7_buf.shrink_to_fit();
    YBFLParaF8_buf.clear();
    YBFLParaF8_buf.shrink_to_fit();

    XiBFLParaF1_buf.clear();
    XiBFLParaF1_buf.shrink_to_fit();
    XiBFLParaF3_buf.clear();
    XiBFLParaF3_buf.shrink_to_fit();
    XiBFLParaF7_buf.clear();
    XiBFLParaF7_buf.shrink_to_fit();
    XiBFLParaF8_buf.clear();
    XiBFLParaF8_buf.shrink_to_fit();

    std::cout << "✓ Memory has been released!" << std::endl;
}



#endif // MEMORYALLOCATOR_H
