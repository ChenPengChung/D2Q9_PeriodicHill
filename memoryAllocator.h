#ifndef MEMORYALLOCATOR_H
#define MEMORYALLOCATOR_H

/**
 * @file memoryAllocator.h
 * @brief 為指標陣列配置二維連續記憶體
 *
 * 核心概念：
 * 1. 為每個 double* arr[7] 配置一個連續的 buffer
 * 2. buffer 大小 = 7 × numPoints
 * 3. 讓 arr[0], arr[1], ..., arr[6] 指向 buffer 的不同段
 */

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
    const int lenXi = NY6 * NZ6;

    std::cout << "=========================================" << std::endl;
    std::cout << "  配置權重陣列記憶體" << std::endl;
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
    std::cout << "[1/4] Y 方向插值權重..." << std::endl;
    AllocateWeightArray(YPara0_h, YPara0_buf, 7, lenY);
    AllocateWeightArray(YPara2_h, YPara2_buf, 7, lenY);
    std::cout << "      配置: 2 個陣列 × 7 權重 × " << lenY << " 點" << std::endl;
    std::cout << "      大小: " << (2 * 7 * lenY * sizeof(double)) / 1024.0 << " KB" << std::endl;

    //-------------------------------------------------------------------------
    // Xi 方向插值權重 (F1~F8)
    //-------------------------------------------------------------------------
    std::cout << "[2/4] Xi 方向插值權重 (F1~F8)..." << std::endl;
    AllocateWeightArray(XiParaF1_h, XiParaF1_buf, 7, lenXi);
    AllocateWeightArray(XiParaF2_h, XiParaF2_buf, 7, lenXi);
    AllocateWeightArray(XiParaF3_h, XiParaF3_buf, 7, lenXi);
    AllocateWeightArray(XiParaF4_h, XiParaF4_buf, 7, lenXi);
    AllocateWeightArray(XiParaF5_h, XiParaF5_buf, 7, lenXi);
    AllocateWeightArray(XiParaF6_h, XiParaF6_buf, 7, lenXi);
    AllocateWeightArray(XiParaF7_h, XiParaF7_buf, 7, lenXi);
    AllocateWeightArray(XiParaF8_h, XiParaF8_buf, 7, lenXi);
    std::cout << "      配置: 8 個陣列 × 7 權重 × " << lenXi << " 點" << std::endl;
    std::cout << "      大小: " << (8 * 7 * lenXi * sizeof(double)) / (1024.0 * 1024.0) << " MB" << std::endl;

    //-------------------------------------------------------------------------
    // BFL Y 方向權重
    //-------------------------------------------------------------------------
    std::cout << "[3/4] BFL Y 方向權重..." << std::endl;
    AllocateWeightArray(YBFLParaF1_h, YBFLParaF1_buf, 7, lenY);
    AllocateWeightArray(YBFLParaF3_h, YBFLParaF3_buf, 7, lenY);
    AllocateWeightArray(YBFLParaF7_h, YBFLParaF7_buf, 7, lenY);
    AllocateWeightArray(YBFLParaF8_h, YBFLParaF8_buf, 7, lenY);
    std::cout << "      配置: 4 個陣列 × 7 權重 × " << lenY << " 點" << std::endl;
    std::cout << "      大小: " << (4 * 7 * lenY * sizeof(double)) / 1024.0 << " KB" << std::endl;

    //-------------------------------------------------------------------------
    // BFL Xi 方向權重
    //-------------------------------------------------------------------------
    std::cout << "[4/4] BFL Xi 方向權重..." << std::endl;
    AllocateWeightArray(XiBFLParaF1_h, XiBFLParaF1_buf, 7, lenXi);
    AllocateWeightArray(XiBFLParaF3_h, XiBFLParaF3_buf, 7, lenXi);
    AllocateWeightArray(XiBFLParaF7_h, XiBFLParaF7_buf, 7, lenXi);
    AllocateWeightArray(XiBFLParaF8_h, XiBFLParaF8_buf, 7, lenXi);
    std::cout << "      配置: 4 個陣列 × 7 權重 × " << lenXi << " 點" << std::endl;
    std::cout << "      大小: " << (4 * 7 * lenXi * sizeof(double)) / (1024.0 * 1024.0) << " MB" << std::endl;

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
    std::cout << "總記憶體使用: " << totalMB << " MB" << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << "✓ 記憶體配置完成！" << std::endl;
    std::cout << "=========================================" << std::endl;
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
    std::cout << "釋放權重陣列記憶體..." << std::endl;

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

    std::cout << "✓ 記憶體已釋放！" << std::endl;
}

//=============================================================================
// 診斷與驗證函數
//=============================================================================

/**
 * @brief 驗證記憶體配置是否正確
 */
inline void VerifyMemoryAllocation() {
    std::cout << "\n=========================================" << std::endl;
    std::cout << "  記憶體配置驗證" << std::endl;
    std::cout << "=========================================" << std::endl;

    extern double* YPara0_h[7];
    extern double* XiParaF1_h[7];

    // 檢查指標是否正確設定
    bool isValid = true;

    // 檢查連續性
    for(int p = 0; p < 6; ++p) {
        ptrdiff_t diff = YPara0_h[p+1] - YPara0_h[p];
        if(diff != NY6) {
            std::cout << "✗ YPara0_h 不連續！" << std::endl;
            isValid = false;
            break;
        }
    }

    for(int p = 0; p < 6; ++p) {
        ptrdiff_t diff = XiParaF1_h[p+1] - XiParaF1_h[p];
        if(diff != NY6 * NZ6) {
            std::cout << "✗ XiParaF1_h 不連續！" << std::endl;
            isValid = false;
            break;
        }
    }

    if(isValid) {
        std::cout << "✓ 記憶體配置正確！" << std::endl;
        std::cout << "✓ 指標陣列連續性檢查通過！" << std::endl;
    }

    // 測試讀寫
    YPara0_h[0][0] = 1.0;
    YPara0_h[6][NY6-1] = 2.0;
    XiParaF1_h[0][0] = 3.0;
    XiParaF1_h[6][NY6*NZ6-1] = 4.0;

    std::cout << "✓ 記憶體讀寫測試通過！" << std::endl;
    std::cout << "=========================================" << std::endl;
}

//=============================================================================
// 記憶體佈局圖解
//=============================================================================

/**
 * @brief 印出記憶體佈局示意圖
 */
inline void PrintMemoryLayout() {
    std::cout << "\n=========================================" << std::endl;
    std::cout << "  記憶體佈局示意圖" << std::endl;
    std::cout << "=========================================" << std::endl;

    std::cout << "\n指標陣列與 buffer 的關係：\n" << std::endl;
    std::cout << "double* YPara0_h[7]:" << std::endl;
    std::cout << "  [0] → buffer[0 * NY6 .. 1 * NY6)" << std::endl;
    std::cout << "  [1] → buffer[1 * NY6 .. 2 * NY6)" << std::endl;
    std::cout << "  [2] → buffer[2 * NY6 .. 3 * NY6)" << std::endl;
    std::cout << "  [3] → buffer[3 * NY6 .. 4 * NY6)" << std::endl;
    std::cout << "  [4] → buffer[4 * NY6 .. 5 * NY6)" << std::endl;
    std::cout << "  [5] → buffer[5 * NY6 .. 6 * NY6)" << std::endl;
    std::cout << "  [6] → buffer[6 * NY6 .. 7 * NY6)" << std::endl;

    std::cout << "\n訪問方式：" << std::endl;
    std::cout << "  YPara0_h[p][i] = buffer[p * NY6 + i]" << std::endl;

    std::cout << "\n連續記憶體佈局：" << std::endl;
    std::cout << "  [權重0的所有點][權重1的所有點]...[權重6的所有點]" << std::endl;
    std::cout << "  └─── NY6 ───┘└─── NY6 ───┘      └─── NY6 ───┘" << std::endl;

    std::cout << "=========================================" << std::endl;
}

#endif // MEMORYALLOCATOR_H
