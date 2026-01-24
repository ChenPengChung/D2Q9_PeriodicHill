#ifndef GLOBALVARIABLES_H
#define GLOBALVARIABLES_H
/**
 * @file globalVariables.h
 * @brief 全域變數宣告（extern）
 *
 * 此檔案包含所有在 initialization.h 中使用的全域變數宣告
 * 實際定義應在 main.cpp 中
 */

#include "variables.h"
//=============================================================================
// 宏觀場變數
//=============================================================================
extern double rho[NY6 * NZ6];           // 密度場
extern double v[NY6 * NZ6];             // Y方向速度
extern double w[NY6 * NZ6];             // Z方向速度
extern double f[9][NY6 * NZ6];          // 分佈函數 (D2Q9)
//=============================================================================
// 外力項
//=============================================================================

extern double Force[2];                 // 外力項 (Fy, Fz)

//=============================================================================
// 網格座標
//=============================================================================

extern double y_global[NY6];            // Y方向物理座標
extern double z_global[NY6 * NZ6];      // Z方向物理座標（含山丘）
extern double xi_h[NZ6];                // 無因次化Z座標（不含山丘）
extern double nonuni_a;                 // 非均勻網格參數 a（預先計算一次）

//=============================================================================
// Y 方向插值權重（一般）
//=============================================================================
/*
extern double* YPara0_h[7];             // F1 使用
extern double* YPara2_h[7];             // F3 使用
*/
//降階版本
extern double* YPara0_h[3];             // F1 使用
extern double* YPara2_h[3];             // F3 使用
//=============================================================================
// Xi 方向插值權重（一般）F1~F8
//=============================================================================

extern double* XiParaF1_h[7];           // F1 (+Y,0)
extern double* XiParaF2_h[7];           // F2 (0,+Z)
extern double* XiParaF3_h[7];           // F3 (-Y,0)
extern double* XiParaF4_h[7];           // F4 (0,-Z)
extern double* XiParaF5_h[7];           // F5 (+Y,+Z)
extern double* XiParaF6_h[7];           // F6 (-Y,+Z)
extern double* XiParaF7_h[7];           // F7 (-Y,-Z)
extern double* XiParaF8_h[7];           // F8 (+Y,-Z)

#endif // GLOBALVARIABLES_H
