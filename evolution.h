#ifndef EVOLUTION_FILE
#define EVOLUTION_FILE

#include "interpolationHillISLBM.h"
#include "MRT_Process.h"
#include "MRT_Matrix.h"
//1.
//面法向量為ｅ_y的單位面積質量流率
double ModifydRho_F1( double F1_in, double f2_old)
{
    double drho = F1_in - f2_old;
    return drho;
}
//2.
//面法向量為－ｅ_y的單位面積質量流率
double ModifydRho_F3( double F3_in, double f4_old)
{
    double drho = F3_in - f4_old;
    return drho;
}
#endif
