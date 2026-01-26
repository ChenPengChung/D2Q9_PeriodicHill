#ifndef INITIALIZATION_FILE
#define INITIALIZATION_FILE
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "globalVariables.h"
#include "initializationTool.h"
using namespace std;

void InitialUsingDftFunc() {
    //正規化離散粒子速度場
    double e[9][2]={{0.0,0.0},{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,1.0},{-1.0,1.0},{-1.0,-1.0},{1.0,-1.0}}; 
    //各個離散速度方向所應得權重、
    double W[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    //宏觀速度場的范數平方 
    double udot;
    for( int k = 0; k < NZ6;  k++ ) {
        for( int j = 0; j < NY6; j++ ) {
            const int index = j*NZ6 + k ;
            //初始條件 
            rho[index] = 1.0;
            v[index] = 0.0;
            w[index] = 0.0;
            udot = v[index]*v[index] + w[index]*w[index];
            //初始化物理空間計算點的插值後一般態分佈函數
            f[0][index] = W[0]*rho[index]*(1.0-1.5*udot);
            for( int dir = 1; dir <= 8; dir++ ) {
                f[dir][index] = W[dir] * rho[index] *( 1.0 +                 //0階項
                                                    3.0 *(e[dir][0] * v[index] + e[dir][1] * w[index])+          //1階項
                                                    4.5 *(e[dir][0] * v[index] + e[dir][1] * w[index] )*(e[dir][0] * v[index] + e[dir][1] * w[index] ) - 1.5*udot );
        }}}
        //離散化\ 宏觀\ 外立\ 場的初始化initilaoization of the discrete macroscopic force term
        Force[0] =  (8.0*niu*Uref)/(LZ*LZ)*5.0; //降低外力係數，原本是 5.0
        Force[1] = 0.0;  // Z 方向無外力
}
//建立Y(主流場方向)方向之均勻網格系統
//計算y_global
void GenerateMesh_Y() {
    double dy;
    int buffr = 3;

    if( Uniform_In_Ydir ){
        dy = LY / (double)(NY6-2*buffr-1); //主流場Steram-Wise方向作為Wet-node boundary 節點佈局
        for( int i = 0; i < NY6; i++ ){
            y_global[i] = dy * ((double)(i-buffr));//配合Hill Function進行座標平移
        }//物理空間計算點在外網格中為節點佈局，換言之，同一個物理空間計算點作為Lattic的中心點，卻作為外網格的節點
    } else {
        cout << "Mesh needs to be uniform in periodic hill problem, exit..." << endl ;
        exit(0);
    }
}

void GenerateMesh_Z() {
    int bufferlayer = 3; //單邊bufferlayer的厚度為3 
    if( Uniform_In_Zdir ){
        cout << "Mesh needs to be non-uniform in z-direction in periodic hill problem, exit..." << endl ;
        exit(0);
    }//Z方向做非均勻網格系統
    // 計算一次非均勻參數 a 並儲存，後續共用
    nonuni_a = GetNonuniParameter(); //計算最合適非均勻參數 
    const double a = nonuni_a;
    
    //計算(不含山丘)離散化無因次化Z座標
    for( int k = bufferlayer; k < NZ6-bufferlayer; k++ ){ //3~NZ6-4 
        xi_h[k] = tanhFunction( LXi, minSize, a, (k-3), (NZ6-7) ) - minSize/2.0;
    }
    //計算(含山丘)離散化全域z座標

    for( int j = 0; j < NY6; j++ ){
        double dy = LY / (double)(NY6-2*bufferlayer-1);
        y_global[j] = dy * ((double)(j-bufferlayer));//配合Hill Function做座標平移
        double total = LZ - HillFunction( y_global[j] ) - minSize;
        for( int k = bufferlayer; k < NZ6-bufferlayer; k++ ){
            z_global[j*NZ6+k] = tanhFunction( total, minSize, a, (k-3), (NZ6-7) ) + 
                                HillFunction( y_global[j] );
        }
        z_global[j*NZ6+2] = HillFunction( y_global[j] );
        z_global[j*NZ6+(NZ6-3)] = (double)LZ;
    }
}



void GetXiParameter(double* XiPara_h[7], double pos_z, double pos_y, int j , int k , int* cell_z){
    //double Extrapolation(double pos , double x1 , double x2 ,  double f0 , double f1 ){
    //double Lagrange_2nd (double pos , double x_i , double x1 , double x2 ){
    //double Lagrange_6th (double pos , double x_i , double x1 , double x2 , double x3 , double x4 , double x5 , double x6){
    //儲存索引使用原始 j，計算使用周期映射後的 j_calc
    const int j_store = j;  // 保留原始 j 用於儲存
    int j_calc = j;         // 用於計算 y_global 的週期映射
    if(j_calc<4) j_calc = j_calc + NY6-7 ;
    if(j_calc>NY6-5) j_calc = j_calc - (NY6-7) ;
    double pos_z0 = pos_z -  0.5*minSize - HillFunction(y_global[j_calc-1]);
    double pos_z1 = pos_z -  0.5*minSize - HillFunction(y_global[j_calc]);
    double pos_z2 = pos_z -  0.5*minSize - HillFunction(y_global[j_calc+1]);
    double L0 = LZ - HillFunction(y_global[j_calc-1]) - minSize;
    double L1 = LZ - HillFunction(y_global[j_calc]) - minSize;
    double L2 = LZ - HillFunction(y_global[j_calc+1]) - minSize;
    double pos_xi0 = (pos_z0) / L0; // 無因次化目標點
    double pos_xi1 = (pos_z1) / L1; // 無因次化目標點
    double pos_xi2 = (pos_z2) / L2; // 無因次化目標點
    //計算該高度在不同y值上的編號
    double index_z0 = Inverse_tanh_index( pos_z0 , L0 , minSize , nonuni_a , (NZ6-7) );
    double index_z1 = Inverse_tanh_index( pos_z1 , L1 , minSize , nonuni_a , (NZ6-7) );
    double index_z2 = Inverse_tanh_index( pos_z2 , L2 , minSize , nonuni_a , (NZ6-7) );//a 為 非均勻網格伸縮參數
    //第一套權重陣列
    double RelationXi_0[7] ; //(j_calc-1)
    if (index_z0 < 3) {
        for(int i = 2 ; i <=6 ; i++) {RelationXi_0[i] = 0.0 ;}
        RelationXi_0[0] = z_global[NZ6*(j_calc-1) + cell_z[NZ6*(j_store)+k + 0 * NY6*NZ6] + 0 ] ;
        RelationXi_0[1] = z_global[NZ6*(j_calc-1) + cell_z[NZ6*(j_store)+k + 0 * NY6*NZ6] + 1 ] ;}
    else if(index_z0 > NZ6-4){
        for(int i = 0 ; i <=4 ; i++) RelationXi_0[i] = 0.0 ;
        RelationXi_0[5] = z_global[NZ6*(j_calc-1) + cell_z[NZ6*(j_store)+k + 0 * NY6*NZ6] + 5 ] ;
        RelationXi_0[6] = z_global[NZ6*(j_calc-1) + cell_z[NZ6*(j_store)+k + 0 * NY6*NZ6] + 6 ] ;
    }else if(index_z0 < 6){
        for(int i = 3 ; i <=6 ; i++) RelationXi_0[i] = 0.0 ;
        for(int i = 0 ; i <=2 ; i++) {
            RelationXi_0[i] = z_global[NZ6*(j_calc-1) + cell_z[NZ6*(j_store)+k + 0 * NY6*NZ6] + i ] ;
        }}
    else if(index_z0 > NZ6-7){
        for(int i = 0 ; i <=3 ; i++) RelationXi_0[i] = 0.0 ;
        for(int i = 4 ; i <=6 ; i++) {
            RelationXi_0[i] = z_global[NZ6*(j_calc-1) + cell_z[NZ6*(j_store)+k + 0 * NY6*NZ6] + i ] ;
        }}
    else{
        for(int i = 0 ; i <7 ; i++){
        RelationXi_0[i] = z_global[NZ6*(j_calc-1) + cell_z[NZ6*(j_store)+k + 0 * NY6*NZ6] + i ] ;
    }}
    for(int i = 0 ; i <=6 ; i++){RelationXi_0[i] = (RelationXi_0[i] - HillFunction( y_global[j_calc-1] ) - 0.5*minSize) / L0 ;} //轉換為無因次化Z座標
    GetParameter_6th2( XiPara_h, pos_xi0 , RelationXi_0 , 0 , j_store , k , index_z0 ); //配置第一套適應性內插權重
    //第二套權重陣列
    double RelationXi_1[7] ; //(j_calc)
    if (index_z1 < 3) {
        for(int i = 2 ; i <=6 ; i++) {RelationXi_1[i] = 0.0 ;}
        RelationXi_1[0] = z_global[NZ6*j_calc + cell_z[NZ6*(j_store)+k + 1 * NY6*NZ6] + 0 ] ;
        RelationXi_1[1] = z_global[NZ6*j_calc + cell_z[NZ6*(j_store)+k + 1 * NY6*NZ6] + 1 ] ;}
    else if(index_z1 > NZ6-4){
        for(int i = 0 ; i <=4 ; i++) RelationXi_1[i] = 0.0 ;
        RelationXi_1[5] = z_global[NZ6*j_calc + cell_z[NZ6*(j_store)+k + 1 * NY6*NZ6] + 5 ] ;
        RelationXi_1[6] = z_global[NZ6*j_calc + cell_z[NZ6*(j_store)+k + 1 * NY6*NZ6] + 6 ] ;
    }else if(index_z1 < 6){
        for(int i = 3 ; i <=6 ; i++) RelationXi_1[i] = 0.0 ;
        for(int i = 0 ; i <=2 ; i++) {
            RelationXi_1[i] = z_global[NZ6*j_calc + cell_z[NZ6*(j_store)+k + 1 * NY6*NZ6] + i ] ;
        }}
    else if(index_z1 > NZ6-7){
        for(int i = 0 ; i <=3 ; i++) RelationXi_1[i] = 0.0 ;
        for(int i = 4 ; i <=6 ; i++) {
            RelationXi_1[i] = z_global[NZ6*j_calc + cell_z[NZ6*(j_store)+k + 1 * NY6*NZ6] + i ] ;
        }}
    else{
        for(int i = 0 ; i <7 ; i++){
        RelationXi_1[i] = z_global[NZ6*j_calc + cell_z[NZ6*(j_store)+k + 1 * NY6*NZ6] + i ] ;
    }}
    for(int i = 0 ; i <=6 ; i++){RelationXi_1[i] = (RelationXi_1[i] - HillFunction( y_global[j_calc] ) - 0.5*minSize) / L1 ;} //轉換為無因次化Z座標
    GetParameter_6th2( XiPara_h, pos_xi1 , RelationXi_1 , 1 , j_store , k , index_z1 ); //配置第二套適應性內插權重
    //第三套權重陣列
    double RelationXi_2[7] ; //(j_calc+1)
    if (index_z2 < 3) {
        for(int i = 2 ; i <=6 ; i++) {RelationXi_2[i] = 0.0 ;}
        RelationXi_2[0] = z_global[NZ6*(j_calc+1) + cell_z[NZ6*(j_store)+k + 2 * NY6*NZ6] + 0 ] ;
        RelationXi_2[1] = z_global[NZ6*(j_calc+1) + cell_z[NZ6*(j_store)+k + 2 * NY6*NZ6] + 1 ] ;}
    else if(index_z2 > NZ6-4){
        for(int i = 0 ; i <=4 ; i++) RelationXi_2[i] = 0.0 ;
        RelationXi_2[5] = z_global[NZ6*(j_calc+1) + cell_z[NZ6*(j_store)+k + 2 * NY6*NZ6] + 5 ] ;
        RelationXi_2[6] = z_global[NZ6*(j_calc+1) + cell_z[NZ6*(j_store)+k + 2 * NY6*NZ6] + 6 ] ;
    }else if(index_z2 < 6){
        for(int i = 3 ; i <=6 ; i++) RelationXi_2[i] = 0.0 ;
        for(int i = 0 ; i <=2 ; i++) {
            RelationXi_2[i] = z_global[NZ6*(j_calc+1) + cell_z[NZ6*(j_store)+k + 2 * NY6*NZ6] + i ] ;
        }}
    else if(index_z2 > NZ6-7){
        for(int i = 0 ; i <=3 ; i++) RelationXi_2[i] = 0.0 ;
        for(int i = 4 ; i <=6 ; i++) {
            RelationXi_2[i] = z_global[NZ6*(j_calc+1) + cell_z[NZ6*(j_store)+k + 2 * NY6*NZ6] + i ] ;
        }}
    else{
        for(int i = 0 ; i <7 ; i++){
        RelationXi_2[i] = z_global[NZ6*(j_calc+1) + cell_z[NZ6*(j_store)+k + 2 * NY6*NZ6] + i ] ;
    }}
    for(int i = 0 ; i <=6 ; i++){RelationXi_2[i] = (RelationXi_2[i] - HillFunction( y_global[j_calc+1] ) - 0.5*minSize) / L2 ;} //轉換為無因次化Z座標
    GetParameter_6th2( XiPara_h, pos_xi2 , RelationXi_2 , 2 , j_store , k , index_z2 ); //配置第三套適應性內插權重
}

//降階版本 
void GetIntrplParameter_Y() {
    for( int j = 3; j < NY6-3; j++ ){
        GetParameter_2nd( YPara0_h, y_global[j]-minSize, y_global, j, j-1 );//使用在F1，作為y方向預配置連乘權重一維連續記憶體 
        GetParameter_2nd( YPara2_h, y_global[j]+minSize, y_global, j, j-1 );//使用在F3，作為y方向預配置連乘權重一維連續記憶體 
    }
}


void GetIntrplParameter_Xi() {
    // 注意：7-point stencil 需要訪問 cellZ 到 cellZ+6
    // 為避免越界，對於 k >= NZ6-16 的點不進行插值參數計算
    // 這些點在 evolution.h 中會使用簡單的 streaming 替代插值
    for( int j = 3; j < NY6-3; j++ ){
        for( int k = 4; k < NZ6-4;  k++ ){
            
            //寫進去 CellZ_F1(3個起點): 
            RelationXi(z_global[j*NZ6+k]         , j , k ,  CellZ_F1 , nonuni_a);
            //寫進去 CellZ_F2 :
            RelationXi(z_global[j*NZ6+k]-minSize , j , k ,  CellZ_F2 , nonuni_a);
             //寫進去 CellZ_F4 :
            RelationXi(z_global[j*NZ6+k]+minSize , j , k ,  CellZ_F4 , nonuni_a);
            // F1 (+Y): 從 (y-Δ, z) 來，z 方向無偏移
            GetXiParameter( XiParaF1_h,  z_global[j*NZ6+k],         y_global[j]-minSize, j , k,  CellZ_F1);
            // F2 (0,+Z): 從 (y, z-Δ) 來
            GetXiParameter( XiParaF2_h,  z_global[j*NZ6+k]-minSize ,y_global[j],         j , k,  CellZ_F2);
            // F3 (-Y,0): 從 (y+Δ, z) 來，z 方向無偏移
            GetXiParameter( XiParaF3_h,  z_global[j*NZ6+k],         y_global[j]+minSize, j , k,  CellZ_F1);
            // F4 (0,-Z): 從 (y, z+Δ) 來
            GetXiParameter( XiParaF4_h,  z_global[j*NZ6+k]+minSize, y_global[j],         j , k,  CellZ_F4);
            // F5 (+Y,+Z): 從 (y-Δ, z-Δ) 來
            GetXiParameter( XiParaF5_h,  z_global[j*NZ6+k]-minSize, y_global[j]-minSize, j , k,  CellZ_F2);
            // F6 (-Y,+Z): 從 (y+Δ, z-Δ) 來
            GetXiParameter( XiParaF6_h,  z_global[j*NZ6+k]-minSize, y_global[j]+minSize, j , k,  CellZ_F2);
            // F7 (-Y,-Z): 從 (y+Δ, z+Δ) 來
            GetXiParameter( XiParaF7_h,  z_global[j*NZ6+k]+minSize, y_global[j]+minSize, j , k,  CellZ_F4);
            // F8 (+Y,-Z): 從 (y-Δ, z+Δ) 來
            GetXiParameter( XiParaF8_h,  z_global[j*NZ6+k]+minSize, y_global[j]-minSize, j , k,  CellZ_F4);
}}
}

void BFLInitialization(double *Q1_h, double *Q3_h, double *Q5_h, double *Q6_h) {
    //此函數是為計算邊界計算點的q值以及邊界計算點的預配置連乘權重一維連續記憶體
    //q: 計算點到壁面的無因次距離
    //delta: BFL 反彈點相對於計算點的偏移量
    //Parameter: 預配置 Lagrange 插值權重

    for(int j = 3; j < NY6-3; j++){
        for(int k = 3; k < NZ6-3; k++){

            //F1 (+Y方向): 左丘邊界
            // 當粒子從 -Y 方向來，可能撞到左丘
            if(IsLeftHill_Boundary_yPlus(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F1的邊界計算點
                double q1 = Left_q_yPlus(y_global[j], z_global[j*NZ6+k]);
                double delta1 = minSize * (1.0 - 2.0*q1);
                // BFL 反彈點在 +Y 方向: y + delta, z 不變
                //GetParameter_6th(YBFLParaF3_h, y_global[j]+delta1, y_global , j , j-3);//F3代表的意思是此權重陣列配合的對象是F3 利用F3來更新F1
                //降階版本
                //GetParameter_2nd(YBFLParaF3_h, y_global[j]+delta1, y_global , j , j-1);//F3代表的意思是此權重陣列配合的對象是F3 利用F3來更新F1
                //GetXiParameter(XiBFLParaF3_h, z_global[j*NZ6+k], y_global[j]+delta1,  j*NZ6+k , j,  k) ;
                //目前先採用 線性內插
                Q1_h[j*NZ6+k] = q1;
            }

            //F3 (-Y方向): 右丘邊界
            // 當粒子從 +Y 方向來，可能撞到右丘
            if(IsRightHill_Boundary_yMinus(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F3的邊界計算點
                double q3 = Right_q_yMinus(y_global[j], z_global[j*NZ6+k]);
                double delta3 = minSize * (1.0 - 2.0*q3);
                // BFL 反彈點在 -Y 方向: y - delta, z 不變
                //GetParameter_6th(YBFLParaF1_h, y_global[j]-delta3, y_global, j , j-3);//F1代表的意思是此權重陣列配合的對象是F1 利用F1來更新F3
                //降階版本
                //GetParameter_2nd(YBFLParaF1_h, y_global[j]-delta3, y_global, j , j-1);//F1代表的意思是此權重陣列配合的對象是F1 利用F1來更新F3
                //GetXiParameter(XiBFLParaF1_h, z_global[j*NZ6+k], y_global[j]-delta3,  j*NZ6+k , j,  k);
                Q3_h[j*NZ6+k] = q3;
            }

            //更新F5 (+Y+Z方向, 45度): 左丘邊界
            // 當粒子從 (-Y,-Z) 方向來，可能撞到左丘
            if(IsLeftHill_Boundary_Diagonal45(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F5的邊界計算點
                double q5 = Left_q_Diagonal45(y_global[j], z_global[j*NZ6+k]);
                double delta5 = minSize * (1.0 - 2.0*q5) / sqrt(2.0);
                // BFL 反彈點在 (+Y,+Z) 方向: y + delta, z + delta
                //GetParameter_6th(YBFLParaF7_h, y_global[j]+delta5, y_global, j, j-3);
                //降階版本
                //GetParameter_2nd(YBFLParaF7_h, y_global[j]+delta5, y_global, j, j-1);//F7代表的意思是此權重陣列配合的對象是F7 利用F7來更新F5
                //GetXiParameter(XiBFLParaF7_h, z_global[j*NZ6+k]+delta5, y_global[j]+delta5,  j*NZ6+k , j,  k);//F7代表的意思是此權重陣列配合的對象是F7 利用F7來更新F5
                Q5_h[j*NZ6+k] = q5;
            }

            //更新F6 (-Y+Z方向, 135度): 右丘邊界
            // 當粒子從 (+Y,-Z) 方向來，可能撞到右丘
            if(IsRightHill_Boundary_Diagonal135(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F6的邊界計算點
                double q6 = Right_q_Diagonal135(y_global[j], z_global[j*NZ6+k]);
                double delta6 = minSize * (1.0 - 2.0*q6) / std::sqrt(2.0);
                // BFL 反彈點在 (-Y,+Z) 方向: y - delta, z + delta
                //GetParameter_6th(YBFLParaF8_h, y_global[j]-delta6, y_global, j, j-3);//F8代表的意思是此權重陣列配合的對象是F8 利用F8來更新F6
                //降階版本
                //GetParameter_2nd(YBFLParaF8_h, y_global[j]-delta6, y_global, j, j-1);//F8代表的意思是此權重陣列配合的對象是F8 利用F8來更新F6
                //GetXiParameter(XiBFLParaF8_h, z_global[j*NZ6+k]+delta6, y_global[j]-delta6,  j*NZ6+k , j,  k);
                Q6_h[j*NZ6+k] = q6;
            }
        }
    }
}
#endif 
