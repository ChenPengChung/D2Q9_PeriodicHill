#ifndef INITIALIZATION_FILE
#define INITIALIZATION_FILE
#include <iostream>
using namespace std;
#include "initializationTool.h"

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
        //離散化宏觀外立場的初始化initilaoization of the discrete macroscopic force term 
        Force[0] =  (8.0*niu*Uref)/(LZ*LZ)*5.0; //0.0001;
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
    double a = GetNonuniParameter(); //計算最合適非均勻參數 
    
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
void GetXiParameter(
    double *XiPara_h[7],    double pos_z,       double pos_y,
    double *Pos_xi,         int IdxToStore,     int k  ) //IdxToStore = now ; k = start 
{
    double L = LZ - HillFunction(pos_y) - minSize;
    //每一個y位置而言每一個計算間都不一樣 
    double pos_xi = LXi * (pos_z - (HillFunction(pos_y)+minSize/2.0)) / L;
    //在Lagrange內插過程中，本應該是無因次化參數的參與，在這邊曲線座標系的處理，不再使用絕對座標
    if( k >= 3 && k <= 6 ){
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, 3 );
    } else if ( k >= NZ6-7 && k <= NZ6-4 ) {
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, NZ6-10 );
    } else {
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, k-3 );
    }    
}

void GetIntrplParameter_Y() {
    for( int i = 3; i < NY6-3; i++ ){
        GetParameter_6th( YPara0_h, y_h[i]-minSize, y_h, i, i-3 );//使用在F1，y方向連乘權重一為
        GetParameter_6th( YPara2_h, y_h[i]+minSize, y_h, i, i-3 );
    }
}


#endif 