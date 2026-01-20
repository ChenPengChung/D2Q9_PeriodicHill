#ifndef INITIALIZATION_FILE
#define INITIALIZATION_FILE
#include <iostream>
#include <cstdlib>
#include <cmath>
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
        //離散化宏觀外立場的初始化initilaoization of the discrete macroscopic force term
        Force[0] =  (8.0*niu*Uref)/(LZ*LZ)*5.0; //0.0001;
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
    double** XiPara_h,    double pos_z,       double pos_y,
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
        GetParameter_6th( YPara0_h, y_global[i]-minSize, y_global, i, i-3 );//使用在F1，作為y方向預配置連乘權重一維連續記憶體 
        GetParameter_6th( YPara2_h, y_global[i]+minSize, y_global, i, i-3 );//使用在F3，作為y方向預配置連乘權重一維連續記憶體 
    }
}


void GetIntrplParameter_Xi() {
    for( int j = 3; j < NY6-3; j++ ){
        for( int k = 3; k < NZ6-3;  k++ ){
            //為什麼需要二維記憶體配置，因為在y_z平面上，
            //每一個物理空間計算點的相應非物理空間計算點的Z方向預配置連乘權重一維連續記憶體都不一樣
            //，會受到座標比例影響內插計算結果。
            //
            // D2Q9 速度方向: F1(+Y), F2(+Z), F3(-Y), F4(-Z), F5(+Y+Z), F6(-Y+Z), F7(-Y-Z), F8(+Y-Z)
            // Streaming: f_i(x,t+dt) = f_i(x - e_i*dt, t) → 從 (y - ey*Δ, z - ez*Δ) 位置取值
            //
            // F1 (+Y,0): 從 (y-Δ, z) 來
            GetXiParameter( XiParaF1_h,  z_global[j*NZ6+k],         y_global[j]-minSize, xi_h, j*NZ6+k, k );
            // F2 (0,+Z): 從 (y, z-Δ) 來
            GetXiParameter( XiParaF2_h,  z_global[j*NZ6+k]-minSize , y_global[j],         xi_h, j*NZ6+k, k );
            // F3 (-Y,0): 從 (y+Δ, z) 來
            GetXiParameter( XiParaF3_h,  z_global[j*NZ6+k],         y_global[j]+minSize, xi_h, j*NZ6+k, k );
            // F4 (0,-Z): 從 (y, z+Δ) 來
            GetXiParameter( XiParaF4_h,  z_global[j*NZ6+k]+minSize, y_global[j],         xi_h, j*NZ6+k, k );
            // F5 (+Y,+Z): 從 (y-Δ, z-Δ) 來
            GetXiParameter( XiParaF5_h,  z_global[j*NZ6+k]-minSize, y_global[j]-minSize, xi_h, j*NZ6+k, k );
            // F6 (-Y,+Z): 從 (y+Δ, z-Δ) 來
            GetXiParameter( XiParaF6_h,  z_global[j*NZ6+k]-minSize, y_global[j]+minSize, xi_h, j*NZ6+k, k );
            // F7 (-Y,-Z): 從 (y+Δ, z+Δ) 來
            GetXiParameter( XiParaF7_h,  z_global[j*NZ6+k]+minSize, y_global[j]+minSize, xi_h, j*NZ6+k, k );
            // F8 (+Y,-Z): 從 (y-Δ, z+Δ) 來
            GetXiParameter( XiParaF8_h,  z_global[j*NZ6+k]+minSize, y_global[j]-minSize, xi_h, j*NZ6+k, k );
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
                GetParameter_6th(YBFLParaF3_h, y_global[j]+delta1, y_global , j , j-3);//F3代表的意思是此權重陣列配合的對象是F3 利用F3來更新F1 
                GetXiParameter(XiBFLParaF3_h, z_global[j*NZ6+k], y_global[j]+delta1, xi_h , j*NZ6+k, k);
                Q1_h[j*NZ6+k] = q1;
            }

            //F3 (-Y方向): 右丘邊界
            // 當粒子從 +Y 方向來，可能撞到右丘
            if(IsRightHill_Boundary_yMinus(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F3的邊界計算點 
                double q3 = Right_q_yMinus(y_global[j], z_global[j*NZ6+k]);
                double delta3 = minSize * (1.0 - 2.0*q3);
                // BFL 反彈點在 -Y 方向: y - delta, z 不變
                GetParameter_6th(YBFLParaF1_h, y_global[j]-delta3, y_global, j , j-3);//F1代表的意思是此權重陣列配合的對象是F1 利用F1來更新F3
                GetXiParameter(XiBFLParaF1_h, z_global[j*NZ6+k], y_global[j]-delta3, xi_h, j*NZ6+k, k);
                Q3_h[j*NZ6+k] = q3;
            }

            //更新F5 (+Y+Z方向, 45度): 左丘邊界
            // 當粒子從 (-Y,-Z) 方向來，可能撞到左丘
            if(IsLeftHill_Boundary_Diagonal45(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F5的邊界計算點
                double q5 = Left_q_Diagonal45(y_global[j], z_global[j*NZ6+k]);
                double delta5 = minSize * (1.0 - 2.0*q5) / sqrt(2.0);
                // BFL 反彈點在 (+Y,+Z) 方向: y + delta, z + delta
                GetParameter_6th(YBFLParaF7_h, y_global[j]+delta5, y_global, j, j-3);
                GetXiParameter(XiBFLParaF7_h, z_global[j*NZ6+k]+delta5, y_global[j]+delta5, xi_h, j*NZ6+k, k);//F7代表的意思是此權重陣列配合的對象是F7 利用F7來更新F5
                Q5_h[j*NZ6+k] = q5;
            }

            //更新F6 (-Y+Z方向, 135度): 右丘邊界
            // 當粒子從 (+Y,-Z) 方向來，可能撞到右丘
            if(IsRightHill_Boundary_Diagonal135(y_global[j], z_global[j*NZ6+k])){//尋找專屬於F6的邊界計算點
                double q6 = Right_q_Diagonal135(y_global[j], z_global[j*NZ6+k]);
                double delta6 = minSize * (1.0 - 2.0*q6) / std::sqrt(2.0);
                // BFL 反彈點在 (-Y,+Z) 方向: y - delta, z + delta
                GetParameter_6th(YBFLParaF8_h, y_global[j]-delta6, y_global, j, j-3);//F8代表的意思是此權重陣列配合的對象是F8 利用F8來更新F6
                GetXiParameter(XiBFLParaF8_h, z_global[j*NZ6+k]+delta6, y_global[j]-delta6, xi_h, j*NZ6+k, k);
                Q6_h[j*NZ6+k] = q6;
            }
        }
    }
}



#endif 
