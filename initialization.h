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
        //離散化\ 宏觀\ 外立\ 場的初始化initilaoization of the discrete macroscopic force term
        Force[0] =  (8.0*niu*Uref)/(LZ*LZ)*1.0; //降低外力係數，原本是 5.0
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



void GetXiParameter(double* XiPara_h[7], double pos_z, double pos_y, int index_xi, int j, int k, int* CellZ_out)
{   //index_xi , j.  k:為當前座標點 
    //CellZ_out[idx_xi] : 為當前座標點的內插成員起始點Z座標
    //越界防呆
    if(j<3) j = j + NY6-7 ;
    if(j>NY6-4) j = j - (NY6-7) ;
    //Z方向係數編號系統(第二格)  = (j+relation)*NZ6 + k  ; relation:1~7
    double L = LZ - HillFunction(pos_y) - minSize;
    //每一個y位置而言每一個計算間都不一樣
    double pos_xi =  (pos_z - (HillFunction(pos_y)+minSize/2.0)) ;
    const double a = nonuni_a; // 已在 GenerateMesh_Z 中計算
    //求解目標點的編號(映射回均勻系統)
    double j_cont = Inverse_tanh_index( pos_xi , L , minSize , a , (NZ6-7) );

    //Z方向成員起始編號
    //在 ISLBM 中，stencil 需要包圍來源點的 k 索引
    //使用 j_cont-3 作為起點，使得來源點落在 stencil 的中間位置 (位置 3)
    //這樣 stencil 覆蓋約 k_src-3 到 k_src+3，共 7 個點
    //
    //重要修正：cell_z1 應該基於 j_cont（映射後的連續座標）而不是傳入的 k
    //因為不同方向的 Fi 有不同的來源位置偏移：
    //  F2 (+Z): 來源在 z-Δ → j_cont 較小
    //  F4 (-Z): 來源在 z+Δ → j_cont 較大，可能接近或超出邊界
    //
    //需要加上 buffer 區域偏移 (+3) 來得到實際的 k 索引
    int k_source = static_cast<int>(std::round(j_cont)) + 3;  // 來源點的 k 索引
    int cell_z1 = k_source - 3;
    if(cell_z1 < 0) cell_z1 = 0;
    if(cell_z1 > NZ6-7) cell_z1 = NZ6-7;
    //預先計算 內插成員起始點座標編號Z座標 
    if(CellZ_out != nullptr) {
        CellZ_out[index_xi] = cell_z1;//儲存Z方向的內插成員座標起始點
        } 
    //y方向為七點版本
    /*/給我一個double 給出相對應七個內插的無因次化座標
    double RelationXi_0[7] ; //(j-3)
    double LT = LZ - HillFunction(y_global[j-3]) - minSize;
    double pos_z2 =  tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_0);//寫入第一套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z2 , RelationXi_0 , 0 , index_xi); //XiPara_h[0][index_xi + 0*NY6*NZ6] ~ XiPara_h[6][index_xi + 0*NY6*NZ6]
    
    double RelationXi_1[7] ; //(j-2)
    LT = LZ - HillFunction(y_global[j-2]) - minSize;
    pos_z2 = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_1);//寫入第二套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z2 , RelationXi_1 , 1 , index_xi); //XiPara_h[0][index_xi + 1*NY6*NZ6] ~ XiPara_h[6][index_xi + 1*NY6*NZ6]
    
    double RelationXi_2[7] ; //(j-1)
    LT = LZ - HillFunction(y_global[j-1]) - minSize;
    pos_z2 = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_2);//寫入第三套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z2 , RelationXi_2 , 2 , index_xi); //XiPara_h[0][index_xi + 2*NY6*NZ6] ~ XiPara_h[6][index_xi + 2*NY6*NZ6]
    
    double RelationXi_3[7] ; //(j)
    LT = LZ - HillFunction(y_global[j]) - minSize;
    pos_z2 = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_3);//寫入第四套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z2 , RelationXi_3 , 3 , index_xi); //XiPara_h[0][index_xi + 3*NY6*NZ6] ~ XiPara_h[6][index_xi + 3*NY6*NZ6]

    double RelationXi_4[7] ; //(j+1)
    LT = LZ - HillFunction(y_global[j+1]) - minSize;
    pos_z2 = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_4);//寫入第五套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z2, RelationXi_4 , 4 , index_xi); //XiPara_h[0][index_xi + 4*NY6*NZ6] ~ XiPara_h[6][index_xi + 4*NY6*NZ6]

    double RelationXi_5[7] ; //(j+2)
    LT = LZ - HillFunction(y_global[j+2]) - minSize;
    pos_z2 = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_5);//寫入第六套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z2 , RelationXi_5 , 5 , index_xi); //XiPara_h[0][index_xi + 5*NY6*NZ6] ~ XiPara_h[6][index_xi + 5*NY6*NZ6]
   
    double RelationXi_6[7] ; //(j+3)
    LT = LZ - HillFunction(y_global[j+3]) - minSize;
    pos_z2 = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_6);//寫入第七套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z2 , RelationXi_6 , 6 , index_xi); //XiPara_h[0][index_xi + 6*NY6*NZ6] ~ XiPara_h[6][index_xi + 6*NY6*NZ6]*/
    
    //y方向為三點版本
    double RelationXi_0[7] ; //(j-1)
    double LT = LZ - HillFunction(y_global[j-1]) - minSize;
    double pos_z2 =  tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_0);//寫入第一套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z2 , RelationXi_0 , 0 , index_xi); //XiPara_h[0][index_xi + 0*NY6*NZ6] ~ XiPara_h[6][index_xi + 0*NY6*NZ6]
    
    double RelationXi_1[7] ; //(j-0)
    LT = LZ - HillFunction(y_global[j-0]) - minSize;
    pos_z2 = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_1);//寫入第二套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z2 , RelationXi_1 , 1 , index_xi); //XiPara_h[0][index_xi + 1*NY6*NZ6] ~ XiPara_h[6][index_xi + 1*NY6*NZ6]
    
    double RelationXi_2[7] ; //(j+1)
    LT = LZ - HillFunction(y_global[j+1]) - minSize;
    pos_z2 = tanhFunction( LT, minSize, a, j_cont , (NZ6-7) ) - minSize/2.0;
    RelationXi(j_cont, LT , minSize , a , (NZ6-7) , RelationXi_2);//寫入第三套垂直方向無因次化Z陣列
    GetParameter_6th2( XiPara_h, pos_z2 , RelationXi_2 , 2 , index_xi); //XiPara_h[0][index_xi + 2*NY6*NZ6] ~ XiPara_h[6][index_xi + 2*NY6*NZ6]
}

//降階版本 
void GetIntrplParameter_Y() {
    for( int j = 3; j < NY6-3; j++ ){
        GetParameter_2nd( YPara0_h, y_global[j]-minSize, y_global, j, j-1 );//使用在F1，作為y方向預配置連乘權重一維連續記憶體 
        GetParameter_2nd( YPara2_h, y_global[j]+minSize, y_global, j, j-1 );//使用在F3，作為y方向預配置連乘權重一維連續記憶體 
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
            // 新增：CellZ_F* 陣列存儲預計算的 stencil 起點，確保與 RelationXi 計算一致
            //
            // F1 (+Y,0): 從 (y-Δ, z) 來，z 方向無偏移
            GetXiParameter( XiParaF1_h,  z_global[j*NZ6+k],         y_global[j]-minSize, j*NZ6+k , j,  k, CellZ_F1);
            // F2 (0,+Z): 從 (y, z-Δ) 來
            GetXiParameter( XiParaF2_h,  z_global[j*NZ6+k]-minSize , y_global[j],        j*NZ6+k , j,  k, CellZ_F2);
            // F3 (-Y,0): 從 (y+Δ, z) 來，z 方向無偏移
            GetXiParameter( XiParaF3_h,  z_global[j*NZ6+k],         y_global[j]+minSize, j*NZ6+k , j,  k, CellZ_F3);
            // F4 (0,-Z): 從 (y, z+Δ) 來
            GetXiParameter( XiParaF4_h,  z_global[j*NZ6+k]+minSize, y_global[j],         j*NZ6+k , j,  k, CellZ_F4);
            // F5 (+Y,+Z): 從 (y-Δ, z-Δ) 來
            GetXiParameter( XiParaF5_h,  z_global[j*NZ6+k]-minSize, y_global[j]-minSize, j*NZ6+k , j,  k, CellZ_F5);
            // F6 (-Y,+Z): 從 (y+Δ, z-Δ) 來
            GetXiParameter( XiParaF6_h,  z_global[j*NZ6+k]-minSize, y_global[j]+minSize, j*NZ6+k , j,  k, CellZ_F6);
            // F7 (-Y,-Z): 從 (y+Δ, z+Δ) 來
            GetXiParameter( XiParaF7_h,  z_global[j*NZ6+k]+minSize, y_global[j]+minSize, j*NZ6+k , j,  k, CellZ_F7);
            // F8 (+Y,-Z): 從 (y-Δ, z+Δ) 來
            GetXiParameter( XiParaF8_h,  z_global[j*NZ6+k]+minSize, y_global[j]-minSize, j*NZ6+k , j,  k, CellZ_F8);
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
                GetParameter_2nd(YBFLParaF3_h, y_global[j]+delta1, y_global , j , j-1);//F3代表的意思是此權重陣列配合的對象是F3 利用F3來更新F1
                GetXiParameter(XiBFLParaF3_h, z_global[j*NZ6+k], y_global[j]+delta1,  j*NZ6+k , j,  k, nullptr) ;
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
                GetParameter_2nd(YBFLParaF1_h, y_global[j]-delta3, y_global, j , j-1);//F1代表的意思是此權重陣列配合的對象是F1 利用F1來更新F3
                GetXiParameter(XiBFLParaF1_h, z_global[j*NZ6+k], y_global[j]-delta3,  j*NZ6+k , j,  k, nullptr);
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
                GetParameter_2nd(YBFLParaF7_h, y_global[j]+delta5, y_global, j, j-1);//F7代表的意思是此權重陣列配合的對象是F7 利用F7來更新F5
                GetXiParameter(XiBFLParaF7_h, z_global[j*NZ6+k]+delta5, y_global[j]+delta5,  j*NZ6+k , j,  k, nullptr);//F7代表的意思是此權重陣列配合的對象是F7 利用F7來更新F5
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
                GetParameter_2nd(YBFLParaF8_h, y_global[j]-delta6, y_global, j, j-1);//F8代表的意思是此權重陣列配合的對象是F8 利用F8來更新F6
                GetXiParameter(XiBFLParaF8_h, z_global[j*NZ6+k]+delta6, y_global[j]-delta6,  j*NZ6+k , j,  k, nullptr);
                Q6_h[j*NZ6+k] = q6;
            }
        }
    }
}
#endif 
