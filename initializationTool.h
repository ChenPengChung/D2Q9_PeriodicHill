#ifndef INITIALIZATIONTOOL_FILE
#define INITIALIZATIONTOOL_FILE
#include <cmath>
#include <variables.h>
#include <model.h>
#define HillHalfWidth 9.0*(54.0/28.0) //半邊山丘寬度

//1.HillFunction 左半丘反函數 (數值法)
// 給定 z，回傳對應的 y
double HillFunction_Inverse_Left(double z) {
    double y_low = 0.0, y_high = HillHalfWidth;
    double y_mid;

    // 左半丘是遞減函數: y 越大，HillFunction(y) 越小
    while (y_high - y_low > 1e-12) {
        y_mid = (y_low + y_high) / 2.0;
        if (HillFunction(y_mid) > z) {
            y_low = y_mid;
        } else {
            y_high = y_mid;
        }
    }
    return y_mid;
}

//1-2.HillFunction 右半丘反函數 (數值法)
// 給定 z，回傳對應的 y
double HillFunction_Inverse_Right(double z) {
    double y_low = LY - HillHalfWidth, y_high = LY;
    double y_mid;

    // 右半丘是遞增函數: y 越大，HillFunction(y) 越大
    while (y_high - y_low > 1e-12) {
        y_mid = (y_low + y_high) / 2.0;
        if (HillFunction(y_mid) < z) {
            y_low = y_mid;
        } else {
            y_high = y_mid;
        }
    }
    return y_mid;
}

//2.
/*L為計算點所包圍的總長度*/
/*MinSize為粒子晶格大小*/
/*a為非均勻網格伸縮參數*/
/*j為計算點編號*/
/*N為計算點為節點之網格數量*/
#define tanhFunction( L , MinSize , a, j, N)\
(           \
    L/2.0 + MinSize/2.0 + ((L/2.0)/a)*tanh((-1.0+2.0*(double)(j)/(double)(N))/2.0*log((1.0+a)/(1.0-a)))\
)
//3.
//伸縮參數a由下函數決定
double GetNonuniParameter() {
    double total = LZ - HillFunction( 0.0 ) - minSize;
    double a_temp[2] = {0.1, 1.0};
    double a_mid;

    double x_temp[2], dx;
    do{
        a_mid = (a_temp[0]+a_temp[1]) / 2.0;
        //dx = Z方向y = 0.0 非均勻網格下的最小間距 (從y=  0,0取)
        //minSize的定義為以LZ-(y = 0.0下的山坡高度)為總長度均勻切割下的網格大小*0.6
        //判斷標準 最小的網格必須要大於minSize
        x_temp[0] = tanhFunction(total, minSize, a_mid, 0, (NZ6-7));
        x_temp[1] = tanhFunction(total, minSize, a_mid, 1, (NZ6-7));
        dx = x_temp[1] - x_temp[0];
        if( dx - minSize >= 0.0 ){
            a_temp[0] = a_mid;
        } else {
            a_temp[1] = a_mid;
        }
    } while (fabs( dx - minSize) > 1e-14 );
    return a_mid;
}
//4.Lagrage Interpolation 公式
double Lagrange_6th(double pos , double x_i , double x1 , double x2 , double x3 , double x4 , double x5 , double x6){
    double Lagrange = (pos - x1)/(x_i - x1)*(pos - x2)/(x_i - x2)*(pos - x3)/(x_i - x3)*(pos - x4)/(x_i - x4)*(pos - x5)/(x_i - x5)*(pos - x6)/(x_i - x6);
    return Lagrange; 
}

//5.產生預配置一維連乘權重陣列
void GetParameter_6th( double** Para , double Position , double* phy , int now , int start){
    Para[0][now] = Lagrange_6th( Position , *(phy+start) ,   *(phy+start+1) ,*(phy+start+2) ,*(phy+start+3) ,*(phy+start+4) ,*(phy+start+5) ,*(phy+start+6) ) ; 
    Para[1][now] = Lagrange_6th( Position , *(phy+start+1) , *(phy+start) ,*(phy+start+2) ,*(phy+start+3) ,*(phy+start+4) ,*(phy+start+5) ,*(phy+start+6) ) ;     
    Para[2][now] = Lagrange_6th( Position , *(phy+start+2) , *(phy+start) ,*(phy+start+1) ,*(phy+start+3) ,*(phy+start+4) ,*(phy+start+5) ,*(phy+start+6) ) ;     
    Para[3][now] = Lagrange_6th( Position , *(phy+start+3) , *(phy+start) ,*(phy+start+1) ,*(phy+start+2) ,*(phy+start+4) ,*(phy+start+5) ,*(phy+start+6) ) ;     
    Para[4][now] = Lagrange_6th( Position , *(phy+start+4) , *(phy+start) ,*(phy+start+1) ,*(phy+start+2) ,*(phy+start+3) ,*(phy+start+5) ,*(phy+start+6) ) ;     
    Para[5][now] = Lagrange_6th( Position , *(phy+start+5) , *(phy+start) ,*(phy+start+1) ,*(phy+start+2) ,*(phy+start+3) ,*(phy+start+4) ,*(phy+start+6) ) ;     
    Para[6][now] = Lagrange_6th( Position , *(phy+start+6) , *(phy+start) ,*(phy+start+1) ,*(phy+start+2) ,*(phy+start+3) ,*(phy+start+4) ,*(phy+start+5) ) ;     
} 
//考察任意物理間計算點是否為曲面邊界的邊界點，分別往y方向 與 45度方想考察．
//由於Z方向必定為Half-Way Bounce-Back邊界條件，且已經確定此方向之邊界計算點為何？
//如何決定邊界計算點？1.以曲面邊界為本體，2.往+y , +z , 45度斜角分別考察半格移動距離：minSize0.5 , minSize0.5 , sqrt(2)*0.5*minSize，就可以知道對於曲面在方向的邊界幾算點，所以可以分累，不同方向的物理間邊界計算點
//6.是否為左山丘的+y方向物理空間邊界計算點(作為之後套用BFL邊界條件的條件)
//考察範圍： 沿著左丘表面往+y方係延伸一個minSize
bool IsLeftHill_Boundary_yPlus(double y , double z){
    //第一步先初步篩選範圍：左下角方形區域為pos_y : [0,54/28.0] ; pos_z : [0,1]
    if( y < 0.0 || y > HillHalfWidth || z < HillFunction(y) || z > 1.0 ){
        return false; //第一步篩選未通過
    }else{//第二步判斷，往左移動會撞倒代表現在高度在移動後座標點之山丘下方
        if( z <= HillFunction(y- minSize) ){
            return true;
        }else{
            return false;
        }
    }
}

//7.是否為右山丘的-y方向物理空間邊界計算點(作為之後套用BFL邊界條件的條件)
//考察範圍： 沿著右丘表面往-y方向延伸一個minSize
bool IsRightHill_Boundary_yMinus(double y , double z){
    //第一步先初步篩選範圍：右下角方形區域為pos_y : [LY-HillHalfWidth, LY] ; pos_z : [0,1]
    if( y < (LY - HillHalfWidth) || y > LY || z < HillFunction(y) || z > 1.0 ){
        return false; //第一步篩選未通過
    }else{//第二步判斷，往右移動會撞倒代表現在高度在移動後座標點之山丘下方 
        if( z <= HillFunction(y + minSize) ){
            return true;
        }else{
            return false;
        }
    }
}
//8.是否為左邊山丘的45度斜向物理空間計算點(作為之後套用BFL邊界條件的判斷)
// 判斷往 (-Y, -Z) 方向移動後是否會進入山丘內部
bool IsLeftHill_Boundary_Diagonal45(double y , double z){
    //第一步先初步篩選範圍：左下角方形區域為pos_y : [0,54/28.0] ; pos_z : [0,1]
    if( y < 0.0 || y > HillHalfWidth || z < HillFunction(y) || z > 1.0 ){
        return false; //第一步篩選未通過
    }else{
        // 移動後的位置
        double y_new = y - minSize;
        double z_new = z - minSize;
        // 判斷移動後是否進入山丘內部 (z_new <= 曲面高度)
        if( z_new <= HillFunction(y_new) ){
            return true;
        }else{
            return false;
        }
    }
}

//9.是否為右邊山丘的135度斜向物理空間計算點(作為之後套用BFL邊界條件的判斷)
// 判斷往 (+Y, -Z) 方向移動後是否會進入山丘內部
bool IsRightHill_Boundary_Diagonal135(double y , double z){
    //第一步先初步篩選範圍：右下角方形區域為pos_y : [LY-HillHalfWidth, LY] ; pos_z : [0,1]
    if( y < (LY - HillHalfWidth) || y > LY || z < HillFunction(y) || z > 1.0 ){
        return false; //第一步篩選未通過
    }else{
        // 移動後的位置
        double y_new = y + minSize;
        double z_new = z - minSize;
        // 判斷移動後是否進入山丘內部 (z_new <= 曲面高度)
        if( z_new <= HillFunction(y_new) ){
            return true;
        }else{
            return false;
        }
    }
}
//10.接續6.bool判斷式，求施加BFL條件之實際與壁面之距離
double Left_q_yPlus(double y , double z){
    if(!IsLeftHill_Boundary_yPlus(y, z)){
        return -1.0; //非邊界點返回-1
    }//返回負值作為錯誤之邊界條件實施判據，因為距離沒有負值
    // 計算 q 值：計算點到曲面的 Y 方向距離 / 格子大小
    return fabs(y - HillFunction_Inverse_Left(z)) / minSize;
}
//11.接續7.bool判斷式，求施加BFL條件之實際與壁面之距離
double Right_q_yMinus(double y , double z){
    if(!IsRightHill_Boundary_yMinus(y, z)){
        return -1.0; //非邊界點返回-1
    }
    // 計算 q 值：計算點到曲面的 Y 方向距離 / 格子大小
    return fabs(y - HillFunction_Inverse_Right(z)) / minSize;
}
//12.接續 8.bool判斷式，求施加BFL條件之右上斜方向邊界計算點與壁面之距離
double Left_q_Diagonal45(double y , double z){
    if(!IsLeftHill_Boundary_Diagonal45(y , z)){
        return -1.0 ; //非邊界點返回-1
    }
    double y_target ;
    //利用區間搜尋法尋找y_target ; 
    double y_temp[2] = {y-minSize ,y} ; 
    //y[0] = y_minSize ; y[1] = y ; 
    double y_middle = 0.0 ;
    double z_middle = 0.0 ; //在45度射線上相對應y_middle的z座標
    double Length_z = 0.0 ;
    do{
        y_middle = (y_temp[0] + y_temp[1]) / 2.0;
        z_middle = y_middle + (z-y) ; 
        Length_z = fabs(z_middle - HillFunction(y_middle));
        if(HillFunction(y_middle) > z_middle){
            //將y_middle往右移動
            y_temp[0] = y_middle;
        }else{
            //將y_middle往左移動
            y_temp[1] = y_middle;
            
        }
    }while(Length_z >= 1e-12) ;
    y_target = y_middle ;
    return fabs(y - y_target) / minSize ;
}
//13.接續 9.bool判斷式，求施加BFL條件之左上斜方向(135度)邊界計算點與壁面之距離
double Right_q_Diagonal135(double y , double z){
    if(!IsRightHill_Boundary_Diagonal135(y , z)){
        return -1.0 ; //非邊界點返回-1
    }
    double y_target ;
    //利用區間搜尋法尋找y_target ;
    double y_temp[2] = {y , y+minSize} ;
    //y[0] = y ; y[1] = y+minSize ;
    double y_middle = 0.0 ;
    double z_middle = 0.0 ; //在135度射線上相對應y_middle的z座標
    double Length_z = 0.0 ;
    do{
        y_middle = (y_temp[0] + y_temp[1]) / 2.0;
        z_middle = -y_middle + (z+y) ; // 135度射線：z = -y + (z0+y0)
        Length_z = fabs(z_middle - HillFunction(y_middle));
        if(HillFunction(y_middle) > z_middle){
            //將y_middle往左移動
            y_temp[1] = y_middle;
        }else{
            //將y_middle往右移動
            y_temp[0] = y_middle;
        }
    }while(Length_z >= 1e-12) ;
    y_target = y_middle ;
    return fabs(y - y_target) / minSize ;
}
#endif