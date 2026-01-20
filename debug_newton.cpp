//=============================================================================
// debug_newton.cpp - Debug Newton-Raphson convergence
//=============================================================================
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "variables.h"
#include "model.h"

// Inline the tanhFunction macro for clarity
inline double tanh_func(double L, double MinSize, double a, int j, int N) {
    return L/2.0 + MinSize/2.0 + ((L/2.0)/a)*tanh((-1.0+2.0*(double)(j)/(double)(N))/2.0*log((1.0+a)/(1.0-a)));
}

int main() {
    double total = LZ - HillFunction(0.0) - minSize;
    const int N = NZ6 - 7;
    
    cout << "Debug Newton-Raphson Method:" << endl;
    cout << "total = " << total << endl;
    cout << "N = " << N << endl;
    cout << endl;
    
    // Newton-Raphson with debug output
    double a = 0.5;
    const int max_iter = 50;
    const double tol = 1e-14;
    const double a_min = 0.1, a_max = 1.0;
    
    for(int iter = 0; iter < max_iter; iter++) {
        double x0 = tanh_func(total, minSize, a, 0, N);
        double x1 = tanh_func(total, minSize, a, 1, N);
        double dx = x1 - x0;
        double f = dx - minSize;
        
        cout << "Iteration " << iter << ":" << endl;
        cout << "  a = " << setprecision(15) << a << endl;
        cout << "  dx = " << dx << endl;
        cout << "  f = " << scientific << setprecision(6) << f << endl;
        
        if (fabs(f) < tol) {
            cout << "  CONVERGED!" << endl;
            cout << endl;
            cout << "Final result: a = " << setprecision(15) << a << endl;
            return 0;
        }
        
        // Compute derivative
        const double da = 1e-8;
        double x0_plus = tanh_func(total, minSize, a + da, 0, N);
        double x1_plus = tanh_func(total, minSize, a + da, 1, N);
        double dx_plus = x1_plus - x0_plus;
        double df_da = (dx_plus - dx) / da;
        
        cout << "  df/da = " << df_da << endl;
        
        if (fabs(df_da) < 1e-15) {
            cout << "  ERROR: Derivative too small!" << endl;
            break;
        }
        
        double a_new = a - f / df_da;
        cout << "  a_new = " << setprecision(15) << a_new << endl;
        
        if (a_new <= a_min || a_new >= a_max) {
            cout << "  ERROR: Out of bounds!" << endl;
            break;
        }
        
        a = a_new;
        cout << endl;
    }
    
    cout << "Newton-Raphson failed to converge" << endl;
    return 1;
}
