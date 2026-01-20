//=============================================================================
// benchmark_convergence.cpp - Benchmark convergence speed comparison
// Compile: g++ -std=c++17 -O3 -o benchmark benchmark_convergence.cpp -lm
// Run: ./benchmark
//=============================================================================
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
using namespace std;

// Include necessary headers
#include "variables.h"
#include "model.h"
#include "initializationTool.h"

// Bisection method (original implementation)
double GetNonuniParameter_Bisection(int& iterations) {
    double total = LZ - HillFunction(0.0) - minSize;
    double a_temp[2] = {0.1, 1.0};
    double a_mid;
    const int N = NZ6 - 7;
    
    double x_temp[2], dx;
    iterations = 0;
    do {
        a_mid = (a_temp[0] + a_temp[1]) / 2.0;
        x_temp[0] = tanhFunction(total, minSize, a_mid, 0, N);
        x_temp[1] = tanhFunction(total, minSize, a_mid, 1, N);
        dx = x_temp[1] - x_temp[0];
        if (dx - minSize >= 0.0) {
            a_temp[0] = a_mid;
        } else {
            a_temp[1] = a_mid;
        }
        iterations++;
    } while (fabs(dx - minSize) > 1e-14);
    return a_mid;
}

// Newton-Raphson method (new implementation)
double GetNonuniParameter_NewtonRaphson(int& iterations) {
    double total = LZ - HillFunction(0.0) - minSize;
    const int N = NZ6 - 7;
    
    double a = 0.5;
    const int max_iter = 50;
    const double tol = 1e-14;
    const double a_min = 0.1, a_max = 1.0;
    
    iterations = 0;
    for (int iter = 0; iter < max_iter; iter++) {
        iterations++;
        
        double x0 = tanhFunction(total, minSize, a, 0, N);
        double x1 = tanhFunction(total, minSize, a, 1, N);
        double dx = x1 - x0;
        double f = dx - minSize;
        
        if (fabs(f) < tol) {
            return a;
        }
        
        // Compute derivative using central difference
        const double da = 1e-8;
        double x0_plus = tanhFunction(total, minSize, a + da, 0, N);
        double x1_plus = tanhFunction(total, minSize, a + da, 1, N);
        double dx_plus = x1_plus - x0_plus;
        double df_da = (dx_plus - dx) / da;
        
        if (fabs(df_da) < 1e-15) {
            break;  // Fallback to bisection
        }
        
        double a_new = a - f / df_da;
        
        if (a_new <= a_min || a_new >= a_max) {
            break;  // Fallback to bisection
        }
        
        a = a_new;
    }
    
    // If Newton-Raphson didn't converge, use bisection
    return GetNonuniParameter_Bisection(iterations);
}

int main() {
    cout << "========================================" << endl;
    cout << "  Convergence Speed Benchmark" << endl;
    cout << "========================================" << endl;
    cout << endl;
    
    // Benchmark Bisection Method
    int iter_bisection = 0;
    auto start1 = chrono::high_resolution_clock::now();
    double a_bisection = GetNonuniParameter_Bisection(iter_bisection);
    auto end1 = chrono::high_resolution_clock::now();
    auto duration1 = chrono::duration_cast<chrono::microseconds>(end1 - start1);
    
    cout << "Bisection Method:" << endl;
    cout << "  Iterations: " << iter_bisection << endl;
    cout << "  Result a = " << setprecision(15) << a_bisection << endl;
    cout << "  Time: " << duration1.count() << " microseconds" << endl;
    cout << endl;
    
    // Benchmark Newton-Raphson Method
    int iter_newton = 0;
    auto start2 = chrono::high_resolution_clock::now();
    double a_newton = GetNonuniParameter_NewtonRaphson(iter_newton);
    auto end2 = chrono::high_resolution_clock::now();
    auto duration2 = chrono::duration_cast<chrono::microseconds>(end2 - start2);
    
    cout << "Newton-Raphson Method:" << endl;
    cout << "  Iterations: " << iter_newton << endl;
    cout << "  Result a = " << setprecision(15) << a_newton << endl;
    cout << "  Time: " << duration2.count() << " microseconds" << endl;
    cout << endl;
    
    // Compare results
    cout << "Comparison:" << endl;
    cout << "  Difference in a: " << scientific << setprecision(6) << fabs(a_bisection - a_newton) << endl;
    cout << "  Iteration reduction: " << iter_bisection - iter_newton << " iterations" << endl;
    cout << "  Speedup: " << fixed << setprecision(2) << (double)iter_bisection / iter_newton << "x" << endl;
    cout << endl;
    
    if (fabs(a_bisection - a_newton) < 1e-12) {
        cout << "✓ Results match within tolerance" << endl;
        return 0;
    } else {
        cout << "✗ Results differ" << endl;
        return 1;
    }
}
