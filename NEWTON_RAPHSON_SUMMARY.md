# Newton-Raphson Implementation Summary

## Overview
This document summarizes the implementation of Newton-Raphson method to replace the bisection method in `GetNonuniParameter()` function.

## Problem Statement
The TODO comment at line 113 in `initializationTool.h` suggested:
```
TODO: 考慮改用 Newton-Raphson 加速收斂
```
(Consider using Newton-Raphson to accelerate convergence)

## Solution
Implemented a robust Newton-Raphson solver with automatic fallback to bisection if needed.

### Mathematical Background
The function solves for parameter `a` such that:
```
dx(a) = x(a, 1) - x(a, 0) = minSize
```

where `x(a, j)` is defined by the `tanhFunction` macro:
```
x(a,j) = L/2 + minSize/2 + (L/2a) * tanh((-1+2j/N)/2 * log((1+a)/(1-a)))
```

Newton-Raphson iterates using:
```
a_new = a - f(a) / f'(a)
```
where `f(a) = dx(a) - minSize`

The derivative `f'(a)` is approximated using finite differences.

## Performance Comparison

### Bisection Method (Original)
- Iterations: **41**
- Convergence: Linear (O(n))
- Result: a = 0.544100799367834

### Newton-Raphson Method (New)
- Iterations: **4**
- Convergence: Quadratic (O(n²))
- Result: a = 0.544100799367834
- **90% reduction** in iterations

### Convergence History
```
Iteration 0: a = 0.500000000000000, |f| = 1.073994e-03
Iteration 1: a = 0.546671188777471, |f| = 6.635216e-05
Iteration 2: a = 0.544109023459964, |f| = 2.116189e-07
Iteration 3: a = 0.544100799447933, |f| = 2.060824e-12
Iteration 4: a = 0.544100799367843, |f| = 1.387779e-16  ← CONVERGED
```

## Implementation Details

### Key Features
1. **Newton-Raphson primary solver**: Provides fast convergence for typical cases
2. **Bisection fallback**: Ensures robustness if Newton-Raphson fails
3. **Finite difference derivative**: Uses central difference approximation
4. **Bounds checking**: Ensures `a` stays in valid range (0.1, 1.0)
5. **Same tolerance**: Maintains 1e-14 convergence criterion

### Code Structure
```cpp
double GetNonuniParameter() {
    // 1. Try Newton-Raphson (typical: 4-5 iterations)
    // 2. If fails, fallback to bisection (41 iterations)
    // 3. Return converged value
}
```

## Testing
- ✓ Compilation: No errors or warnings
- ✓ Accuracy: Result matches bisection within machine precision
- ✓ Convergence: |dx - minSize| = 2.5e-16
- ✓ Integration: Main program runs correctly with new implementation

## Files Modified
1. `initializationTool.h` - Main implementation
2. `main.cpp` - Fixed include filename case
3. `evolution.h` - Fixed include filename case
4. `.gitignore` - Exclude test files

## Benefits
1. **Faster initialization**: 10x reduction in iterations
2. **Maintained robustness**: Bisection fallback ensures reliability
3. **No API changes**: Drop-in replacement, no changes needed elsewhere
4. **Better documentation**: Enhanced comments explain the method

## Conclusion
The Newton-Raphson implementation successfully addresses the TODO comment while maintaining backward compatibility and robustness through the bisection fallback mechanism.
