import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange

# 設定字體為 Times New Roman
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'stix'  # 數學符號也用類似 Times 的字體

# Runge's function
def runge(x):
    return 1 / (1 + 25 * x**2)

# Generate interpolation points
def evenly_spaced(a, b, n):
    return np.linspace(a, b, n)

# Chebyshev nodes
def chebyshev_nodes(a, b, n):
    k = np.arange(n)
    nodes = np.cos((2*k + 1) / (2*n) * np.pi)  # [-1, 1]
    return (a + b) / 2 + (b - a) / 2 * nodes   # 轉換到 [a, b]

# Hyperbolic tangent mesh (類似 ISLBM 使用的 tanh 網格)
def tanh_nodes(a, b, n, stretch_factor=0.9):
    """
    生成 tanh 分佈的節點 (邊界加密)
    stretch_factor: 控制邊界加密程度 (0 < a < 1)，越大越密集在邊界
    公式: x = L/2 + (L/2)/a * tanh((-1 + 2*j/N) / 2 * ln((1+a)/(1-a)))
    """
    L = b - a
    center = (a + b) / 2
    j = np.arange(n)
    N = n - 1
    a_param = stretch_factor
    
    # tanh 轉換
    arg = (-1.0 + 2.0 * j / N) / 2.0 * np.log((1 + a_param) / (1 - a_param))
    nodes = center + (L / 2) / a_param * np.tanh(arg)
    return nodes

# Set up arrays x and y ( with n = 11 )
a, b, n = -1, 1, 11
x_true = np.linspace(a, b, 500)
y_true = runge(x_true)

# Evenly spaced nodes
x_even = evenly_spaced(a, b, n)
y_even = runge(x_even)

# Chebyshev nodes
x_cheb = chebyshev_nodes(a, b, n)
y_cheb = runge(x_cheb)

# Hyperbolic tangent nodes (邊界加密)
x_tanh = tanh_nodes(a, b, n, stretch_factor=0.85)
y_tanh = runge(x_tanh)

# Lagrange interpolation
poly_even = lagrange(x_even, y_even)
poly_cheb = lagrange(x_cheb, y_cheb)
poly_tanh = lagrange(x_tanh, y_tanh)
y_interp_even = poly_even(x_true)
y_interp_cheb = poly_cheb(x_true)
y_interp_tanh = poly_tanh(x_true)

# 計算最大誤差
max_err_even = np.max(np.abs(y_interp_even - y_true))
max_err_cheb = np.max(np.abs(y_interp_cheb - y_true))
max_err_tanh = np.max(np.abs(y_interp_tanh - y_true))

print("=" * 60)
print("Lagrange Interpolation Error Comparison (n = {})".format(n))
print("=" * 60)
print(f"Evenly Spaced Nodes:    Max Error = {max_err_even:.6f}")
print(f"Chebyshev Nodes:        Max Error = {max_err_cheb:.6f}")
print(f"Hyperbolic Tanh Nodes:  Max Error = {max_err_tanh:.6f}")
print("=" * 60)

# Plot comparison (2x2 grid)
fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)

# Plot 1: Evenly spaced (Runge phenomenon)
ax1 = axes[0, 0]
ax1.plot(x_true, y_true, 'k-', label="True Function", linewidth=2)
ax1.plot(x_true, y_interp_even, 'b--', label="Lagrange Interpolation", linewidth=2)
ax1.scatter(x_even, y_even, color='red', s=50, zorder=5, label="Evenly Spaced Nodes")
ax1.set_xlabel("x"), ax1.set_ylabel("y")
ax1.set_title(f"Evenly Spaced Nodes (n={n})\n→ Runge Phenomenon! (Max Err = {max_err_even:.4f})")
ax1.legend(fontsize=9), ax1.grid(True)
ax1.set_ylim(-0.5, 1.5)

# Plot 2: Chebyshev nodes (no oscillation)
ax2 = axes[0, 1]
ax2.plot(x_true, y_true, 'k-', label="True Function", linewidth=2)
ax2.plot(x_true, y_interp_cheb, 'g--', label="Lagrange Interpolation", linewidth=2)
ax2.scatter(x_cheb, y_cheb, color='orange', s=50, zorder=5, label="Chebyshev Nodes")
ax2.set_xlabel("x"), ax2.set_ylabel("y")
ax2.set_title(f"Chebyshev Nodes (n={n})\n→ No Oscillation! (Max Err = {max_err_cheb:.4f})")
ax2.legend(fontsize=9), ax2.grid(True)
ax2.set_ylim(-0.5, 1.5)

# Plot 3: Hyperbolic tanh nodes
ax3 = axes[1, 0]
ax3.plot(x_true, y_true, 'k-', label="True Function", linewidth=2)
ax3.plot(x_true, y_interp_tanh, 'm--', label="Lagrange Interpolation", linewidth=2)
ax3.scatter(x_tanh, y_tanh, color='purple', s=50, zorder=5, label="Tanh Nodes (a=0.85)")
ax3.set_xlabel("x"), ax3.set_ylabel("y")
ax3.set_title(f"Hyperbolic Tanh Nodes (n={n})\n→ Reduced Oscillation (Max Err = {max_err_tanh:.4f})")
ax3.legend(fontsize=9), ax3.grid(True)
ax3.set_ylim(-0.5, 1.5)

# Plot 4: Node distribution comparison
ax4 = axes[1, 1]
ax4.scatter(x_even, np.ones_like(x_even)*3, color='red', s=80, label='Evenly Spaced', marker='o')
ax4.scatter(x_cheb, np.ones_like(x_cheb)*2, color='orange', s=80, label='Chebyshev', marker='s')
ax4.scatter(x_tanh, np.ones_like(x_tanh)*1, color='purple', s=80, label='Tanh (a=0.85)', marker='^')
ax4.set_xlabel("x")
ax4.set_yticks([1, 2, 3])
ax4.set_yticklabels(['Tanh', 'Chebyshev', 'Evenly'])
ax4.set_title("Node Distribution Comparison\n(Edge clustering reduces Runge phenomenon)")
ax4.legend(fontsize=9, loc='upper right')
ax4.grid(True, axis='x')
ax4.set_xlim(-1.1, 1.1)

plt.tight_layout()
plt.savefig("runge_phenomenon_comparison.png", dpi=300)
plt.show()

# 額外測試不同 stretch factor 的 tanh 網格
print("\n" + "=" * 60)
print("Tanh Mesh: Effect of Stretch Factor (a)")
print("=" * 60)
for a_param in [0.5, 0.7, 0.85, 0.9, 0.95]:
    x_test = tanh_nodes(-1, 1, n, stretch_factor=a_param)
    y_test = runge(x_test)
    poly_test = lagrange(x_test, y_test)
    y_interp_test = poly_test(x_true)
    max_err_test = np.max(np.abs(y_interp_test - y_true))
    print(f"  a = {a_param:.2f}:  Max Error = {max_err_test:.6f}")