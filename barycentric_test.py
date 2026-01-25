"""
Barycentric Interpolation vs Lagrange Interpolation
證明兩者產生相同的結果，因此 Barycentric 無法解決 Runge 振盪
"""
import numpy as np
import matplotlib.pyplot as plt

# 設定字體為 Times New Roman
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'stix'

# Runge's function
def runge(x):
    return 1 / (1 + 25 * x**2)

# =============================================================================
# Standard Lagrange Interpolation (O(n²))
# =============================================================================
def lagrange_interpolation(x_nodes, y_nodes, x_eval):
    """標準 Lagrange 插值 - 直接計算"""
    n = len(x_nodes)
    result = 0.0
    for j in range(n):
        # 計算 L_j(x)
        Lj = 1.0
        for k in range(n):
            if k != j:
                Lj *= (x_eval - x_nodes[k]) / (x_nodes[j] - x_nodes[k])
        result += y_nodes[j] * Lj
    return result

# =============================================================================
# Barycentric Lagrange Interpolation (O(n) after precomputation)
# =============================================================================
def compute_barycentric_weights(x_nodes):
    """預計算重心權重 w_j = 1 / prod_{k≠j}(x_j - x_k)"""
    n = len(x_nodes)
    weights = np.ones(n)
    for j in range(n):
        for k in range(n):
            if k != j:
                weights[j] /= (x_nodes[j] - x_nodes[k])
    return weights

def barycentric_interpolation(x_nodes, y_nodes, weights, x_eval, eps=1e-15):
    """
    重心插值公式:
    L(x) = sum_j (w_j / (x - x_j)) * y_j / sum_j (w_j / (x - x_j))
    
    優點: 數值更穩定，計算效率 O(n)
    """
    n = len(x_nodes)
    
    # 檢查是否正好在節點上
    for j in range(n):
        if abs(x_eval - x_nodes[j]) < eps:
            return y_nodes[j]
    
    numerator = 0.0
    denominator = 0.0
    for j in range(n):
        temp = weights[j] / (x_eval - x_nodes[j])
        numerator += temp * y_nodes[j]
        denominator += temp
    
    return numerator / denominator

# =============================================================================
# Test: Compare Lagrange vs Barycentric
# =============================================================================
def main():
    # 參數設定
    n_nodes = 15  # 插值節點數（使用奇數以確保 Runge 振盪明顯）
    
    # 均勻分佈節點（最容易產生 Runge 振盪）
    x_nodes = np.linspace(-1, 1, n_nodes)
    y_nodes = runge(x_nodes)
    
    # 預計算 Barycentric 權重
    weights = compute_barycentric_weights(x_nodes)
    
    # 評估點
    x_fine = np.linspace(-1, 1, 500)
    y_true = runge(x_fine)
    
    # 計算兩種方法的插值結果
    y_lagrange = np.array([lagrange_interpolation(x_nodes, y_nodes, x) for x in x_fine])
    y_barycentric = np.array([barycentric_interpolation(x_nodes, y_nodes, weights, x) for x in x_fine])
    
    # 計算差異
    diff = np.abs(y_lagrange - y_barycentric)
    max_diff = np.max(diff)
    
    # 計算誤差
    error_lagrange = np.max(np.abs(y_lagrange - y_true))
    error_barycentric = np.max(np.abs(y_barycentric - y_true))
    
    print("="*60)
    print("Barycentric vs Lagrange Interpolation Comparison")
    print("="*60)
    print(f"Number of nodes: {n_nodes}")
    print(f"\n--- Results ---")
    print(f"Lagrange Max Error:     {error_lagrange:.6f}")
    print(f"Barycentric Max Error:  {error_barycentric:.6f}")
    print(f"Difference between methods: {max_diff:.2e}")
    print(f"\n>>> 兩種方法的結果相同！Barycentric 只是計算方式不同。")
    print(f">>> 因此 Barycentric 【無法解決】Runge 振盪問題。")
    print("="*60)
    
    # 繪圖
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Lagrange 結果
    ax1 = axes[0, 0]
    ax1.plot(x_fine, y_true, 'k-', linewidth=2, label='True: $f(x) = 1/(1+25x^2)$')
    ax1.plot(x_fine, y_lagrange, 'b--', linewidth=1.5, label=f'Lagrange (n={n_nodes})')
    ax1.scatter(x_nodes, y_nodes, c='red', s=50, zorder=5, label='Interpolation nodes')
    ax1.set_xlabel('x')
    ax1.set_ylabel('f(x)')
    ax1.set_title(f'Standard Lagrange Interpolation\nMax Error = {error_lagrange:.4f}')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim([-0.5, 1.5])
    
    # 2. Barycentric 結果
    ax2 = axes[0, 1]
    ax2.plot(x_fine, y_true, 'k-', linewidth=2, label='True: $f(x) = 1/(1+25x^2)$')
    ax2.plot(x_fine, y_barycentric, 'g--', linewidth=1.5, label=f'Barycentric (n={n_nodes})')
    ax2.scatter(x_nodes, y_nodes, c='red', s=50, zorder=5, label='Interpolation nodes')
    ax2.set_xlabel('x')
    ax2.set_ylabel('f(x)')
    ax2.set_title(f'Barycentric Interpolation\nMax Error = {error_barycentric:.4f}')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([-0.5, 1.5])
    
    # 3. 兩種方法的差異（應該極小）
    ax3 = axes[1, 0]
    ax3.plot(x_fine, diff, 'r-', linewidth=1.5)
    ax3.set_xlabel('x')
    ax3.set_ylabel('|Lagrange - Barycentric|')
    ax3.set_title(f'Difference Between Methods\nMax Diff = {max_diff:.2e}')
    ax3.grid(True, alpha=0.3)
    ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    # 4. 誤差分佈
    ax4 = axes[1, 1]
    error_dist = np.abs(y_lagrange - y_true)
    ax4.plot(x_fine, error_dist, 'b-', linewidth=1.5, label='Error distribution')
    ax4.axhline(y=error_lagrange, color='r', linestyle='--', label=f'Max error = {error_lagrange:.4f}')
    ax4.set_xlabel('x')
    ax4.set_ylabel('Absolute Error')
    ax4.set_title('Error Distribution (Runge Phenomenon)')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.suptitle('Proof: Barycentric Formula Cannot Solve Runge Phenomenon', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('barycentric_comparison.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print("\n圖片已儲存為 barycentric_comparison.png")
    
    # ==========================================================================
    # 額外測試：Chebyshev 節點 + Barycentric
    # ==========================================================================
    print("\n" + "="*60)
    print("Bonus: Chebyshev Nodes + Barycentric")
    print("="*60)
    
    # Chebyshev 節點
    k = np.arange(n_nodes)
    x_cheby = np.cos((2*k + 1) / (2*n_nodes) * np.pi)
    y_cheby = runge(x_cheby)
    weights_cheby = compute_barycentric_weights(x_cheby)
    
    y_cheby_interp = np.array([barycentric_interpolation(x_cheby, y_cheby, weights_cheby, x) 
                               for x in x_fine])
    error_cheby = np.max(np.abs(y_cheby_interp - y_true))
    
    print(f"Uniform nodes (Barycentric):   Max Error = {error_barycentric:.6f}")
    print(f"Chebyshev nodes (Barycentric): Max Error = {error_cheby:.6f}")
    print(f"\n>>> 結論：解決 Runge 振盪需要【Chebyshev 節點】，不是 Barycentric 公式！")
    print("="*60)

if __name__ == "__main__":
    main()
