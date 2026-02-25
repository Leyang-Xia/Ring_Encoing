# Ring-Based XOR Optimization for Erasure Coding

A C++ implementation that explores ring mapping techniques to optimize XOR operations in erasure coding systems, specifically targeting Galois Field arithmetic optimizations through quotient ring mappings.

## Overview

This project implements a novel approach to reduce XOR operations in erasure coding by mapping elements from GF(2^8) to a quotient ring defined by a reducible polynomial. The system generates encoding matrices, converts them to bit matrices, and applies the Uber optimization tool to minimize XOR operations.

## Mathematical Foundation

### Base Field
- **Field**: GF(2^8) with irreducible polynomial `x^8+x^4+x^3+x^2+1` (0x11D)

### Ring Mapping
- **Target Ring**: GF(2)[x] / ((x^8+x^4+x^3+x^2+1)(x^2+x+1))
- **Ring Modulus**: 0x753 (reducible polynomial)
- **Mapping Strategy**: f(x) → f(x) + g(x)(x^8+x^4+x^3+x^2+1) where g(x) ∈ {0, 1, x, x+1}
- **Optimization Goal**: Minimize Hamming weight through optimal g(x) selection

## 面向初学者的完整数学方案介绍（建议先读）

这一节按“零基础有限域 + 编码初学者”来写，回答三个核心问题：

1. 原来的编码到底在做什么？
2. 这个 Ring-Based 方法本质在做什么？
3. 为什么它可以这样做而不改变编码正确性？

---

### 1) 原始编码在做什么：本质是“矩阵乘法”

先不谈 Ring 映射，任何线性纠删码编码都可以写成：

- 数据向量：`x = (x_0, x_1, ..., x_{k-1})^T`
- 编码矩阵：`A`（元素在 `GF(2^8)`）
- 输出向量：`y = A x`

第 `i` 行就是：

`y_i = Σ_j A_{i,j} * x_j`（加法和乘法都在 `GF(2^8)` 中）

在这个仓库中，大体流程是：

- 先构造 Vandermonde 矩阵并求逆，得到编码矩阵 `A`；
- 再把 `GF(2^8)` 的矩阵转成二进制 bitmatrix（供 Uber/X-Sets 统计和优化 XOR 调度）。

所以“原始编码”并不神秘：就是在 `GF(2^8)` 上做线性代数。

---

### 2) 为什么会关心 XOR 次数

`GF(2^8)` 元素可以看成 8 位二进制多项式（系数在 `GF(2)`）。

- 域加法 = 按位异或（XOR）
- 域乘法 = 多项式乘法后再模不可约多项式约简

在软件实现里，编码耗时常被“需要执行多少 XOR”主导，因此项目会把矩阵展开成 bit-level 线性变换，并用 Uber/X-Sets 去找更省 XOR 的计算顺序。

---

### 3) `GF(2^8)` 的数学模型（最关键的基础）

这里固定：

- `P(x) = x^8 + x^4 + x^3 + x^2 + 1`（十六进制 `0x11D`）
- `GF(2^8) = GF(2)[x] / (P(x))`

意思是：两个多项式只要相差 `P(x)` 的倍数，就代表同一个域元素。

---

### 4) Ring-Based 方法本质：给每个元素选一个“等价提升”

项目里定义：

- `Q(x) = x^2 + x + 1`
- `M(x) = P(x)Q(x) = x^10 + x^9 + x^8 + x^6 + x^4 + x + 1`（`0x753`）
- 目标环：`R = GF(2)[x] / (M(x))`

对任意原始系数 `a`（它是 `GF(2^8)` 元素），构造：

`a' = a + g(x)P(x)`，其中 `g(x) ∈ {0, 1, x, x+1}`

这四个 `g` 恰好对应 4 种“提升版本”（因为这里本质上只需要看 `g mod Q`，而 `deg(Q)=2`，所以 `GF(2)[x]/(Q)` 一共只有 `2^2=4` 个元素）。

直觉上：

- `a` 不变其数学意义；
- 但二进制形态会变（位分布可能更稀疏）；
- 稀疏度改变会影响后续 XOR 成本。

所以这个方法不是“改了码的定义”，而是“在等价类里换表示”。

---

### 5) 为什么这样做是对的：严格正确性说明

定义投影（回映射）：

`π: R -> GF(2^8),   π([f] mod M) = [f] mod P`

因为 `P | M`，这个投影是良定义的环同态。

现在看一个输出分量：

`y'_i = Σ_j (A_{i,j} + g_{i,j}P) x_j`（在 `R` 中计算）

把它投影回 `GF(2^8)`：

`π(y'_i) = Σ_j π((A_{i,j} + g_{i,j}P) x_j)`
`       = Σ_j (π(A_{i,j}x_j) + π(g_{i,j}Px_j))`
`       = Σ_j (A_{i,j}x_j + 0)`（因为 `P` 的倍数对 `mod P` 来说都是 0）
`       = Σ_j A_{i,j} x_j`
`       = y_i`

因为 `g_{i,j}P mod P = 0`，额外项全部消失。

**结论**：不管你怎么选 `g(x)`，只要最后按 `mod P` 回来，结果和原始编码完全一致。

---

### 6) 从“代数正确”到“工程可用”的直观理解

可以把它理解为：

- 原始编码在一个“8 位坐标系”里做；
- Ring 方法先临时换到一个“10 位坐标系”里做中间计算；
- 但这个 10 位空间中有专门设计的投影能无损回到原坐标系；
- 你利用了这段“中间坐标自由度”去压低 XOR 成本。

这和“换一组基做同一个线性变换”有类似味道：变的是实现路径，不是数学目标。

---

### 7) 本仓库对应实现流程（把数学落到代码）

1. 生成 `GF(2^8)` 编码矩阵 `A`（Vandermonde + invert）。
2. 将 `A` 转为 bitmatrix，得到原始 XOR 成本。
3. 对每个 `A_{i,j}` 在四个候选 `a + gP` 中挑一个（当前实现是偏向局部稀疏度）。
4. 生成 ring-mapped 矩阵 `A'` 的 mixed bitmatrix（`10x8` 子块拼接）。
5. 对 `A'` 路径做 XOR 优化。
6. 叠加“从 10 位回 8 位”的逆映射开销，和原始路径比较总 XOR。

因此最终比较的是：

`Total_Ring_XOR = XOR(A' 路径) + XOR(回映射)`

vs

`Total_Original_XOR = XOR(A 路径)`

---

### 8) 这个方法的优化空间在哪里

当前策略通常是“逐元素选最小汉明重量”的启发式，它容易实现，但不保证全局最优。

原因是：真实总代价取决于整个矩阵的结构共享、调度复用和逆映射耦合，不是单个元素的位数相加这么简单。

可以继续优化的方向：

- 把目标从“单元素 popcount”升级为“矩阵级总代价估计”；
- 对并列最优候选引入二级目标（例如结构共享友好性）；
- 做局部搜索/坐标下降，而不是一次性独立决定每个元素的 `g`。

---

### 9) 给初学者的最短总结

- 原始编码：`GF(2^8)` 上做 `y = A x`。
- Ring 方法：把每个系数临时替换为 `a + gP`，利用表示自由度优化 XOR。
- 正确性来源：最后 `mod P` 投影，所有 `gP` 项都会消失，结果严格等于原编码输出。
- 所以它是“同一个数学映射的不同实现路径”，不是改码本身。

## Features

- **Matrix Generation**: Vandermonde matrix creation and inversion
- **Bitmatrix Conversion**: GF(2^8) to binary matrix transformation
- **Ring Mapping**: Element-wise mapping with sparsity optimization
- **Mixed Bitmatrix**: GF(2^10) × GF(2^8) multiplication matrices
- **Uber Integration**: Automatic XOR count optimization
- **Inverse Mapping**: Proper reduction from ring back to base field

## Dependencies

### Required Libraries
```bash
# GF-Complete (Galois Field arithmetic)
sudo apt-get install libgf-complete-dev

# Jerasure (Erasure coding library)
sudo apt-get install libjerasure-dev
```

### Build Requirements
- C++17 compatible compiler (GCC 7+ or Clang 5+)
- Make build system
- Uber optimization tool (should be in project directory)

## Building

```bash
# Clone the repository
git clone https://github.com/Leyang-Xia/Ring_Encoing.git
cd Ring_Encoing

# Compile the project
g++ -std=c++17 -O3 ring_based_XOR.cpp -lgf_complete -ljerasure -o ring_xor

# Make sure Uber tool is executable
chmod +x Uber
```

## Usage

### Basic Execution
```bash
./ring_xor
```

### Interactive Mode
The program offers several optimization levels:
1. Choose optimization level (1-5) when prompted
2. Review original matrix optimization results
3. Examine ring mapping statistics
4. Compare optimization improvements

### Output Files
- `matrix_vandermonde_gf8.txt` - Original Vandermonde matrix
- `matrix_A_gf8.txt` - Encoding matrix (inverted Vandermonde)
- `matrix_A_bitmatrix.txt` - Binary representation for Uber
- `matrix_A_ring_mapped.txt` - Ring-mapped encoding matrix
- `matrix_A_ring_bitmatrix.txt` - Mixed bitmatrix (100×80)
- `inverse_mapping_matrix.txt` - Reduction matrix

## Configuration

Key parameters in `Config` namespace:
```cpp
constexpr int FIELD_SIZE = 8;           // GF(2^8)
constexpr int RING_SIZE = 10;           // 10-bit ring
constexpr int MATRIX_SIZE = 10;         // 10×10 matrices
constexpr uint16_t PRIMITIVE_POLY = 0x11D;  // Base field polynomial
constexpr uint16_t RING_MODULUS = 0x753;    // Ring modulus
```

## Architecture

### Core Classes
- **`GaloisField`**: RAII wrapper for GF-Complete operations
- **`MatrixHandler`**: Matrix generation, inversion, and bitmatrix conversion
- **`RingMapper`**: Ring mapping logic with sparsity optimization
- **`UberOptimizer`**: Interface to Uber optimization tool
- **`FileManager`**: File I/O operations for matrices
- **`ErasureCodingApp`**: Main application orchestrator

### Key Algorithms
1. **Vandermonde Matrix Generation**: α^(i×j) pattern
2. **Matrix Inversion**: Jerasure-based inversion
3. **Custom Bitmatrix Conversion**: Respects specified primitive polynomial
4. **Optimal g(x) Selection**: Hamming weight minimization
5. **Mixed Bitmatrix Generation**: Ring×Field multiplication

## Performance Results

The system compares XOR counts between:
- Original GF(2^8) matrix optimization
- Ring-mapped matrix optimization + inverse mapping overhead
- Reports percentage improvement/degradation

Typical output:
```
--- Ring Mapping Comparison Results ---
Original matrix XORs: 1250
Ring-mapped matrix XORs: 980
Inverse mapping XORs: 180
Total ring mapping XORs: 1160
XOR improvement: 90 (BETTER)
Percentage improvement: 7.2%
```

## Research Applications

This implementation supports research in:
- Galois Field arithmetic optimization
- Sparse matrix techniques for erasure coding
- Ring-theoretic approaches to linear algebra
- XOR minimization in coding theory

## Contributing

1. Fork the repository
2. Create a feature branch
3. Follow C++17 best practices and existing code style
4. Add tests for new functionality
5. Submit a pull request

## License

This project is provided for research and educational purposes. See individual library licenses for GF-Complete and Jerasure.

## References

- GF-Complete: Fast Galois Field arithmetic library
- Jerasure: Library for erasure coding in storage systems
- Uber: XOR optimization tool for binary matrices

## Authors

Leyang Xia - Initial implementation and ring mapping theory

## Contact

For questions or collaboration: spicycurrykk@gmail.com 
