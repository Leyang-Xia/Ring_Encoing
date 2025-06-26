# Ring-Based XOR Optimization for Erasure Coding

A C++ implementation that explores ring mapping techniques to optimize XOR operations in erasure coding systems, specifically targeting Galois Field arithmetic optimizations through quotient ring mappings.

## Overview

This project implements a novel approach to reduce XOR operations in erasure coding by mapping elements from GF(2^8) to a quotient ring defined by a reducible polynomial. The system generates encoding matrices, converts them to bit matrices, and applies the Uber optimization tool to minimize XOR operations.

## Mathematical Foundation

### Base Field
- **Field**: GF(2^8) with irreducible polynomial `x^8+x^4+x^3+x^2+1` (0x11D)
- **Primitive element**: α = 2

### Ring Mapping
- **Target Ring**: GF(2)[x] / ((x^8+x^4+x^3+x^2+1)(x^2+x+1))
- **Ring Modulus**: 0x753 (reducible polynomial)
- **Mapping Strategy**: f(x) → f(x) + g(x)(x^8+x^4+x^3+x^2+1) where g(x) ∈ {0, 1, x, x+1}
- **Optimization Goal**: Minimize Hamming weight through optimal g(x) selection

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