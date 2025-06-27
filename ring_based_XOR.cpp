#ifndef RING_BASED_XOR_HPP
#define RING_BASED_XOR_HPP

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstdio>    // For popen, pclose, stderr
#include <array>     // For std::array
#include <algorithm> // For std::copy
#include <cstring>   // For strerror
#include <errno.h>   // For errno
#include <iomanip>   // For std::setw
#include <memory>
#include <sstream>   // For std::istringstream

// Windows compatibility for process status macros
#ifdef _WIN32
#define WIFEXITED(status) (((status) & 0x7f) == 0)
#define WEXITSTATUS(status) (((status) & 0xff00) >> 8)
#else
#include <sys/wait.h> // For WIFEXITED, WEXITSTATUS
#endif

// GF-Complete and Jerasure headers
extern "C" {
#include <gf_complete.h>
#include <jerasure.h>
}

// Configuration constants
namespace Config {
    constexpr int FIELD_SIZE = 8;            // GF(2^8) - base field
    constexpr int RING_SIZE = 10;            // 10 bits for ring polynomial degree
    constexpr int MATRIX_SIZE = 10;          // 10x10 matrix - smaller for testing
    constexpr uint16_t PRIMITIVE_POLY = 0x11D; // x^8+x^4+x^3+x^2+1 for GF(2^8)
    constexpr uint16_t RING_MODULUS = 0x753; // (x^8+x^4+x^3+x^2+1)(x^2+x+1) GF(2) multiplication
    constexpr gf_val_32_t ALPHA = 2;         // Primitive element
    
    const std::string POLY_DESCRIPTION = "x^8+x^4+x^3+x^2+1";
    const std::string RING_POLY_DESCRIPTION = "(x^8+x^4+x^3+x^2+1)(x^2+x+1)";
    const std::string VANDERMONDE_FILE = "matrix_vandermonde_gf8.txt";
    const std::string ENCODING_MATRIX_FILE = "matrix_A_gf8.txt";
    const std::string BITMATRIX_FILE = "matrix_A_bitmatrix.txt";
    const std::string BITMATRIX_READABLE_FILE = "matrix_A_bitmatrix_readable.txt";
    
    // Ring mapping files
    const std::string RING_MAPPED_MATRIX_FILE = "matrix_A_ring_mapped.txt";
    const std::string RING_BITMATRIX_FILE = "matrix_A_ring_bitmatrix.txt";
    const std::string RING_BITMATRIX_READABLE_FILE = "matrix_A_ring_bitmatrix_readable.txt";
    const std::string INVERSE_MAPPING_MATRIX_FILE = "inverse_mapping_matrix.txt";
}

// Forward declarations
class GaloisField;
class MatrixHandler;
class UberOptimizer;
class FileManager;
class RingMapper;

/**
 * @brief RAII wrapper for Galois Field operations
 */
class GaloisField {
public:
    explicit GaloisField(int field_size);
    ~GaloisField();
    
    // Delete copy constructor and assignment operator
    GaloisField(const GaloisField&) = delete;
    GaloisField& operator=(const GaloisField&) = delete;
    
    // Move constructor and assignment operator
    GaloisField(GaloisField&& other) noexcept;
    GaloisField& operator=(GaloisField&& other) noexcept;
    
    gf_t* getField() { return &m_gf_field; }
    bool isInitialized() const { return m_initialized; }
    
    gf_val_32_t multiply(gf_val_32_t a, gf_val_32_t b);
    gf_val_32_t power(gf_val_32_t base, int exponent);

private:
    gf_t m_gf_field;
    bool m_initialized;
};

/**
 * @brief Handles matrix operations and conversions
 */
class MatrixHandler {
public:
    explicit MatrixHandler(GaloisField& gf);
    
    std::vector<gf_val_32_t> generateVandermondeMatrix(int rows, int cols, gf_val_32_t alpha);
    std::vector<gf_val_32_t> invertMatrix(const std::vector<gf_val_32_t>& matrix, int size);
    std::vector<int> convertToBitmatrix(const std::vector<gf_val_32_t>& matrix, int rows, int cols);
    
    void printMatrixPreview(const std::vector<gf_val_32_t>& matrix, int rows, int cols, const std::string& name);
    
    // Custom matrix to bitmatrix conversion with specified polynomial
    std::vector<int> customMatrixToBitmatrix(const std::vector<gf_val_32_t>& matrix, int k, int m, int w, uint16_t primitive_poly);

private:
    GaloisField& m_gf;
    
    // Custom Galois field operations with specified polynomial
    gf_val_32_t galoisMultiply(gf_val_32_t a, gf_val_32_t b, int w, uint16_t primitive_poly);
};

/**
 * @brief Handles Uber optimization operations (Implementation mode only)
 */
class UberOptimizer {
public:
    struct OptimizationResult {
        int xor_count;
        int direct_xor_count;
        int overhead;
        double ratio;
        bool success;
        int level_used;
    };
    
    // Default optimization with level 2
    OptimizationResult optimize(const std::string& matrix_filename);
    
    // Optimization with specific level
    OptimizationResult optimize(const std::string& matrix_filename, int level);
    
    // Get available optimization levels (typically 1-5)
    static std::vector<int> getAvailableLevels() { return {1, 2, 3, 4, 5}; }
    
private:
    int parseXorCount(const std::string& output);
    int calculateDirectXorCount(int bit_rows, int bit_cols);
    std::string buildUberCommand(const std::string& matrix_filename, int level);
};

/**
 * @brief Implements ring mapping from GF(2^8) to a quotient ring
 * 
 * MATHEMATICAL BASIS:
 * - Base field: GF(2^8) with irreducible polynomial x^8+x^4+x^3+x^2+1 (0x11D)
 * - Target ring: GF(2)[x] / ((x^8+x^4+x^3+x^2+1)(x^2+x+1)) = GF(2)[x] / (x¹⁰+x⁹+x⁸+x⁶+x⁴+x+1) = GF(2)[x] / (0x753)
 * 
 * RING MAPPING THEORY:
 * - The target modulus is REDUCIBLE: (x^8+x^4+x^3+x^2+1)(x^2+x+1)
 * - This allows proper inverse mapping via modular reduction
 * - Maps f(x) ∈ GF(2^8) to f(x) + g(x)(x^8+x^4+x^3+x^2+1) in the quotient ring
 * - Inverse mapping: reduce modulo (x^8+x^4+x^3+x^2+1) to get back to GF(2^8)
 * - g(x) ∈ {0, 1, x, x+1} is chosen to minimize Hamming weight
 */
class RingMapper {
public:
    struct MappingResult {
        std::vector<gf_val_32_t> mapped_matrix;
        std::vector<gf_val_32_t> inverse_mapping_matrix;
        std::vector<int> g_choices;  // Store g(x) choice for each element
        int hamming_weight_original;
        int hamming_weight_mapped;
        bool success;
    };
    
    explicit RingMapper(GaloisField& base_field);
    ~RingMapper();
    
    // Map a single element from GF(2^8) to GF(2^10)
    gf_val_32_t mapElement(gf_val_32_t element, int g_coeff = 0);
    
    // Map back from GF(2^10) to GF(2^8)
    gf_val_32_t mapElementBack(gf_val_32_t ring_element);
    
    // Find optimal g coefficient for sparsity
    int findOptimalG(gf_val_32_t element);
    
    // Map entire matrix with optimal g choices
    MappingResult mapMatrix(const std::vector<gf_val_32_t>& matrix, int rows, int cols);
    
    // Generate mixed bitmatrix for A' (GF2^10) × x (GF2^8) multiplication
    std::vector<int> generateMixedBitmatrix(const std::vector<gf_val_32_t>& mapped_matrix, 
                                           int rows, int cols);
    
    // Generate inverse mapping matrix
    std::vector<gf_val_32_t> generateInverseMappingMatrix(const std::vector<int>& g_choices);
    
    // Calculate Hamming weight of matrix when converted to bitmatrix
    int calculateHammingWeight(const std::vector<gf_val_32_t>& matrix, int rows, int cols, int field_size);
    
    // Print mapping statistics
    void printMappingStatistics(const MappingResult& result);

private:
    GaloisField& m_base_field;
    std::unique_ptr<GaloisField> m_ring_field;
    
    // Cache for element mappings
    std::vector<gf_val_32_t> m_element_map_g0;  // g(x) = 0
    std::vector<gf_val_32_t> m_element_map_g1;  // g(x) = 1
    std::vector<gf_val_32_t> m_element_map_gx;  // g(x) = x
    std::vector<gf_val_32_t> m_element_map_g1x; // g(x) = x + 1
    std::vector<int> m_optimal_g;               // Optimal g for each element (0, 1, 2, 3)
    
    void initializeMappingCache();
    gf_val_32_t polynomialMultiply(gf_val_32_t a, gf_val_32_t b, gf_val_32_t modulus);
    gf_val_32_t polynomialMod(gf_val_32_t dividend, gf_val_32_t divisor);
    gf_val_32_t multiplyRingWithField(gf_val_32_t ring_element, gf_val_32_t field_element);
    gf_val_32_t generateReductionCoefficient(int position, const std::vector<int>& g_choices);
};

/**
 * @brief Handles file I/O operations
 */
class FileManager {
public:
    static void writeMatrixToFile(const std::vector<gf_val_32_t>& matrix, int rows, int cols, 
                                  const std::string& filename, const std::string& description);
    
    static void writeBitmatrixToUberFile(const std::vector<int>& bitmatrix, int rows, int cols, 
                                         const std::string& filename);
    
    static void writeBitmatrixReadable(const std::vector<int>& bitmatrix, int rows, int cols, 
                                       const std::string& filename, int field_size);
};

/**
 * @brief Main application class that orchestrates the entire process
 */
class ErasureCodingApp {
public:
    ErasureCodingApp();
    int run();

private:
    void printConfiguration();
    void processMatrices();
    void runOptimizationWithOptions();
    void processRingMapping();
    void runRingMappingComparison();

    void printOptimizationMenu();
    bool getUserOptimizationChoice(int& level);
    void displayOptimizationResult(const UberOptimizer::OptimizationResult& result);
    void displayRingMappingComparison(const UberOptimizer::OptimizationResult& original, 
                                    const UberOptimizer::OptimizationResult& mapped,
                                    const UberOptimizer::OptimizationResult& inverse);
    void displayRingMappingComparisonNew(const UberOptimizer::OptimizationResult& original, 
                                       const UberOptimizer::OptimizationResult& mapped,
                                       int inverse_xor_count);
    
    std::unique_ptr<GaloisField> m_gf;
    std::unique_ptr<MatrixHandler> m_matrix_handler;
    std::unique_ptr<UberOptimizer> m_optimizer;
    std::unique_ptr<RingMapper> m_ring_mapper;
    
    std::vector<gf_val_32_t> m_vandermonde_matrix;
    std::vector<gf_val_32_t> m_encoding_matrix;
    std::vector<int> m_bitmatrix;
    
    // Ring mapping results
    RingMapper::MappingResult m_ring_mapping_result;
    std::vector<int> m_ring_bitmatrix;
    int m_inverse_xor_count;  // XOR count for inverse mapping
    
    // Optimization results
    UberOptimizer::OptimizationResult m_original_optimization_result;
};

#endif // RING_BASED_XOR_HPP



GaloisField::GaloisField(int field_size) : m_initialized(false) {
    if (gf_init_easy(&m_gf_field, field_size) != 0) {
        m_initialized = true;
    }
}

GaloisField::~GaloisField() {
    if (m_initialized) {
        gf_free(&m_gf_field, 0);
    }
}

GaloisField::GaloisField(GaloisField&& other) noexcept 
    : m_gf_field(other.m_gf_field), m_initialized(other.m_initialized) {
    other.m_initialized = false;
}

GaloisField& GaloisField::operator=(GaloisField&& other) noexcept {
    if (this != &other) {
        if (m_initialized) {
            gf_free(&m_gf_field, 0);
        }
        m_gf_field = other.m_gf_field;
        m_initialized = other.m_initialized;
        other.m_initialized = false;
    }
    return *this;
}

gf_val_32_t GaloisField::multiply(gf_val_32_t a, gf_val_32_t b) {
    return m_gf_field.multiply.w32(&m_gf_field, a, b);
}

gf_val_32_t GaloisField::power(gf_val_32_t base, int exponent) {
    if (exponent == 0) return 1;
    
    gf_val_32_t result = 1;
    for (int i = 0; i < exponent; ++i) {
        result = multiply(result, base);
    }
    return result;
}

MatrixHandler::MatrixHandler(GaloisField& gf) : m_gf(gf) {}

std::vector<gf_val_32_t> MatrixHandler::generateVandermondeMatrix(int rows, int cols, gf_val_32_t alpha) {
    std::vector<gf_val_32_t> matrix(rows * cols);
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (i == 0 || j == 0) {
                matrix[i * cols + j] = 1;  // alpha^0 = 1
            } else {
                matrix[i * cols + j] = m_gf.power(alpha, i * j);
            }
        }
    }
    
    return matrix;
}

std::vector<gf_val_32_t> MatrixHandler::invertMatrix(const std::vector<gf_val_32_t>& matrix, int size) {
    // Convert to int array for Jerasure
    std::vector<int> matrix_int(size * size);
    for (size_t i = 0; i < matrix.size(); ++i) {
        matrix_int[i] = static_cast<int>(matrix[i]);
    }
    
    std::vector<int> inv_int(size * size);
    int result = jerasure_invert_matrix(matrix_int.data(), inv_int.data(), size, Config::FIELD_SIZE);
    if (result != 0) {
        throw std::runtime_error("Failed to invert matrix!");
    }
    
    // Convert back to gf_val_32_t
    std::vector<gf_val_32_t> inv_matrix(size * size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            inv_matrix[i * size + j] = static_cast<gf_val_32_t>(inv_int[i * size + j]);
        }
    }
    
    return inv_matrix;
}

std::vector<int> MatrixHandler::convertToBitmatrix(const std::vector<gf_val_32_t>& matrix, int rows, int cols) {
    // Use our custom implementation with the specified primitive polynomial
    return customMatrixToBitmatrix(matrix, cols, rows, Config::FIELD_SIZE, Config::PRIMITIVE_POLY);
}

void MatrixHandler::printMatrixPreview(const std::vector<gf_val_32_t>& matrix, int rows, int cols, const std::string& name) {
    constexpr int PREVIEW_SIZE = 5;
    
    std::cout << "\nFirst " << PREVIEW_SIZE << "x" << PREVIEW_SIZE << " elements of " << name << ":" << std::endl;
    for (int i = 0; i < std::min(PREVIEW_SIZE, rows); ++i) {
        for (int j = 0; j < std::min(PREVIEW_SIZE, cols); ++j) {
            std::cout << std::hex << std::setw(4) << matrix[i * cols + j] << " ";
        }
        std::cout << (cols > PREVIEW_SIZE ? "..." : "") << std::endl;
    }
    if (rows > PREVIEW_SIZE) std::cout << "..." << std::endl;
    std::cout << std::dec << std::endl;
}

// Custom Galois field multiplication with specified primitive polynomial
gf_val_32_t MatrixHandler::galoisMultiply(gf_val_32_t a, gf_val_32_t b, int w, uint16_t primitive_poly) {
    gf_val_32_t result = 0;
    gf_val_32_t temp_a = a;
    gf_val_32_t temp_b = b;
    
    // Polynomial multiplication in GF(2^w)
    while (temp_b) {
        if (temp_b & 1) {
            result ^= temp_a;
        }
        temp_a <<= 1;
        temp_b >>= 1;
        
        // Reduce if degree exceeds w
        if (temp_a & (1 << w)) {
            temp_a ^= primitive_poly;
        }
    }
    return result;
}

// Custom matrix to bitmatrix conversion based on Jerasure's implementation
// but using our specified primitive polynomial
std::vector<int> MatrixHandler::customMatrixToBitmatrix(const std::vector<gf_val_32_t>& matrix, int k, int m, int w, uint16_t primitive_poly) {
    std::vector<int> bitmatrix(k * m * w * w);
    int rowelts = k * w;
    int rowindex = 0;
    
    for (int i = 0; i < m; i++) {
        int colindex = rowindex;
        for (int j = 0; j < k; j++) {
            gf_val_32_t elt = matrix[i * k + j];
            
            for (int x = 0; x < w; x++) {
                for (int l = 0; l < w; l++) {
                    bitmatrix[colindex + x + l * rowelts] = ((elt & (1 << l)) ? 1 : 0);
                }
                // Multiply by the primitive element (2) using our custom multiplication
                // Note: We use our custom galoisMultiply to ensure we use the specified primitive_poly
                // galois_single_multiply from Jerasure uses predefined polynomials that may differ
                elt = galoisMultiply(elt, 2, w, primitive_poly);
            }
            colindex += w;
        }
        rowindex += rowelts * w;
    }
    
    return bitmatrix;
}

// RingMapper implementation
RingMapper::RingMapper(GaloisField& base_field) : m_base_field(base_field) {
    // Note: We're now working with a ring (not a field) defined by the reducible polynomial
    // So we don't initialize a GaloisField for the ring
    initializeMappingCache();
}

RingMapper::~RingMapper() = default;

void RingMapper::initializeMappingCache() {
    int field_size = 1 << Config::FIELD_SIZE;  // 2^8 = 256
    m_element_map_g0.resize(field_size);
    m_element_map_g1.resize(field_size);
    m_element_map_gx.resize(field_size);
    m_element_map_g1x.resize(field_size);
    m_optimal_g.resize(field_size);
    
    // GF(2^8) modulus: x^8+x^4+x^3+x^2+1 = 0x11D
    gf_val_32_t base_modulus = Config::PRIMITIVE_POLY;
    
    for (int elem = 0; elem < field_size; ++elem) {
        // Map with g(x) = 0: f(x) -> f(x)
        m_element_map_g0[elem] = elem;
        
        // Map with g(x) = 1: f(x) -> f(x) + 1*(x^8+x^4+x^3+x^2+1)
        m_element_map_g1[elem] = elem ^ base_modulus;
        
        // Map with g(x) = x: f(x) -> f(x) + x*(x^8+x^4+x^3+x^2+1)
        // x*(x^8+x^4+x^3+x^2+1) = x^9+x^5+x^4+x^3+x = 0x238
        gf_val_32_t x_times_base = base_modulus << 1;  // Shift left by 1 for multiplication by x
        m_element_map_gx[elem] = elem ^ x_times_base;
        
        // Map with g(x) = x+1: f(x) -> f(x) + (x+1)*(x^8+x^4+x^3+x^2+1)
        // (x+1)*(x^8+x^4+x^3+x^2+1) = x^9+x^5+x^4+x^3+x + x^8+x^4+x^3+x^2+1
        //                           = x^9+x^8+x^5+x^2+x+1
        gf_val_32_t x1_times_base = x_times_base ^ base_modulus;
        m_element_map_g1x[elem] = elem ^ x1_times_base;
        
        // Reduce modulo the ring polynomial if necessary
        if (m_element_map_g1[elem] >= (1 << Config::RING_SIZE)) {
            m_element_map_g1[elem] = polynomialMod(m_element_map_g1[elem], Config::RING_MODULUS);
        }
        if (m_element_map_gx[elem] >= (1 << Config::RING_SIZE)) {
            m_element_map_gx[elem] = polynomialMod(m_element_map_gx[elem], Config::RING_MODULUS);
        }
        if (m_element_map_g1x[elem] >= (1 << Config::RING_SIZE)) {
            m_element_map_g1x[elem] = polynomialMod(m_element_map_g1x[elem], Config::RING_MODULUS);
        }
        
        // Find optimal g by comparing Hamming weights of all four options
        int weight_g0 = __builtin_popcount(m_element_map_g0[elem]);
        int weight_g1 = __builtin_popcount(m_element_map_g1[elem]);
        int weight_gx = __builtin_popcount(m_element_map_gx[elem]);
        int weight_g1x = __builtin_popcount(m_element_map_g1x[elem]);
        
        // Find the option with minimum Hamming weight
        int min_weight = weight_g0;
        int optimal_choice = 0;
        
        if (weight_g1 < min_weight) {
            min_weight = weight_g1;
            optimal_choice = 1;
        }
        if (weight_gx < min_weight) {
            min_weight = weight_gx;
            optimal_choice = 2;
        }
        if (weight_g1x < min_weight) {
            min_weight = weight_g1x;
            optimal_choice = 3;
        }
        
        m_optimal_g[elem] = optimal_choice;
    }
}

gf_val_32_t RingMapper::mapElement(gf_val_32_t element, int g_coeff) {
    if (element >= m_element_map_g0.size()) {
        return element;  // Invalid element, return as-is
    }
    
    switch (g_coeff) {
        case 0: return m_element_map_g0[element];   // g(x) = 0
        case 1: return m_element_map_g1[element];   // g(x) = 1
        case 2: return m_element_map_gx[element];   // g(x) = x
        case 3: return m_element_map_g1x[element];  // g(x) = x + 1
        default: return m_element_map_g0[element];  // Fallback to g(x) = 0
    }
}

gf_val_32_t RingMapper::mapElementBack(gf_val_32_t ring_element) {
    // Map back by taking modulo the base field modulus
    return polynomialMod(ring_element, Config::PRIMITIVE_POLY);
}

int RingMapper::findOptimalG(gf_val_32_t element) {
    if (element >= m_optimal_g.size()) {
        return 0;
    }
    return m_optimal_g[element];
}

RingMapper::MappingResult RingMapper::mapMatrix(const std::vector<gf_val_32_t>& matrix, int rows, int cols) {
    MappingResult result;
    result.mapped_matrix.resize(rows * cols);
    result.g_choices.resize(rows * cols);
    result.success = true;
    
    // Calculate original Hamming weight
    result.hamming_weight_original = calculateHammingWeight(matrix, rows, cols, Config::FIELD_SIZE);
    
    // Map each element with optimal g
    for (int i = 0; i < rows * cols; ++i) {
        gf_val_32_t original_elem = matrix[i];
        int optimal_g = findOptimalG(original_elem);
        
        result.mapped_matrix[i] = mapElement(original_elem, optimal_g);
        result.g_choices[i] = optimal_g;
    }
    
    // Calculate mapped Hamming weight
    result.hamming_weight_mapped = calculateHammingWeight(result.mapped_matrix, rows, cols, Config::RING_SIZE);
    
    // Generate inverse mapping matrix
    result.inverse_mapping_matrix = generateInverseMappingMatrix(result.g_choices);
    
    return result;
}

std::vector<int> RingMapper::generateMixedBitmatrix(const std::vector<gf_val_32_t>& mapped_matrix, 
                                                   int rows, int cols) {
    // Create mixed bitmatrix: 100×80
    // Each element of A' (GF2^10) creates a 10×8 sub-block for multiplication with GF2^8 elements
    
    int result_rows = Config::RING_SIZE * rows;      // 100 rows
    int result_cols = Config::FIELD_SIZE * cols;     // 80 cols
    std::vector<int> mixed_bitmatrix(result_rows * result_cols, 0);
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            gf_val_32_t ring_element = mapped_matrix[i * cols + j];
            
            // For each ring element, create a 10×8 sub-block
            // This represents multiplication of ring_element with all possible GF2^8 values
            for (int bit_row = 0; bit_row < Config::RING_SIZE; ++bit_row) {
                for (int bit_col = 0; bit_col < Config::FIELD_SIZE; ++bit_col) {
                    // Calculate the contribution of this bit position
                    // ring_element × (2^bit_col) and check bit_row of result
                    
                    gf_val_32_t test_value = 1 << bit_col;  // 2^bit_col in GF2^8
                    gf_val_32_t product = multiplyRingWithField(ring_element, test_value);
                    
                    // Check if bit_row is set in the product
                    int bit_value = (product >> bit_row) & 1;
                    
                    int global_row = i * Config::RING_SIZE + bit_row;
                    int global_col = j * Config::FIELD_SIZE + bit_col;
                    mixed_bitmatrix[global_row * result_cols + global_col] = bit_value;
                }
            }
        }
    }
    
    return mixed_bitmatrix;
}

std::vector<gf_val_32_t> RingMapper::generateInverseMappingMatrix(const std::vector<int>& g_choices) {
    // Create the inverse mapping matrix that projects GF(2^10) results back to GF(2^8)
    // This should be an 8x10 matrix based on the polynomial x^10 + x^3 + 1
    
    int output_size = Config::FIELD_SIZE;  // 8 bits for GF(2^8)
    int input_size = Config::RING_SIZE;    // 10 bits for GF(2^10)
    std::vector<gf_val_32_t> inverse_matrix(output_size * input_size, 0);
    
    // Build the reduction matrix based on x^10 + x^3 + 1 = 0
    // This means x^10 = x^3 + 1, so x^9 = (x^3 + 1)/x = x^2 + x^(-1)
    // But we need to work in the polynomial ring directly
    
    // Identity mapping for x^7 down to x^0
    for (int i = 0; i < output_size; ++i) {
        inverse_matrix[i * input_size + i] = 1;  // x^i -> x^i directly
    }
    
    // Reduction from GF(2^10) to GF(2^8) based on the base field polynomial
    // Base field polynomial: x^8+x^4+x^3+x^2+1 (0x11D)
    // We need to reduce x^8 and x^9 modulo this polynomial
    
    // From x^8+x^4+x^3+x^2+1 = 0, we get: x^8 = x^4+x^3+x^2+1
    // And: x^9 = x*(x^8) = x*(x^4+x^3+x^2+1) = x^5+x^4+x^3+x
    
    // x^8 column (index 8): x^8 = x^4+x^3+x^2+1
    inverse_matrix[0 * input_size + 8] = 1;  // x^0 (constant term)
    inverse_matrix[2 * input_size + 8] = 1;  // x^2
    inverse_matrix[3 * input_size + 8] = 1;  // x^3
    inverse_matrix[4 * input_size + 8] = 1;  // x^4
    
    // x^9 column (index 9): x^9 = x^5+x^4+x^3+x
    inverse_matrix[1 * input_size + 9] = 1;  // x^1
    inverse_matrix[3 * input_size + 9] = 1;  // x^3
    inverse_matrix[4 * input_size + 9] = 1;  // x^4
    inverse_matrix[5 * input_size + 9] = 1;  // x^5
    
    return inverse_matrix;
}

// Helper method to generate reduction coefficients
gf_val_32_t RingMapper::generateReductionCoefficient(int position, const std::vector<int>& g_choices) {
    // Generate a non-trivial coefficient for the inverse mapping
    // This simulates the complexity of reducing from GF(2^10) to GF(2^8)
    
    gf_val_32_t base_coeff = 1;
    
    // Add complexity based on the g choice at this position
    if (position < static_cast<int>(g_choices.size())) {
        if (g_choices[position] == 1) {
            // For elements that used g(x)=1, we need more complex reduction
            base_coeff = this->polynomialMod(Config::PRIMITIVE_POLY >> 1, Config::PRIMITIVE_POLY);
            if (base_coeff == 0) base_coeff = 3;  // Ensure non-zero
        } else {
            // For elements that used g(x)=0, simpler reduction
            base_coeff = 1 + (position % 7);  // Small variation
        }
    }
    
    // Ensure the coefficient is valid in GF(2^8)
    return base_coeff & 0xFF;
}

int RingMapper::calculateHammingWeight(const std::vector<gf_val_32_t>& matrix, int rows, int cols, int field_size) {
    int total_weight = 0;
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            gf_val_32_t elem = matrix[i * cols + j];
            // Count 1s in the binary representation
            total_weight += __builtin_popcount(elem);
        }
    }
    
    return total_weight;
}

void RingMapper::printMappingStatistics(const MappingResult& result) {
    std::cout << "\n--- Ring Mapping Statistics ---" << std::endl;
    std::cout << "Original Hamming weight: " << result.hamming_weight_original << std::endl;
    std::cout << "Mapped Hamming weight: " << result.hamming_weight_mapped << std::endl;
    
    if (result.hamming_weight_original > 0) {
        double improvement = (double)(result.hamming_weight_original - result.hamming_weight_mapped) / result.hamming_weight_original * 100.0;
        std::cout << "Sparsity improvement: " << std::fixed << std::setprecision(1) << improvement << "%" << std::endl;
    }
}

gf_val_32_t RingMapper::polynomialMultiply(gf_val_32_t a, gf_val_32_t b, gf_val_32_t modulus) {
    gf_val_32_t result = 0;
    while (b) {
        if (b & 1) {
            result ^= a;
        }
        a <<= 1;
        b >>= 1;
    }
    return polynomialMod(result, modulus);
}

gf_val_32_t RingMapper::polynomialMod(gf_val_32_t dividend, gf_val_32_t divisor) {
    if (divisor == 0) return dividend;
    
    int divisor_degree = 31 - __builtin_clz(divisor);  // Degree of divisor
    int dividend_degree = (dividend == 0) ? -1 : 31 - __builtin_clz(dividend);
    
    while (dividend_degree >= divisor_degree && dividend != 0) {
        int shift = dividend_degree - divisor_degree;
        dividend ^= (divisor << shift);
        dividend_degree = (dividend == 0) ? -1 : 31 - __builtin_clz(dividend);
    }
    
    return dividend;
}

gf_val_32_t RingMapper::multiplyRingWithField(gf_val_32_t ring_element, gf_val_32_t field_element) {
    // Multiply ring_element (in GF2^10) with field_element (treated as element in GF2^10)
    // Since field_element < 2^8, it's automatically valid in GF2^10
    
    gf_val_32_t result = 0;
    gf_val_32_t temp_ring = ring_element;
    gf_val_32_t temp_field = field_element;
    
    // Polynomial multiplication in GF2^10
    while (temp_field) {
        if (temp_field & 1) {
            result ^= temp_ring;
        }
        temp_ring <<= 1;
        temp_field >>= 1;
        
        // Reduce if degree exceeds ring size
        if (temp_ring >= (1 << Config::RING_SIZE)) {
            temp_ring = polynomialMod(temp_ring, Config::RING_MODULUS);
        }
    }
    
    return polynomialMod(result, Config::RING_MODULUS);
}

UberOptimizer::OptimizationResult UberOptimizer::optimize(const std::string& matrix_filename) {
    // Default: use level 2
    return optimize(matrix_filename, 2);
}

UberOptimizer::OptimizationResult UberOptimizer::optimize(const std::string& matrix_filename, int level) {
    OptimizationResult result{};
    result.level_used = level;
    
    std::string command = buildUberCommand(matrix_filename, level);
    std::string full_output;
    
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Failed to run Uber command: " << command << std::endl;
        return result;
    }
    
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        full_output += buffer;
    }
    
    int status = pclose(pipe);
    if (status == -1) {
        std::cerr << "Warning: pclose() failed for command: " << command << std::endl;
        return result;
    }
    
    if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {
        std::cerr << "Warning: Uber command exited with status " << WEXITSTATUS(status) << std::endl;
    }
    
    result.xor_count = parseXorCount(full_output);
    if (result.xor_count != -1) {
        int bit_size = Config::FIELD_SIZE * Config::MATRIX_SIZE;
        result.direct_xor_count = calculateDirectXorCount(bit_size, bit_size);
        result.overhead = result.xor_count - result.direct_xor_count;
        result.ratio = static_cast<double>(result.xor_count) / result.direct_xor_count;
        result.success = true;
    }
    
    return result;
}

int UberOptimizer::parseXorCount(const std::string& output) {
    // Try to find "Total XOR operations: " in the output (this is the actual format from Uber)
    size_t pos = output.rfind("Total XOR operations: ");
    if (pos != std::string::npos) {
        try {
            // Extract the number after "Total XOR operations: "
            std::string substr = output.substr(pos + 22);  // 22 is length of "Total XOR operations: "
            // Find the first newline to get just the number
            size_t newline = substr.find('\n');
            if (newline != std::string::npos) {
                substr = substr.substr(0, newline);
            }
            return std::stoi(substr);
        } catch (...) {
            // Fall through to other parsing methods
        }
    }
    
    // Fallback: Try to find "Total XORs: " (alternative format)
    pos = output.rfind("Total XORs: ");
    if (pos != std::string::npos) {
        try {
            std::string substr = output.substr(pos + 12);
            size_t newline = substr.find('\n');
            if (newline != std::string::npos) {
                substr = substr.substr(0, newline);
            }
            return std::stoi(substr);
        } catch (...) {
            // Fall through to other parsing methods
        }
    }
    
    // Final fallback: Try to parse the last line as a number
    size_t last_newline = output.find_last_of('\n');
    if (last_newline != std::string::npos && last_newline > 0) {
        std::string last_line = output.substr(last_newline + 1);
        try {
            return std::stoi(last_line);
        } catch (...) {
            // Parsing failed
        }
    }
    
    std::cerr << "Failed to parse XOR count from Uber output." << std::endl;
    return -1;
}

std::string UberOptimizer::buildUberCommand(const std::string& matrix_filename, int level) {
    return "./Uber I " + std::to_string(level) + " < " + matrix_filename;
}

int UberOptimizer::calculateDirectXorCount(int bit_rows, int bit_cols) {    
    std::ifstream file(Config::BITMATRIX_FILE);
    if (!file) {
        return bit_rows * bit_cols - bit_rows; 
    }
    
    int total_xors = 0;
    std::string line;
    
    while (std::getline(file, line)) {
        if (!line.empty()) {
            int ones_count = 0;
            for (char c : line) {
                if (c == '1') ones_count++;
            }
            int xors_needed = std::max(0, ones_count - 1);
            total_xors += xors_needed;
        }
    }
    
    file.close();
    return total_xors;
}

void FileManager::writeMatrixToFile(const std::vector<gf_val_32_t>& matrix, int rows, int cols, 
                                   const std::string& filename, const std::string& description) {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    file << "# " << description << " (" << rows << "x" << cols << ")\n";
    file << "# Elements in GF(2^" << Config::FIELD_SIZE << ") hexadecimal representation\n\n";
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file << std::hex << std::setw(4) << matrix[i * cols + j];
            if (j < cols - 1) file << " ";
        }
        file << std::endl;
    }
}

void FileManager::writeBitmatrixToUberFile(const std::vector<int>& bitmatrix, int rows, int cols, 
                                          const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file << bitmatrix[i * cols + j];
        }
        file << std::endl;
    }
}

void FileManager::writeBitmatrixReadable(const std::vector<int>& bitmatrix, int rows, int cols, 
                                        const std::string& filename, int field_size) {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    // Determine the block structure based on matrix dimensions
    if (rows == 100 && cols == 80) {
        // Mixed bitmatrix: 100x80 for GF(2^10) × GF(2^8)
        file << "# Mixed bitmatrix representation (" << rows << "x" << cols << ")\n";
        file << "# Each 10x8 block corresponds to one GF(2^10) element × GF(2^8) element\n";
        file << "# Rows grouped by 10 (GF(2^10) bits), columns grouped by 8 (GF(2^8) bits)\n\n";
        
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                file << bitmatrix[i * cols + j];
                if ((j + 1) % 8 == 0 && j < cols - 1) file << " ";  // Group by 8 for GF(2^8)
            }
            file << std::endl;
            if ((i + 1) % 10 == 0 && i < rows - 1) file << std::endl;  // Group by 10 for GF(2^10)
        }
    } else {
        // Regular square bitmatrix
        file << "# Bitmatrix representation (" << rows << "x" << cols << ")\n";
        file << "# Each " << field_size << "x" << field_size << " block corresponds to one GF(2^" << field_size << ") element\n\n";
        
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                file << bitmatrix[i * cols + j];
                if ((j + 1) % field_size == 0 && j < cols - 1) file << " ";
            }
            file << std::endl;
            if ((i + 1) % field_size == 0 && i < rows - 1) file << std::endl;
        }
    }
}

ErasureCodingApp::ErasureCodingApp() {
    m_gf = std::make_unique<GaloisField>(Config::FIELD_SIZE);
    if (!m_gf->isInitialized()) {
        throw std::runtime_error("Failed to initialize GF(2^" + std::to_string(Config::FIELD_SIZE) + ")");
    }
    
    m_matrix_handler = std::make_unique<MatrixHandler>(*m_gf);
    m_optimizer = std::make_unique<UberOptimizer>();
    m_ring_mapper = std::make_unique<RingMapper>(*m_gf);
}

int ErasureCodingApp::run() {
    try {
        printConfiguration();
        processMatrices();
        runOptimizationWithOptions();
        processRingMapping();
        runRingMappingComparison();
        
        std::cout << "\nCompleted successfully." << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

void ErasureCodingApp::printConfiguration() {
    std::cout << "=== Configuration ===" << std::endl;
    std::cout << "Base field: GF(2^" << Config::FIELD_SIZE << ") with irreducible polynomial: " << Config::POLY_DESCRIPTION << std::endl;
    std::cout << "Base polynomial representation: 0x" << std::hex << Config::PRIMITIVE_POLY << std::dec << std::endl;
    std::cout << "Target ring: GF(2)[x] / " << Config::RING_POLY_DESCRIPTION << std::endl;
    std::cout << "Ring modulus representation: 0x" << std::hex << Config::RING_MODULUS << std::dec << std::endl;
    std::cout << "Matrix size: " << Config::MATRIX_SIZE << "x" << Config::MATRIX_SIZE << std::endl;
    std::cout << "Ring mapping: GF(2^" << Config::FIELD_SIZE << ") -> quotient ring (reducible modulus)" << std::endl;
}

void ErasureCodingApp::processMatrices() {
    // Generate Vandermonde matrix
    m_vandermonde_matrix = m_matrix_handler->generateVandermondeMatrix(Config::MATRIX_SIZE, Config::MATRIX_SIZE, Config::ALPHA);
    
    // Invert to get encoding matrix
    m_encoding_matrix = m_matrix_handler->invertMatrix(m_vandermonde_matrix, Config::MATRIX_SIZE);
    
    // Convert to bitmatrix
    m_bitmatrix = m_matrix_handler->convertToBitmatrix(m_encoding_matrix, Config::MATRIX_SIZE, Config::MATRIX_SIZE);
    
    // Save matrices to files
    FileManager::writeMatrixToFile(m_vandermonde_matrix, Config::MATRIX_SIZE, Config::MATRIX_SIZE, 
                                   Config::VANDERMONDE_FILE, "Original Vandermonde Matrix");
    
    FileManager::writeMatrixToFile(m_encoding_matrix, Config::MATRIX_SIZE, Config::MATRIX_SIZE, 
                                   Config::ENCODING_MATRIX_FILE, "Encoding Matrix A (Inverted Vandermonde)");
    
    int bit_size = Config::FIELD_SIZE * Config::MATRIX_SIZE;
    FileManager::writeBitmatrixToUberFile(m_bitmatrix, bit_size, bit_size, Config::BITMATRIX_FILE);
    FileManager::writeBitmatrixReadable(m_bitmatrix, bit_size, bit_size, Config::BITMATRIX_READABLE_FILE, Config::FIELD_SIZE);
    
    std::cout << "Matrices generated and saved to files." << std::endl;
}

void ErasureCodingApp::runOptimizationWithOptions() {
    std::cout << "\n--- Uber Optimization ---" << std::endl;
    
    int level;
    
    if (!getUserOptimizationChoice(level)) {
        level = 2;  // Default
    }
    
    m_original_optimization_result = m_optimizer->optimize(Config::BITMATRIX_FILE, level);
    displayOptimizationResult(m_original_optimization_result);
}

void ErasureCodingApp::printOptimizationMenu() {
    std::cout << "\nChoose optimization level (1-5) or press Enter for default (2): ";
}

bool ErasureCodingApp::getUserOptimizationChoice(int& level) {
    printOptimizationMenu();
    
    std::string input;
    std::getline(std::cin, input);
    
    // If user just pressed Enter, use defaults
    if (input.empty()) {
        return false;
    }
    
    // Parse input
    std::istringstream iss(input);
    
    if (!(iss >> level)) {
        std::cout << "Invalid input format." << std::endl;
        return false;
    }
    
    // Validate level
    if (level < 1 || level > 5) {
        std::cout << "Invalid level. Please choose 1-5." << std::endl;
        return false;
    }
    
    return true;
}

void ErasureCodingApp::displayOptimizationResult(const UberOptimizer::OptimizationResult& result) {
    if (result.success) {
        std::cout << "\n--- Optimization Results ---" << std::endl;
        std::cout << "Level: " << result.level_used << std::endl;
        std::cout << "Optimized XOR operations: " << result.xor_count << std::endl;
        std::cout << "Direct calculation XORs: " << result.direct_xor_count << std::endl;
        std::cout << "Optimization savings: " << (result.direct_xor_count - result.xor_count) << " XORs" << std::endl;
        std::cout << "Efficiency ratio: " << std::fixed << std::setprecision(2) << result.ratio << "x" << std::endl;
        
        double percentage_saved = ((double)(result.direct_xor_count - result.xor_count) / result.direct_xor_count) * 100.0;
        std::cout << "Percentage saved: " << std::fixed << std::setprecision(1) << percentage_saved << "%" << std::endl;
    } else {
        std::cout << "Failed to calculate XOR count." << std::endl;
    }
}

void ErasureCodingApp::processRingMapping() {
    std::cout << "\n=== Ring Mapping ===" << std::endl;
    
    // Perform ring mapping
    m_ring_mapping_result = m_ring_mapper->mapMatrix(m_encoding_matrix, Config::MATRIX_SIZE, Config::MATRIX_SIZE);
    
    if (m_ring_mapping_result.success) {
        // Print mapping statistics
        m_ring_mapper->printMappingStatistics(m_ring_mapping_result);
        
        // Generate mixed bitmatrix for A' (GF2^10) × x (GF2^8) multiplication
        m_ring_bitmatrix = m_ring_mapper->generateMixedBitmatrix(m_ring_mapping_result.mapped_matrix, 
                                                               Config::MATRIX_SIZE, Config::MATRIX_SIZE);
        
        // Calculate XOR count for inverse mapping using Uber optimization
        int inverse_rows = Config::FIELD_SIZE;     // 8 rows 
        int inverse_cols = Config::RING_SIZE;      // 10 cols
        
        // Convert inverse mapping matrix to int vector for Uber (it's already 0/1)
        std::vector<int> inverse_bitmatrix(inverse_rows * inverse_cols);
        for (int i = 0; i < inverse_rows; ++i) {
            for (int j = 0; j < inverse_cols; ++j) {
                inverse_bitmatrix[i * inverse_cols + j] = 
                    static_cast<int>(m_ring_mapping_result.inverse_mapping_matrix[i * inverse_cols + j]);
            }
        }
        
        // Save as bitmatrix format for Uber (reuse existing filename with .txt extension)
        std::string inverse_bitmatrix_file = "inverse_mapping_matrix_uber.txt";
        FileManager::writeBitmatrixToUberFile(inverse_bitmatrix, inverse_rows, inverse_cols, inverse_bitmatrix_file);
        
        // Use Uber to optimize inverse mapping with level 3 (as requested)
        auto inverse_optimization_result = m_optimizer->optimize(inverse_bitmatrix_file, 3);
        
        int inverse_xor_count = 0;
        if (inverse_optimization_result.success) {
            inverse_xor_count = inverse_optimization_result.xor_count;
            std::cout << "Inverse mapping optimization (level 3): " << inverse_xor_count << " XORs" << std::endl;
        } else {
            // Fallback to direct calculation if Uber fails
            std::cout << "Uber optimization failed for inverse mapping, using direct calculation" << std::endl;
            for (int i = 0; i < inverse_rows; ++i) {
                int ones_in_row = 0;
                for (int j = 0; j < inverse_cols; ++j) {
                    if (m_ring_mapping_result.inverse_mapping_matrix[i * inverse_cols + j] == 1) {
                        ones_in_row++;
                    }
                }
                if (ones_in_row > 1) {
                    inverse_xor_count += (ones_in_row - 1);
                }
            }
        }
        
        int total_inverse_xors = inverse_xor_count * Config::MATRIX_SIZE;
        
        // Save ring-mapped matrices to files
        FileManager::writeMatrixToFile(m_ring_mapping_result.mapped_matrix, Config::MATRIX_SIZE, Config::MATRIX_SIZE, 
                                      Config::RING_MAPPED_MATRIX_FILE, "Ring-Mapped Encoding Matrix");
        
        int ring_bit_rows = Config::RING_SIZE * Config::MATRIX_SIZE;
        int ring_bit_cols = Config::FIELD_SIZE * Config::MATRIX_SIZE;
        FileManager::writeBitmatrixToUberFile(m_ring_bitmatrix, ring_bit_rows, ring_bit_cols, Config::RING_BITMATRIX_FILE);
        FileManager::writeBitmatrixReadable(m_ring_bitmatrix, ring_bit_rows, ring_bit_cols, 
                                           Config::RING_BITMATRIX_READABLE_FILE, Config::RING_SIZE);
        
        FileManager::writeMatrixToFile(m_ring_mapping_result.inverse_mapping_matrix, inverse_rows, inverse_cols, 
                                      Config::INVERSE_MAPPING_MATRIX_FILE, "Inverse Mapping Matrix");
        
        m_inverse_xor_count = total_inverse_xors;
        std::cout << "Ring mapping completed. Files saved." << std::endl;
    } else {
        std::cout << "Ring mapping failed!" << std::endl;
    }
}

void ErasureCodingApp::runRingMappingComparison() {
    if (!m_ring_mapping_result.success) {
        std::cout << "\nSkipping ring mapping comparison due to mapping failure." << std::endl;
        return;
    }
    
    std::cout << "\n=== Ring Mapping Comparison ===" << std::endl;
    
    // Optimize ring-mapped matrix  
    auto mapped_result = m_optimizer->optimize(Config::RING_BITMATRIX_FILE, 2);
    
    // Display comparison
    displayRingMappingComparisonNew(m_original_optimization_result, mapped_result, m_inverse_xor_count);
}

void ErasureCodingApp::displayRingMappingComparison(const UberOptimizer::OptimizationResult& original, 
                                                   const UberOptimizer::OptimizationResult& mapped,
                                                   const UberOptimizer::OptimizationResult& inverse) {
    std::cout << "\n--- Ring Mapping Comparison Results ---" << std::endl;
    
    if (original.success && mapped.success && inverse.success) {
        std::cout << "Original matrix XORs: " << original.xor_count << std::endl;
        std::cout << "Ring-mapped matrix XORs: " << mapped.xor_count << std::endl;
        std::cout << "Inverse mapping XORs: " << inverse.xor_count << std::endl;
        
        int total_ring_xors = mapped.xor_count + inverse.xor_count;
        std::cout << "Total ring mapping XORs: " << total_ring_xors << std::endl;
        
        int xor_improvement = original.xor_count - total_ring_xors;
        std::cout << "XOR improvement: " << xor_improvement << " (" << 
                     (xor_improvement > 0 ? "BETTER" : "WORSE") << ")" << std::endl;
        
        if (original.xor_count > 0) {
            double percentage_improvement = (double)xor_improvement / original.xor_count * 100.0;
            std::cout << "Percentage improvement: " << std::fixed << std::setprecision(1) << 
                         percentage_improvement << "%" << std::endl;
        }
        
        if (xor_improvement > 0) {
            std::cout << "✓ Ring mapping provides XOR optimization!" << std::endl;
        } else {
            std::cout << "✗ Ring mapping does not provide XOR optimization." << std::endl;
        }
    } else {
        std::cout << "Some optimizations failed." << std::endl;
    }
}

void ErasureCodingApp::displayRingMappingComparisonNew(const UberOptimizer::OptimizationResult& original, 
                                                      const UberOptimizer::OptimizationResult& mapped,
                                                      int inverse_xor_count) {
    std::cout << "\n--- Ring Mapping Comparison Results ---" << std::endl;
    
    if (original.success && mapped.success) {
        std::cout << "Original matrix XORs: " << original.xor_count << std::endl;
        std::cout << "Ring-mapped matrix XORs: " << mapped.xor_count << std::endl;
        std::cout << "Inverse mapping XORs: " << inverse_xor_count << std::endl;
        
        int total_ring_xors = mapped.xor_count + inverse_xor_count;
        std::cout << "Total ring mapping XORs: " << total_ring_xors << std::endl;
        
        int xor_improvement = original.xor_count - total_ring_xors;
        std::cout << "XOR improvement: " << xor_improvement << " (" << 
                     (xor_improvement > 0 ? "BETTER" : "WORSE") << ")" << std::endl;
        
        if (original.xor_count > 0) {
            double percentage_improvement = (double)xor_improvement / original.xor_count * 100.0;
            std::cout << "Percentage improvement: " << std::fixed << std::setprecision(1) << 
                         percentage_improvement << "%" << std::endl;
        }
        
        if (xor_improvement > 0) {
            std::cout << "✓ Ring mapping provides XOR optimization!" << std::endl;
        } else {
            std::cout << "✗ Ring mapping does not provide XOR optimization." << std::endl;
        }
    } else {
        std::cout << "Some optimizations failed." << std::endl;
    }
}

int main() {
    try {
        ErasureCodingApp app;
        return app.run();
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    }
}

