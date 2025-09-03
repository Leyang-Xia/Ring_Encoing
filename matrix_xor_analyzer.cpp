#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstdio>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <sstream>
#include <random>
#include <ctime>
#include <cstring>

// Windows compatibility for process status macros
#ifdef _WIN32
#define WIFEXITED(status) (((status) & 0x7f) == 0)
#define WEXITSTATUS(status) (((status) & 0xff00) >> 8)
#else
#include <sys/wait.h>
#endif

// Configuration constants for 300x200 matrix analysis
namespace Config {
    constexpr int MATRIX_ROWS = 300;
    constexpr int MATRIX_COLS = 200;
    
    const std::string BINARY_MATRIX_FILE = "matrix_300x200.txt";
    const std::string XSETS_RESULT_FILE = "xsets_optimization_result.txt";
}

/**
 * @brief Handles binary matrix generation and basic operations
 */
class BinaryMatrixHandler {
public:
    // Generate random binary matrix
    std::vector<std::vector<int>> generateRandomMatrix(int rows, int cols, double sparsity = 0.5);
    
    // Generate structured binary matrix (more realistic for coding theory)
    std::vector<std::vector<int>> generateStructuredMatrix(int rows, int cols);
    
    // Save matrix to file in X-Sets compatible format
    void saveMatrixToFile(const std::vector<std::vector<int>>& matrix, const std::string& filename);
    
    // Calculate direct XOR count
    int calculateDirectXORs(const std::vector<std::vector<int>>& matrix);
    
    // Print matrix statistics
    void printMatrixStatistics(const std::vector<std::vector<int>>& matrix);

private:
    int countOnesInRow(const std::vector<int>& row);
};

/**
 * @brief Handles X-Sets optimization operations
 */
class XSetsOptimizer {
public:
    enum class Technique {
        MW,
        MW_SS,
        UBER_XSET
    };
    
    struct OptimizationResult {
        int xor_count;
        int direct_xor_count;
        int savings;
        double efficiency_ratio;
        double percentage_saved;
        bool success;
        Technique technique_used;
        std::string technique_name;
        std::string output;
    };
    
    // Run X-Sets optimization with default technique (MW_SS)
    OptimizationResult optimize(const std::string& matrix_filename);
    
    // Run X-Sets optimization with specific technique
    OptimizationResult optimize(const std::string& matrix_filename, Technique technique);
    
    // Run multiple techniques and compare
    std::vector<OptimizationResult> optimizeAllTechniques(const std::string& matrix_filename);
    
    static std::string techniqueToString(Technique technique);

private:
    int parseXorCount(const std::string& output);
    std::string buildXSetsCommand(const std::string& matrix_filename, Technique technique);
};

/**
 * @brief Main application class for matrix XOR analysis
 */
class MatrixXORAnalyzer {
public:
    MatrixXORAnalyzer();
    int run();

private:
    void printConfiguration();
    void generateTestMatrix();
    void analyzeDirectXORs();
    void runXSetsAnalysis();
    void displayResults();
    void displayComparison(const XSetsOptimizer::OptimizationResult& result);
    void displayMultipleComparison(const std::vector<XSetsOptimizer::OptimizationResult>& results);

    std::unique_ptr<BinaryMatrixHandler> m_matrix_handler;
    std::unique_ptr<XSetsOptimizer> m_optimizer;
    
    std::vector<std::vector<int>> m_binary_matrix;
    int m_direct_xor_count;
    std::vector<XSetsOptimizer::OptimizationResult> m_optimization_results;
};

// BinaryMatrixHandler implementation
std::vector<std::vector<int>> BinaryMatrixHandler::generateRandomMatrix(int rows, int cols, double sparsity) {
    std::vector<std::vector<int>> matrix(rows, std::vector<int>(cols, 0));
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = (dis(gen) < sparsity) ? 1 : 0;
        }
    }
    
    return matrix;
}

std::vector<std::vector<int>> BinaryMatrixHandler::generateStructuredMatrix(int rows, int cols) {
    std::vector<std::vector<int>> matrix(rows, std::vector<int>(cols, 0));
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> bit_dis(0, 1);
    std::uniform_int_distribution<> ones_dis(2, 8); // Each row will have 2-8 ones
    
    for (int i = 0; i < rows; ++i) {
        int num_ones = ones_dis(gen);
        std::vector<int> positions(cols);
        std::iota(positions.begin(), positions.end(), 0);
        std::shuffle(positions.begin(), positions.end(), gen);
        
        for (int j = 0; j < num_ones && j < cols; ++j) {
            matrix[i][positions[j]] = 1;
        }
    }
    
    return matrix;
}

void BinaryMatrixHandler::saveMatrixToFile(const std::vector<std::vector<int>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    for (const auto& row : matrix) {
        for (int bit : row) {
            file << bit;
        }
        file << std::endl;
    }
}

int BinaryMatrixHandler::calculateDirectXORs(const std::vector<std::vector<int>>& matrix) {
    int total_xors = 0;
    
    for (const auto& row : matrix) {
        int ones_count = countOnesInRow(row);
        if (ones_count > 1) {
            total_xors += (ones_count - 1);
        }
    }
    
    return total_xors;
}

void BinaryMatrixHandler::printMatrixStatistics(const std::vector<std::vector<int>>& matrix) {
    int total_ones = 0;
    int total_elements = 0;
    int min_ones_per_row = Config::MATRIX_COLS;
    int max_ones_per_row = 0;
    
    for (const auto& row : matrix) {
        int ones_in_row = countOnesInRow(row);
        total_ones += ones_in_row;
        total_elements += row.size();
        min_ones_per_row = std::min(min_ones_per_row, ones_in_row);
        max_ones_per_row = std::max(max_ones_per_row, ones_in_row);
    }
    
    double sparsity = (double)total_ones / total_elements;
    
    std::cout << "\n=== 矩阵统计信息 ===" << std::endl;
    std::cout << "矩阵大小: " << matrix.size() << " x " << matrix[0].size() << std::endl;
    std::cout << "总元素数: " << total_elements << std::endl;
    std::cout << "1的总数: " << total_ones << std::endl;
    std::cout << "稀疏度: " << std::fixed << std::setprecision(3) << sparsity << std::endl;
    std::cout << "每行最少1的个数: " << min_ones_per_row << std::endl;
    std::cout << "每行最多1的个数: " << max_ones_per_row << std::endl;
    std::cout << "平均每行1的个数: " << std::fixed << std::setprecision(1) 
              << (double)total_ones / matrix.size() << std::endl;
}

int BinaryMatrixHandler::countOnesInRow(const std::vector<int>& row) {
    int count = 0;
    for (int bit : row) {
        if (bit == 1) count++;
    }
    return count;
}

// XSetsOptimizer implementation
XSetsOptimizer::OptimizationResult XSetsOptimizer::optimize(const std::string& matrix_filename) {
    return optimize(matrix_filename, Technique::MW_SS);
}

XSetsOptimizer::OptimizationResult XSetsOptimizer::optimize(const std::string& matrix_filename, Technique technique) {
    OptimizationResult result{};
    result.technique_used = technique;
    result.technique_name = techniqueToString(technique);
    
    std::string command = buildXSetsCommand(matrix_filename, technique);
    
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
        std::cerr << "错误: 无法运行X-Sets命令: " << command << std::endl;
        return result;
    }
    
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result.output += buffer;
    }
    
    int status = pclose(pipe);
    if (status == -1) {
        std::cerr << "警告: pclose()失败，命令: " << command << std::endl;
        return result;
    }
    
    if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {
        std::cerr << "警告: X-Sets命令退出状态 " << WEXITSTATUS(status) << std::endl;
    }
    
    result.xor_count = parseXorCount(result.output);
    if (result.xor_count != -1) {
        result.success = true;
    }
    
    return result;
}

std::vector<XSetsOptimizer::OptimizationResult> XSetsOptimizer::optimizeAllTechniques(const std::string& matrix_filename) {
    std::vector<OptimizationResult> results;
    std::vector<Technique> techniques = {Technique::MW, Technique::MW_SS, Technique::UBER_XSET};
    
    for (Technique tech : techniques) {
        std::cout << "正在运行 " << techniqueToString(tech) << " 技术..." << std::endl;
        OptimizationResult result = optimize(matrix_filename, tech);
        if (result.success) {
            results.push_back(result);
        } else {
            std::cout << "技术 " << techniqueToString(tech) << " 失败" << std::endl;
        }
    }
    
    return results;
}

std::string XSetsOptimizer::techniqueToString(Technique technique) {
    switch (technique) {
        case Technique::MW: return "MW";
        case Technique::MW_SS: return "MW_SS";
        case Technique::UBER_XSET: return "UBER_XSET";
        default: return "MW_SS";
    }
}

int XSetsOptimizer::parseXorCount(const std::string& output) {
    // 寻找 "Total XOR operations: " 
    size_t pos = output.rfind("Total XOR operations: ");
    if (pos != std::string::npos) {
        try {
            std::string substr = output.substr(pos + 22);
            size_t newline = substr.find('\n');
            if (newline != std::string::npos) {
                substr = substr.substr(0, newline);
            }
            return std::stoi(substr);
        } catch (...) {
            // 继续尝试其他解析方法
        }
    }
    
    // 尝试解析最后一行作为数字
    size_t last_newline = output.find_last_of('\n');
    if (last_newline != std::string::npos && last_newline > 0) {
        std::string last_line = output.substr(last_newline + 1);
        try {
            return std::stoi(last_line);
        } catch (...) {
            // 解析失败
        }
    }
    
    std::cerr << "无法从X-Sets输出中解析XOR次数" << std::endl;
    return -1;
}

std::string XSetsOptimizer::buildXSetsCommand(const std::string& matrix_filename, Technique technique) {
    // X-Sets 参数: thresh=10, nstart=10, give-up-when=100, technique
    return "./X-Sets 10 10 100 " + techniqueToString(technique) + " < " + matrix_filename;
}

// MatrixXORAnalyzer implementation
MatrixXORAnalyzer::MatrixXORAnalyzer() {
    m_matrix_handler = std::make_unique<BinaryMatrixHandler>();
    m_optimizer = std::make_unique<XSetsOptimizer>();
    m_direct_xor_count = 0;
}

int MatrixXORAnalyzer::run() {
    try {
        printConfiguration();
        generateTestMatrix();
        analyzeDirectXORs();
        runXSetsAnalysis();
        displayResults();
        
        std::cout << "\n分析完成。" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }
}

void MatrixXORAnalyzer::printConfiguration() {
    std::cout << "=== 300x200 矩阵XOR分析配置 ===" << std::endl;
    std::cout << "矩阵大小: " << Config::MATRIX_ROWS << " x " << Config::MATRIX_COLS << std::endl;
    std::cout << "分析目标: 比较直接计算与X-Sets优化的XOR次数" << std::endl;
    std::cout << "输出文件: " << Config::BINARY_MATRIX_FILE << std::endl;
}

void MatrixXORAnalyzer::generateTestMatrix() {
    std::cout << "\n=== 生成测试矩阵 ===" << std::endl;
    
    // 生成结构化矩阵（更符合实际编码理论应用）
    m_binary_matrix = m_matrix_handler->generateStructuredMatrix(Config::MATRIX_ROWS, Config::MATRIX_COLS);
    
    // 保存矩阵到文件
    m_matrix_handler->saveMatrixToFile(m_binary_matrix, Config::BINARY_MATRIX_FILE);
    
    // 显示矩阵统计信息
    m_matrix_handler->printMatrixStatistics(m_binary_matrix);
    
    std::cout << "矩阵已生成并保存到 " << Config::BINARY_MATRIX_FILE << std::endl;
}

void MatrixXORAnalyzer::analyzeDirectXORs() {
    std::cout << "\n=== 直接XOR次数分析 ===" << std::endl;
    
    m_direct_xor_count = m_matrix_handler->calculateDirectXORs(m_binary_matrix);
    
    std::cout << "直接计算XOR次数: " << m_direct_xor_count << std::endl;
    std::cout << "计算方法: 每行中1的个数减1，然后求和" << std::endl;
}

void MatrixXORAnalyzer::runXSetsAnalysis() {
    std::cout << "\n=== X-Sets优化分析 ===" << std::endl;
    
    // 运行所有可用的X-Sets技术
    m_optimization_results = m_optimizer->optimizeAllTechniques(Config::BINARY_MATRIX_FILE);
    
    // 为每个成功的结果计算统计信息
    for (auto& result : m_optimization_results) {
        if (result.success) {
            result.direct_xor_count = m_direct_xor_count;
            result.savings = result.direct_xor_count - result.xor_count;
            result.efficiency_ratio = (double)result.xor_count / result.direct_xor_count;
            result.percentage_saved = ((double)result.savings / result.direct_xor_count) * 100.0;
        }
    }
}

void MatrixXORAnalyzer::displayResults() {
    std::cout << "\n=== 分析结果 ===" << std::endl;
    
    std::cout << "\n基准结果:" << std::endl;
    std::cout << "矩阵大小: " << Config::MATRIX_ROWS << " x " << Config::MATRIX_COLS << std::endl;
    std::cout << "直接计算XOR次数: " << m_direct_xor_count << std::endl;
    
    if (m_optimization_results.empty()) {
        std::cout << "\n没有成功的优化结果。" << std::endl;
        return;
    }
    
    // 显示每个技术的结果
    for (const auto& result : m_optimization_results) {
        displayComparison(result);
    }
    
    // 显示所有技术的比较
    displayMultipleComparison(m_optimization_results);
}

void MatrixXORAnalyzer::displayComparison(const XSetsOptimizer::OptimizationResult& result) {
    std::cout << "\n--- " << result.technique_name << " 优化结果 ---" << std::endl;
    std::cout << "优化后XOR次数: " << result.xor_count << std::endl;
    std::cout << "节省的XOR次数: " << result.savings << std::endl;
    std::cout << "效率比例: " << std::fixed << std::setprecision(2) << result.efficiency_ratio << "x" << std::endl;
    std::cout << "节省百分比: " << std::fixed << std::setprecision(1) << result.percentage_saved << "%" << std::endl;
    
    if (result.savings > 0) {
        std::cout << "✓ " << result.technique_name << " 提供了优化效果!" << std::endl;
    } else {
        std::cout << "✗ " << result.technique_name << " 没有提供优化效果。" << std::endl;
    }
}

void MatrixXORAnalyzer::displayMultipleComparison(const std::vector<XSetsOptimizer::OptimizationResult>& results) {
    if (results.size() < 2) return;
    
    std::cout << "\n=== 所有技术比较 ===" << std::endl;
    
    // 找出最佳结果
    auto best_result = std::min_element(results.begin(), results.end(),
        [](const XSetsOptimizer::OptimizationResult& a, const XSetsOptimizer::OptimizationResult& b) {
            return a.xor_count < b.xor_count;
        });
    
    std::cout << "\n技术排名（按XOR次数从少到多）:" << std::endl;
    
    // 排序并显示
    std::vector<XSetsOptimizer::OptimizationResult> sorted_results = results;
    std::sort(sorted_results.begin(), sorted_results.end(),
        [](const XSetsOptimizer::OptimizationResult& a, const XSetsOptimizer::OptimizationResult& b) {
            return a.xor_count < b.xor_count;
        });
    
    for (int i = 0; i < sorted_results.size(); ++i) {
        const auto& result = sorted_results[i];
        std::cout << std::setw(2) << (i+1) << ". " 
                  << std::setw(12) << result.technique_name << ": " 
                  << std::setw(8) << result.xor_count << " XORs "
                  << "(" << std::fixed << std::setprecision(1) << result.percentage_saved << "% 节省)" << std::endl;
    }
    
    std::cout << "\n最佳技术: " << best_result->technique_name 
              << " (节省 " << best_result->savings << " XORs, " 
              << std::fixed << std::setprecision(1) << best_result->percentage_saved << "%)" << std::endl;
}

int main() {
    try {
        MatrixXORAnalyzer analyzer;
        return analyzer.run();
    } catch (const std::exception& e) {
        std::cerr << "致命错误: " << e.what() << std::endl;
        return 1;
    }
}
