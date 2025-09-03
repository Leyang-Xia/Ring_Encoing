#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstdio>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <random>
#include <cstring>

#ifdef _WIN32
#define WIFEXITED(status) (((status) & 0x7f) == 0)
#define WEXITSTATUS(status) (((status) & 0xff00) >> 8)
#else
#include <sys/wait.h>
#endif

// 优化结果结构
struct OptimizationResult {
    std::string technique;
    int xor_count;
    int savings;
    double percentage_saved;
    bool success;
};

// 加载矩阵文件
std::vector<std::vector<int>> loadMatrix(const std::string& filename) {
    std::vector<std::vector<int>> matrix;
    std::ifstream file(filename);
    
    if (!file) {
        throw std::runtime_error("无法打开文件: " + filename);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::vector<int> row;
        
        // 解析行：支持空格分隔或连续格式
        if (line.find(' ') != std::string::npos) {
            std::istringstream iss(line);
            std::string token;
            while (iss >> token) {
                if (token == "0" || token == "1") {
                    row.push_back(std::stoi(token));
                }
            }
        } else {
            for (char c : line) {
                if (c == '0' || c == '1') {
                    row.push_back(c - '0');
                }
            }
        }
        
        if (!row.empty()) {
            matrix.push_back(row);
        }
    }
    
    if (matrix.empty()) {
        throw std::runtime_error("矩阵文件为空");
    }
    
    return matrix;
}

// 生成300x200测试矩阵
std::vector<std::vector<int>> generateTestMatrix() {
    std::vector<std::vector<int>> matrix(300, std::vector<int>(200, 0));
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> ones_dis(2, 8);
    
    for (int i = 0; i < 300; ++i) {
        int num_ones = ones_dis(gen);
        std::vector<int> positions(200);
        std::iota(positions.begin(), positions.end(), 0);
        std::shuffle(positions.begin(), positions.end(), gen);
        
        for (int j = 0; j < num_ones && j < 200; ++j) {
            matrix[i][positions[j]] = 1;
        }
    }
    
    return matrix;
}

// 保存矩阵到文件
void saveMatrix(const std::vector<std::vector<int>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("无法写入文件: " + filename);
    }
    
    for (const auto& row : matrix) {
        for (int bit : row) {
            file << bit;
        }
        file << std::endl;
    }
}

// 计算直接XOR次数
int calculateDirectXORs(const std::vector<std::vector<int>>& matrix) {
    int total_xors = 0;
    
    for (const auto& row : matrix) {
        int ones_count = 0;
        for (int bit : row) {
            if (bit == 1) ones_count++;
        }
        if (ones_count > 1) {
            total_xors += (ones_count - 1);
        }
    }
    
    return total_xors;
}

// 运行X-Sets优化
OptimizationResult runXSets(const std::string& matrix_file, const std::string& technique) {
    OptimizationResult result;
    result.technique = technique;
    result.success = false;
    
    std::string command = "./X-Sets 10 10 100 " + technique + " < " + matrix_file;
    
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
        return result;
    }
    
    std::string output;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        output += buffer;
    }
    
    pclose(pipe);
    
    // 解析XOR次数
    size_t pos = output.rfind("Total XOR operations: ");
    if (pos != std::string::npos) {
        try {
            std::string substr = output.substr(pos + 22);
            size_t newline = substr.find('\n');
            if (newline != std::string::npos) {
                substr = substr.substr(0, newline);
            }
            result.xor_count = std::stoi(substr);
            result.success = true;
        } catch (...) {
            // 解析失败
        }
    }
    
    return result;
}

// 主分析函数
void analyzeMatrix(const std::string& matrix_file) {
    try {
        // 加载矩阵
        auto matrix = loadMatrix(matrix_file);
        std::cout << "矩阵大小: " << matrix.size() << " x " << matrix[0].size() << std::endl;
        
        // 计算直接XOR次数
        int direct_xors = calculateDirectXORs(matrix);
        std::cout << "直接计算: " << direct_xors << " XORs" << std::endl;
        
        // 运行X-Sets优化
        std::vector<std::string> techniques = {"MW", "MW_SS", "UBER_XSET"};
        std::vector<OptimizationResult> results;
        
        for (const auto& tech : techniques) {
            std::cout << "运行 " << tech << "..." << std::endl;
            auto result = runXSets(matrix_file, tech);
            if (result.success) {
                result.savings = direct_xors - result.xor_count;
                result.percentage_saved = ((double)result.savings / direct_xors) * 100.0;
                results.push_back(result);
            }
        }
        
        // 显示结果
        std::cout << "\n--- 优化结果 ---" << std::endl;
        for (const auto& result : results) {
            std::cout << result.technique << ": " << result.xor_count 
                      << " XORs (节省 " << result.savings << ", " 
                      << std::fixed << std::setprecision(1) << result.percentage_saved << "%)" << std::endl;
        }
        
        if (!results.empty()) {
            auto best = *std::min_element(results.begin(), results.end(),
                [](const OptimizationResult& a, const OptimizationResult& b) {
                    return a.xor_count < b.xor_count;
                });
            std::cout << "\n最佳: " << best.technique << " (节省 " << best.savings << " XORs)" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        // 生成300x200矩阵并分析
        std::cout << "=== 矩阵XOR分析 ===" << std::endl;
        std::cout << "生成300x200测试矩阵..." << std::endl;
        
        auto matrix = generateTestMatrix();
        std::string filename = "matrix_300x200.txt";
        saveMatrix(matrix, filename);
        std::cout << "矩阵大小: 300 x 200" << std::endl;
        
        analyzeMatrix(filename);
        
    } else if (argc == 2) {
        // 分析指定文件
        std::string matrix_file = argv[1];
        
        std::ifstream file(matrix_file);
        if (!file) {
            std::cerr << "错误: 文件不存在 '" << matrix_file << "'" << std::endl;
            return 1;
        }
        file.close();
        
        std::cout << "=== 矩阵XOR分析 ===" << std::endl;
        analyzeMatrix(matrix_file);
        
    } else {
        std::cout << "用法:" << std::endl;
        std::cout << "  " << argv[0] << "                    # 生成300x200矩阵并分析" << std::endl;
        std::cout << "  " << argv[0] << " <matrix_file>      # 分析指定矩阵文件" << std::endl;
        return 1;
    }
    
    return 0;
}
