#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

// Jerasure headers
#include "jerasure.h"
#include "reed_sol.h"
#include "galois.h"

typedef struct {
    int k;                  // Number of data symbols
    int m;                  // Number of parity symbols
    int w;                  // Galois field size (2^w)
    int* matrix;            // Encoding matrix
    int* bitmatrix;         // Binary matrix
} RSAnalyzer;

// Initialize RS analyzer
RSAnalyzer* create_analyzer(int k, int m, int w) {
    RSAnalyzer* analyzer = malloc(sizeof(RSAnalyzer));
    if (!analyzer) return NULL;
    
    analyzer->k = k;
    analyzer->m = m;
    analyzer->w = w;
    
    // Create custom parity check matrix P (m×k)
    // For RS(6,4): P = [1 1 1 1; 1 2 3 4]
    analyzer->matrix = malloc(k * m * sizeof(int));
    if (!analyzer->matrix) {
        free(analyzer);
        return NULL;
    }
    
    if (k == 4 && m == 2) {
        // Custom P matrix: [1 1 1 1; 1 2 3 4]
        // Row 0: [1, 1, 1, 1]
        analyzer->matrix[0] = 1; analyzer->matrix[1] = 1; 
        analyzer->matrix[2] = 1; analyzer->matrix[3] = 1;
        // Row 1: [1, 2, 3, 4]
        analyzer->matrix[4] = 1; analyzer->matrix[5] = 2; 
        analyzer->matrix[6] = 3; analyzer->matrix[7] = 4;
    } else {
        // For other configurations, use default pattern
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < k; j++) {
                if (i == 0) {
                    analyzer->matrix[i * k + j] = 1;  // First row: all 1s
                } else {
                    analyzer->matrix[i * k + j] = j + 1;  // Pattern: 1,2,3,4,...
                }
            }
        }
    }
    
    // Convert to bitmatrix using Jerasure function
    analyzer->bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, analyzer->matrix);
    if (!analyzer->bitmatrix) {
        free(analyzer->matrix);
        free(analyzer);
        return NULL;
    }
    
    printf("Created RS(%d, %d) with custom parity matrix P over GF(2^%d)\n", k + m, k, w);
    return analyzer;
}

// Free analyzer memory
void free_analyzer(RSAnalyzer* analyzer) {
    if (analyzer) {
        if (analyzer->matrix) free(analyzer->matrix);
        if (analyzer->bitmatrix) free(analyzer->bitmatrix);
        free(analyzer);
    }
}

// Print the complete generator matrix G = [I | P]
void print_matrix(RSAnalyzer* analyzer) {
    printf("\nParity Check Matrix P (GF(2^%d), %d×%d):\n", analyzer->w, analyzer->m, analyzer->k);
    for (int i = 0; i < analyzer->m; i++) {
        for (int j = 0; j < analyzer->k; j++) {
            printf("%3x ", analyzer->matrix[i * analyzer->k + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Save bitmatrix for Uber tool
int save_bitmatrix(RSAnalyzer* analyzer, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Error: Cannot open file %s for writing\n", filename);
        return 0;
    }
    
    int rows = analyzer->m * analyzer->w;
    int cols = analyzer->k * analyzer->w;
    
    // Write bitmatrix in format expected by Uber (rows of 0s and 1s)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(file, "%d", analyzer->bitmatrix[i * cols + j]);
            if (j < cols - 1) fprintf(file, " ");
        }
        fprintf(file, "\n");
    }
    
    fclose(file);
    printf("Bitmatrix saved to %s\n", filename);
    return 1;
}

// Parse Uber output to get XOR count
int parse_uber_output(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error: Cannot open Uber output file %s\n", filename);
        return -1;
    }
    
    char line[256];
    int xor_count = -1;
    
    // Look for "Total XOR operations:" line
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "Total XOR operations:")) {
            char* colon = strchr(line, ':');
            if (colon) {
                xor_count = atoi(colon + 1);
                break;
            }
        }
    }
    
    fclose(file);
    return xor_count;
}

// Run Uber optimization
int run_uber_optimization(RSAnalyzer* analyzer, const char* matrix_file) {
    char command[512];
    char output_file[256];
    int best_xors = INT_MAX;
    int modes[] = {3};  // I 3
    int best_mode = 3;
    
    for (int i = 0; i < 1; i++) {
        int L = modes[i];
        
        // Create command and output filename
        snprintf(output_file, sizeof(output_file), "uber_output_%d_%d_I%d.txt", 
                analyzer->k, analyzer->m, L);
        snprintf(command, sizeof(command), "./Uber I %d < %s > %s 2>/dev/null", L, matrix_file, output_file);
        
        // Run Uber
        int result = system(command);
        if (result != 0) continue;
        
        // Parse output
        int xors = parse_uber_output(output_file);
        printf("Uber optimized XOR count: %d (I %d)\n", xors, L);
        return xors;
    }
    
    if (best_xors < INT_MAX) {
        printf("Uber optimized XOR count: %d (I %d)\n", best_xors, best_mode);
        return best_xors;
    }
    
    return -1;
}



// Analyze one configuration
void analyze_configuration(int k, int m) {
    printf("\n=== Analysis for RS(%d, %d) ===\n", k + m, k);
    
    RSAnalyzer* analyzer = create_analyzer(k, m, 8);
    if (!analyzer) {
        printf("Error: Failed to create analyzer\n");
        return;
    }
    
    print_matrix(analyzer);
    
    // Save bitmatrix and run Uber optimization
    char matrix_filename[64];
    snprintf(matrix_filename, sizeof(matrix_filename), "matrix_%d_%d.txt", k, m);
    
    if (save_bitmatrix(analyzer, matrix_filename)) {
        int uber_xors = run_uber_optimization(analyzer, matrix_filename);
        
        if (uber_xors >= 0) {
            // Calculate XORs per bit using formula: XOR次数/8/(k-1)
            double uber_per_bit = (double)uber_xors / (8.0 * (k-1));
            printf("Uber XORs per bit: %.3f\n", uber_per_bit);
        }
    }
    
    free_analyzer(analyzer);
}

int main() {
    printf("Reed-Solomon XOR Analysis with Custom Generator Matrix\n");
    printf("Using custom G = [I | P] where P = [1 1 1 1; 1 2 3 4] for RS(6,4)\n");
    printf("GF(2^8) with irreducible polynomial x^8 + x^4 + x^3 + x^2 + 1 (0x11D)\n");
    printf("Optimized using Uber I 3\n");
    printf("=================================================================\n");
    
    // Test configurations: k+2 for k = 4,5,6,7
    int k_values[] = {4, 5, 6, 7};
    int num_configs = sizeof(k_values) / sizeof(k_values[0]);
    int m = 2; // Always 2 parity symbols
    
    // Analyze each configuration
    for (int i = 0; i < num_configs; i++) {
        analyze_configuration(k_values[i], m);
    }
    
    printf("\nAnalysis complete. All results shown above.\n");
    
    return 0;
} 
