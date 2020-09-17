#include "dgemm.h"
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <immintrin.h>

void dgemm(float alpha, const float *a, const float *b, float beta, float *c) {
    float alpha_[8] = {alpha,alpha,alpha,alpha,alpha,alpha,alpha,alpha};
    float c_array[8] = {0,0,0,0,0,0,0,0};
    float temp[8] = {0,0,0,0,0,0,0,0};
    float temp_sum = 0.0;
    __m256 alpha_array = _mm256_load_ps(alpha_);
    int avx_matrix_size = MATRIX_SIZE - MATRIX_SIZE%8;

    for (int32_t i = 0; i < MATRIX_SIZE; i++) {
        for (int32_t j = 0; j < MATRIX_SIZE; j++) {
            c[i * MATRIX_SIZE + j] *= beta;
            __m256 c_array_ = _mm256_load_ps(c_array);
            for (int k = 0; k < avx_matrix_size; k+=8) {
                __m256 a_array = _mm256_loadu_ps(a + i * MATRIX_SIZE + k);
                __m256 b_array = _mm256_loadu_ps(b + j * MATRIX_SIZE + k);
                __m256 alpha_a = _mm256_mul_ps(alpha_array, a_array);
                __m256 alpha_a_b = _mm256_mul_ps(alpha_a, b_array);
                c_array_ = _mm256_add_ps (c_array_, alpha_a_b);
            }
            _mm256_storeu_ps(temp, c_array_);
            for(auto& num : temp)
                temp_sum += num;
            for (int k = avx_matrix_size; k < MATRIX_SIZE; k++) {
                temp_sum += alpha * a[i * MATRIX_SIZE + k] * b[j * MATRIX_SIZE + k];
            }
            c[i * MATRIX_SIZE + j] += temp_sum;
            temp_sum=0.0f;
        }
    }
}

int main(int, char **) {
    float alpha, beta;

    // mem allocations
    int mem_size = MATRIX_SIZE * MATRIX_SIZE * sizeof(float);
//    auto a = (float *) std::aligned_alloc(32,mem_size);
//    auto b = (float *) std::aligned_alloc(32,mem_size);
//    auto c = (float *) std::aligned_alloc(32,mem_size);
    auto a = (float *) malloc(mem_size);
    auto b = (float *) malloc(mem_size);
    auto c = (float *) malloc(mem_size);

    // check if allocated
    if (nullptr == a || nullptr == b || nullptr == c) {
        printf("Memory allocation failed\n");
        if (nullptr != a) free(a);
        if (nullptr != b) free(b);
        if (nullptr != c) free(c);
        return 0;
    }

    generateProblemFromInput(alpha, a, b, beta, c);

    std::cerr << "Launching dgemm step." << std::endl;
    dgemm(alpha, a, b, beta, c);

    outputSolution(c);

    free(a);
    free(b);
    free(c);
    return 0;
}
