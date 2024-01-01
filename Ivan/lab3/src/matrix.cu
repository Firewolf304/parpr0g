#include "../include/includes.h"
class matrix_gpu;
__device__ static void getFromThread(int degree, int * counter, int n, int stringSize, int size) {
    int index_x = threadIdx.x + blockIdx.x * blockDim.x;
    int index_y = threadIdx.y + blockIdx.y * blockDim.y;
    int index_z = threadIdx.z + blockIdx.z * blockDim.z;
    int value = index_x + index_y * blockDim.y + index_z * blockDim.z;
    if(value > size ) return;
    int output = 0;
    int schet = 1;
    int count = 0;
    while(value) {
        output += schet * (value % degree);
        value /= degree;
        schet *= 10;
        count++;
    }
    if(output > 9) {
        int previousDigit = output % 10;
        output /= 10;
        while (output > 0) {
            int currentDigit = output % 10;
            output /= 10;
            if (currentDigit > previousDigit) {
                return;
            }
            previousDigit = currentDigit;
        }
    }
    atomicAdd(counter, 1);
}
__global__ void kernel(int str_size, int sub_size, int * count, int sizer) {
    getFromThread(str_size - sub_size + 1, count, sub_size, str_size, sizer);
}

class matrix_gpu {
public:
    matrix_gpu(int n, string str = "xyxxz") {
        this->str = str;
        this->n = n;
    }
    string str;
    int n = 3;
    thrust::device_vector<string> map_gpu;
    void get_number() {
        if(this->str.length() < this->n) {return;}
        if(this->n == 0) {cout << "Result GPU: "<< this->n << endl; return;}
        int * pointer;
        unsigned long long count = 0;
        cudaMalloc(&pointer, sizeof(count));
        cudaMemcpy(pointer, &count, sizeof(count), cudaMemcpyHostToDevice);
        int size = pow(this->str.length() - this->n + 1, this->n);
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, 0);
        int threads_per_block = std::pow(2, str.length() - this->n );
        int blocks_per_grid = (size + threads_per_block - 1) / threads_per_block;
        if (blocks_per_grid > deviceProp.maxThreadsPerBlock) {
            blocks_per_grid = deviceProp.maxThreadsPerBlock;
            threads_per_block = (size + blocks_per_grid - 1) / blocks_per_grid;
        }
        kernel  <<<dim3(blocks_per_grid, 1, 1), dim3( threads_per_block, 1, 1 )>>> (this->str.length(), this->n, pointer, size);
        cudaDeviceSynchronize();
        cudaMemcpy(&count, pointer, sizeof(count), cudaMemcpyDeviceToHost);
        cudaFree(pointer);
        if(count == 0) { // sorry, but for large numbers
            count = tgamma( (double)(this->str.length() + 1) ) / (tgamma( (double)(n + 1) ) * tgamma(  (double)(this->str.length() - n + 1) ) );
        }
        cout << "Result GPU: "<< count << endl;
    }
};
class matrix_cpu {
public:
    matrix_cpu(int n, string str = "xyxxz") {
        this->n = n;
        this->str=str;
    }
    std::vector<string> map;
    int n = 3;
    string str;
    void get_number() {
        bool skep = false;
        if(this->str.length() == this->n) { cout << "Result CPU: "<< 1 << endl; skep = true;}
        if(this->str.length() < this->n) {skep = true;}

        int count = 0;
        if(this->str.length() == this->n) {}
        std::vector<int> offset (this->n,0);
        std::vector<int> iter (this->n,0);
        int maxSIZE = this->str.length() - this->n;
        string line;
        for(int i = 0; i < this->n; i++) { line.push_back(' '); }
        for(int i = 0; i < this->n; i++) {offset[i] = i; line[i] = str[i]; }
        int c = 0;
        this->map.push_back("");
        for(int i = 0; i < this->n; i++ ) {
            this->map.back() += this->str[offset[i] + iter[i]];
        }
        count++;
        for(;c != this->n && !skep;) {
            if(iter.back() + 1 <= maxSIZE) {
                iter.back()++;
            }
            else {
                for (int i = 0; i < n; i++) {
                    if (iter[i] >= maxSIZE) {
                        if(iter[i - 1] + 1 <= maxSIZE) {
                            iter[i - 1]++;
                        }
                        if (iter[i - 1] <= maxSIZE) {
                            iter[i] = iter[i - 1];
                        }
                        for (int j = i; j < n; j++) {
                            iter[j] = iter[i];
                        }
                    }
                }
            }
            this->map.push_back("");
            for(int i = 0; i < this->n; i++ )
                this->map.back() += this->str[offset[i] + iter[i]];
            c = 0;
            for(auto d : iter) {
                if(d == maxSIZE) {
                    c++;
                }
            }
            count++;
        }
        cout << "Result CPU: "<< count << endl;
    }
};
