#include <iostream>
#include <string>
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include <stack>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
//#include <cmath>
#include <math.h>
#include <unistd.h>

#define NS_PER_SECOND 100000000000000
/*
        Задание 17
 Найти количество всех различных последовательностей символов заданного размера в
 заданной строке с учетом порядка следования символов (например, строка xyxxz содержит
 подстроки:
        x (3), y (1), z (1) размера 1,
        xy (1), xx (3), xz (3), yx (2), yz (1) размера 2,
        xyx (2) , xyz (1), xxx (1), xxz (3), yxx (1), yxz(2) размера 3,
        xyxx (1), xyxz (2), xxxz (1), yxxz (1) размера 4,
        xyxxz (1) размера 5;
 в скобках указано количество подстрок; других подстрок нет; таким образом, подстрок размера 1 ровно 5, размера 2 – 10, размера 3 – 10 и т.д.)
 */

using std::string;
using std::cout;
using std::cin;
using std::endl;
void sub_timespec(struct timespec t1, struct timespec t2, struct timespec *td)
{
    td->tv_nsec = t2.tv_nsec - t1.tv_nsec;
    td->tv_sec  = t2.tv_sec - t1.tv_sec;
    if (td->tv_sec > 0 && td->tv_nsec < 0)
    {
        td->tv_nsec += NS_PER_SECOND;
        td->tv_sec--;
    }
    else if (td->tv_sec < 0 && td->tv_nsec > 0)
    {
        td->tv_nsec -= NS_PER_SECOND;
        td->tv_sec++;
    }
}
//string str = "abcde";

std::vector<string> map_cpu;
thrust::device_vector<string> map_gpu;

/*int cycle_cpu(string piece, int a){
    string save = piece;
    for (int i = a + 1; i < str.length(); i++) {
        piece += str[i];
        if (piece.length() >= n ){
            cout<<piece<<endl;
            count +=1;
        }
        else cycle_cpu(piece, i);
        piece = save;
    }
}*/
void cycle_cpu(string str = "xyxxz", int n = 3) {
    bool skep = false;
    if(str.length() == n) { cout << "CPU: "<< 1 << endl; skep = true;}
    if(str.length() < n) {skep = true;}

    int count = 0;
    if(str.length() == n) {}
    std::vector<int> offset (n,0);
    std::vector<int> iter (n,0);
    int maxSIZE = str.length() - n;
    string line;
    for(int i = 0; i < n; i++) { line.push_back(' '); }
    for(int i = 0; i < n; i++) {offset[i] = i; line[i] = str[i]; }
    int c = 0;
    map_cpu.push_back("");
    for(int i = 0; i < n; i++ ) {
        map_cpu.back() += str[offset[i] + iter[i]];
    }
    count++;
    for(;c != n && !skep;) {
        if(iter.back() + 1 <= maxSIZE) {
            iter.back()++;
        }
        else { // остальной участок просто разбить на потоки и длины так, чтобы каждый элемент двигался от 0 д maxSIZE
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
        map_cpu.push_back("");
        for(int i = 0; i < n; i++ ) {
            map_cpu.back() += str[offset[i] + iter[i]];
        }
        c = 0;
        for(auto d : iter) {
            if(d == maxSIZE) {
                c++;
            }
        }
        count++;
    }
    /*for(auto d : map_cpu) {
        cout << d << endl;
    }*/
    cout << "CPU: "<< count << endl;
}
__device__ class debugger {
public:
    __device__ void debug(int value, const char* message) {
        printf("%d:%d:%d => %d => %s\n", threadIdx.x,threadIdx.y, threadIdx.z, value, message);
    }
    __device__ void debug(int value) {
        printf("%d:%d:%d => %d\n", threadIdx.x,threadIdx.y, threadIdx.z, value);
    }
    __device__ void debug(const char* message) {
        printf("%d:%d:%d => %s\n", threadIdx.x,threadIdx.y, threadIdx.z, message);
    }


    __device__ debugger operator<<(const char* message) {
        debug(message);
    }
    __device__ debugger operator<<(const int value) {
        debug(value);
    }
};
__device__ debugger debug;
__device__ __host__ int factorial(int n)
{
    return (n==1 || n==0) ? 1: n * factorial(n - 1);
}

__device__ void get_values(int stepen, int * counter, int n, int str_size, int sizer) {
    //int value = threadIdx.x + blockIdx.x * blockDim.x; // very small
    int index_x = threadIdx.x + blockIdx.x * blockDim.x;
    int index_y = threadIdx.y + blockIdx.y * blockDim.y;
    int index_z = threadIdx.z + blockIdx.z * blockDim.z;
    //int value = index_x + index_y * gridDim.x * blockDim.x + index_z * gridDim.x * blockDim.x * blockDim.y * gridDim.y;
    //int value = index_x + index_y * blockDim.x + index_z * blockDim.x * blockDim.y;
    int value = index_x + index_y * blockDim.y + index_z * blockDim.z;
    int save = value;
    if(value > sizer ) return;
    int output = 0;
    int schet = 1;
    int count = 0;

    while(value) {
        output += schet * (value % stepen);
        value /= stepen;
        schet *= 10;
        count++;
    }
    int copy = output;
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
    //printf("X:%d + Y:%d * Z:%d = %d | stepen=%d out=%d value=%d %d\n", threadIdx.x,threadIdx.y, threadIdx.z, value, stepen, copy, save, sizer );
    atomicAdd(counter, 1);
}
__global__ void kernel(int str_size, int sub_size, int * count, int sizer) { //
    //get_values()
    get_values(str_size - sub_size + 1, count, sub_size, str_size, sizer);
    //int globalThreadId = threadIdx.x + blockIdx.x * blockDim.x;
    //printf("%d %d %d %d\n", threadIdx.x,threadIdx.y, threadIdx.z, globalThreadId );
}
void cycle_gpu(string str = "xyxxz", int n = 3) {
    if(str.length() < n) {return;}
    if(n == 0) {cout << "GPU: "<< n << endl; return;}
    int count = 0;

    int * ptr;
    cudaMalloc(&ptr, sizeof(int));
    cudaMemcpy(ptr, &count, sizeof(count), cudaMemcpyHostToDevice);
    int size = pow(str.length() - n + 1, n);
    int sizer = factorial( str.length() ) / (factorial( n ) * factorial(  str.length() - n ));
    /*dim3 block_dim(8, 8, 8); // Задаем размер блока (x, y, z)
    dim3 grid_dim((size + block_dim.x - 1) / block_dim.x, (n + block_dim.y - 1) / block_dim.y, (n + block_dim.z - 1) / block_dim.z); // Вычисляем размер сетки (x, y, z)
    kernel<<<grid_dim, block_dim>>>( str.length(), n, ptr);*/
    //small size
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    int threads_per_block = std::pow(2, str.length() - n );
    int blocks_per_grid = (size + threads_per_block - 1) / threads_per_block;

    // Ограничение на максимальное количество блоков
    // cout << blocks_per_grid << endl;
    if (blocks_per_grid > deviceProp.maxThreadsPerBlock) {
        blocks_per_grid = deviceProp.maxThreadsPerBlock;
        threads_per_block = (size + blocks_per_grid - 1) / blocks_per_grid;
    }
    dim3 block_size(
            threads_per_block,
            1,
            1
            );
    dim3 grid_size(blocks_per_grid, 1, 1);
    kernel  <<<grid_size, block_size>>> ( str.length(), n, ptr, size);
    cudaDeviceSynchronize();
    cudaMemcpy( &count,ptr, sizeof(int), cudaMemcpyDeviceToHost);
    cout << "GPU: "<< count << endl;
    cudaFree(ptr);
}


int main() {
    timespec start, finish, delta;
    clock_gettime(CLOCK_REALTIME, &start);
    cycle_gpu("xyxxzxxdaww", 3);
    clock_gettime(CLOCK_REALTIME, &finish);
    sub_timespec(start, finish, &delta);
    cout << "Execute time = " << delta.tv_sec << "," << delta.tv_nsec << " took seconds\n";
    clock_gettime(CLOCK_REALTIME, &start);
    cycle_cpu("xyxxzxxdaww", 3);
    clock_gettime(CLOCK_REALTIME, &finish);
    sub_timespec(start, finish, &delta);
    cout << "Execute time = " << delta.tv_sec << "," << delta.tv_nsec << " took seconds\n";
    //cycle_gpu("xyxxzxxdawwawda", 5);
    //cycle_cpu("xyxxzxxdawwawda", 5);
    /*
        string piece = "";
        for(int i = 0; i <= str.length() - n; i ++){
        piece = str[i];
        if (piece.length() >= n && int(piece[piece.length() - 1]) != 0){
            cout<<piece<<endl;
            count +=1;
        }
        else cycle(piece, i);
    }*/

}
