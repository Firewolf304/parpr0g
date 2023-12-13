#include <iostream>
#include <string>
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include <stack>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <cmath>
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
;
//string str = "abcde";

int count = 0;
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
    std::vector<int> offset (n,0);
    std::vector<int> iter (n,0);
    int maxSIZE = str.length() - n;
    string line(n,' ');
    for(int i = 0; i < n; i++) {offset[i] = i; line[i] = str[i]; }
    int c = 0;
    map_cpu.push_back("");
    for(int i = 0; i < n; i++ ) {
        map_cpu.back() += str[offset[i] + iter[i]];
    }
    count++;
    while(c != n) {
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
}

__device__ void get_values(int stepen, int n, int str_size) {

    /*int value = 2;
    int schet = 1;
    thrust::device_vector<int> & H;
    H.push_back()
    while(value) {
        output += schet * (value % stepen);
        value /= stepen;
        schet *= 10;
    }*/
}
__global__ void kernel(thrust::device_ptr<int> count, int str_size, int sub_size) {
    //get_values()
    printf("HELLO\n");
}
void cycle_gpu(string str = "xyxxz", int n = 3) {
    int count = 0;
    thrust::device_ptr<int> ptr = thrust::device_pointer_cast<int>(&count);
    cudaMemcpy(ptr.get(), &count, sizeof(count), cudaMemcpyHostToDevice);
    kernel  <<<
        dim3(1,1,1),
        dim3(pow(n - str.length(), n),1, 1 )
    >>>( ptr , str.length(), n);
    //cudaMemcpy(ptr, thrust::raw_pointer_cast(map_gpu.data()), map_gpu.size() * sizeof(string), cudaMemcpyHostToDevice);
    cudaMemcpy(&count, ptr.get(), sizeof(count), cudaMemcpyDeviceToHost);
    cudaFree(ptr.get());
    cudaDeviceSynchronize();
}


int main() {
    string piece = "";
    //cycle_cpu("",-1);
    cycle_gpu();
    for(auto d : map_cpu) {
        cout << d << endl;
    }
    /*for(int i = 0; i <= str.length() - n; i ++){
        piece = str[i];
        if (piece.length() >= n && int(piece[piece.length() - 1]) != 0){
            cout<<piece<<endl;
            count +=1;
        }
        else cycle(piece, i);
    }*/
    cout<<count;
}
