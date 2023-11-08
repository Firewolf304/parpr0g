#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_device_runtime_api.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <functional>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <omp.h>
//#include "matrixer.h"
using std::cout;
using std::cin;
using std::endl;
using std::vector;
class make_cube;
__global__ void MakeMat_byKernel(make_cube * obj , thrust::device_ptr<float> data);
__global__ void collect_byKernel(make_cube * obj , thrust::device_ptr<float> data, float* output, int x, int y, int z);
struct dataId {
    struct th {
        unsigned int x;
        unsigned int y;
        unsigned int z;
    } thread;
    struct bl {
        unsigned int x;
        unsigned int y;
        unsigned int z;
    } block;
    struct blDim {
        unsigned int x;
        unsigned int y;
        unsigned int z;
    } blockDim;
};
/*__global__ void test(thrust::device_ptr<int> H) {

    int a = 1 + 1;
    (H[0][0]) = a;
}
int main() {

    //std::cout << "Hello, World!" << std::endl;
    int N = 5;
    thrust::device_vector<int> H(N);

    test<<<1, 1>>>( thrust::device_pointer_cast(H.data()) );
    cout << H[0] << endl;
    return 0;
}
*/



/*
Нынешнее представление:
                    +--------+                   представление координат
                   /|       /|                    +------------------+
                  / |      / |                    |  y               |
Координата y -> y+--------+  |                    |  ^    z          |
                 |  |     |  |                    |  |   ^           |
  Координата z ->| z+-----|--+                    |  |  /            |
                 | /      | /                     |  | /             |
                 |/       |/                      |  |/              |
                 0--------+x    <- координата x   |  0----------> x  |
              (0;0;0)                             +------------------+
              радиус цилиндра считается по стороне x z от координат центра
              T1 находится на стороне YZ при Z=0
              T2 находится на обратной стороне YZ при Z=len(Z)-1
     */

    class cpu {
        typedef vector<vector<vector<double>>> matrix_type;
    public:

        cpu(int x, int y, int z){
            this->x = x;
            this->y = y;
            this->z = z;
        }
        cpu(int x, int y, int z, double cylinder_radius) : cpu(x,y,z) {
            /*if(cylinder_radius > (float)(x/2) && cylinder_radius > (float)(y/2) && cylinder_radius > (float)(z/2))  {
                throw std::runtime_error("Cylinder is out of cube");
            }*/
            this->radius = cylinder_radius;
            matrix = matrix_type(x, vector<vector<double>>(y,vector<double>(z, 0)));
        }
        ~cpu() {
            this->matrix.clear();
        }
        std::string file_name = "cpumatrix";
        int x=2,y=2,z=2, radius = 2;
        double T1 = 20.0f;              // температура 1 стороны
        double T2 = 10.0f;              // температура противоположной стороны
        double T_bottom = 0.0f;         // температура нижней грани
        double alpha = 0.05351f;           // коэф теплопроводности
        bool show_iter = false;
        /*-----matrix-----*/
        matrix_type matrix;
        void computeMat() {
            for(int i = 0; i < this->x; i++) {
                for(int j = 0; j < this->y; j++) {
                    for(int k = 0; k < this->z; k++) {
                        if(k == 0) {
                            this->matrix[i][j][k] = this->T1;
                        }
                        else if(k == this->z-1) {
                            this->matrix[i][j][k] = this->T2;
                        }
                        /*else {
                            if(std::sqrt( std::pow( i - (this->x / 2) , 2 ) + std::pow( j - (this->y / 2) , 2 ) + std::pow( k - (this->z / 2) , 2 ) ) )
                        }*/
                    }
                }
            }
        }
        void genValue() {
            srand(time(NULL));
            this->T1 = 100 - (rand() % 50);
            this->T2 = 100 - (rand() % 50);
        }
        void iteration (double eps = 0.0001f, int maxIteration=10000) { // using basic matrix iterator

            double norm = 1;
            int count = 0;
            while (norm > eps && count < maxIteration) {
                // перебор каждой точки
                for(int i = 0; i < this->x; i++) {
                    for (int j = 0; j < this->y; j++) {
                        for (int k = 0; k < this->z; k++) {

                            // просмотр соседей и сумма новой температуры
                            double oldTemp = this->matrix[i][j][k];
                            double newTemp = 0.0f;
                            if(this->show_iter){ cout << "newtemp = "; }
                            for(int x1 = ((i>0) ? i-1 : i) ; x1 <= i+1 && x1 < this->x; x1++) {
                                for(int y1 = ((j>0) ? j-1 : j); y1 <= j+1 && y1 < this->y; y1++) {
                                    for(int z1 = ((k>0) ? k-1 : k); z1 <= k+1 && z1 < this->z; z1++) {

                                        if (this->radius < (float) std::sqrt(std::pow(x1 - (float) (this->x / 2), 2) +
                                                                             std::pow(z1 - (float) (this->z / 2),
                                                                                      2)) /*&& std::abs(y1 - (this->y / 2)) <= (this->y / 2)*/) {
                                            if(this->show_iter){cout << " + " << this->matrix[x1][y1][z1];}
                                            newTemp += this->matrix[x1][y1][z1];
                                        }

                                    }
                                }
                            }
                            //newTemp *= this->alpha; // формула распределения: dT/dt = alpha * ( (d^2T/dx^2) + (d^2T/dy^2) + (d^2T/dz^2) )
                            // т.е. нужно взять от зависимой точки сумму всех соседей и умножить на коэф распределения

                            if(this->show_iter) {cout << endl;}
                            this->matrix[i][j][k] = newTemp * this->alpha;
                            if(this->show_iter) {cout << "newtemp=" << newTemp << " matrix=" << this->matrix[i][j][k] << " (" << i << "," << j << "," << k << ")" << endl;}
                            norm = std::abs( newTemp - oldTemp );
                        }
                    }
                }
                count++;
            }

        }
        void show() {
            for (int i = 0; i < this->x; i++) {
                for (int k = 0; k < z; k++) {
                    cout << this->matrix [i][y-1][k] << " ";
                }
                cout << endl;
            }
            cout << "-------------" << endl;
            for (int j = 0; j < y; j++) {
                for (int i = 0; i < x; i++) {
                    cout << this->matrix[i][j][(int)(z/2)] << " ";
                }
                cout << endl;
            }
            cout << "-------------" << endl;
            for (int j = 0; j < y; j++) {
                for (int k = 0; k < z; k++) {
                    cout << this->matrix[(int)(x/2)][j][k] << " ";
                }
                cout << endl;
            }
            cout << "-------------" << endl;
            for (int j = 0; j < y; j++) {
                for (int i = 0; i < x; i++) {
                    cout << this->matrix[i][j][(int)(z-1)] << " ";
                }
                cout << endl;
            }
            cout << "-------------" << endl;
            for (int j = 0; j < y; j++) {
                for (int i = 0; i < x; i++) {
                    cout << this->matrix[i][j][(int)(0)] << " ";
                }
                cout << endl;
            }
        }

        void visualize() {
            /*
             Пример отрисовки:
             Матрицы 3 на 3:
             1 0 1
             0 0 0
             1 0 1 все подобны

             отрисовка результат:

                 |   +---------+
                 |  / 1  0  1 /|
                 | / 0  0  0 /1|
                 |/ 1  0  1 /00|
               y +---------+101+
                 | 1  0  1 |00/
                 | 0  0  0 |1/
                 | 1  0  1 |/
        ---------+---------+x----------
                /|
               / |
              z  |




             метод итерации по оси xz :
             for (int i = 0, x_pos = z-1; i < x; i++, x_pos--) { // z-1 т.к. мы рассматриваем высоту строки по X в зависимости от координаты Z, т.е. кол-во строк равно Z, как и их высота
                for (int k = 0; k < z; k++) {
                    if(x_pos>=0 || i == 0 || i = x-1) {      // добавить пробелы, пока есть строки по xz или это начало строки
                        cout << " ";
                    }
                    else{
                        cout << answer[i][y-1][k] << " "; // y-1 т.к. смотрим только вверх
                    }
                }
                cout << endl;
            }
             */

        }
        matrix_type operator() () { // output matrix
            return this->matrix;
        }
    private:

    };

     __global__ class make_cube {
    private:
        __host__ __device__ int index(int x, int y, int z) {
            return x + this->x * (y + this->y * z);
            //return x * (this->y * this->z) + y * this->x + z;
        }

    public:
         std::string file_name = "gpumatrix";
        thrust::device_vector<float> matrix;
        int x=2,y=2,z=2, radius = 2;
        float T1 = 20.0f;              // температура 1 стороны
        float T2 = 10.0f;              // температура противоположной стороны
        float T_bottom = 0.0f;         // температура нижней грани
        float alpha = 0.05351f;           // коэф теплопроводности
        float eps = 0.01f;
        int maxIteration=50;
        bool run_test = false;                  // Ручной откладчик потоков и блоков
        bool show_iter = false;                 // Ручной откладчик итерации
        __host__ void genValue() {
            srand(time(NULL));
            this->T1 = 100 - (rand() % 50);
            this->T2 = 100 - (rand() % 50);
        }
        __host__ thrust::device_vector<float> operator() () { // output matrix
            return this->matrix;
        }
        __host__ __device__ make_cube() { }
        __host__ __device__ make_cube(int x, int y, int z){
            this->x = x;
            this->y = y;
            this->z = z;
        }
        __host__ __device__ ~make_cube() {
            this->matrix.shrink_to_fit();
            //this->matrix.clear();

        }

        __host__ void machine_info() {
            int deviceCount;
            cudaGetDeviceCount(&deviceCount);
            cout << "Detected " << deviceCount << " devices:" << endl;
            for (int device = 0; device < deviceCount; ++device) {
                cudaDeviceProp deviceProp;
                cudaGetDeviceProperties(&deviceProp, device);
                cout << "Device " << device << ": " << deviceProp.name << endl;
                cout << "  Compute Capability: " << deviceProp.major << "." << deviceProp.minor << endl;
                cout << "  Total Global Memory: " << deviceProp.totalGlobalMem / (1024 * 1024) << " MB" << endl;
                cout << "  Multiprocessors: " << deviceProp.multiProcessorCount << endl;
                cout << "  Max Threads per Block: " << deviceProp.maxThreadsPerBlock << endl;
                cout << "  Max Threads per Multiprocessor: " << deviceProp.maxThreadsPerMultiProcessor << endl;
                cout << "  Warp Size: " << deviceProp.warpSize << endl;
                cout << "  Max Blocks per Grid: " << deviceProp.maxGridSize[0] << " x " << deviceProp.maxGridSize[1] << " x " << deviceProp.maxGridSize[2] << endl;
                cout << "  Max Threads per Dim: " << deviceProp.maxThreadsDim[0] << " x " << deviceProp.maxThreadsDim[1] << " x " << deviceProp.maxThreadsDim[2] << endl;
                cout << "  Memory Clock Rate: " << deviceProp.memoryClockRate / 1e3 << " MHz" << endl;
                cout << "  Memory Bus Width: " << deviceProp.memoryBusWidth << " bits" << endl;
                cout << "  L2 Cache Size: " << deviceProp.l2CacheSize / 1024 << " KB" << endl;
                cout << "  Max Shared Memory per Block: " << deviceProp.sharedMemPerBlock / 1024 << " KB" << endl;
                cout << "  Max Registers per Block: " << deviceProp.regsPerBlock << endl;
                cout << "  Clock Rate: " << deviceProp.clockRate / 1e3 << " MHz" << endl;
                cout << "  Texture Alignment: " << deviceProp.textureAlignment << " bytes" << endl;
                cout << "  GPU Overlap: " << (deviceProp.deviceOverlap ? "Yes" : "No") << endl;
                cout << "  Kernel Execution Timeout: " << (deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No") << endl;
                cout << "  Concurrent Kernels: " << (deviceProp.concurrentKernels ? "Yes" : "No") << endl;
                cout << "  ECC Memory: " << (deviceProp.ECCEnabled ? "Yes" : "No") << std::endl;
                cout << "  Unified Addressing: " << (deviceProp.unifiedAddressing ? "Yes" : "No") << endl;
            }

        }

        __host__ make_cube(int x, int y, int z, float cylinder_radius) : make_cube(x, y, z) {
            if(cylinder_radius > (float)(x/2) || cylinder_radius > (float)(z/2))  {
                throw std::runtime_error("Cylinder is out of cube");
            }
            this->radius = cylinder_radius;
            matrix = thrust::device_vector<float>(x*y*z, 0);
            //matrix.resize(x*y*z);
            //std::fill(matrix.begin(), matrix.end(),0);
            //float* pd_vec = thrust::raw_pointer_cast(matrix.data());
            //pd_vec[0] = 20;
        }
        /*
         Паралельное заполнение элементов, логика такова:
            - 1 поток в строке(x) и столбце(y), а 256 на каждую долготу(z)
            - один блок в столбце(y) и строке(x), несколько блоков в глубине(z)
        */
        void computeMat() {
            thrust::device_ptr<float> ptr = thrust::device_pointer_cast<float>(this->matrix.data());
            make_cube * obj;
            cudaMalloc(&obj, sizeof(make_cube));
            cudaMemcpy(obj, this, sizeof(make_cube), cudaMemcpyHostToDevice);

            //MakeMat_byKernel<<< dim3(16,16,16), dim3((this->x+ 15) / 16, (this->y+ 15) / 16, (this->z+ 15) / 16) >>>(obj, ptr );
            dim3 blocks(2, 2, 2);
            MakeMat_byKernel<<<
            dim3(
                    (this->x+ blocks.x-1) / blocks.x,
                    (this->y+ blocks.y-1) / blocks.y,
                    (this->z+ blocks.z-1) / blocks.z),
            blocks
            >>>(obj, ptr );

            if(this->run_test) {
                std::cerr << "----------------BLOCK TEST----------------" << std::endl;
                cudaDeviceSynchronize();
                sleep(1);
                std::cerr << "VERIFY: BLOCKS=" << ((this->x+ blocks.x-1) / blocks.x) * ((this->y+ blocks.y-1) / blocks.y) * ((this->z+ blocks.z-1) / blocks.z)
                          << " THREADS=" << blocks.x*blocks.y*blocks.z << endl;
                std::cerr << "----------------END BLOCK TEST----------------" << std::endl;
            }
            cudaMemcpy(this, obj, sizeof(make_cube), cudaMemcpyDeviceToHost);
            cudaFree(obj);
            cudaDeviceSynchronize();

        }
        void iterateMat() {
            //int N = this->y*this->x*this->z;
            auto collector = [this](int i, int j, int k ) -> float {  // вызов kernel из host
                float value = 0;
                float * value_ptr;
                cudaMalloc(&value_ptr, sizeof(float));
                cudaMemcpy(value_ptr, &value, sizeof(float), cudaMemcpyHostToDevice);

                thrust::device_ptr<float> ptr = thrust::device_pointer_cast<float>(this->matrix.data());
                make_cube * obj;
                cudaMalloc(&obj, sizeof(make_cube));
                cudaMemcpy(obj, this, sizeof(make_cube), cudaMemcpyHostToDevice);
                /*
                  for(int x1 = ((i>0) ? i-1 : i) ; x1 <= i+1 && x1 < this->x; x1++) {
                                for(int y1 = ((j>0) ? j-1 : j); y1 <= j+1 && y1 < this->y; y1++) {
                                    for(int z1 = ((k>0) ? k-1 : k); z1 <= k+1 && z1 < this->z; z1++) {
                 */
                //MakeMat_byKernel<<< dim3(this->x, this->y-1, this->z), dim3((this->x+ 15) / 16, (this->y+ 15) / 16, (this->z+ 15) / 16) >>>(obj, ptr );
                dim3 blocks(1,1,1);

                collect_byKernel <<<
                dim3(
                        ((i > 0) ? 3 : 2), // Исправлено на 3, чтобы учесть случай, когда i равно 0
                        ((j > 0) ? 3 : 2), // Исправлено на 3, чтобы учесть случай, когда j равно 0
                        ((k > 0) ? 3 : 2)  // Исправлено на 3, чтобы учесть случай, когда k равно 0
                ),blocks

                >>>(obj, ptr, value_ptr, ((i>0)? i-1 : i), ((j>0)? j-1 : j), ((k>0)? k-1 : k) );
                /*dim3 blocks(2,2,3);
                collect_byKernel <<<blocks,
                dim3(
                        1,
                        1,
                        1
                )
                >>>(obj, ptr, value_ptr, ((i>0)? i-1 : i), ((j>0)? j-1 : j), ((k>0)? k-1 : k) );*/
                cudaDeviceSynchronize();
                cudaMemcpy(this, obj, sizeof(make_cube), cudaMemcpyDeviceToHost);
                cudaMemcpy(&value, value_ptr, sizeof(float), cudaMemcpyDeviceToHost);                   // еблан забывать float?
                cudaFree(obj);
                cudaFree(value_ptr);

                return value;
            };
            if(this->run_test) {
                cudaDeviceSynchronize();
                std::cerr << "----------------BLOCK TEST----------------" << std::endl;
                cudaDeviceSynchronize();
                //sleep(1);
                int i=0,j=3,k=0;    // test input
                int test = collector(i, j, k);
                cudaDeviceSynchronize();
                sleep(1.2);
                std::cerr << "VERIFY BLOCKS: " << (i>0 ? 3 : 2) * (j>0 ? 3 : 2) * (k>0 ? 3 : 2) << endl;
                std::cerr << "OUT: " << test << endl;
                std::cerr << "----------------END BLOCK TEST----------------" << std::endl;
                return;
            }
            float norm = 1;
            int count = 0;
            while (norm > eps && count < maxIteration) {
                // перебор каждой точки
                for(int i = 0; i < this->x; i++) {
                    for (int j = 0; j < this->y; j++) {
                        for (int k = 0; k < this->z; k++) {
                            // просмотр соседей и сумма новой температуры
                            float oldTemp = this->matrix[index(i, j, k)];
                            if(this->show_iter) { cout << "newtemp="; }
                            float newTemp = collector(i, j, k);
                            if(this->show_iter) { cout << endl; }
                            //newTemp *= this->alpha; // формула распределения: dT/dt = alpha * ( (d^2T/dx^2) + (d^2T/dy^2) + (d^2T/dz^2) )
                            // т.е. нужно взять от зависимой точки сумму всех соседей и умножить на коэф распределения
                            this->matrix[index(i, j, k)] = newTemp * this->alpha;
                            if(this->show_iter) {
                                cout << "newtemp=" << newTemp << " matrix=" << matrix[index(i, j, k)] << " (" << i << "," << j << "," << k << ")" << endl;
                            }
                            norm = std::abs( newTemp - oldTemp );
                            //printf("----------\n");
                        }
                    }
                }
                count++;
            }

        }
        __device__ void MakeMat(thrust::device_ptr<float> mat, make_cube * obj, dataId config ) {
            int x1 = config.thread.x + config.block.x * config.blockDim.x;
            int y1 = config.thread.y + config.block.y * config.blockDim.y;
            int z1 = config.thread.z + config.block.z * config.blockDim.z;
            //printf("Getted: (%d,%d,%d)\n", x1, y1, z1);
            if(z1 == 0  ) {
                /*if(obj->run_test) {
                    printf("Getted: (%d,%d,%d)\n", x1, y1, z1);
                }*/
                mat[ index(x1, y1, z1) ] = this-> T1;
                //this->matrix[ index_device(i,j,k) ] = this->T1;
                //this->matrix[i][j][k] = this->T1;
            }
            else if(z1 == (this->z)-1) {
                /*if(obj->run_test) {
                    printf("Getted: (%d,%d,%d)\n", x1, y1, z1);
                }*/
                mat[ index(x1, y1, z1) ] = this->T2;
                //this->matrix[i][j][k] = this->T2;
            }
        }
        __device__ void collectNeighbors (thrust::device_ptr<float> mat, dataId config, float* output, int xpos, int ypos, int zpos, make_cube * obj) { // using basic matrix iterator

            int x1 = config.thread.x + config.block.x * config.blockDim.x; // допустим, коорда подана: (3,0,0), решает как 3+0, 3+1, 3+2. Если коорда (0,0,0), выдаст 0+0, 0+1
            int y1 = config.thread.y + config.block.y * config.blockDim.y;
            int z1 = config.thread.z + config.block.z * config.blockDim.z;
            if(x1+xpos >= this->x || y1+ypos >= this->y || z1+zpos >= this->z) return ;
            if(obj->run_test) {
                printf("Getted: (%d,%d,%d) + (%d,%d,%d) => (%d,%d,%d) output=%.4f checked=%d\n", x1, y1, z1, xpos, ypos, zpos, x1+xpos, y1+ypos, z1+zpos, *output,
                       this->radius < (float)std::sqrt(std::pow(x1 - (float)(this->x / 2), 2) + std::pow(z1 - (float)(this->z / 2), 2)));
            }
            if(obj->radius < (float)std::sqrt(std::pow(xpos+ x1 - (float)(obj->x / 2), 2) + std::pow(zpos + z1 - (float)(obj->z / 2), 2)))
            {
                atomicAdd(output, mat[index(xpos+x1, ypos+y1, zpos+z1)]);      // на каждом блоке отдельные адреса, нельзя просто суммировать адреса
                if(obj->show_iter) {
                    float temp = mat[index(xpos + x1, ypos + y1, zpos + z1)];
                    printf("+%.1f", temp);
                }
                //(*output) += mat[index(x1, y1, z1)];

            }
            //printf("Getted: (%d,%d,%d), out=%d\n",x1,y1,z1, output);

        }
        __host__ void show(int precision = 3 ) {
            cout.setf(std::ios::fixed);
            cout.precision(precision);
            cudaDeviceSynchronize();
            for (int i = 0; i < this->x; i++) {
                for (int k = 0; k < this->z; k++) {
                    cout << this->matrix[index(i, (this->y)-1, k)] << " ";
                }
                cout << endl;
            }
            cout << "-------------" << endl;
            for (int j = 0; j < this->y; j++) {
                for (int i = 0; i < this->x; i++) {
                    cout << this->matrix[index(i, j, (int)((this->z)/2))] << " ";
                }
                cout << endl;
            }
            cout << "-------------" << endl;
            for (int j = 0; j < this->y; j++) {
                for (int k = 0; k < this->z; k++) {
                    cout << this->matrix[index((int)((this->x)/2), j, k)] << " ";
                }
                cout << endl;
            }
            cout << "-------------" << endl;
            for (int j = 0; j < this->y; j++) {
                for (int i = 0; i < this->x; i++) {
                    cout << this->matrix[index(i, j, (int)(this->z-1))] << " ";
                }
                cout << endl;
            }
            cout << "-------------" << endl;
            for (int j = 0; j < this->y; j++) {
                for (int i = 0; i < this->x; i++) {
                    cout << this->matrix[index(i, j, (int)(0))] << " ";
                }
                cout << endl;
            }
        }
        __host__ void write_file(std::string prefix = "file_", int precision = 3 ) {
            std::ofstream answer( prefix + this->file_name , std::ios_base::trunc);
            answer.setf(std::ios::fixed);
            answer.precision(precision);
            cudaDeviceSynchronize();
            answer << "------Сторона XZ сверху-------" << endl;
            for (int i = 0; i < this->x; i++) {
                for (int k = 0; k < this->z; k++) {
                    answer << this->matrix[index(i, (this->y)-1, k)] << " ";
                }
                answer << endl;
            }
            answer << "-------Сторона XY в разрезе------" << endl;
            for (int j = 0; j < this->y; j++) {
                for (int i = 0; i < this->x; i++) {
                    answer << this->matrix[index(i, j, (int)((this->z)/2))] << " ";
                }
                answer << endl;
            }
            answer << "-------Сторона YZ при в разрезе------" << endl;
            for (int j = 0; j < this->y; j++) {
                for (int k = 0; k < this->z; k++) {
                    answer << this->matrix[index((int)((this->x)/2), j, k)] << " ";
                }
                answer << endl;
            }
            answer << "-------Сторона YZ при Z=0------" << endl;
            for (int j = 0; j < this->y; j++) {
                for (int i = 0; i < this->x; i++) {
                    answer << this->matrix[index(i, j, (int)(this->z-1))] << " ";
                }
                answer << endl;
            }
            answer << "-------Сторона YZ при Z=Z-1------" << endl;
            for (int j = 0; j < this->y; j++) {
                for (int i = 0; i < this->x; i++) {
                    answer << this->matrix[index(i, j, (int)(0))] << " ";
                }
                answer << endl;
            }
            answer.close();
        }
    };



__global__ void MakeMat_byKernel(make_cube * obj , thrust::device_ptr<float> data) {
    //printf("Config: blocks (%d,%d,%d), threads (%d,%d,%d), blockDim(%d,%d,%d)\n",blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y,threadIdx.z,blockDim.x,blockDim.y,blockDim.z);
    dataId config = { {threadIdx.x, threadIdx.y, threadIdx.z}, {blockIdx.x, blockIdx.y, blockIdx.z}, {blockDim.x, blockDim.y, blockDim.z} };
    //printf("Verify: blocks (%d,%d,%d), threads (%d,%d,%d), blockDim(%d,%d,%d)\n", config.block.x, config.block.y, config.block.z, config.thread.x, config.thread.y,config.thread.z,config.blockDim.x,config.blockDim.y,config.blockDim.z);
    obj->MakeMat( data, obj, config );
    if(obj->run_test) {
        int threadId = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
        int blockId = blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z);
        int numThreads = blockDim.x * blockDim.y * blockDim.z;
        int numBlocks = gridDim.x * gridDim.y * gridDim.z;
        if (threadId == 0 && blockId == 0) {
            // Вывод информации только из одного потока в блоке 0
            printf("\nMakeMat_byKernel:\n");
            printf("\tNumber of Threads: %d\n", numThreads);
            printf("\tNumber of Blocks: %d\n", numBlocks);
            printf("\tblockDim: (%d, %d, %d)\n", blockDim.x, blockDim.y, blockDim.z);
            printf("\tgridDim: (%d, %d, %d)\n", gridDim.x, gridDim.y, gridDim.z);
            printf("\tDebug values:\n");
            printf("\t\tX=%d, Y=%d, Z=%d, radius=%d\n", obj->x, obj->y, obj->z, obj->radius);
            printf("\t\trun_test=%d, maxIteration=%d, eps=%f \n", obj->run_test, obj->maxIteration, obj->eps);
            printf("\t\talpha=%f, T1=%f, T2=%f \n", obj->alpha, obj->T1, obj->T2);
        }
    }
}
__global__ void collect_byKernel(make_cube * obj , thrust::device_ptr<float> data, float* output, int x, int y, int z ) {
    //printf("Config: blocks (%d,%d,%d), threads (%d,%d,%d)\n",blockDim.x, blockDim.y, blockDim.z, threadIdx.x, threadIdx.y,threadIdx.z);
    dataId config = { {threadIdx.x, threadIdx.y, threadIdx.z}, {blockIdx.x, blockIdx.y, blockIdx.z}, {blockDim.x, blockDim.y, blockDim.z} };
    obj->collectNeighbors( data, config, output, x,y,z, obj );
    if(obj->run_test) {
        int threadId = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
        int blockId = blockIdx.x + gridDim.x * (blockIdx.y + gridDim.y * blockIdx.z);
        int numThreads = blockDim.x * blockDim.y * blockDim.z;
        int numBlocks = gridDim.x * gridDim.y * gridDim.z;
        if (threadId == 0 && blockId == 0) {
            printf("\ncollect_byKernel:\n");
            // Вывод информации только из одного потока в блоке 0
            printf("\tNumber of Threads: %d\n", numThreads);
            printf("\tNumber of Blocks: %d\n", numBlocks);
            printf("\tblockDim: (%d, %d, %d)\n", blockDim.x, blockDim.y, blockDim.z);
            printf("\tgridDim: (%d, %d, %d)\n", gridDim.x, gridDim.y, gridDim.z);
            printf("\tDebug values:\n");
            printf("\t\tX=%d, Y=%d, Z=%d, radius=%d\n", obj->x, obj->y, obj->z, obj->radius);
            printf("\t\trun_test=%d, maxIteration=%d, eps=%f \n", obj->run_test, obj->maxIteration, obj->eps);
            printf("\t\talpha=%f, T1=%f, T2=%f \n", obj->alpha, obj->T1, obj->T2);
        }
    }
    __syncthreads();

}



/*__device__ void oo() {

}
__global__ void kernel() {
    oo();
}*/



int main() {
    cout.setf(std::ios::fixed); cout.precision(3);
    std::cout << "Hello, World!" << std::endl;
    int x = 5, y = 3, z = 5;
    float radius = 2.0f;
    float alpha = 0.001f;
    float eps = 0.01f;
    int maxIter = 50;
    bool run_test = false;
    bool show_iter = false;
    if(x%2 != 0 || y%2 != 0 || z%2 != 0) {
        std::cerr << "WARNING: scary input (multiples of 2 are required)! " << endl;
    }

    //gpu_main<<<1,1>>>( thrust::device_pointer_cast<float>( matrix.data() ), x,y,z );
    //show();
    //cube_gpu.show();
    //kernel<<<1,1>>>();
    // =================cpu zone
    cpu cube_cpu(x, y, z, radius);
    cube_cpu.alpha = alpha;
    cube_cpu.computeMat();
    cube_cpu.iteration(eps, maxIter);
    // =================cpu zone
    // =================make_cube zone
    make_cube cube_gpu(x, y, z, radius);
    //cube_gpu.machine_info();
    cube_gpu.run_test = run_test;
    cube_gpu.show_iter = show_iter;
    cube_gpu.file_name = "gpumatrix";
    if(cube_gpu.run_test) {
        std::cerr << "******************************" << std::endl;
        std::cerr << "*       RUNNED RUN_TEST      *" << std::endl;
        std::cerr << "*  (p.s. iterator disabled)  *" << std::endl;
        std::cerr << "******************************" << std::endl;
    }
    cube_gpu.alpha = alpha; //0.05351f;
    cube_gpu.eps = eps;
    cube_gpu.maxIteration = maxIter;
    //make_cube * obj;
    cube_gpu.computeMat();
    cube_gpu.iterateMat();
    //cube_gpu.write_file("output_");
    cube_gpu.show();
    // =================make_cube zone
    return 0;
}

/*__device__ __host__ int x=5,y=5,z=5, radius = 2;
__device__ __host__ float T1 = 20.0f;              // температура 1 стороны
__device__ __host__ float T2 = 10.0f;              // температура противоположной стороны
__device__ __host__ float T_bottom = 0.0f;         // температура нижней грани
__device__ __host__ float alpha = 0.001;           // коэф теплопроводности
thrust::device_vector<float> matrix(x*y*z, 0);
__host__ void genValue() {
    srand(time(NULL));
    T1 = 100 - (rand() % 50);
    T2 = 100 - (rand() % 50);
}
__host__ __device__ int index(int PointX, int PointY, int PointZ) {
    return PointX + x * (PointY + y * z);
}
__device__ void MakeMat( thrust::device_ptr<float> mat, int x_dev, int y_dev, int z_dev) {
    for(int i = 0; i < x_dev; i++) {
        for(int j = 0; j < y_dev; j++) {
            for(int k = 0; k < z_dev; k++) {
                if(k == 0) {
                    mat[ index(i, j, k) ] = T1;
                    //this->matrix[ index_device(i,j,k) ] = this->T1;
                    //this->matrix[i][j][k] = this->T1;
                }
                else if(k == z_dev-1) {
                    mat[ index(i, j, k) ] = T2;
                    //this->matrix[i][j][k] = this->T2;
                }

            }
        }
    }
}
__host__ void show() {
    cudaDeviceSynchronize();
    for (int i = 0; i < x; i++) {
        for (int k = 0; k < z; k++) {
            cout << matrix[index(i, (y)-1, k)] << " ";
        }
        cout << endl;
    }
    cout << "-------------" << endl;
    for (int j = 0; j < y; j++) {
        for (int i = 0; i < x; i++) {
            cout << matrix[index(i, j, (int)((z)/2))] << " ";
        }
        cout << endl;
    }
    cout << "-------------" << endl;
    for (int j = 0; j < y; j++) {
        for (int k = 0; k < z; k++) {
            cout << matrix[index((int)((x)/2), j, k)] << " ";
        }
        cout << endl;
    }
}

__global__ void gpu_main( thrust::device_ptr<float> data, int i, int j, int k ) {

    MakeMat(data, x,y,z);
    //show();

}*/