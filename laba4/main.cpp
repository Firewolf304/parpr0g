#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_TARGET_OPENCL_VERSION 300

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "math.h"
#include <CL/cl2.hpp>
/*
 Вычислить методом последовательных приближений распределение значений температуры в
 точках верхней границы плоской прямоугольной области размером m * n, имеющей внутри
 круглый вырез. Теплопроводность материала области конечна и не равна нулю. Левая и
 нижняя граница области имеют постоянную температуру 0. Правый верхний угол имеет
 постоянную температуру T, значения температур точек правой границы постоянны и
 равномерно убывают от T до 0. Все остальные точки в начальный момент имеют температуру Т3.
*/
using std::vector;


using std::vector;
using std::cout;
using std::cin;
using std::endl;


class make_matrix_gpu {
    int index(int x, int y) {
        return this->m * y + x;
    }
/*matrix example:
 n
 --------------
 0 x x x T    |
 0 x x x (T-1)|
 0 x x x (T-2)|
 0 x x x (T-3)|
 0 0 0 0 0    |  m*/
public:
    //=================variables=================
    vector<double> matrix;
    int n = 3;      // like y
    int m = 3;      // like x
    double T = 20;
    int t3 = 3;     // ???
    double radius = 1.0f;
    double alpha = 0.05351f; // коэф распределения (лучше коэф газа Sulfur dioxide)
    //=================variables=================

    //=================GPU=================
    std::vector<cl::Platform> platforms;
    cl::Platform platform;
    std::vector<cl::Device> devices;
    cl::Device device;
    cl::Context context;
    cl::CommandQueue commandList;
    /*__global double* output, int m, int n, float step, double radius, int t3*/
    std::string makeMatrix =   "__kernel void makeMatrix(__global double* output, int m, int n, float step, double radius, int t3) {"
                               "    __private int x = get_global_id(0);\n"
                               "    __private int y = get_global_id(1);\n"
                               //"    output[m * y + x] = sqrt((float)y);"
                               "    if(radius < sqrt((float) (pow((float) (x - (float) (n / 2)), 2) + pow((float) (y - (float) (m / 2)),2)))) {"
                               "        if (x == 0 && y >= 0 || x >= 0 && y == n - 1) {"
                               "            output[m * y + x] = 0;"
                               "        }else if (x == m - 1 && y > 0) {"
                               "            "
                               "        } else {"
                               "            if (y == 0 && x == m - 1) return;"
                               "            output[m * (y) + (x)] = t3;"
                               "        }"
                               "    }"
                               "}";
    /*0__global double* input, 1__global float* output, 2__private int xpos, 3__private int ypos, 4int m, 5int n, 6double radius*/
    std::string programRun =  "__kernel void collectNeighbors(__global double* input, __global double* output, int xpos, int ypos, int m, int n, double radius) {"
                              " int x = get_global_id(0);"
                              " int y = get_global_id(1);"
                              //" *output = xpos;"
                              " if(x+xpos >= n || y+ypos >= m) return ;"
                              " if(radius < sqrt((float) (pow((float) (xpos+x - (float) (n / 2)), 2) + pow((float) (ypos+y - (float) (m / 2)),2)))) {"
                              //"       *output += input[m * (ypos+ y) + xpos + x];"
                              //"       atomic_add(output, input[m * (ypos+ y) + xpos + x]);"
                              //"         *output += input[m * (ypos+ y) + xpos + x];"
                                "       *output += x;"
                              " }"
                              "}";
    double valueGet(int x, int y) {
        double temp = 0;
        this->Readyoffset = cl::NDRange(0, 0);
        this->ReadylocalSize = cl::NDRange(1, 1);
        this->ReadyglobalSizeNDRange = cl::NDRange(3, 3);

        cl::Buffer inputArray(this->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, this->matrix.size() * sizeof(double), this->matrix.data());
        cl::Buffer outputValue(this->context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(temp), &temp);
        this->ReadyKernelRun.setArg(0, inputArray);
        this->ReadyKernelRun.setArg(1, outputValue);
        this->ReadyKernelRun.setArg(2, (x>0) ? x-1 : x);
        this->ReadyKernelRun.setArg(3, (y>0) ? y-1 : y);
        this->commandList.enqueueNDRangeKernel(ReadyKernelRun, this->Readyoffset, this->ReadyglobalSizeNDRange, this->ReadylocalSize);
        this->commandList.enqueueReadBuffer(outputValue, CL_TRUE, 0, sizeof(temp), &temp);
        return temp;
        std::string remove = programRun;
    }

    cl::Program ReadyProgramRun;
    cl::Kernel ReadyKernelRun;
    cl::NDRange Readyoffset;
    cl::NDRange ReadylocalSize;
    cl::NDRange ReadyglobalSizeNDRange;
    //=================GPU=================

    static const void checkError(cl_int err, const std::string& message) {
        if (err != CL_SUCCESS) {
            std::cerr << "Ошибка: " << message << " (" << err << ")" << std::endl;
            exit(err);
        }
    }
    void init_gpu(int IDplatform = 0, int IDdevice = 0) {
        this->matrix = vector<double>(this->n * this->m, 0);
        cl::Platform::get(&(this->platforms));
        checkError(this->platforms.size() > 0 ? CL_SUCCESS : -1, "No platforms");
        this->platform = this->platforms[IDplatform];

        platform.getDevices(CL_DEVICE_TYPE_ALL, &(this->devices));
        checkError(this->devices.size() > 0 ? CL_SUCCESS : -1, "No devices");
        this->device = this->devices[IDdevice];
        this->context = cl::Context(this->device);
        this->commandList = cl::CommandQueue(this->context, this->device);
        //cl::Program::Sources sources(1, std::make_pair(programRun.c_str(), programRun.length() + 1));

        //init gpu program
        this->ReadyProgramRun = cl::Program(this->context, programRun);
        checkError(this->ReadyProgramRun.build({ this->device }, "-cl-std=CL2.0 -cl-fast-relaxed-math -cl-mad-enable"), "Ошибка сборки программы");
        this->ReadyKernelRun = cl::Kernel(this->ReadyProgramRun, "collectNeighbors");
        this->ReadyKernelRun.setArg(4, this->m);
        this->ReadyKernelRun.setArg(5, this->n);
        this->ReadyKernelRun.setArg(6, this->radius);
    }
    void make_gpu() {
        float step = (float)(T / (n-1));
        matrix[ index(this->m-1,0) ] = this->T;
        for (int i = 1; i < this->n; i++) {
            this->matrix [ index(this->m-1,i) ] = (double) (this->matrix[index(this->m-1,i-1)] - step);
        }
        cl::Buffer outputBuffer(this->context, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, this->matrix.size() * sizeof(double), this->matrix.data());
        vector<size_t> global = { (size_t)this->m, (size_t)this->n }; // ПОЗИЦИЯ ДОЛЖНА БЫТЬ РАВНА РАСЧЕТУ index!

        // BUILD
        cl::Program program(this->context, makeMatrix);
        checkError(program.build({ this->device }, "-cl-std=CL2.0 -Werror -cl-fast-relaxed-math -cl-mad-enable"), "Ошибка сборки программы");

        // ARGS
        cl::Kernel kernel(program, "makeMatrix");
        kernel.setArg(0, outputBuffer);
        kernel.setArg(1, this->m);
        kernel.setArg(2, this->n);
        kernel.setArg(3, step);
        kernel.setArg(4, this->radius);
        kernel.setArg(5, this->t3);

        // RUN
        cl::NDRange offset(0, 0);
        cl::NDRange localSize(1, 1);
        cl::NDRange globalSizeNDRange(global[0], global[1]);
        this->commandList.enqueueNDRangeKernel(kernel, offset, globalSizeNDRange, localSize);
        this->commandList.enqueueReadBuffer(outputBuffer, CL_TRUE, 0, this->matrix.size() * sizeof(double), this->matrix.data());
    }
    void device_info() {
        if(this->platforms.empty()) return;
        int platform_id = 0;
        int device_id = 0;
        for(std::vector<cl::Platform>::iterator it = this->platforms.begin(); it != this->platforms.end(); ++it) {
            cl::Platform platform(*it);

            std::cout << "Platform ID: " << platform_id++ << std::endl;
            std::cout << "Platform Name: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
            std::cout << "Platform Vendor: " << platform.getInfo<CL_PLATFORM_VENDOR>() << std::endl;

            std::vector<cl::Device> devices;
            platform.getDevices(CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_CPU, &devices);

            for(std::vector<cl::Device>::iterator it2 = devices.begin(); it2 != devices.end(); ++it2) {
                cl::Device device(*it2);
                std::cout  << "\tDevice " << device_id++ << ": " << std::endl;
                std::cout << "\t\tDevice Name: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
                std::cout << "\t\tDevice Vendor: " << device.getInfo<CL_DEVICE_VENDOR>() << std::endl;
                std::cout << "\t\tDevice Version: " << device.getInfo<CL_DEVICE_VERSION>() << std::endl;
                switch (device.getInfo<CL_DEVICE_TYPE>()) {
                    case 4:
                        std::cout << "\t\tDevice Type: GPU" << std::endl;
                        break;
                    case 2:
                        std::cout << "\t\tDevice Type: CPU" << std::endl;
                        break;
                    default:
                        std::cout << "\t\tDevice Type: unknown" << std::endl;
                }
                std::cout << "\t\tDevice Max Compute Units: " << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl;
                std::cout << "\t\tDevice Global Memory: " << device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << " bytes" << std::endl;
                std::cout << "\t\tDevice Max Clock Frequency: " << device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
                std::cout << "\t\tDevice Max Memory Allocation: " << device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() << std::endl;
                std::cout << "\t\tDevice Local Memory: " << device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;
                std::cout << "\t\tDevice Available: " << device.getInfo<CL_DEVICE_AVAILABLE>() << std::endl;
            }
            std::cout << std::endl;
        }
    }

    vector<double> get_top_matrix(vector<vector<double>> mat) {
        vector<double> out;
        for(auto d : mat[0]) {
            out.push_back(d);
        }
        return out;
    }
    void show_matrix(int precision = 3) {
        cout.setf(std::ios::fixed);
        cout.precision(precision);
        for (int y = 0; y < this->n; y++) {
            for (int x = 0; x < this->m; x++) {
                cout << this->matrix[index(x, y)] << " ";
            }
            cout << endl;
        }
    }
    vector<double> operator()() {
        return matrix;
    }



    void iteration(double eps = 0.0001f, int maxIteration=10000) {


        double norm = 1;
        int count = 0;
        while (norm > eps && count < maxIteration) {
            for(int i = 0; i < this->n; i++) {
                for (int j = 0; j < this->m; j++) {
                    double oldTemp = this->matrix[ index(i, j) ];
                    float newTemp = valueGet(j,i);
                    cout << newTemp << " " << i << ":" << j << endl;
                    this->matrix[ index(i, j) ] = newTemp * this->alpha;
                    norm = std::abs(newTemp - oldTemp);
                }
            }
            count++;
        }
    }
};

int main() {
    std::cout << "Hello, World!" << std::endl;
    make_matrix_gpu matrix;
    matrix.m = 9;
    matrix.n = 9;
    matrix.radius = 1.5f;
    matrix.T = 100;
    matrix.t3 = 1;
    matrix.init_gpu();
    matrix.device_info();
    matrix.make_gpu();
    matrix.show_matrix();
    cout << endl;
    matrix.iteration(0.01f, 1);
    matrix.show_matrix();
    /*
    matrix.n = 9; matrix.m = 9;
    matrix.T = 100;
    matrix.t3 = 100;
    matrix.radius = 1.5f;
    matrix.alpha = 0.05351f;

    matrix.make();
    matrix.show_matrix();

    matrix.iteration(0.01f, 100);
    matrix.show_matrix();*/
    return 0;
}
