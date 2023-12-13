//#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_TARGET_OPENCL_VERSION 120

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "math.h"
#include <CL/cl2.hpp>
#include <csignal>
#include <time.h>
#include <unistd.h>
#include <fstream>

#define NS_PER_SECOND 1000000000
/*
 Вычислить методом последовательных приближений распределение значений температуры в
 точках верхней границы плоской прямоугольной области размером m * n, имеющей внутри
 круглый вырез. Теплопроводность материала области конечна и не равна нулю. Левая и
 нижняя граница области имеют постоянную температуру 0. Правый верхний угол имеет
 постоянную температуру T, значения температур точек правой границы постоянны и
 равномерно убывают от T до 0. Все остальные точки в начальный момент имеют температуру Т3.
*/
using std::vector;
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


class make_matrix_cpu {
/*matrix example:
 n
 --------------
 0 x x x T    |
 0 x x x (T-1)|
 0 x x x (T-2)|
 0 x x x (T-3)|
 0 0 0 0 0    |  m*/
public:

    vector<vector<double>> matrix;
    int n = 3;      // like y
    int m = 3;      // like x
    double T = 20;
    int t3 = 3;     // ???
    double radius = 1.0f;
    double alpha = 0.05351f; // коэф распределения (лучше коэф газа Sulfur dioxide)
    void show_matrix(int accuracy = 3) {
        cout.setf(std::ios::fixed); cout.precision(accuracy);
        for(int ypos = 0; ypos < this->matrix.size(); ypos++) {
            for(int xpos = 0; xpos < this->matrix[ypos].size(); xpos++) {
                cout << matrix[ypos][xpos] << " ";
            }
            cout << endl;
        }
        /*for(auto y : this->matrix) {
            for(auto x : y) {
                cout << x << " ";
            }
            cout << endl;
        }*/
    }
    void write_file(std::string filename = "output_cpu", int precision = 3 ) {
        std::ofstream answer( filename.c_str() , std::ios_base::trunc);
        answer.setf(std::ios::fixed);
        answer.precision(precision);
        answer << "TOP:" << endl;
        for (int x = 0; x < this->m; x++) {
            answer << this->matrix[0][x] << " ";
        }
        answer << "\nOTHER:" << endl;
        for (int y = 0; y < this->n; y++) {
            for (int x = 0; x < this->m; x++) {
                answer << matrix[y][x] << " ";
            }
            answer << endl;
        }

        answer.close();
    }
    void show_matrix(vector<vector<double>> mat, int accuracy = 3) {
        cout.setf(std::ios::fixed); cout.precision(accuracy);
        for(auto y : mat) {
            for(auto x : y) {
                cout << x << " ";
            }
            cout << endl;
        }
    }
    void show_matrix(vector<double> mat) {
        for(auto x : mat) {
            cout << x << " ";
        }
        cout << endl;
    }
    vector<double> get_top_matrix(vector<vector<double>> mat) {
        vector<double> out;
        for(auto d : mat[0]) {
            out.push_back(d);
        }
        return out;
    }
    vector<vector<double>> get_another_matrix(vector<vector<double>> mat) {
        vector<vector<double>> out;
        out.push_back({});
        for(int y = 1; y < this->n; y++) {
            for(int x = 0; x < this-> m; x++) {
                out[y].push_back(mat[y][x]);
            }
            out.push_back({});
        }
        return out;
    }
    void make() {
        float step = (float)(T / (n-1));
        this->matrix = vector<vector<double>>(this->n, vector<double>(this->m, 0));
        matrix[0][m-1] = this->T;

        for(int y = 0 ; y < this->n; y++ ) {
            for(int x = 0; x < this->m; x++) {
                if(this->radius < (float) std::sqrt(std::pow(x - (float) (this->n / 2), 2) +
                                                    std::pow(y - (float) (this->m / 2),2))) {
                    if (x == 0 && y >= 0 || x >= 0 && y == this->n - 1) { // left/down place
                        this->matrix[y][x] = 0;
                    } else if (x == this->m - 1 && y > 0) { // right place
                        this->matrix[y][x] = (double) (this->matrix[y - 1][x] - step);
                    } else { // another
                        if (y == 0 && x == this->m - 1) continue;
                        //this->matrix[y][x] = (double)(rand() % this->t3) * ((rand() % 1000000) > 500000 ? -1 : 1);
                        this->matrix[y][x] = this->t3;
                    }
                }
            }
        }
    }
    vector<vector<double>> operator()() {
        return matrix;
    }
    void iteration(double eps = 0.0001f, int maxIteration=10000) {
        double norm = 1;
        int count = 0;
        while (norm > eps && count < maxIteration) {
            for(int i = 0; i < this->n; i++) {
                for (int j = 0; j < this->m; j++) {
                    double oldTemp = this->matrix[i][j];
                    double newTemp = 0.0f;
                    if(this->radius < (float) std::sqrt(std::pow(i - (float)(this->n / 2), 2) +
                                                        std::pow(j - (float)(this->m / 2), 2))) {
                        for (int y1 = ((i > 0) ? i - 1 : i); y1 <= i + 1 && y1 < this->n; y1++) {
                            for (int x1 = ((j > 0) ? j - 1 : j); x1 <= j + 1 && x1 < this->m; x1++) {
                                if (this->radius < (float) std::sqrt(std::pow(i - (float) (this->n / 2), 2) +
                                                                     std::pow(j - (float) (this->m / 2), 2))) {
                                    sleep(0.0f);
                                    newTemp += this->matrix[y1][x1];
                                }
                            }
                        }
                        this->matrix[i][j] = newTemp * this->alpha;
                        sleep(0.0f);
                        norm = std::abs(newTemp - oldTemp);
                    }
                }
            }
            count++;
        }
    }
};

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
    bool test = false;
    //service
    double koefX = 0;
    double koefY = 0;
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
                               "    if(radius < sqrt((float) (pow((float) (y - (float) (n / 2)), 2) + pow((float) (x - (float) (m / 2)),2)))) {"
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
    std::string programRun =  " __kernel void collectNeighbors(__global double* input, __global double* output, int xpos, int ypos, int m, int n, double radius, __global double * testX, __global double * testY) {"
                              " __private int x = get_global_id(0);"
                              " __private int y = get_global_id(1);"
                              //"barrier(CLK_GLOBAL_MEM_FENCE);"
                              " if(x+xpos >= m || y+ypos >= n){ return ;}"
                              " if(radius < sqrt((float) (pow((float)(xpos + x - (m / 2)), 2) + "
                              "                           pow((float)(ypos + y - (n / 2)), 2)))) {"
                              "       if(testX != NULL && testY != NULL) {"
                              "       testX[x] = x + xpos;"
                              "       testY[y] = y + ypos;}"
                              "         "
                              "       double save = input[m * (ypos+ y) + (xpos + x)];"
                              "       output[get_global_size(0) * y + x] = save;"
                              "       "
                              //"       testX[x] = save;"
                              //"         *output += input[m * (ypos+ y) + xpos + x];"
                              //"       *output = xpos;"
                              "}"
                              ""
                              "}";
    vector<double> testX;
    vector<double> testY;
    cl::Buffer testBuffX;
    cl::Buffer testBuffY;
    double valueGet(int x, int y) {
        double temp = 0;
        this->ReadylocalSize = cl::NDRange(1, 1);
        this->ReadyglobalSizeNDRange = cl::NDRange((x > 0 ? 3 : 2), (y > 0 ? 3 : 2));
        cl::Buffer inputArray(this->context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, this->matrix.size() * sizeof(double), this->matrix.data());
        vector<double> answer((x > 0 ? 3 : 2) * (y > 0 ? 3 : 2), 0);
        cl::Buffer outputValue(this->context, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, answer.size() * sizeof(double), answer.data());
        if(test) {
            testBuffX = cl::Buffer(this->context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, testX.size() * sizeof(double), testX.data());
            testBuffY = cl::Buffer(this->context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, testY.size() * sizeof(double), testY.data());
        }
        this->ReadyKernelRun.setArg(0, inputArray);
        this->ReadyKernelRun.setArg(1, outputValue);
        this->ReadyKernelRun.setArg(2, (x>0) ? x-1 : x);
        this->ReadyKernelRun.setArg(3, (y>0) ? y-1 : y);
        (this->test) ? this->ReadyKernelRun.setArg(7, testBuffX) : this->ReadyKernelRun.setArg(7, NULL);
        (this->test) ? this->ReadyKernelRun.setArg(8, testBuffY) : this->ReadyKernelRun.setArg(8, NULL);
        this->commandList.enqueueNDRangeKernel(ReadyKernelRun, cl::NullRange, this->ReadyglobalSizeNDRange, this->ReadylocalSize);
        this->commandList.enqueueReadBuffer(outputValue, CL_TRUE, 0, answer.size() * sizeof(double), answer.data());
        if(test)
            this->commandList.enqueueReadBuffer(testBuffX, CL_TRUE, 0, testX.size() * sizeof(double), testX.data());
        this->commandList.enqueueReadBuffer(testBuffY, CL_TRUE, 0, testY.size() * sizeof(double), testY.data());
        //this->commandList.finish();
        for(auto d : answer) {
            temp += d;
        }
        if (test) {
            cout << x << ":" << y << " temp=" << temp << " ready=" << temp * this->alpha << " test=" << endl << "X:\t";

            for (auto i: testX) {
                cout << i << " ";
            }
            cout << endl << "Y:\t";
            for (auto i: testY) {
                cout << i << " ";
            }
            cout << endl << "A:\t";
            for (auto i: answer) {
                cout << i << " ";
            }
            cout << endl;
        }
        return temp;
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
        this->koefX = this->m / 2;
        this->koefY = this->n / 2;
        this->testX = vector<double>(12, -1);
        this->testY = vector<double>(12, -1);
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
        cl_int err = this->ReadyProgramRun.build(this->device, "-cl-fast-relaxed-math -cl-mad-enable");
        checkError(err, "Ошибка сборки программы"); //-cl-fast-relaxed-math -cl-mad-enable -cl-std=CL2.0
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
        checkError(program.build(this->device, ""), "Ошибка сборки программы"); //-cl-std=CL2.0 -Werror -cl-fast-relaxed-math -cl-mad-enable

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

    void write_file(std::string filename = "output_gpu", int precision = 3 ) {
        std::ofstream answer( filename.c_str() , std::ios_base::trunc);
        answer.setf(std::ios::fixed);
        answer.precision(precision);
        answer << "TOP:" << endl;
        for (int x = 0; x < this->m; x++) {
            answer << this->matrix[index(x, 0)] << " ";
        }
        answer << "\nOTHER:" << endl;
        for (int y = 0; y < this->n; y++) {
            for (int x = 0; x < this->m; x++) {
                answer << this->matrix[index(x, y)] << " ";
            }
            answer << endl;
        }

        answer.close();
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
                    if(this->radius < std::sqrt(std::pow(i - this->koefY, 2) +
                                                std::pow(j - this->koefX, 2))) {
                        float newTemp = valueGet(j, i);
                        norm = std::abs(newTemp - this->matrix[index(j, i)]);
                        //cout << newTemp << " " << i << ":" << j << endl;
                        this->matrix[index(j, i)] = newTemp * this->alpha;
                        if(this->test) {
                            show_matrix();
                            cout << std::string(60, '-') << endl;
                        }
                        //
                    }
                }
            }
            count++;
        }
    }
};

int main() {
    std::cout << "Hello, World!" << std::endl;
    timespec start, finish, delta;
    make_matrix_gpu matrix_gpu;
    make_matrix_cpu matrix_cpu;



    int n = 9, m = 9;
    int T = 100, t3 = 100;
    double alpha = 0.05351f, radius = 1.5f;
    double eps = 0.01f;
    int iter_count = 1;
    bool gpu_test = false;
    int precision = 3;
    bool show_gpu = false;
    //  INPUT
    std::string yes;
    cout << "Run test? (y/n): ";
    if (!getline(cin>>std::ws,yes,'\n')) return true;
    gpu_test = (yes == "y");
    cout << "Show GPU information? (y/n): ";
    if (!getline(cin>>std::ws,yes,'\n')) return true;
    show_gpu = (yes == "y");
    cout << "Enter size (x,y,radius): "; cin >> m >> n >> radius;
    cout << "Enter temperature (T, t3): "; cin >> T >> t3;
    cout << "Enter alpha: "; cin >> alpha;
    cout << "Enter eps: "; cin >> eps;
    cout << "Enter maximum of iterations: "; cin >> iter_count;

    if(m%2 != 0 || n%2 != 0) {
        std::cerr << "WARNING: scary input (multiples of 2 are required)!" << endl;
    }

    cout << std::string(30,'=') << "GPU" << std::string(30,'=') << endl;
    matrix_gpu.n = n; matrix_gpu.m = m;
    matrix_gpu.T = T;
    matrix_gpu.t3 = t3;
    matrix_gpu.radius = radius;
    matrix_gpu.alpha = alpha;
    matrix_gpu.test = gpu_test;
    matrix_gpu.init_gpu();
    if(show_gpu) {
        matrix_gpu.device_info();
    }
    clock_gettime(CLOCK_REALTIME, &start);
    matrix_gpu.make_gpu();
    //matrix_gpu.show_matrix(precision);
    matrix_gpu.iteration(eps, iter_count);
    clock_gettime(CLOCK_REALTIME, &finish);
    sub_timespec(start, finish, &delta);
    matrix_gpu.write_file();
    cout << "Execute time = " << delta.tv_sec << "," << delta.tv_nsec << " took seconds\n\n";

    cout << std::string(30,'=') << "CPU" << std::string(30,'=') << endl;
    matrix_cpu.n = n; matrix_cpu.m = m;
    matrix_cpu.T = T;
    matrix_cpu.t3 = t3;
    matrix_cpu.radius = radius;
    matrix_cpu.alpha = alpha;
    clock_gettime(CLOCK_REALTIME, &start);
    matrix_cpu.make();
    //matrix_cpu.show_matrix();
    //cout << "---------" << endl;
    matrix_cpu.iteration(eps, iter_count);
    clock_gettime(CLOCK_REALTIME, &finish);
    sub_timespec(start, finish, &delta);
    matrix_cpu.write_file();
    cout << "Execute time = " << delta.tv_sec << "," << delta.tv_nsec << " took seconds\n\n";
    return 0;
}
