/*
 refactoring code to:
    #include <stdio.h>
    #include <CL/cl_platform.h>
    #include <CL/cl.h>
    #include <malloc.h>
 */
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <csignal>
#include <time.h>
#include <unistd.h>
#include <fstream>
#include <CL/cl_platform.h>
#include <CL/cl.h>
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
    cl_int clerr;
    cl_platform_id* platforms;
    cl_uint qty_platforms = 0;
    cl_platform_id platform;
    cl_uint *qty_devices;
    cl_device_id ** devices;
    cl_device_id device;
    cl_context context;
    cl_command_queue commandList;
    /*__global double* output, int m, int n, float step, double radius, int t3*/
    char* makeMatrix =   "__kernel void makeMatrix(__global double* output, int m, int n, float step, double radius, int t3) {"
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
    /*0__global double* input, 1__global float* output, 2__private int xpos, 3__private int ypos, 4int m, 5int n, 6double radius, __global 7double * testX, __global 8double * testY*/
    char* programRun =  " __kernel void collectNeighbors(__global double* input, __global double* output, int xpos, int ypos, int m, int n, double radius, __global double * testX, __global double * testY) {"
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
    //cl::Buffer testBuffX;
    //cl::Buffer testBuffY;
    cl_program ReadyProgramRun;
    cl_kernel ReadyKernelRun;
    //cl::NDRange Readyoffset;
    //cl::NDRange ReadylocalSize;
    //cl::NDRange ReadyglobalSizeNDRange;
    //=================GPU=================

    const void checkError(int err, const std::string& message) {
        if (err != CL_SUCCESS) {
            std::cerr << message << ": " << errCodeToString(err) << " (" << err << ")" << std::endl;
            exit(err);
        }
    }
    char* errCodeToString(int err){
        switch (err) {
            case CL_SUCCESS:                            return "Success!";
            case CL_DEVICE_NOT_FOUND:                   return "Device not found.";
            case CL_DEVICE_NOT_AVAILABLE:               return "Device not available";
            case CL_COMPILER_NOT_AVAILABLE:             return "Compiler not available";
            case CL_MEM_OBJECT_ALLOCATION_FAILURE:      return "Memory object allocation failure";
            case CL_OUT_OF_RESOURCES:                   return "Out of resources";
            case CL_OUT_OF_HOST_MEMORY:                 return "Out of host memory";
            case CL_PROFILING_INFO_NOT_AVAILABLE:       return "Profiling information not available";
            case CL_MEM_COPY_OVERLAP:                   return "Memory copy overlap";
            case CL_IMAGE_FORMAT_MISMATCH:              return "Image format mismatch";
            case CL_IMAGE_FORMAT_NOT_SUPPORTED:         return "Image format not supported";
            case CL_BUILD_PROGRAM_FAILURE:              return "Program build failure";
            case CL_MAP_FAILURE:                        return "Map failure";
            case CL_INVALID_VALUE:                      return "Invalid value";
            case CL_INVALID_DEVICE_TYPE:                return "Invalid device type";
            case CL_INVALID_PLATFORM:                   return "Invalid platform";
            case CL_INVALID_DEVICE:                     return "Invalid device";
            case CL_INVALID_CONTEXT:                    return "Invalid context";
            case CL_INVALID_QUEUE_PROPERTIES:           return "Invalid queue properties";
            case CL_INVALID_COMMAND_QUEUE:              return "Invalid command queue";
            case CL_INVALID_HOST_PTR:                   return "Invalid host pointer";
            case CL_INVALID_MEM_OBJECT:                 return "Invalid memory object";
            case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    return "Invalid image format descriptor";
            case CL_INVALID_IMAGE_SIZE:                 return "Invalid image size";
            case CL_INVALID_SAMPLER:                    return "Invalid sampler";
            case CL_INVALID_BINARY:                     return "Invalid binary";
            case CL_INVALID_BUILD_OPTIONS:              return "Invalid build options";
            case CL_INVALID_PROGRAM:                    return "Invalid program";
            case CL_INVALID_PROGRAM_EXECUTABLE:         return "Invalid program executable";
            case CL_INVALID_KERNEL_NAME:                return "Invalid kernel name";
            case CL_INVALID_KERNEL_DEFINITION:          return "Invalid kernel definition";
            case CL_INVALID_KERNEL:                     return "Invalid kernel";
            case CL_INVALID_ARG_INDEX:                  return "Invalid argument index";
            case CL_INVALID_ARG_VALUE:                  return "Invalid argument value";
            case CL_INVALID_ARG_SIZE:                   return "Invalid argument size";
            case CL_INVALID_KERNEL_ARGS:                return "Invalid kernel arguments";
            case CL_INVALID_WORK_DIMENSION:             return "Invalid work dimension";
            case CL_INVALID_WORK_GROUP_SIZE:            return "Invalid work group size";
            case CL_INVALID_WORK_ITEM_SIZE:             return "Invalid work item size";
            case CL_INVALID_GLOBAL_OFFSET:              return "Invalid global offset";
            case CL_INVALID_EVENT_WAIT_LIST:            return "Invalid event wait list";
            case CL_INVALID_EVENT:                      return "Invalid event";
            case CL_INVALID_OPERATION:                  return "Invalid operation";
            case CL_INVALID_GL_OBJECT:                  return "Invalid OpenGL object";
            case CL_INVALID_BUFFER_SIZE:                return "Invalid buffer size";
            case CL_INVALID_MIP_LEVEL:                  return "Invalid mip-map level";
            default: return "Unknown";
        }
    }
    void init_gpu(int IDplatform = 0, int IDdevice = 0) {
        clerr = CL_SUCCESS;
        cl_uint ui;
        this->koefX = this->m / 2;
        this->koefY = this->n / 2;
        this->testX = vector<double>(12, -1);
        this->testY = vector<double>(12, -1);
        this->matrix = vector<double>(this->n * this->m, 0);
        clerr = clGetPlatformIDs(0, NULL, &qty_platforms);
        platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id)*qty_platforms);
        devices = (cl_device_id**)malloc(sizeof(cl_device_id*)*qty_platforms);
        qty_devices = (cl_uint*)malloc(sizeof(cl_uint)*qty_platforms);
        checkError( clGetPlatformIDs(qty_platforms, platforms, NULL), "No platforms" );
        for (ui=0; ui < qty_platforms; ui++){
            checkError( clGetDeviceIDs(platforms[ui], CL_DEVICE_TYPE_ALL,0, NULL, &qty_devices[ui]), "No devices" );
            if(qty_devices[ui]){
                devices[ui] = (cl_device_id*)malloc(qty_devices[ui] * sizeof(cl_device_id));
                checkError(clGetDeviceIDs(platforms[ui], CL_DEVICE_TYPE_ALL, qty_devices[ui], devices[ui], NULL), "Error bind device");
            }
        }
        //init gpu program
        this->context = clCreateContext(0, qty_devices[0], devices[0], NULL, NULL, &clerr);
        this->commandList = clCreateCommandQueue(this->context, devices[0][0], CL_QUEUE_PROFILING_ENABLE, &clerr);
        this->ReadyProgramRun = clCreateProgramWithSource(this->context, 1, (const char**)&(programRun), NULL, &clerr);
        checkError(clBuildProgram(this->ReadyProgramRun, 0, NULL,  NULL /*clcompileflags*/, NULL, NULL), "Error compile");
        size_t log_size;
        checkError(clGetProgramBuildInfo(this->ReadyProgramRun, devices[0][0], CL_PROGRAM_BUILD_LOG, 0, NULL, NULL), "Error clGetProgramBuildInfo");
        this->ReadyKernelRun = clCreateKernel(this->ReadyProgramRun, "collectNeighbors", &clerr);
        checkError(clSetKernelArg(this->ReadyKernelRun, 4, sizeof ( int ), &( this->m)), "Error bind arg" );
        checkError(clSetKernelArg(this->ReadyKernelRun, 5, sizeof ( int ), &( this->n)), "Error bind arg" );
        checkError(clSetKernelArg(this->ReadyKernelRun, 6, sizeof ( double ), &( this->radius)), "Error bind arg" );
    }
    void make_gpu() {
        float step = (float)(T / (n-1));
        matrix[ index(this->m-1,0) ] = this->T;
        for (int i = 1; i < this->n; i++) {
            this->matrix [ index(this->m-1,i) ] = (double) (this->matrix[index(this->m-1,i-1)] - step);
        }
        cl_mem outputBuffer = clCreateBuffer(this->context, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, this->matrix.size() * sizeof(double), this->matrix.data(), NULL );
        vector<size_t> global = { (size_t)this->m, (size_t)this->n }; // ПОЗИЦИЯ ДОЛЖНА БЫТЬ РАВНА РАСЧЕТУ index!
        cl_program program = clCreateProgramWithSource(this->context, 1, (const char**)&makeMatrix, NULL,  &clerr );
        checkError(clBuildProgram(program, 0, NULL,  NULL /*clcompileflags*/, NULL, NULL), "Error compile make_gpu");
        cl_kernel kernel = clCreateKernel(program, "makeMatrix", &clerr);

        checkError(clSetKernelArg(kernel, 0, sizeof ( cl_mem ), &( outputBuffer )), "Error bind arg" );
        checkError(clSetKernelArg(kernel, 1, sizeof ( int ), &( this->m )), "Error bind arg" );
        checkError(clSetKernelArg(kernel, 2, sizeof ( int ), &( this->n )), "Error bind arg" );
        checkError(clSetKernelArg(kernel, 3, sizeof ( float ), &( step )), "Error bind arg" );
        checkError(clSetKernelArg(kernel, 4, sizeof ( double ), &( this->radius )), "Error bind arg" );
        checkError(clSetKernelArg(kernel, 5, sizeof ( int ), &( this->t3 )), "Error bind arg" );
        


        /*float step = (float)(T / (n-1));
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
    */
    }
    /*double valueGet(int x, int y) {
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
    }*/
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
    /*std::string yes;
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
    cout << "Enter maximum of iterations: "; cin >> iter_count;*/

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
    clock_gettime(CLOCK_REALTIME, &start);
    matrix_gpu.make_gpu();

}