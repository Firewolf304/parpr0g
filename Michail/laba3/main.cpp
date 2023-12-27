#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <iostream>
#include <vector>
#include <CL/cl.h>
using std::cout;
using std::cin;
using std::endl;
using std::vector;

void sieve_of_Eratosthenes(int n){
    //unsigned long long *numbers = (unsigned long long*) calloc(n, sizeof(unsigned long long));
    //unsigned long long numbers[n];
    //auto * numbers = new unsigned long long[n];
    vector<unsigned long long> numbers(n,0);
    int countPrimal = 0;
    for(int i = 0; i < n; i++){
        numbers[i] = i;}

    for(int i = 2; i <= sqrt(n); i++){ //Решето Эратосфена
        for (int j = i * 2; j < n; j += i){
            numbers[j] = 0;
        }
    }
    for(int i = 3; i < n; i += 2){ //Считаем кол-во простых чисел для инициализации массива простых чисел
        countPrimal += (numbers[i] != 0);}

    vector<int> primal(countPrimal+1, 0);
    countPrimal = 0;
    primal[countPrimal++] = 2;

    for(int i = 3; i < n; i += 2){  //Заполняем массив простыми числами
        if (numbers[i] != 0){
            primal[countPrimal] = numbers[i];
            //cout<<primal[countPrimal]<<" ";
            countPrimal++;
        }
    }

//22.Найти числа a>b>c, большие заданного числа A,
//имеющие одинаковые количества различных простых делителей такие, что a-b=b-c.

    int A = 1300;
    A += (A%2 == 0);
    int C = A + 1;
    int B;
    int i = 0;
    int countA = 0, countC = 0, countB = 0;
    int N = 5;
    while(1){
        B = (A + C) / 2;
        while(primal[i] * 2 < (B)){
            if (A % primal[i] == 0 && B % primal[i] == 0 || C % primal[i] == 0 && B % primal[i] == 0 || C % primal[i] == 0 && A % primal[i] == 0)
                break;
            countA += (A % primal[i]) == 0; countC += (C % primal[i]) == 0; countB += ((B) % primal[i]) == 0;
            i++;
        }
        if (countA == countC and countA == countB and countA != 0 && B != A && B != C && A-B==B-C){
            cout << "A = "<< A << "; B = " << B << "; C = " << C << endl;
            break;
        }

        i = 0;
        countA = 0, countC = 0, countB = 0;
        A += 2 * (C % 1000 == 0);
        C += 1 + (A - C) * (C % 1000 == 0); //Оптимизация ифа сверху
    }
    cout<<"Count div = "<<countA<<endl; //Изменения от сюда(68-89) строка
    cout<< A <<" ";
    int j = 0, count = 0;
    while (j < countA){
        if ((A % primal[count]) == 0) {cout << primal[count]<<" ";
            j++;}
        count++;
    }
    j = 0, count = 0;
    cout<<endl << (B) <<" ";
    while (j < countA){
        if ((B % primal[count]) == 0) {cout << primal[count]<<" ";
            j++;}
        count++;
    }
    j = 0, count = 0;
    cout<<endl << C <<" ";
    while (j < countA){
        if ((C % primal[count]) == 0) {cout << primal[count]<<" ";
            j++;}
        count++;
    } //и до сюда

}
char* errCodeToString(int err){
    switch (err) {
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
const void checkError(cl_int err, const std::string& message) {
    if (err != CL_SUCCESS) {
        std::cerr << "Ошибка '" << message << "' (" << errCodeToString(err) << " = " << err << ")" << std::endl;
        exit(err);
    }
}

void sieve_of_Eratosthenes_gpu(int n){
    cl_int clerr;
    cl_uint qty_platforms = 0;
    cl_platform_id* platforms;
    cl_uint ui;
    cl_uint * qty_devices;
    cl_device_id **devices;
    cl_context cntx;
    cl_event kEvent;
    cl_command_queue cq;
    cl_device_id * ds; // devices
    cl_program p;
    cl_kernel k;
    const char * kernel = "__kernel void kernelfunc(__global int* Aout, __global int* primal, __global int* countAout,  __global int* Bout, __global int* Cout, __global int* min, __global int* offset) {"
                          "int A = get_global_id(0) + *offset;\n"
                          "int C = A + 1;\n"
                          "int B;\n"
                          "int i = 0;\n"
                          "int countA = 0, countC = 0, countB = 0;\n"
                          "int flag = 0;"
                          "for (C = A+1; C < A - (A % 1000) + 1000; C++){\n"
                          "            B = (A + C) / 2;\n"
                          "            while(primal[i] * 2 <= (C)){"
                          "                if (A % primal[i] == 0 && B % primal[i] == 0 || C % primal[i] == 0 && B % primal[i] == 0 || C % primal[i] == 0 && A % primal[i] == 0)\n"
                          "                    break;\n"
                          "                countA += (A % primal[i]) == 0; countC += (C % primal[i]) == 0; countB += ((B) % primal[i]) == 0;\n"
                          "                i++;\n"
                          "            }\n"
                          "            i = 0;\n"
                          "            if (countA == countC && countA == countB && countA != 0 && B != A && B != C && A-B==B-C){\n"
                          "                flag = 1;"
                          "                break;\n"
                          "            }\n"
                          "            countA = 0, countC = 0, countB = 0;\n"
                          "\n"
                          "} "
                          "if(flag) {"
                          "Aout[get_global_id(0)] = A; "
                          "Bout[get_global_id(0)] = countA; "
                          "Cout[get_global_id(0)] = C;"
                          "*countAout = countA;"
                          ""
                          "}"
                          //"Aout[0] = 9999;"
                          "}";

    // INIT GPU
    checkError(clGetPlatformIDs(0, NULL, &qty_platforms), "Error get platforms");
    qty_devices = (cl_uint*)malloc(sizeof(cl_uint)*qty_platforms);
    platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id)*qty_platforms);
    devices = (cl_device_id**)malloc(sizeof(cl_device_id*)*qty_platforms);
    checkError(clGetPlatformIDs(qty_platforms, platforms, NULL), "Error get platforms");
    for (ui=0; ui < qty_platforms; ui++){
        checkError( clGetDeviceIDs(platforms[ui], CL_DEVICE_TYPE_ALL,0, NULL, &qty_devices[ui]), "Error get device ID");
        if(qty_devices[ui]){
            devices[ui] = (cl_device_id*)malloc(qty_devices[ui] * sizeof(cl_device_id));
            checkError(clGetDeviceIDs(platforms[ui], CL_DEVICE_TYPE_ALL, qty_devices[ui], devices[ui], NULL), "Error get device ID");
        }
    }
    cntx = clCreateContext(0, qty_devices[0], devices[0], NULL, NULL, &clerr);
    cq = clCreateCommandQueue(cntx, devices[0][0], CL_QUEUE_PROFILING_ENABLE, &clerr);
    p = clCreateProgramWithSource(cntx, 1, (const char**)&kernel, NULL, &clerr);
    clerr = clBuildProgram(p, 0, NULL,  NULL, NULL, NULL);
    if (clerr == CL_BUILD_PROGRAM_FAILURE) {
        std::cerr << "Build error: " << errCodeToString(clerr) << endl;
        size_t log_size;
        clGetProgramBuildInfo(p, devices[0][0], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
        char *log = (char *) malloc(log_size);
        clGetProgramBuildInfo(p, devices[0][0], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
        cout << log << endl;
        exit(-1);
    }
    k = clCreateKernel(p, "kernelfunc", &clerr);
    checkError(clerr, "Error kernel");


    vector<unsigned long long> numbers(n,0);
    int countPrimal = 0;
    for(int i = 0; i < n; i++){
        numbers[i] = i;}

    for(int i = 2; i <= sqrt(n); i++){ //Решето Эратосфена
        for (int j = i * 2; j < n; j += i){
            numbers[j] = 0;
        }
    }
    for(int i = 3; i < n; i += 2){ //Считаем кол-во простых чисел для инициализации массива простых чисел
        countPrimal += (numbers[i] != 0);}

    //int primal[countPrimal+1];
    int * primal = new int[countPrimal+1];
    countPrimal = 0;
    primal[countPrimal++] = 2;

    for(int i = 3; i < n; i += 2){  //Заполняем массив простыми числами
        if (numbers[i] != 0){
            primal[countPrimal] = numbers[i];
            //cout<<primal[countPrimal]<<" ";
            countPrimal++;
        }
    }

//22.Найти числа a>b>c, большие заданного числа A,
//имеющие одинаковые количества различных простых делителей такие, что a-b=b-c.

    int A = 1300;
    A += (A%2 == 0);
    int B = 0;
    int countA = 0, countC = 0, countB = 0;
    int C = 0;

    //vector<int> valuesA(1000, 0);
    //vector<int> valuesB(1000, 0);
    //vector<int> valuesC(1000, 0);
    int * valuesA = (int*)calloc( 1000, sizeof(int));
    int * valuesB = (int*)calloc( 1000, sizeof(int));
    int * valuesC = (int*)calloc( 1000, sizeof(int));
    /*cl_mem Ainp = clCreateBuffer (cntx , CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR ,
        sizeof ( A ), &A , NULL );
    cl_mem Binp = clCreateBuffer (cntx , CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR ,
        sizeof ( B ), &B , NULL );
    cl_mem Cinp = clCreateBuffer (cntx , CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR ,
        sizeof ( C ), &C , NULL );*/
    cl_mem Ainp = clCreateBuffer (cntx , CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR ,
                                  1000 * sizeof ( int ), valuesA , NULL );
    cl_mem Binp = clCreateBuffer (cntx , CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR ,
                                  1000 * sizeof ( int ), valuesB , NULL );
    cl_mem Cinp = clCreateBuffer (cntx , CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR ,
                                  1000 * sizeof ( int ), valuesC , NULL );
    cl_mem primaleInp = clCreateBuffer (cntx , CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR ,
                                        (countPrimal + 1) * sizeof(int), primal , &clerr );
    checkError(clerr, "Invalid input");
    cl_mem countAinp = clCreateBuffer (cntx , CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR ,
                                  sizeof ( int ), &countA , NULL );
    int min = 999999;
    cl_mem minn = clCreateBuffer (cntx , CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR ,
                                       sizeof ( int ), &min , NULL );
    cl_mem offset = clCreateBuffer (cntx , CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR ,
                                    sizeof ( int ), &A , NULL );

    //__global int* Aout, __local int* primal, __global int* countAout, __global int* Bout, __global int* Cout
    checkError(clSetKernelArg (k, 0, sizeof ( cl_mem ), &Ainp), "Invalid arg");
    checkError(clSetKernelArg (k, 1, sizeof ( cl_mem ), &primaleInp), "Invalid arg");
    checkError(clSetKernelArg (k, 2, sizeof ( cl_mem ), &countAinp), "Invalid arg");
    checkError(clSetKernelArg (k, 3, sizeof ( cl_mem ), &Binp), "Invalid arg");
    checkError(clSetKernelArg (k, 4, sizeof ( cl_mem ), &Cinp), "Invalid arg");
    checkError(clSetKernelArg (k, 5, sizeof ( cl_mem ), &minn), "Invalid arg");
    checkError(clSetKernelArg (k, 6, sizeof ( cl_mem ), &offset), "Invalid arg");

    kEvent = clCreateUserEvent(cntx, &clerr);
    size_t global_work_size [1] = {1000};
    //size_t offset [1] = {static_cast<size_t>(A)};

    checkError(clEnqueueNDRangeKernel (cq, k, 1, NULL, global_work_size , NULL /*local_work_size*/ , 0, NULL , &kEvent ), "Error start");

    checkError(clEnqueueReadBuffer (cq , Ainp, CL_TRUE , 0,	 1000 * sizeof ( int ), valuesA , 0, NULL , NULL ), "error get value");
    checkError(clEnqueueReadBuffer (cq , Binp, CL_TRUE , 0,	 1000 * sizeof ( int ), valuesB , 0, NULL , NULL ), "error get value");
    checkError(clEnqueueReadBuffer (cq , Cinp, CL_TRUE , 0,	 1000 * sizeof ( int ), valuesC , 0, NULL , NULL ), "error get value");
    checkError(clEnqueueReadBuffer (cq , countAinp, CL_TRUE , 0,	 sizeof ( countA ), &countA , 0, NULL , NULL ), "error get value");
    cout <<"A = "<< A << "; B = " << B << "; C = " << C << "; test=" << valuesA[0] <<endl;
    cout <<"A = "<< A << "; B = " << B << "; C = " << C<<endl;
    clReleaseMemObject(Ainp);
    clReleaseMemObject(Binp);
    clReleaseMemObject(Cinp);
    clReleaseMemObject(countAinp);
    clReleaseKernel (k);
    clReleaseProgram (p);
    clReleaseCommandQueue (cq );
    clReleaseContext (cntx );
    int index = 0;
    for(; index < 1000; index++) {
        if(valuesA[index] != 0) {
            break;
        }
    }
    A = valuesA[index]; countA = valuesB[index]; C = valuesC[index];
    B = (A+C)/2;
    cout<<"Count div = "<<countA<<endl; //Изменения от сюда(68-89) строка
    cout<< A <<" ";
    int j = 0, count = 0;
    while (j < countA && count < countPrimal){
        if ((A % primal[count]) == 0) {cout << primal[count]<<" ";
            j++;}
        count++;
    }
    j = 0, count = 0;
    cout<<endl << (B) <<" ";
    while (j < countA){
        if ((B % primal[count]) == 0) {cout << primal[count]<<" ";
            j++;}
        count++;
    }
    j = 0, count = 0;
    cout<<endl << C <<" ";
    while (j < countA){
        if ((C % primal[count]) == 0) {cout << primal[count]<<" ";
            j++;}
        count++;
    } //и до сюда
    cout<<endl;

}

int main() {
    sieve_of_Eratosthenes_gpu(10000);
    sieve_of_Eratosthenes(10000);
}
