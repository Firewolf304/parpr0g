#include "../include/includes.h"

class calculate_gpu {
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
    bool test_polindrom(int a){
        int countDigit = 0, t = 0;
        countDigit -= counter_digit(a) % 2 == 1;
        while (a>=9){
            t += first_digit(a) != a%10;
            a -= first_digit(a) * pow(10, counter_digit(a)-1);
            a /= 10;
        }
        return (!t);
    }
    int sum_square(int a, int b){
        int sum = 0;
        for(a; a<=b; a++){
            sum += pow(a,2);
        }
        return sum;
    }
    int counter_digit(int a){
        int c = 0;
        while (a){
            c++;
            a /= 10;
        }
        return c;
    }

    int first_digit(int a){
        int c = counter_digit(a)-1;
        for(int i = 0; i < c; i++){
            a/=10;
        }
        return a;
    }
public:
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
    string kernel = "__kernel void kernelfunc(__global int * N, __global int * retI, __global int * retJ) {"
                    "   int i = get_global_id(0);"
                    "   int j = get_global_id(1);\n"
                    "   __private int sum = 0;\n"
                    "   if(i < 2 || j < (i + 1)) {return;}\n"
                    "   for(int a = i; a<=j; a++){ "
                    "       sum += pow((float)a,2);"
                    "   }"

                    "   if (sum > *N) return;\n"
                    "   if (sum == *N) {\n"
                    //"       printf(\"Getted %d %d %d | %d\\n\", *N, i, j, sum); "
                    "       *retI = i; *retJ = j;"
                    "       return;\n"
                    "   }"
                    "}";

    void init() {
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
    }

    void calculate(int N ) {
        int sum;
        cl_mem NBuff;
        int i = -1, j = -1, Nr;
        cl_mem retI = clCreateBuffer (this->cntx , CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR ,sizeof ( int ), &i , NULL );
        cl_mem retJ = clCreateBuffer (this->cntx , CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR ,sizeof ( int ), &j , NULL );
        checkError(clSetKernelArg (k, 1, sizeof ( cl_mem ), &retI), "Invalid arg");
        checkError(clSetKernelArg (k, 2, sizeof ( cl_mem ), &retJ), "Invalid arg");
        while (N > -1){
            if (test_polindrom(N)){
                size_t global_work_size [2] = {(size_t)sqrt(N), (size_t)sqrt(N)};
                //cout << (size_t)sqrt(N) << " " << N << endl;
                NBuff = clCreateBuffer (this->cntx , CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof ( int ), &N , NULL );
                checkError(clSetKernelArg (k, 0, sizeof ( cl_mem ), &NBuff), "Invalid arg");
                kEvent = clCreateUserEvent(cntx, &clerr);
                checkError(clEnqueueNDRangeKernel (cq, k, 2, NULL, global_work_size , NULL /*local_work_size*/ , 0, NULL , &kEvent ), "Error start");
                checkError(clEnqueueReadBuffer (cq , retI, CL_TRUE , 0,	 sizeof ( i ), &i , 0, NULL , NULL ), "error get value");
                checkError(clEnqueueReadBuffer (cq , retJ, CL_TRUE , 0,	 sizeof ( j ), &j , 0, NULL , NULL ), "error get value");
                if(i > -1 && j >-1 ) break;
            }
            N--;
        }
        cout << N << " " << i << " " << j << endl;
    }
};