#include "../include/includes.h"
int main() {
    int N = 500000; // 4334 37 39 // 818 2 13
    cout << "Insert N: ";
    cin >> N;
    timer time;
    calculate_gpu gpu;
    gpu.init();
    time << timer::START;
    cpu_func(N);
    time << timer::FINISH;
    cout << "Execute time CPU: " << time.sub_timespec() << endl;
    time << timer::START;
    gpu.calculate(N);
    time << timer::FINISH;
    cout << "Execute time GPU: " << time.sub_timespec() << endl;
    return 0;
}