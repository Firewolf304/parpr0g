#include "../include/includes.h"

int main() {
    timer time;
    string text = "xyxxxxxzxyx";
    int n = 3;
    std::cout << "Hello, World!" << std::endl;
    cout << "Insert text: ";
    getline(cin>>std::ws,text,'\n');
    cout << "Insert substring: ";
    cin >> n;
    matrix_gpu gpu(n, text);
    time << timer::START;
    gpu.get_number();
    time << timer::FINISH;
    cout << "Execute time: " << time.sub_timespec() << endl;
    matrix_cpu cpu(n, text);
    time << timer::START;
    cpu.get_number();
    time << timer::FINISH;
    cout << "Execute time: " << time.sub_timespec() << endl;
    return 0;
}
