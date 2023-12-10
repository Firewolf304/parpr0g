#include <iostream>
#include <string>
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include <stack>

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
string str = "xyxxz";
//string str = "abcde";

int n = 3, count = 0;

int cycle_cpu(string piece, int a){
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
}
void process_stack() {
    std::vector<int> offset (n,0);
    std::vector<int> iter (n,0);
    int maxSIZE = str.length() - n;
    string line(n,' ');
    for(int i = 0; i < n; i++) {offset[i] = i; line[i] = str[i]; }
    int c = 0;
    for(auto dd : iter ) {
        cout << dd << " ";
    }
    cout << endl;
    while(c != n) {
        if(iter.back() + 1 <= maxSIZE) {
            iter.back()++;
        }
        else {
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
        for(auto dd : iter ) {
            cout << dd << " ";
        }
        cout << endl;
        c = 0;
        for(auto d : iter) {
            if(d == maxSIZE) {
                c++;
            }
        }

    }
}
__global__ void cycle_gpu(string piece, int a) {

}

int main() {
    string piece = "";
    cycle_cpu("",-1);
    process_stack();
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
