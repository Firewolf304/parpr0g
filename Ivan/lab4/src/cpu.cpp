#include "../include/includes.h"
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
void cpu_func (int N) {
    int sum;
    while (N > -1){
        if (test_polindrom(N)){
            for (int i = 2; i < sqrt(N); i++) {
                for (int j = i + 1; j < sqrt(N); j++) {
                    //sum = sum_square(i, j);
                    int sum = 0;
                    for(int a = i; a<=j; a++){
                        sum += pow(a,2);
                    }
                    //if (sum > N) break;
                    if (sum == N) {
                        cout << N << " " << i << " " << j << endl;
                        return;
                    }
                }
            }
        }
        N--;
    }
}
