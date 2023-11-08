#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <omp.h>

using namespace std;


void gen_matrixA(double **a, int n){
    srand(time(NULL));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if (i == j)
                a[i][j] = double(rand() % 10 + 100) / 10;
            else
                a[i][j] = double(rand() % 10) / 10;
        }
    }
}

void gen_matrixB(double *b, int n){
    srand(time(NULL)+1);
    for(int i = 0; i < n; i++){
        b[i] = double(rand() % 100) / 10;
    }
}

int main() {
    int n, powEps;
    double eps,w, beginTime;

    ifstream f("matrix.txt");

    cout << "Matrix scale:";
    cin >> n;
    cout << "eps = e * 10 ^ -";
    cin >> powEps;
    eps = pow(10, -1*(powEps+1));
    cout << "Omega:";
    cin >> w;
    /*n=3;
    eps = 0.00001;
    w = 0.9;*/

    double **a = new double *[n+1]; //Выделяем память под строки
    for (int i = 0; i < n; i++)
        a[i] = new double[n+1];
    double b[n + 1];
    double x[n + 1];
    gen_matrixA(a, n);
    gen_matrixB(b, n);
    ofstream out1 ("matrix.txt");
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            out1 << a[i][j] << ' ';
            x[i] = b[i]/ a[i][i];
        }
        out1 << b[i] << endl;
    }
    out1.close();



    beginTime = omp_get_wtime();
    double norm;
    do {
        norm=0;
        for (int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 1; j < n; ++j) {
                if (i != j) {
                    sum = sum + a[i][j] * x[j];
                }
            }
            double Next = (1 - w) * x[i] + w * (b[i] - sum) / a[i][i];

            if (fabs(Next - x[i]) > norm) {
                norm = fabs(Next - x[i]);
            }
            x[i] = Next;
        }
    } while (norm > eps);
    cout<<"Time for linear method: "<<beginTime - omp_get_wtime();

    ofstream out2 ("answerLinerial.txt");
    for (int i = 0; i < n; i++) {
        out2 << setprecision(powEps) << x[i] << ' ';
        x[i] = b[i]/ a[i][i];
    }
    out2.close();


    beginTime = omp_get_wtime();
    do {
        norm=0;
        int i = 0;
#pragma omp parallel for shared(a, x, norm) private(i)
        for (i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 1; j < n; ++j) {
                if (i != j) {
                    sum = sum + a[i][j] * x[j];
                }
            }
            double Next = (1 - w) * x[i] + w * (b[i] - sum) / a[i][i];

            if (fabs(Next - x[i]) > norm) {
                norm = fabs(Next - x[i]);
            }
            x[i] = Next;
        }
    } while (norm > eps);
    cout<<"Time for parallel method: "<<beginTime - omp_get_wtime();

    ofstream out ("answerParallel.txt");
    for (int i = 0; i < n; i++) {
        out << setprecision(powEps) << x[i] << ' ';
        x[i] = 0;
    }
    out.close();

}

