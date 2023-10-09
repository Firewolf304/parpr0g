#include <iostream>
#include <omp.h>
#include <vector>
#include <string>
#include <math.h>
#include <chrono>
#include <thread>


using std::vector;
using std::cout;
using std::cin;
using std::endl;
int Det(vector<vector<double>> matr, int n)
{
    int temp = 0;   //временная переменная для хранения определителя
    int k = 1;      //степень
    if(n < 1){
        cout<<"Not true size of matrix" << endl;
        return 0;
    }
    else if (n == 1)
        temp = matr[0][0];
    else if (n == 2)
        temp = matr[0][0] * matr[1][1] - matr[1][0] * matr[0][1];
    else{
        for(int i = 0; i < n; i++){
            int m = n - 1;
            vector<vector<double>>temp_matr(m, vector<double>(m, 0));
            temp = temp + k * matr[0][i] * Det(temp_matr, m);
            k = -k;
            //FreeMem(temp_matr, m);
        }
    }
    return temp;
}
vector<vector<double>> matrix_genA(int N, int max = 5) {
    srand(time(NULL));
    vector<vector<double>> r;
    for(int i = 0; i < N; i++) {
        r.push_back({});
        for(int j = 0; j < N; j++) {
            //std::this_thread::sleep_for(std::chrono::milliseconds(1));
            //int neg = ((std::rand() % 1) ? -1 : 1 );
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            r[i].push_back(   ((double)(std::rand() % max) / 10) );

        }
    }
    for(int y = 0; y < N; y++) {
        double summ = 0;
        for(int x = 0; x < N; x++) {
            if(x != y) {
                summ += r[y][x];
            }
        }
        if(r[y][y] < summ ) {
            srand(time(NULL));
            r[y][y] = summ + ((double)(std::rand() % max) / 10) + 0.1f; // 1 чтобы избавиться от 0, все равно рандом
        }
    }
    return r;
}
vector<double> matrix_genB(int N, int max = 100) {
    srand(time(NULL));
    vector<double> r;
    for(int j = 0; j < N; j++) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        r.push_back((double)(std::rand() % max)/ 10);
    }
    return r;
}
bool normal(vector<vector<double>> mass) {
    for(int y = 0; y < mass.size(); y++) {
        double summ = 0;
        for(int x = 0; x < mass[y].size(); x++) {
            if(x != y) {
                summ += mass[y][x];
            }
        }
        if(mass[y][y] < summ ) {
            return false;
        }
    }
    return true;
}
vector<double> jacobi (vector<vector<double>> A, vector<double> F, const double& eps = 0.001)
/*
 * A - переменные слау
 * F - переменные после равно
 */
{
    double itime, ftime, exec_time;
    itime = omp_get_wtime();
    vector<double> out(F.size(), 0);   // X
    vector<double> TempX(F.size());    // временный x
    double norm = 1;
    while(norm > eps) {
        for (int i = 0; i < F.size(); i++) {
            TempX[i] = (double)F[i];
            for (int j = 0; j < F.size(); j++) {
                if (i != j) {
                    TempX[i] -= (double)(A[i][j] * out[j]);
                }
            }
            TempX[i] /= (double)A[i][i];
        }
        norm = fabs(out[0] - TempX[0]);
        for (int i = 0; i < F.size(); i++) {
            if (fabs(out[i] - TempX[i]) > norm) {
                norm = fabs(out[i] - TempX[i]);
            }
            out[i] = (double)TempX[i];
        }
    }
    ftime = omp_get_wtime();
    exec_time = ftime - itime;
    cout << "TIME linel: " << exec_time << endl;
    return out;
}
vector<double> jacobi_omp (vector<vector<double>> A, vector<double> F, const double& eps = 0.001)
/*
 * A - переменные слау
 * F - переменные после равно
 */
{
    double itime, ftime, exec_time;
    itime = omp_get_wtime();
    vector<double> out(F.size(), 0);   // X
    vector<double> TempX(F.size());    // временный x
    double norm = 1;
    while(norm > eps) {
        int i = 0;
//#pragma omp parallel shared(norm, out, TempX) private(i)
    #pragma omp parallel for shared(TempX) private(i)
            for (i = 0; i < F.size(); i++) {
                TempX[i] = (double)F[i];
                for (int j = 0; j < F.size(); j++) {
                    if (i != j) {
                        TempX[i] -= (double)(A[i][j] * out[j]);
                    }
                }
                TempX[i] /= (double)A[i][i];
            }
            norm = fabs(out[0] - TempX[0]);
    #pragma omp parallel for shared(norm, out) private(i)
            for (i = 0; i < F.size(); i++) {
                if (fabs(out[i] - TempX[i]) > norm) {
                    norm = fabs(out[i] - TempX[i]);
                }
                out[i] = (double)TempX[i];
            }
    }
    ftime = omp_get_wtime();
    exec_time = ftime - itime;
    cout << "TIME parallel: " << exec_time << endl;
    return out;
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    //vector<vector<double>> a = {{9,0,-12},{0,51,12},{-12,5,24}};
    //vector<vector<double>> a = {{9.2f,2.5f,-3.7f},{0.9f,9.0f,0.2f},{4.5f,-1.6f,-10.3f}};
    //vector<double> b = {-17.5f,4.4f,-22.1f};
    int i = 2, maxval;
    cout << "Insert size of square matrix: ";
    cin >> i;
    cout << "Insert max values: ";
    cin >> maxval;
    vector<vector<double>> a = matrix_genA(i,maxval);
    vector<double> b = matrix_genB(i,maxval);
    bool check = normal(a);
    while(!normal(a)) {
        for(int y = 0 ; y < i; y++) {
            for(int x = 0 ; x < i; x++) {
                if(x == i-1) {
                    cout << a[y][x] << " ";
                    cout << b[y] << " ";
                }
                else {
                    cout << a[y][x] << " ";
                }
            }
            cout << endl;
        }
        cout << "Recheck" << endl;
        a = matrix_genA(i,maxval);
        b = matrix_genB(i,maxval);
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }

    cout << "Generated " << i << " size its " << (normal(a) ? true : false) << endl;
    for(int y = 0 ; y < i; y++) {
        for(int x = 0 ; x < i; x++) {
            if(x == i-1) {
                cout << a[y][x] << " ";
                cout << b[y] << " ";
            }
            else {
                cout << a[y][x] << " ";
            }
        }
        cout << endl;
    }
    vector<double> out = jacobi(a, b);
    vector<double> outomp = jacobi_omp(a, b);
    std::cout << "ANSWER:\nlinel" << std::endl;
    for(auto d : out) {
        cout << d << endl;
    }
    std::cout << "parallel" << std::endl;
    for(auto d : outomp) {
        cout << d << endl;
    }
    /*for(int i = 2; i <= 300; i++) {
        vector<vector<double>> a = matrix_genA(i);
        vector<double> b = matrix_genB(i);
        cout << "Generated " << i << endl;
        vector<double> out = jacobi(a, b);
        vector<double> outomp = jacobi_omp(a, b);
    }*/
    /*for(auto d : out) {
        cout << d << endl;
    }
    std::cout << "parallel" << std::endl;
    for(auto d : outomp) {
        cout << d << endl;
    }*/
    return 0;
}
