#include <iostream>
#include <omp.h>
#include <vector>
#include <string>
#include <math.h>
#include <chrono>
#include <thread>
#include <fstream>


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
vector<double> jacobi (vector<vector<double>> A, vector<double> F, double eps = 0.001f)
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
vector<double> jacobi_omp (vector<vector<double>> A, vector<double> F, double eps = 0.001f)
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
std::vector<std::string> split(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}
int main() {
    std::cout << "Hello, World!" << std::endl;
    cout << "Insert random matrix (y/n): ";
    std::string yes;
    if (!getline(cin,yes,'\n')) return true;
    vector<vector<double>> a;
    vector<double> b;
    if (yes == "y") {
        int i = 2, maxval;
        cout << "Insert size of square matrix: ";
        cin >> i;
        cout << "Insert max values: ";
        cin >> maxval;
        a = matrix_genA(i, maxval);
        b = matrix_genB(i, maxval);
        std::ofstream matrix("matrix.txt", std::ios_base::trunc);
        cout << "Generated " << i << " size" << endl;
        for(int y = 0 ; y < b.size(); y++) {
            for(int x = 0 ; x < b.size(); x++) {
                if(x == b.size()-1) {
                    matrix << a[y][x] << " ";
                    matrix << b[y] << " ";
                }
                else {
                    matrix << a[y][x] << " ";
                }
            }
            matrix << "\n";
        }
        matrix.close();
    }
    else {
        int i = 2;
        cout << "Insert size of square matrix: ";
        cin >> i;
        a = vector<vector<double>>(i, vector<double>(i, 0));
        b = vector<double>(i, 0);
        for(int y =0; y < i; y++ ) {
            std::string line;
            getline(cin>>std::ws,line);
            vector<std::string> d = split(line, " ");
            for(int x = 0; x <= i; x++) { // т.к. учитываем, что в матрице есть вектор b
                if(x == i) {
                    b[y] = std::stoi(d[x]);
                }
                else {
                    a[y][x] = (double)std::stoi(d[x]);
                }
            }
            line.clear();
        }
        bool norm = normal(a);
        cout << "Inserted " << b.size() << " size its " << (norm ? "true" : "false") << endl;
        if(!norm) {
            perror("Infinite answer");
            exit(-1);
        }
    }
    double eps;
    cout << "Insert eps: ";
    cin >> eps;
    vector<double> out = jacobi(a, b, eps);
    vector<double> outomp = jacobi_omp(a, b, eps);
    a.clear();
    b.clear();
    std::ofstream answerlinel("answerlinel.txt", std::ios_base::trunc);
    std::ofstream answerparallel("answerparallel.txt", std::ios_base::trunc);
    int count=outomp.size();
    for(int i = 0; i < outomp.size(); i++) {
        answerlinel << out[i] << "\n";
        answerparallel << outomp[i] << "\n";
        if(out[i] == outomp[i] ) {
            count--;
        }
    }
    cout << "Matching answer: " << ((count == 0) ? "true" : "false" ) << endl;
    answerlinel.close();
    answerparallel.close();
    return 0;
}
