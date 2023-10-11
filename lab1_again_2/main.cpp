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
vector<vector<double>> generatematrixA(int N, int max = 5) {
    vector<vector<double>> r;

    for(int i = 0; i < N; i++) {
        r.push_back({});
        for(int j = 0; j < N; j++) {
            srand(time(NULL));
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
vector<double> generatematrixB(int N, int max = 100) {
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
vector<double> jacobi (vector<vector<double>> A, vector<double> F, double eps = 0.001f)
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
    while(norm > eps)
    {
        for (int i = 0; i < F.size(); i++)
        {
            TempX[i] = (double)F[i];
            for (int j = 0; j < F.size(); j++)
            {
                if (i != j)
                {
                    TempX[i] -= A[i][j] * out[j];
                }
            }
            TempX[i] /= A[i][i];
        }
        norm = fabs(out[0] - TempX[0]);
        for (int i = 0; i < F.size(); i++) {
            if (fabs(out[i] - TempX[i]) > norm) {
                norm = fabs(out[i] - TempX[i]);
            }
            out[i] = TempX[i];
        }
    }
    ftime = omp_get_wtime();
    exec_time = ftime - itime;
    cout << "Seconds of parallel: " << exec_time << endl;
    return out;
}
vector<double> jacobi_omp (vector<vector<double>> A, vector<double> F, double eps = 0.001f)
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
    cout << "Seconds of parallel: " << exec_time << endl;
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
    //vector<vector<double>> matrixA = {{9,0,-12},{0,51,12},{-12,5,24}};
    //vector<vector<double>> matrixA = {{9.2f,2.5f,-3.7f},{0.9f,9.0f,0.2f},{4.5f,-1.6f,-10.3f}};
    //vector<double> matrixB = {-17.5f,4.4f,-22.1f};

    cout << "Insert random matrix (y/n): ";
    std::string yes;
    if (!getline(cin,yes,'\n')) return true;

    vector<vector<double>> matrixA;
    vector<double> matrixB;
    if (yes == "y") {
        int i = 2, maxval;
        cout << "Insert size of square matrix: ";
        cin >> i;
        cout << "Insert max values: ";
        cin >> maxval;
        matrixA = generatematrixA(i, maxval);
        matrixB = generatematrixB(i, maxval);
        /*bool check = normal(matrixA);
        while (!normal(matrixA)) {
            for (int y = 0; y < i; y++) {
                for (int x = 0; x < i; x++) {
                    if (x == i - 1) {
                        cout << matrixA[y][x] << " ";
                        cout << matrixB[y] << " ";
                    } else {
                        cout << matrixA[y][x] << " ";
                    }
                }
                cout << endl;
            }
            cout << "Recheck generate" << endl;

            matrixA = generatematrixA(i, maxval);
            matrixB = generatematrixB(i, maxval);
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        }*/
        std::ofstream matrix("matrix.txt", std::ios_base::trunc);
        cout << "Generated " << i << " size" << endl;
        for(int y = 0 ; y < matrixB.size(); y++) {
            for(int x = 0 ; x < matrixB.size(); x++) {
                if(x == matrixB.size() - 1) {
                    matrix << matrixA[y][x] << " ";
                    matrix << matrixB[y] << " ";
                }
                else {
                    matrix << matrixA[y][x] << " ";
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
        matrixA = vector<vector<double>>(i, vector<double>(i, 0));
        matrixB = vector<double>(i, 0);
        for(int y =0; y < i; y++ ) {
            std::string line;
            getline(cin>>std::ws,line);
            vector<std::string> d = split(line, " ");
            for(int x = 0; x <= i; x++) { // т.к. учитываем, что в матрице есть вектор matrixB
                if(x == i) {
                    matrixB[y] = std::stoi(d[x]);
                }
                else {
                    matrixA[y][x] = (double)std::stoi(d[x]);
                }
            }
            line.clear();
        }
        bool norm = normal(matrixA);
        cout << "Inserted " << matrixB.size() << " size its " << (norm ? "true" : "false") << endl;
        if(!norm) {
            perror("Infinite answer");
            exit(-1);
        }
    }
    double eps;
    cout << "Insert eps: ";
    cin >> eps;
    vector<double> out = jacobi(matrixA, matrixB, eps);
    vector<double> outomp = jacobi_omp(matrixA, matrixB, eps);
    matrixA.clear();
    matrixB.clear();
    std::ofstream answerlinel("answerlinel.txt", std::ios_base::trunc);
    std::ofstream answerparallel("answerparallel.txt", std::ios_base::trunc);
    //std::cout << "ANSWER:\nlinel" << std::endl;
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
