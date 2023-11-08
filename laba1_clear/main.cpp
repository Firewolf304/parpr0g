#include <omp.h>
#include <string>
#include <math.h>
#include <chrono>
#include <thread>
#include <fstream>
#include <iostream>
using std::vector;
using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::rand;
#define null NULL


void answer (vector<vector<double>> A, vector<double> F, vector<double> *answer_ptr, double eps = 0.001f)
// w - параметр релаксации
// A - основная матрица
// F - матрица свободных членов
// answer_ptr - вывод
{
    double itime = omp_get_wtime();
    int size = F.size();
    //vector<double> Xk(size);
    vector<double> Zk(size);
    vector<double> Rk(size);
    vector<double> Sz(size);
    double alpha, beta, mf;
    double Spr, Spr1, Spz;
    int max = 100000;
    int i,j,kl=1;
/* Вычисляем сумму квадратов элементов вектора F*/
    for(mf = 0,i = 0; i < size; i++){
        mf += F[i] * F[i];
    }


/* Задаем начальное приближение корней. В Хk хранятся значения корней
 * к-й итерации. */
    for(i = 0; i < size; i++){
        (*answer_ptr)[i] = 0.2;
    }

/* Задаем начальное значение r0 и z0. */
    for(i = 0; i < size; i++){
        for(Sz[i]=0,j = 0; j < size; j++)
            Sz[i] += A[i][j] * (*answer_ptr)[j];
        Rk[i] = F[i] - Sz[i];
        Zk[i] = Rk[i];
    }
    int Iteration = 0;
    do{
        Iteration++;
        /* Вычисляем числитель и знаменатель для коэффициента
         * alpha = (rk-1,rk-1)/(Azk-1,zk-1) */
        Spz = 0;
        Spr = 0;
        for(i = 0; i < size; i++){
            for(Sz[i] = 0, j = 0; j < size; j++){
                Sz[i] += A[i][j] * Zk[j];
            }
            Spz += Sz[i] * Zk[i];
            Spr += Rk[i] * Rk[i];
        }
        alpha = Spr / Spz;             /*  alpha    */


/* Вычисляем вектор решения: xk = xk-1+ alpha * zk-1,
    вектор невязки: rk = rk-1 - alpha * A * zk-1 и числитель для betaa равный (rk,rk) */
        Spr1 = 0;
        for(i = 0; i < size; i++){
            (*answer_ptr)[i] += alpha * Zk[i];
            Rk[i] -= alpha * Sz[i];
            Spr1 += Rk[i] * Rk[i];
            //cout << "Iter #" << kl;
            //cout << " " << "X[" << i << "] = " << Xk[i] << endl;
        }
        //cout << endl;
        kl++;

/* Вычисляем  beta  */
        beta = Spr1 / Spr;

/* Вычисляем вектор спуска: zk = rk+ beta * zk-1 */
        for(i = 0; i < size; i++)
            Zk[i] = Rk[i] + beta * Zk[i];
    }
/* Проверяем условие выхода из итерационного цикла  */
    while(Spr1 / mf > eps * eps && Iteration <= max);

    double ftime = omp_get_wtime();
    double exec_time = ftime - itime;
    cout << "Seconds of parallel: " << exec_time << endl;
}
std::vector<std::string> splitter(std::string s, std::string simb) {
    size_t pos_start = 0, pos_end, delim_len = simb.length();
    std::string token;
    std::vector<std::string> res;
    while ((pos_end = s.find(simb, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }
    res.push_back (s.substr (pos_start));
    return res;
}
vector<double> answer_omp (vector<vector<double>> A, vector<double> F,  double eps = 0.001f)
// w - параметр релаксации
// A - основная матрица
// F - матрица свободных членов
// answer_ptr - вывод
{
    double itime = omp_get_wtime();
    vector<double> answer_ptr(F.size(), 0);
    int size = F.size();
    //vector<double> Xk(size);
    vector<double> Zk(size);
    vector<double> Rk(size);
    vector<double> Sz(size);
    double alpha, beta, mf;
    double Spr, Spr1, Spz;
    int max = 100000;
    int i,j,kl=1;
/* Вычисляем сумму квадратов элементов вектора F*/
    mf = 0;
#pragma omp parallel for private(i) reduction(+ : mf)
    for(i = 0; i < size; i++){
        mf = F[i] * F[i];
    }
/* Задаем начальное приближение корней. В Хk хранятся значения корней
 * к-й итерации. */
#pragma omp parallel for shared(answer_ptr) private(i)
    for(i = 0; i < size; i++){
        answer_ptr[i] = 0.2;
    }

/* Задаем начальное значение r0 и z0. */
#pragma omp parallel for shared(Sz, Rk, Zk) private(i)
    for(i = 0; i < size; i++){
        double check = 0;
#pragma omp parallel for private(j) reduction(+ : check)
        for(j = 0; j < size; j++)
        {
            check = A[i][j] * answer_ptr[j];
        }
        Sz[i] = check;
        Rk[i] = F[i] - Sz[i];
        Zk[i] = Rk[i];
    }
    int Iteration = 0;
    do{
        Iteration++;
        /* Вычисляем числитель и знаменатель для коэффициента
         * alpha = (rk-1,rk-1)/(Azk-1,zk-1) */
        Spz = 0;
        Spr = 0;
#pragma omp parallel for shared(Sz, Spz) private(i)
        for(i = 0; i < size; i++){
            double check = 0;
#pragma omp parallel for private(j) reduction(+ : check)
            for( j = 0; j < size; j++){
                check += A[i][j] * Zk[j];
            }
            Sz[i] = check;
            Spz += Sz[i] * Zk[i];
            Spr += Rk[i] * Rk[i];
        }
        alpha = Spr / Spz;
        Spr1 = 0;
        #pragma omp parallel for shared(answer_ptr, Rk, Spr1) private(i)
        for(i = 0; i < size; i++){
            answer_ptr[i] += alpha * Zk[i];
            Rk[i] -= alpha * Sz[i];
            Spr1 += Rk[i] * Rk[i];
        }
        kl++;
        beta = Spr1 / Spr;
        #pragma omp parallel for shared(Zk) private(i)
        for(i = 0; i < size; i++)
            Zk[i] = Rk[i] + beta * Zk[i];
    }
    while(Spr1 / mf > eps * eps && Iteration <= max);

    double ftime = omp_get_wtime();
    double exec_time = ftime - itime;
    cout << "Seconds of omp: " << exec_time << endl;
    return answer_ptr;
}
vector<vector<double>> generatematrixA(int n, int maxValue = 100) {
    srand(time(null));
    vector<vector<double>> returned(n, vector<double>(n, 0));
    for(int y = 0; y < n; y++) {
        for(int x = y; x < n; x++) {
            //std::this_thread::sleep_for(std::chrono::milliseconds(5));
            returned[y][x] =  (double)(rand() % maxValue);
        }
    }
    for(int y = 0; y < n; y++) {
        for(int x = 0; x < y; x++) {
            returned[y][x] = returned[x][y];
        }
    }
    // проверка неправильной матрицы
    for(int y = 0; y < n; y++) {
        double summ = 0;
        for(int x = 0; x < n; x++) {
            if(x != y) {
                summ += returned[y][x];
            }
        }
        if(returned[y][y] < summ ) {
            srand(time(NULL));
            returned[y][y] = summ + ((double)(rand() % maxValue) / 10) + 0.1f; // 1 чтобы избавиться от 0, все равно рандом
        }
    }
    return returned;
}
vector<double> generatematrixB(int n, int maxValue = 100) {
    vector<double> returned;
    for(int x = 0; x < n; x++) {
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        returned.push_back((double) (rand() % maxValue));
    }
    return returned;
}
bool checkMatrix(vector<vector<double>> returned) {
    for(int y = 0; y < returned.size(); y++) {
        double summ = 0;
        for(int x = 0; x < returned.size(); x++) {
            if(x != y) {
                summ += returned[y][x];
            }
        }
        if(returned[y][y] < summ ) {
            return false;
        }
    }
    return true;
}
int main() {

    vector<vector<double>> matrixA;
    vector<double> matrixB;
    int i = 2;
    cout << "Enter size: ";
    cin >> i;
    string check;
    int maxValue;
    cout << "Enter max value: ";
    cin >> maxValue;
    matrixA = generatematrixA(i, maxValue);
    matrixB = generatematrixB(i, maxValue);
    std::ofstream matrix("matrix_AB.txt", std::ios_base::trunc);
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
    double eps = 0.01f;
    cout << "Enter eps: ";
    cin >> eps;
    cin.ignore(INT_MAX, '\n');
    vector<double> answermat(i, 0);
    answer(matrixA, matrixB, &answermat, eps);
    vector<double> answermat_omp = answer_omp(matrixA, matrixB, eps);
    matrixA.clear();
    matrixB.clear();
    std::ofstream answerlinel("alinel.txt", std::ios_base::trunc);
    std::ofstream answerparallel("aparallel.txt", std::ios_base::trunc);
    int count=answermat_omp.size();
    for(int i = 0; i < answermat_omp.size(); i++) {
        answerlinel << answermat_omp[i] << "\n";

        answerparallel << answermat_omp[i] << "\n";
        if(answermat_omp[i] == answermat_omp[i] ) {
            count--;
        }
    }

    answerlinel.close();
    answerparallel.close();
    return 0;
}
