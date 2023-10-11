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


void answer (vector<vector<double>> A, vector<double> F, vector<double> *answer_ptr, double w = 1,  double eps = 0.001f)
// w - параметр релаксации
// A - основная матрица
// F - матрица свободных членов
// answer_ptr - вывод
{
    double itime = omp_get_wtime();

    int k = 0;
    double norm;

    do {
        norm=0;
        for (int i = 0; i < F.size(); ++i) {
            double sum = 0;
            for (int j = 0; j < F.size(); ++j) {
                if (i != j) {
                    sum = sum + A[i][j] * (*answer_ptr)[j];
                }
            }
            double Next = (1-w)*(*answer_ptr)[i] + w*(F[i] - sum) / A[i][i];
            if (fabs(Next - (*answer_ptr)[i]) > norm) {
                norm = fabs(Next - (*answer_ptr)[i]);
            }
            (*answer_ptr)[i] = Next;
        }
        k++;
    } while (norm > eps);

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
vector<double> answer_omp (vector<vector<double>> A, vector<double> F, double w = 1,  double eps = 0.001f)
// w - параметр релаксации
// A - основная матрица
// F - матрица свободных членов
// answer_ptr - вывод
{
    double itime = omp_get_wtime();
    vector<double> answer_ptr(F.size(), 0);
    int k = 0;
    double norm;
    do {
        norm=0;
        int i = 0;
        #pragma omp parallel for shared(A,F, answer_ptr, norm) private(i)
        for (i = 0; i < F.size(); i++) {
            double sum = 0;
            for (int j = 0; j < F.size(); ++j) {
                if (i != j) {
                    sum = sum + A[i][j] * answer_ptr[j];
                }
            }
            double Next = (1-w)*answer_ptr[i] + w*(F[i] - sum) / A[i][i];
            if (fabs(Next - answer_ptr[i]) > norm) {
                norm = fabs(Next - answer_ptr[i]);
            }
            answer_ptr[i] = Next;
        }
        k++;
    } while (norm > eps);

    double ftime = omp_get_wtime();
    double exec_time = ftime - itime;
    cout << "Seconds of omp: " << exec_time << endl;
    return answer_ptr;
}
vector<vector<double>> generatematrixA(int n, int maxValue = 100) {
    srand(time(null));
    vector<vector<double>> returned;
    for(int y = 0; y < n; y++) {
        returned.push_back({});
        for(int x = 0; x < n; x++) {
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
            returned[y].push_back((double) (rand() % maxValue));
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
    cout << "Generate matrix? (y/n) : ";
    if (!getline(cin>>std::ws,check,'\n')) {
        return true;
    }
    if(check=="y") {
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
    }
    else {
        matrixA = vector<vector<double>>(i, vector<double>(i, 0));
        matrixB = vector<double>(i, 0);
        for(int y =0; y < i; y++ ) {
            std::string line;
            getline(cin>>std::ws,line);
            vector<std::string> d = splitter(line, " ");
            for(int x = 0; x <= i; x++) {
                if(x == i) {
                    matrixB[y] = std::stoi(d[x]);
                }
                else {
                    matrixA[y][x] = (double)std::stoi(d[x]);
                }
            }
            line.clear();
        }
    }
    double eps = 0.01f;
    cout << "Enter eps: ";
    cin >> eps;
    cin.ignore(INT_MAX, '\n');
    double w = 1;
    cout << "Enter numeric parameter (1-2): ";
    cin >> w;

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
    vector<double> answermat(i, 0);
    answer(matrixA, matrixB, &answermat, w, eps);
    vector<double> answermat_omp = answer_omp(matrixA, matrixB, w, eps);
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
    /*for(int i = 0; i < answermat.size(); i++) {
        cout << answermat[i] << "\t\t" << answermat_omp[i] << endl;
    }*/
    answerlinel.close();
    answerparallel.close();
    return 0;
}
