#include <iostream>
#include <omp.h>
#include <vector>
#include <math.h>

using std::vector;
using std::cout;
using std::cin;
using std::endl;
class make_matrix {
/*matrix example:
 n
 --------------
 0 x x x T    |
 0 x x x (T-1)|
 0 x x x (T-2)|
 0 x x x (T-3)|
 0 0 0 0 0    |  m*/
public:
    void clear_matrix() {
        this->matrix.clear();
        this->matrix.push_back({});
        for(int y = 0 ; y < this->n; y++ ) {
            for(int x = 0; x < this->m; x++) {
                this->matrix[y].push_back(0);
            }
            this->matrix.push_back({});
        }
    }
    vector<vector<double>> matrix;
    int n = 3;      // like y
    int m = 3;      // like x
    double T = 20;
    int t3 = 3;     // ???
    void show_matrix() {
        for(auto y : this->matrix) {
            for(auto x : y) {
                cout << x << " ";
            }
            cout << endl;
        }
    }
    void show_matrix(vector<vector<double>> mat) {
        for(auto y : mat) {
            for(auto x : y) {
                cout << x << " ";
            }
            cout << endl;
        }
    }
    void show_matrix(vector<double> mat) {
        for(auto x : mat) {
            cout << x << " ";
        }
        cout << endl;
    }
    vector<double> get_top_matrix(vector<vector<double>> mat) {
        vector<double> out;
        for(auto d : mat[0]) {
            out.push_back(d);
        }
        return out;
    }
    vector<vector<double>> get_another_matrix(vector<vector<double>> mat) {
        vector<vector<double>> out;
        out.push_back({});
        for(int y = 1; y < this->n; y++) {
            for(int x = 0; x < this-> m; x++) {
                out[y].push_back(mat[y][x]);
            }
            out.push_back({});
        }
        return out;
    }
    void make() {
        float step = (float)(T / (n-1));
        clear_matrix();
        matrix[0][m-1] = this->T;

        for(int y = 0 ; y < this->n; y++ ) {
            for(int x = 0; x < this->m; x++) {
                if(x==0 && y >= 0 || x>=0 && y == this->n-1 ) { // left/down place
                    this->matrix[y][x] = 0;
                }
                else if(x == this->m-1 && y > 0) { // right place
                    this->matrix[y][x] = (double)(this->matrix[y-1][x] - step);
                }

                else { // another
                    if(y== 0 && x == this->m-1) continue;
                    this->matrix[y][x] = (double)(rand() % this->t3) * ((rand() % 1000000) > 500000 ? -1 : 1);
                }
            }
        }
    }
    vector<vector<double>> operator()() {
        return matrix;
    }
};


vector<double> jacobi (vector<vector<double>> A, vector<double> F, const double& eps = 0.001)
/*
* A - переменные слау
* F - переменные после равно
*/
{
    vector<double> out(F.size(), 0); // X
    vector<double> TempX(F.size()); // временный x
    double norm = 1;
    while(norm > eps) {
        for (int i = 0; i < F.size(); i++) {
            TempX[i] = F[i];
            for (int j = 0; j < F.size(); j++) {
                if (i != j) {
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
    return out;
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    make_matrix matrix;
    matrix.n = 3; matrix.m = 4;
    matrix.T = 10;
    matrix.t3 = 10;
    matrix.make();
    matrix.show_matrix();
    cout << "------" << endl;
    vector<double> top = matrix.get_top_matrix(matrix());
    matrix.show_matrix(top);
    cout << "------" << endl;
    int c = 0;
    #pragma omp parallel sections num_threads(100)
    {
        #pragma omp section
        {
            c++;
            cout << omp_get_thread_num() << endl;
        }
        #pragma omp section
        {
            c++;
            cout << omp_get_thread_num() << endl;
        }
    }
    cout << c << endl;
    //matrix.show_matrix();
    /*int c = 0;
    #pragma omp parallel num_threads(3)
        {
            cout << omp_get_thread_num() << endl;
            for (int i = 0; i < 100; i++)
                c++;
        }

    cout << c << endl;
    */

    return 0;
}

