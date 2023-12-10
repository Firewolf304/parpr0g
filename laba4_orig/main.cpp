#include <iostream>
#include <vector>
#include <string>
#include <cmath>
/*
 Вычислить методом последовательных приближений распределение значений температуры в
 точках верхней границы плоской прямоугольной области размером m * n, имеющей внутри
 круглый вырез. Теплопроводность материала области конечна и не равна нулю. Левая и
 нижняя граница области имеют постоянную температуру 0. Правый верхний угол имеет
 постоянную температуру T, значения температур точек правой границы постоянны и
 равномерно убывают от T до 0. Все остальные точки в начальный момент имеют температуру Т3.
*/
using std::vector;


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
    double radius = 1.0f;
    double alpha = 0.05351f; // коэф распределения (лучше коэф газа Sulfur dioxide)
    void show_matrix(int accuracy = 3) {
        cout.setf(std::ios::fixed); cout.precision(accuracy);
        for(auto y : this->matrix) {
            for(auto x : y) {
                cout << x << " ";
            }
            cout << endl;
        }
    }
    void show_matrix(vector<vector<double>> mat, int accuracy = 3) {
        cout.setf(std::ios::fixed); cout.precision(accuracy);
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
        this->matrix = vector<vector<double>>(this->m, vector<double>(this->n, 0));
        matrix[0][m-1] = this->T;

        for(int y = 0 ; y < this->n; y++ ) {
            for(int x = 0; x < this->m; x++) {
                if(this->radius < (float) std::sqrt(std::pow(x - (float) (this->n / 2), 2) +
                                                    std::pow(y - (float) (this->m / 2),2))) {
                    if (x == 0 && y >= 0 || x >= 0 && y == this->n - 1) { // left/down place
                        this->matrix[y][x] = 0;
                    } else if (x == this->m - 1 && y > 0) { // right place
                        this->matrix[y][x] = (double) (this->matrix[y - 1][x] - step);
                    } else { // another
                        if (y == 0 && x == this->m - 1) continue;
                        //this->matrix[y][x] = (double)(rand() % this->t3) * ((rand() % 1000000) > 500000 ? -1 : 1);
                        this->matrix[y][x] = this->t3;
                    }
                }
            }
        }
    }
    vector<vector<double>> operator()() {
        return matrix;
    }
    void iteration(double eps = 0.0001f, int maxIteration=10000) {
        double norm = 1;
        int count = 0;
        while (norm > eps && count < maxIteration) {
            for(int i = 0; i < this->n; i++) {
                for (int j = 0; j < this->m; j++) {

                        double oldTemp = this->matrix[i][j];
                        double newTemp = 0.0f;
                        for (int y1 = ((j > 0) ? j - 1 : j); y1 <= j + 1 && y1 < this->m; y1++) {
                            for (int x1 = ((i > 0) ? i - 1 : i); x1 <= i + 1 && x1 < this->n; x1++) {
                                if(this->radius < (double) std::sqrt(std::pow(i - (double) (this->n / 2), 2) +
                                                                     std::pow(j - (float) (this->m / 2),2))) {
                                    newTemp += this->matrix[y1][x1];
                                }
                            }
                        }
                        this->matrix[i][j] = newTemp * this->alpha;
                        norm = std::abs(newTemp - oldTemp);

                }
            }
            count++;
        }
    }
};

int main() {
    std::cout << "Hello, World!" << std::endl;
    make_matrix matrix;
    matrix.n = 9; matrix.m = 9;
    matrix.T = 100;
    matrix.t3 = 100;
    matrix.radius = 1.5f;
    matrix.alpha = 0.05351f;

    matrix.make();
    matrix.show_matrix();

    matrix.iteration(0.01f, 100);
    matrix.show_matrix();
    return 0;
}
