#include <iostream>
#include <vector>
#include <cmath>

using std::cout;
using std::cin;
using std::endl;
using std::vector;
/*     Нынешнее представление:
                    +--------+                   представление координат
                   /|       /|                    +------------------+
                  / |      / |                    |  y               |
Координата y -> y+--------+  |                    |  ^    z          |
                 |  |     |  |                    |  |   ^           |
  Координата z ->| z+-----|--+                    |  |  /            |
                 | /      | /                     |  | /             |
                 |/       |/                      |  |/              |
                 0--------+x    <- координата x   |  0----------> x  |
              (0;0;0)                             +------------------+
              радиус цилиндра считается по стороне x z от координат центра
              T1 находится на стороне YZ при Z=0
              T2 находится на обратной стороне YZ при Z=len(Z)-1*/
class make_cube {
    typedef vector<vector<vector<double>>> matrix_type;
public:

    make_cube(int x, int y, int z){
        this->x = x;
        this->y = y;
        this->z = z;
    }
    make_cube(int x, int y, int z, double cylinder_radius) : make_cube(x,y,z) {
        /*if(cylinder_radius > (float)(x/2) && cylinder_radius > (float)(y/2) && cylinder_radius > (float)(z/2))  {
            throw std::runtime_error("Cylinder is out of cube");
        }*/
        this->radius = cylinder_radius;
        matrix = matrix_type(x, vector<vector<double>>(y,vector<double>(z, 0)));
    }

    int x=2,y=2,z=2, radius = 2;
    double T1 = 20.0f;              // температура 1 стороны
    double T2 = 10.0f;              // температура противоположной стороны
    double T_bottom = 0.0f;         // температура нижней грани
    double alpha = 0.05351f;           // коэф теплопроводности

    /*-----matrix-----*/
    matrix_type matrix;
    void computeMat() {
        for(int i = 0; i < this->x; i++) {
            for(int j = 0; j < this->y; j++) {
                for(int k = 0; k < this->z; k++) {
                    if(k == 0) {
                        this->matrix[i][j][k] = this->T1;
                    }
                    else if(k == this->z-1) {
                        this->matrix[i][j][k] = this->T2;
                    }
                    /*else {
                        if(std::sqrt( std::pow( i - (this->x / 2) , 2 ) + std::pow( j - (this->y / 2) , 2 ) + std::pow( k - (this->z / 2) , 2 ) ) )
                    }*/
                }
            }
        }
    }
    void genValue() {
        srand(time(NULL));
        this->T1 = 100 - (rand() % 50);
        this->T2 = 100 - (rand() % 50);
    }
    void iteration (double eps = 0.0001f, int maxIteration=10000) { // using basic matrix iterator

        double norm = 1;
        int count = 0;
        while (norm > eps && count < maxIteration) {
            // перебор каждой точки
            for(int i = 0; i < this->x; i++) {
                for (int j = 0; j < this->y; j++) {
                    for (int k = 0; k < this->z; k++) {

                        // просмотр соседей и сумма новой температуры
                        double oldTemp = this->matrix[i][j][k];
                        double newTemp = 0.0f;
                        cout << "newtemp = ";
                        for(int x1 = ((i>0) ? i-1 : i) ; x1 <= i+1 && x1 < this->x; x1++) {
                            for(int y1 = ((j>0) ? j-1 : j); y1 <= j+1 && y1 < this->y; y1++) {
                                for(int z1 = ((k>0) ? k-1 : k); z1 <= k+1 && z1 < this->z; z1++) {

                                    if (this->radius < (float) std::sqrt(std::pow(x1 - (float) (this->x / 2), 2) +
                                                                         std::pow(z1 - (float) (this->z / 2),
                                                                                  2)) /*&& std::abs(y1 - (this->y / 2)) <= (this->y / 2)*/) {
                                        cout << " + " << this->matrix[x1][y1][z1];
                                        newTemp += this->matrix[x1][y1][z1];
                                    }

                                }
                            }
                        }
                        //newTemp *= this->alpha; // формула распределения: dT/dt = alpha * ( (d^2T/dx^2) + (d^2T/dy^2) + (d^2T/dz^2) )
                                                  // т.е. нужно взять от зависимой точки сумму всех соседей и умножить на коэф распределения

                        cout << endl;
                        this->matrix[i][j][k] = newTemp * this->alpha;
                        cout << "newtemp=" << newTemp << " matrix=" << this->matrix[i][j][k] << " (" << i << "," << j << "," << k << ")" << endl;
                        norm = std::abs( newTemp - oldTemp );
                    }
                }
            }
            count++;
        }

    }
    void show() {
        for (int i = 0; i < x; i++) {
            for (int k = 0; k < z; k++) {
                cout << this->matrix [i][y-1][k] << " ";
            }
            cout << endl;
        }
        cout << "-------------" << endl;
        for (int j = 0; j < y; j++) {
            for (int i = 0; i < x; i++) {
                cout << this->matrix[i][j][(int)(z/2)] << " ";
            }
            cout << endl;
        }
        cout << "-------------" << endl;
        for (int j = 0; j < y; j++) {
            for (int k = 0; k < z; k++) {
                cout << this->matrix[(int)(x/2)][j][k] << " ";
            }
            cout << endl;
        }
        cout << "-------------" << endl;
        for (int j = 0; j < y; j++) {
            for (int i = 0; i < x; i++) {
                cout << this->matrix[i][j][(int)(z-1)] << " ";
            }
            cout << endl;
        }
        cout << "-------------" << endl;
        for (int j = 0; j < y; j++) {
            for (int i = 0; i < x; i++) {
                cout << this->matrix[i][j][(int)(0)] << " ";
            }
            cout << endl;
        }
    }

    void visualize() {
        /*
         Пример отрисовки:
         Матрицы 3 на 3:
         1 0 1
         0 0 0
         1 0 1 все подобны

         отрисовка результат:

             |   +---------+
             |  / 1  0  1 /|
             | / 0  0  0 /1|
             |/ 1  0  1 /00|
           y +---------+101+
             | 1  0  1 |00/
             | 0  0  0 |1/
             | 1  0  1 |/
    ---------+---------+x----------
            /|
           / |
          z  |




         метод итерации по оси xz :
         for (int i = 0, x_pos = z-1; i < x; i++, x_pos--) { // z-1 т.к. мы рассматриваем высоту строки по X в зависимости от координаты Z, т.е. кол-во строк равно Z, как и их высота
            for (int k = 0; k < z; k++) {
                if(x_pos>=0 || i == 0 || i = x-1) {      // добавить пробелы, пока есть строки по xz или это начало строки
                    cout << " ";
                }
                else{
                    cout << answer[i][y-1][k] << " "; // y-1 т.к. смотрим только вверх
                }
            }
            cout << endl;
        }
         */

    }
    matrix_type operator() () { // output matrix
        return this->matrix;
    }
private:

};

int main() {
    cout << cout.setf(std::ios::fixed) << cout.precision(3) << endl;
    std::cout << "Hello, World!" << std::endl;
    int x = 10, y = 4, z = 10;

    make_cube cube(x,y,z, 1.0f);
    cube.alpha= 0.05351f; //0.05351f;
    cube.computeMat();
    cube.iteration(0.00001f, 1);
    vector<vector<vector<double>>> answer = cube();

    return 0;
}
