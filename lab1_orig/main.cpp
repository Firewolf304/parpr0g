#include <iostream>
#include <vector>
#include <string>
#include <math.h>
using std::vector;
using std::cout;
using std::endl;
vector<double> jacobi (vector<vector<double>> A, vector<double> F, const double& eps = 0.001)
/*
 * A - переменные слау
 * F - переменные после равно
 */
{
    vector<double> out(F.size(), 0);   // X
    vector<double> TempX(F.size());    // временный x
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
    vector<vector<double>> a = {{8,0,-12},{0,51,12},{-12,12,24}};
    vector<double> b = {58,-41,-88};
    vector<double> out = jacobi(a,b);
    for(auto d : out) {
        cout << d << endl;
    }
    return 0;
}
