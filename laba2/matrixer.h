//
// Created by firewolf304 on 05.11.23.
//

#ifndef LABA2_MATRIXER_H
#define LABA2_MATRIXER_H
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
using std::cout;
using std::cin;
using std::endl;

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
              T2 находится на обратной стороне YZ при Z=len(Z)-1
     */
__device__ class make_cube1 {
private:
    __device__ int index(int x, int y, int z) {
        x + this->x * (y + this->y * z);
    }
    __device__ void MakeMat( thrust::device_ptr<double> mat ) {
        for(int i = 0; i < this->x; i++) {
            for(int j = 0; j < this->y; j++) {
                for(int k = 0; k < this->z; k++) {
                    if(k == 0) {
                        mat[ index(i,j,k) ] = this->T1;
                        //this->matrix[ index(i,j,k) ] = this->T1;
                        //this->matrix[i][j][k] = this->T1;
                    }
                    else if(k == this->z-1) {
                        mat[ index(i,j,k) ] = this->T2;
                        //this->matrix[i][j][k] = this->T2;
                    }
                }
            }
        }
    }
public:
    thrust::device_vector<double> matrix;
    int x=5,y=5,z=5, radius = 2;
    double T1 = 20.0f;              // температура 1 стороны
    double T2 = 10.0f;              // температура противоположной стороны
    double T_bottom = 0.0f;         // температура нижней грани
    double alpha = 0.001;           // коэф теплопроводности
    void genValue();
    thrust::device_vector<double> operator() ();
    __host__ make_cube() { }
    __host__ make_cube(int x, int y, int z){
        this->x = x;
        this->y = y;
        this->z = z;
    }
    __host__ make_cube(int x, int y, int z, double cylinder_radius) : make_cube(x,y,z) {
        //if(cylinder_radius > (float)(x/2) && cylinder_radius > (float)(y/2) && cylinder_radius > (float)(z/2))  {
        //    throw std::runtime_error("Cylinder is out of cube");
        //}
        this->radius = cylinder_radius;
        matrix = thrust::device_vector<double>(x*y*z, 0);
    }

    __host__ void computeMat() {
        dim3 blocks(1,1,256);
        make_cube::MakeMat<<< blocks, dim3(1,1, (this->y + blocks.z - 1) / blocks.z) >>>( thrust::device_pointer_cast<double>(this->matrix.data()) );
    }
};
#endif //LABA2_MATRIXER_H
