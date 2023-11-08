//
// Created by firewolf304 on 05.11.23.
//
#include "matrixer.h"

/*__device__ void make_cube::MakeMat( thrust::device_ptr<double> mat ) {
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
}*/

__host__ void make_cube::genValue() {
    srand(time(NULL));
    this->T1 = 100 - (rand() % 50);
    this->T2 = 100 - (rand() % 50);
}
__host__ thrust::device_vector<double> make_cube::operator() () { // output matrix
    return this->matrix;
}
/*__host__ void make_cube::computeMat() {
    dim3 blocks(1,1,256);
    make_cube::MakeMat<<< blocks, dim3(1,1, (this->y + blocks.z - 1) / blocks.z) >>>( thrust::device_pointer_cast<double>(this->matrix.data()) );
}*/
