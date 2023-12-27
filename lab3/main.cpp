#include <iostream>
#include <mpi.h>
#include <cmath>
#include <vector>
#define NS_PER_SECOND 1000000000000
void sub_timespec(struct timespec t1, struct timespec t2, struct timespec *td)
{
    td->tv_nsec = t2.tv_nsec - t1.tv_nsec;
    td->tv_sec  = t2.tv_sec - t1.tv_sec;
    if (td->tv_sec > 0 && td->tv_nsec < 0)
    {
        td->tv_nsec += NS_PER_SECOND;
        td->tv_sec--;
    }
    else if (td->tv_sec < 0 && td->tv_nsec > 0)
    {
        td->tv_nsec -= NS_PER_SECOND;
        td->tv_sec++;
    }
}
class par_iter {
public:
    int ranker, sizer;
    par_iter(int rank, int size)  { this->ranker= rank; this->sizer =size; };
    bool CheckBasis(int n) {
        if (n <= 1)
            return false;
        for (int i = 2; i <= sqrt(n); ++i) {
            if (n % i == 0) {
                return false;
            }
        }
        return true;
    }
    int run(int N) {


        return findMaxNumbers(N);
    }
    int findMaxNumbers(int N) {
        std::vector<int> primes;
        for (int i = 2; i <= N; ++i) {
            if (CheckBasis(i)) {
                primes.push_back(i);
            }
        }
        int primesPerProcess = primes.size() / this->sizer;
        int start = this->ranker * primesPerProcess;
        int end = (this->ranker == this->sizer - 1) ? primes.size() : start + primesPerProcess;
        for (; N >= 28; N--) {
            for (int i = start; i < end; i++) {
                int value = pow(primes[i], 2) + pow(primes[i], 3) + pow(primes[i], 4);
                if (N == value) {
                    return pow(primes[i], 2) + pow(primes[i], 3) + pow(primes[i], 4);
                }
                if (N < value)
                    break;
            }
        }


        return -1;
    }
};

class lin_iter {
public:
    int ranker, sizer;
    bool CheckBasis(int num) {
        if (num < 2) {
            return false;
        }
        for (int i = 2; i <= sqrt(num); ++i) {
            if (num % i == 0) {
                return false;
            }
        }
        return true;
    }


    int findMaxNumber(int N) {
        std::vector<int> primes;
        for (int i = 2; i <= N; ++i) {
            if (CheckBasis(i)) {
                primes.push_back(i);
            }
        }

        int maxNumber = 0;
        for (; N >= 28; N--) {
            for (int i = 0; i < primes.size(); ++i) {
                int value = pow(primes[i], 2) + pow(primes[i], 3) + pow(primes[i], 4);
                if (N == value) {
                    return pow(primes[i], 2) + pow(primes[i], 3) + pow(primes[i], 4);
                }
                if (N < value)
                    break;
            }
        }
        return maxNumber;
    }
};

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    par_iter iter_par(rank, size);
    lin_iter iter_lin;



    int N;
    if(rank == 0) {
        std::cout << "Введите число N: ";
        std::cin >> N;
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    timespec start, finish, delta;
    clock_gettime(CLOCK_REALTIME, &start);
    int result_par = iter_par.findMaxNumbers(N);


    int globalMax;
    MPI_Reduce(&result_par, &globalMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        clock_gettime(CLOCK_REALTIME, &finish);
        sub_timespec(start, finish, &delta);
        std::cout << "Выполнена параллель за " << delta.tv_sec << "," << delta.tv_nsec << " секунд" << std::endl;
        if (globalMax != -1) {
            std::cout << "Максимальное число, меньшее " << N << ": " << globalMax << std::endl;
        } else {
            std::cout << "Нет подходящего числа" << std::endl;
        }

        clock_gettime(CLOCK_REALTIME, &start);
        int result_lin = iter_lin.findMaxNumber(N);
        clock_gettime(CLOCK_REALTIME, &finish);
        sub_timespec(start, finish, &delta);
        std::cout << "Выполнена линель за " << delta.tv_sec << "," << delta.tv_nsec << " секунд" << std::endl;
        if (result_lin != -1) {
            std::cout << "Максимальное число, меньшее " << N << ": " << globalMax << std::endl;
        } else {
            std::cout << "Нет подходящего числа" << std::endl;
        }


    }




    MPI_Finalize();

    return 0;
}