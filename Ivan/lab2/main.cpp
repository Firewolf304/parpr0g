#include <iostream>
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <vector>
#define NS_PER_SECOND 10000000
timespec one, two, delta;

using namespace std;

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
bool istrue(int n) {
    if (n <= 1)
    {
        return false;
    }
    for (int i = 2; i <= sqrt(n); ++i)
    {
        if (n % i == 0)
        {
            return false;
        }
    }
    return true;
}
int algorithm(int N, int rank, int size) {
    vector<int> primes;
    for (int i = 2; i <= N; ++i)
    {
        if (istrue(i))
        {
            primes.push_back(i);
        }
    }
    int coef = primes.size() / size;
    int from = rank * coef;
    int to = (rank == size - 1) ? primes.size() : from + coef;
    while ( N >= 28 )
    {
        for (int i = from; i < to; i++)
        {
            int value = pow(primes[i], 2) + pow(primes[i], 3) + pow(primes[i], 4);
            if (N == value)
            {
                return pow(primes[i], 2) + pow(primes[i], 3) + pow(primes[i], 4);
            }
            if (N < value)
            {
                break;
            }
        }
        N--;
    }
    return -1;
}

int N;
int main(int argc, char *argv[]) {
    int rank;
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(rank == 0)
    {
        cout << "Input N: "; cin >> N;
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    clock_gettime(CLOCK_REALTIME, &one);
    int result_par = algorithm(N, rank, size);
    int globalMax;
    MPI_Reduce(&result_par, &globalMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        clock_gettime(CLOCK_REALTIME, &two);
        sub_timespec(one, two, &delta);
        cout << "Parallel:" << delta.tv_sec << "," << delta.tv_nsec << "sec" << endl;
        if (globalMax != -1)
        {
            cout << "Maximum less of N = " << globalMax << endl;
        }
        else
        {
            cout << "No number" << endl;
        }
        clock_gettime(CLOCK_REALTIME, &one);
        int result_lin = algorithm(N, 0, 1);
        clock_gettime(CLOCK_REALTIME, &two);
        sub_timespec(one, two, &delta);
        cout << "Linel:" << delta.tv_sec << "," << delta.tv_nsec << "sec" << endl;
        if (result_lin != -1)
        {
            cout << "Maximum less of N = " << globalMax << endl;
        }
        else
        {
            cout << "No number" << endl;
        }
    }
    MPI_Finalize();
    return 0;
}
