#include <iostream>
#include <vector>
#include <cmath>

// Функция для проверки, является ли число простым
bool isPrime(int num) {
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

// Функция для генерации списка простых чисел до заданного значения
std::vector<int> generatePrimes(int limit) {
    std::vector<int> primes;
    for (int i = 2; i <= limit; ++i) {
        if (isPrime(i)) {
            primes.push_back(i);
        }
    }
    return primes;
}

// Функция для нахождения максимального числа
int findMaxNumber(int N) {
    // Генерируем простые числа до N
    std::vector<int> primes = generatePrimes(N);

    int maxNumber = 0;

    // Проходим по всем возможным комбинациям степеней 2, 3 и 4
    for (int i = 0; i < primes.size(); ++i) {
        for (int j = 0; j < primes.size(); ++j) {
            for (int k = 0; k < primes.size(); ++k) {
                int currentNumber = pow(2, i) + pow(3, j) + pow(4, k);
                // Проверяем, чтобы число было меньше N
                if (currentNumber < N && currentNumber > maxNumber) {
                    maxNumber = currentNumber;
                }
            }
        }
    }
    return maxNumber;
}

int main() {
    int N = 28; // Замените N на ваше значение
    int result = findMaxNumber(N);

    std::cout << "Максимальное число, меньшее " << N
              << ", представленное суммой степеней 2, 3 и 4 простых чисел: " << result << std::endl;

    return 0;
}