/*
* Napisz program w języku „C/C++”, realizujący metody:
* (a) Picarda
* (b) bisekcji
* (c) Newtona
* (d) siecznych
* rozwiązywania pojedynczych algebraicznych równań nieliniowych. Zastosuj program do
* przykładów z zadania 1. Zastosuj trzy niezależne kryteria zakończenia iteracji. Zadbaj o to, aby
* wyprowadzać na konsolę wyniki pośrednie obliczeń dla każdej iteracji, tak aby możliwe było
* obserwowanie zbieżności kolejnych przybliżeń pierwiastków i porównanie liczby iteracji
* niezbędnych do uzyskania rozwiązania o zadanej dokładności przez każdą z metod. W szczególności
* oblicz jak zmienia się estymator błędu rozwiązania oraz residuum równania w trakcie iteracji.
*/

#include <iostream>
#include <cmath>
#include <functional>

#define ITERATIONS 100
#define TOLX 1e-12
#define TOLF 1e-12

double function_tanh(double);

double function_sinh(double);

double phi_tanh(double);

double phi_sinh(double);

void picard(const std::function<double(double)> &, const std::function<double(double)> &);

void bisection(double, double, const std::function<double(double)> &);

void newton(double, const std::function<double(double)> &, const std::function<double(double)> &);

void secant(double, double, const std::function<double(double)> &);

int main() {
    picard(function_tanh, phi_tanh);
    picard(function_sinh, phi_sinh);
    return 0;
}

double function_tanh(double x) {
    return std::tanh(x) + 2.0 * (x - 1.0);
}

double function_sinh(double x) {
    return std::sinh(x) + x / 4.0 - 1.0;
}

double phi_tanh(double x) {
    return -1.0 / (2.0 * std::cosh(x) * std::cosh(x));
}

double phi_sinh(double x) {
    return -4.0 * std::cosh(x);
}

void picard(const std::function<double(double)> &f, const std::function<double(double)> &phi) {
    std::cout << "----------------------- Picard -----------------------" << std::endl;
    std::printf("%-5s%-17s%-17s%-17s\n", "i", "Wartosc", "Estymator", "Residuum");

    double x_0 = 1.0;
    for (int i = 1; i <= ITERATIONS; i++) {
        const double x_1 = phi(x_0);
        const double estimator = std::fabs(x_1 - x_0);
        const double residuum = std::fabs(f(x_1));

        std::printf("%-3d %-10.12f   %-10.12f   %-10.12f\n", i, x_1, estimator, residuum);

        if (estimator < TOLX || residuum < TOLF || std::fabs(x_1) > 1.0) break;

        x_0 = x_1;
    }
    std::cout << std::endl;
}

void bisection(double a, double b, const std::function<double(double)> &f) {
}
