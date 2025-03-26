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

#define N_MAX 100
#define TOLX 1e-12
#define TOLF 1e-12

double function_tanh(double);

double function_sinh(double);

double derivative_tanh(double);

double derivative_sinh(double);

double phi_tanh(double);

double phi_sinh(double);

void picard(const std::function<double(double)> &, const std::function<double(double)> &);

void bisection(double, double, const std::function<double(double)> &);

void newton(double, const std::function<double(double)> &, const std::function<double(double)> &);

void secant(double, double, const std::function<double(double)> &);

int main() {
    std::cout << "FUNKCJA: f(x) = tanh(x) + 2(x - 1)" << std::endl;
    picard(function_tanh, phi_tanh);
    bisection(0.0, 2.0, function_tanh);
    newton(5.0, function_tanh, derivative_tanh);
    secant(-1.0, 5.0, function_tanh);

    std::cout << "FUNKCJA: f(x) = sinh(x) + x/4 - 1" << std::endl;
    bisection(0.0, 2.0, function_sinh);
    newton(5.0, function_sinh, derivative_sinh);
    secant(-1.0, 5.0, function_sinh);
    return 0;
}

double function_tanh(double x) {
    return std::tanh(x) + 2.0 * (x - 1.0);
}

double function_sinh(double x) {
    return std::sinh(x) + x / 4.0 - 1.0;
}

double derivative_tanh(double x) {
    return 1.0 / (std::cosh(x) * std::cosh(x)) + 2.0;
}

double derivative_sinh(double x) {
    return std::cosh(x) + 0.25;
}

double phi_tanh(double x) {
    return 0.5 * (2.0 - std::tanh(x));
}

double phi_sinh(double x) {
    return 4.0 * (1.0 - std::sinh(x));
}

void picard(const std::function<double(double)> &f, const std::function<double(double)> &phi) {
    std::cout << "---------------------- Picard ----------------------" << std::endl;
    std::printf("%-5s%-17s%-17s%-17s\n", "i", "Wartosc", "Estymator", "Residuum");

    double x_0 = 1.0;
    for (int i = 1; i <= N_MAX; i++) {
        const double x_1 = phi(x_0);
        const double estimator = std::fabs(x_1 - x_0);
        const double residuum = std::fabs(f(x_1));

        std::printf("%-3d %-10.12f   %-10.12f   %-10.12f\n", i, x_1, estimator, residuum);

        if (estimator < TOLX && residuum < TOLF) break;

        x_0 = x_1;
    }
    std::cout << std::endl;
}

void bisection(double a, double b, const std::function<double(double)> &f) {
    std::cout << "--------------------- Bisekcja ---------------------" << std::endl;
    std::printf("%-5s%-17s%-17s%-17s\n", "i", "Wartosc", "Estymator", "Residuum");
    double x_n;
    for (int i = 1; i <= N_MAX; i++) {
        x_n = (a + b) / 2.0;
        const double estimator = std::fabs(b - a) / 2.0;
        const double residuum = std::fabs(f(x_n));

        std::printf("%-3d %-10.12f   %-10.12f   %-10.12f\n", i, x_n, estimator, residuum);

        if (estimator < TOLX && residuum < TOLF) break;

        if ((f(a) < 0 && f(x_n) > 0) || (f(a) > 0 && f(x_n) < 0))
            b = x_n;
        else
            a = x_n;
    }
    std::cout << std::endl;
}

void newton(double x_1, const std::function<double(double)> &f, const std::function<double(double)> &fd) {
    std::cout << "---------------------- Newton ----------------------" << std::endl;
    std::printf("%-5s%-17s%-17s%-17s\n", "i", "Wartosc", "Estymator", "Residuum");

    for (int i = 1; i <= N_MAX; i++) {
        const double x_n = x_1 - f(x_1) / fd(x_1);
        const double estimator = std::fabs(x_n - x_1);
        const double residuum = std::fabs(f(x_n));

        std::printf("%-3d %-10.12f   %-10.12f   %-10.12f\n", i, x_n, estimator, residuum);

        if (estimator < TOLX && residuum < TOLF) break;

        x_1 = x_n;
    }
    std::cout << std::endl;
}

void secant(double x_1, double x_2, const std::function<double(double)> &f) {
    std::cout << "---------------------- Sieczne ---------------------" << std::endl;
    std::printf("%-5s%-17s%-17s%-17s\n", "i", "Wartosc", "Estymator", "Residuum");

    for (int i = 1; i <= N_MAX; i++) {
        const double x_n = x_2 - f(x_2) * (x_2 - x_1) / (f(x_2) - f(x_1));
        const double estimator = std::fabs(x_n - x_2);
        const double residuum = std::fabs(f(x_n));

        std::printf("%-3d %-10.12f   %-10.12f   %-10.12f\n", i, x_n, estimator, residuum);

        if (estimator < TOLX && residuum < TOLF) break;

        x_1 = x_2;
        x_2 = x_n;
    }
    std::cout << std::endl;
}
