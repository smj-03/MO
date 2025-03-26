/*
* Napisz program w języku „C/C++”, umożliwiający „doświadczalne” wyznaczenie liczby bitów
* mantysy oraz tzw. epsylona maszynowego, tj. najmniejszej liczby E takiej, że fl(E + 1) > 1. Wyznacz
* epsylon maszynowy dla zmiennych typu float i double i sprawdź czy da się go wyznaczyć dla
* zmiennych typu long double. Sprawdź też ile dokładnych cyfr znaczących posiada epsylon maszynowy.
* Aby znaleźć odpowiedź na pytanie jak napisać taki program, zacznij od wyjaśnienia kwestii jaki jest
* związek E z precyzją arytmetyki.
*/

// Notatki
// Epsylon maszynowy - Rożnica pomiędzy 1, a następną dostępną liczbą.
// E = x+ - x-, Zakładając, że x- = 1
// E = 2 * v = 2 * 2^[-1 * (t+1)] = 2^(-t)

#include <iostream>
#include <iomanip>

int main() {
    // FLOAT
    int t_float = 0;
    float E_float = 1.0f;
    float temp_float;

    while (1.0f + E_float > 1.0f) {
        t_float++;
        temp_float = E_float;
        E_float /= 2.0f;
    }

    std::cout << "Liczba bitow mantysy dla float: " << t_float - 1 << std::endl;
    std::cout << "Epsylon maszynowy dla float: " << std::fixed << std::setprecision(t_float) << temp_float << std::endl
            << std::endl;

    // DOUBLE
    int t_double = 0;
    double E_double = 1.0;
    double temp_double;

    while (1.0 + E_double > 1.0) {
        t_double++;
        temp_double = E_double;
        E_double /= 2.0;
    }

    std::cout << "Liczba bitow mantysy dla double: " << t_double - 1 << std::endl;
    std::cout << "Epsylon maszynowy dla double: " << std::fixed << std::setprecision(t_double) << temp_double <<
            std::endl << std::endl;

    // LONG DOUBLE
    int t_long_double = 0;
    long double E_long_double = 1.0L;
    long double temp_long_double;

    while (1.0L + E_long_double > 1.0L) {
        t_long_double++;
        temp_long_double = E_long_double;
        E_long_double /= 2.0L;
    }

    std::cout << "Liczba bitow mantysy dla long double: " << t_long_double - 1 << std::endl;
    std::cout << "Epsylon maszynowy dla long double: " << std::fixed << std::setprecision(t_long_double) <<
            temp_long_double << std::endl;

    return 0;
}
