/*
* Zaimplementuj w języku „C/C++” algorytm obliczający przybliżone wartości funkcji
* f(x) = x^3/{6[sinh(x)-x]} dla x ∈ [10^-10, 10^3], korzystając z funkcji standardowej sinh().
* W oparciu o zbiór dokładnych wartości tej funkcji, udostępniony przez prowadzącego zajęcia,
* zbadaj jak zmieniają się błędy względne przybliżenia funkcji w tym algorytmie, w zależności od x.
* W tym celu wykonaj rysunek przedstawiający zależność logarytmu dziesiętnego z bezwzględnej wartości
* błędu względnego od logarytmu dziesiętnego z x. Z wykresu odczytaj zakres zmiennej x, w którym błąd
* względny pozostaje na poziomie błędu reprezentacji, oraz zakresy zmiennej x, w których błąd
* względny jest większy. Wyjaśnij przyczyny obserwowanych zmian błędów. Na tej podstawie
* zaproponuj alternatywny sposób obliczania wartości funkcji f(x) w sytuacjach gdy obserwowany
* błąd jest duży. Dokonaj stosownej modyfikacji programu, tak aby uzyskać błąd względny na
* poziomie błędu reprezentacji (czyli tzw. dokładność maszynową) dla dowolnego x ∈ [10^-10, 10^3].
* Jaki typ zmiennych należy zastosować i dlaczego?
* Do wykonania rysunku w tym ćwiczeniu (a także w niektórych dalszych ćwiczeniach) najlepiej użyć
* programu GNUPLOT (dostępnego za darmo z Internetu).
*/

#include <iomanip>
#include <iostream>
#include <valarray>
#include <fstream>
#include <string>
#include <cmath>

long double function(long double x) {
    return x * x * x / (6.0L * (std::sinhl(x) - x));
}

long double relative_error(long double accurate, long double approximate) {
    return std::fabsl(approximate - accurate) / accurate;
}

long double sin_h(long double x) {
    long double sum = 0.0L;
    long double term = x;

    for (int n = 1; n < 1000; n++) {
        sum += term;
        term *= (x * x) / ((2.0L * n) * (2.0L * n + 1));
    }
    return sum;
}

long double alternative_function(long double x) {
    return x * x * x / (6.0L * (sin_h(x) - x));
}

int main() {
    std::ifstream input("../data/Lab2/dane.txt");
    std::ofstream output_1("../data/Lab2/wyniki_1.txt");
    std::ofstream output_2("../data/Lab2/wyniki_2.txt");
    if (!input || !output_1 || !output_2) {
        std::cerr << "Error: Blad przy otwieraniu pliku!" << std::endl;
        return 1;
    }

    std::string line;
    for (int i = 0; i < 3; i++) std::getline(input, line);

    long double log_x, x, result;
    while (input >> log_x >> x >> result) {
        // std::printf("%.5f %.20e %.20f\n", log_x, x, result);
        long double log_error = std::log10l(relative_error(result, function(x)));
        output_1 << std::fixed << std::setprecision(5) << log_x << " ";
        output_1 << std::scientific << std::setprecision(20) << log_error << std::endl;

        long double log_error_alt = std::log10l(relative_error(result, alternative_function(x)));
        output_2 << std::fixed << std::setprecision(5) << log_x << " ";
        output_2 << std::scientific << std::setprecision(20) << log_error_alt << std::endl;
    }

    input.close();
    output_1.close();
    output_2.close();

    return 0;
}
