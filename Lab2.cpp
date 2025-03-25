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

double function(double x) {
    return x * x * x / (6.0 * (std::sinh(x) - x));
}

double relative_error(double accurate, double approximate) {
    return std::fabs(approximate - accurate) / accurate;
}

double sin_h(double x) {
    double result = 0.0;
    double temp = x;

    for (int n = 1; n < 1000; n++) {
        result += temp;
        temp *= x / (double) (2 * n);
    }

    return result;
}

double alternative_function(double x) {
    return x * x * x / (6.0 * (sin_h(x) - x));
}

int main() {
    std::ifstream input("../data/dane_do_laboratorium_2.txt");
    std::ofstream output_1("../results/wyniki_1_lab_2.txt");
    std::ofstream output_2("../results/wyniki_2_lab_2.txt");
    if (!input || !output_1 || !output_2) {
        std::cerr << "Error: Blad przy otwieraniu pliku!" << std::endl;
        return 1;
    }

    std::string line;
    for (int i = 0; i < 3; i++) std::getline(input, line);

    double log_x, x, result;
    while (input >> log_x >> x >> result) {
        // std::printf("%.5f %.20e %.20f\n", log_x, x, result);
        double log_error = std::log10(relative_error(result, function(x)));
        output_1 << std::fixed << std::setprecision(5) << log_x << " ";
        output_1 << std::scientific << std::setprecision(20) << log_error << std::endl;

        double log_error_alt = std::log10(relative_error(result, alternative_function(x)));
        output_2 << std::fixed << std::setprecision(5) << log_x << " ";
        output_2 << std::scientific << std::setprecision(20) << log_error_alt << std::endl;
    }

    input.close();
    output_1.close();
    output_2.close();

    return 0;
}

/*
plot \
"C:/Users/szymon.chwastek/Documents/C++/PK/MO/Zad_1/results/wyniki_1_lab_2.txt" using 1:2 with lines title "log10(x) od log10(error)", \
"C:/Users/szymon.chwastek/Documents/C++/PK/MO/Zad_1/results/wyniki_2_lab_2.txt" using 1:2 with lines title "log10(x) od log10(alt_error)", \
-16 with lines lc rgb "red" title "Blad reprezentacji"
*/
