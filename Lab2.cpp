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

#include <iostream>

int main() {
    std::cout << "Hello World!" << std::endl;
    return 0;
}
