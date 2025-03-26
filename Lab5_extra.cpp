/*
* Dany jest układ równań z macierzą
*
*     [ 1+e  1    1           ]                   [ 6+e  ]
* A = [  1  1+e  -3           ] oraz wektorem b = [ 6+2e ]
*     [  1   1   300   5      ]                   [ 6+2e ]
*     [  1   1   -6   200  -7 ]                   [ 6+e  ]
*
* Znajdź najpierw rozwiązanie analityczne, a następnie numeryczne, przyjmując coraz mniejsze wartości
* e = 10^-5, 10^-6, 10^-7, itd.
* Porównaj rozwiązania numeryczne z analitycznym, i wyjaśnij ewentualne obserwowane zmiany
* błędu rozwiązań numerycznych.
*/

// Notatki
// Np. metoda wyznaczników