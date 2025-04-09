/*
* Napisz program w języku „C/C++”, rozwiązujący układ czterech równań liniowych metodami
* iteracyjnymi: (a) Jacobiego, (b) Gaussa-Seidela, (c) SOR z parametrem ω = 1/2, a następnie zastosuj
* ten program do rozwiązania układu równań liniowych Ax = b, gdzie
*
*     [ 50  5  4  3  2 ]      [ 140 ]                                        [ 6 ]
*     [  1 40  1  2  3 ]      [  67 ]                                        [ 6 ]
* A = [  4  5 30 -5 -4 ], b = [  62 ]. Przyjmij przybliżenie początkowe x0 = [ 6 ]
*     [ -3 -2 -1 20  0 ]      [  89 ]                                        [ 6 ]
*     [  1  2  3  4 30 ]      [ 153 ]                                        [ 6 ]
*
* Zastosuj trzy niezależne kryteria zakończenia iteracji. Zadbaj o to, aby wyprowadzać na konsolę
* wyniki pośrednie obliczeń dla każdej iteracji, tak aby możliwe było obserwowanie zbieżności
* kolejnych przybliżeń pierwiastków i porównanie liczby iteracji niezbędnych do uzyskania (za
* pomocą różnych metod) rozwiązania o zadanej dokładności bezwzględnej. W szczególności oblicz
* jak zmienia się estymator błędu rozwiązania oraz residuum układu w trakcie kolejnych iteracji.
*/

#include <iostream>

int main() {
    std::cout << "Hello World!" << std::endl;
    return 0;
}