/*
* Napisz program w języku „C/C++”, realizujący algorytm Thomasa dla macierzy trój-diagonalnej o
* dowolnych rozmiarach N x N, a następnie zastosuj ten program do rozwiązania układu równań
* Ax = b, w którym
*
*     [ 100 -1                ]                 [ 199 ]
*     [  2  200  -3           ]                 [ 195 ]
* A = [      4   300   5      ] oraz wektor b = [ 929 ]
*     [          -6   200  -7 ]                 [ 954 ]
*     [               -8  100 ]                 [ 360 ]
*
* Program należy zrealizować w postaci dwóch odrębnych procedur: jednej, która operuje wyłącznie
* na macierzy A, i drugiej, która operuje wyłącznie na wektorze b, korzystając z wyników działania
* procedury pierwszej.
* Uwaga: ponieważ macierz trój-diagonalna jest macierzą rzadką, więc w programie
* NIE NALEŻY używać tablic kwadratowych do reprezentacji macierzy A.
*/

// Notatki
// Macierz A w postaci 3 wektorów z przekątnymi