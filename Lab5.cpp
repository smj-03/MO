/*
* Dana jest macierz
*     [ 5    4    3    2    1 ]                 [ 37 ]
*     [ 10   8    7    6    5 ]                 [ 99 ]
* A = [-1    2   -3    4   -5 ] oraz wektor b = [ -9 ]
*     [ 6    5   -4    3   -2 ]                 [ 12 ]
*     [ 1    2    3    4    5 ]                 [ 53 ]
*
* Napisz program w języku "C/C++", realizujący dekompozycję LU macierzy A, przy zastosowaniu
* eliminacji Gaussa z częściowym wyborem elementu podstawowego, a następnie rozwiązujący układ
* równań Ax = b.
*
* Uwaga: należy zaimplementować wariant dekompozycji omawiany na wykładzie.
*
* Program należy zrealizować w postaci dwóch odrębnych procedur: jednej, która operuje wyłącznie
* na macierzy A, i drugiej, która operuje wyłącznie na wektorze b, korzystając z wyników działania
* procedury pierwszej.
*/

// Nie używać A[i][j], zastosować pomocniczy wektor index[i]=i
// A[index[i][j], b[index[i]][j]

// Created by szymon.chwastek on 3/26/2025.
// 5.2 Metoda Wyznaczników
// 6 3 tablice przechowujace przekatne
