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

// Notatki
// Nie używać A[i][j], zastosować pomocniczy wektor index[i]=i
// A[index[i]][j], b[index[i]][j]

#include <iostream>
#include <math.h>

#define N 5

double A[N][N] =
{
    {5.0, 4.0, 3.0, 2.0, 1.0},
    {10.0, 8.0, 7.0, 6.0, 5.0},
    {-1.0, 2.0, -3.0, 4.0, -5.0},
    {6.0, 5.0, -4.0, 3.0, -2.0},
    {1.0, 2.0, 3.0, 4.0, 5.0}
};

double b[N] = {37.0, 99.0, -9.0, 12.0, 53.0};

int index[N] = {0, 1, 2, 3, 4};

double L[N][N] = {
    {0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0, 0.0}
};

void lu_decomposition();

void change_pivot(int);

void print_matrix(double [N][N]);

int main() {
    lu_decomposition();
    std::cout << "L MATRIX:" << std::endl;
    print_matrix(A);
    std::cout << "U MATRIX:" << std::endl;
    print_matrix(L);
    return 0;
}

void lu_decomposition() {
    for (int i = 0; i < N - 1; i++) {
        if (A[index[i]][i] == 0.0) change_pivot(i);
        for (int k = i + 1; k < N; k++) {
            double quotient = A[index[k]][i] / A[index[i]][i];
            L[k][i] = quotient;
            for (int j = 0; j < N; j++)
                A[index[k]][j] -= A[index[i]][j] * quotient;
        }
    }

    for (int i = 0; i < N; i++) L[index[i]][i] = 1;
}

void change_pivot(int i) {
    double max = -INFINITY;
    int max_i, temp_i;
    for (int j = i + 1; j < N; j++)
        if (A[index[j]][i] > max) {
            max = A[index[j]][i];
            max_i = j;
        }
    temp_i = index[i];
    index[i] = max_i;
    index[max_i] = temp_i;
}

void print_matrix(double M[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            printf("%5.2f ", M[index[i]][j]);
        std::cout << std::endl;
    }
}
