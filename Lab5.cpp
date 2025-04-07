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
#include <cmath>

#define N 5

int index[N] = {0, 1, 2, 3, 4};

void lu_decomposition(double [N][N], double [N][N]);

void change_pivot(double [N][N], int);

void solve_equation(double [N][N], double [N][N], const double [N], double [N]);

void print_matrix(double [N][N]);

int main() {
    double A[N][N] =
    {
        {5.0, 4.0, 3.0, 2.0, 1.0},
        {10.0, 8.0, 7.0, 6.0, 5.0},
        {-1.0, 2.0, -3.0, 4.0, -5.0},
        {6.0, 5.0, -4.0, 3.0, -2.0},
        {1.0, 2.0, 3.0, 4.0, 5.0}
    };

    double b[N] = {37.0, 99.0, -9.0, 12.0, 53.0};

    double L[N][N], x[N];

    lu_decomposition(A, L);
    std::cout << "L MATRIX:" << std::endl;
    print_matrix(A);
    std::cout << std::endl;

    std::cout << "U MATRIX:" << std::endl;
    print_matrix(L);
    solve_equation(A, L, b, x);
    std::cout << std::endl;

    std::cout << "SOLUTION:" << std::endl;
    for (int i = 0; i < N; i++) printf("x%d: %.2f\n", i, x[i]);
    return 0;
}

void lu_decomposition(double M[N][N], double L[N][N]) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            L[i][j] = 0;

    for (int i = 0; i < N - 1; i++) {
        if (M[index[i]][i] == 0.0) change_pivot(M, i);
        for (int k = i + 1; k < N; k++) {
            double quotient = M[index[k]][i] / M[index[i]][i];
            L[k][i] = quotient;
            for (int j = 0; j < N; j++)
                M[index[k]][j] -= M[index[i]][j] * quotient;
        }
    }

    for (int i = 0; i < N; i++) L[index[i]][i] = 1;
}

void change_pivot(double M[N][N], const int i) {
    double max = -INFINITY;
    int max_i = i;
    for (int j = i + 1; j < N; j++)
        if (M[index[j]][i] > max) {
            max = M[index[j]][i];
            max_i = j;
        }
    const int temp_i = index[i];
    index[i] = max_i;
    index[max_i] = temp_i;
}

void solve_equation(double U[N][N], double L[N][N], const double b[N], double x[N]) {
    double y[N];

    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) sum += L[index[i]][j] * y[j];
        y[i] = b[index[i]] - sum;
    }

    for (int i = N - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < N; j++) sum += U[index[i]][j] * x[j];
        x[i] = (y[i] - sum) / U[index[i]][i];
    }
}

void print_matrix(double M[N][N]) {
    for (const int i : index) {
        for (int j = 0; j < N; j++)
            printf("%5.2f ", M[i][j]);
        std::cout << std::endl;
    }
}
