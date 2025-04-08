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

#include <iostream>

#define N 5

void lu_decomposition_n(double [N], double [N], double [N], double [N]);

void lu_decomposition_r(double [N], double [N], double [N], double [N]);

void solve_equation(double [N], double [N], double [N], double [N]);

void print_vector(double [N]);

int main() {
    double Ad[N] = {100.0, 200.0, 300.0, 200.0, 100.0};
    double Al[N - 1] = {2.0, 4.0, -6.0, -8.0};
    double Au[N - 1] = {-1.0, -3.0, 5.0, -7.0};

    double b[N] = {199.0, 195.0, 929.0, 954.0, 360.0};

    double n[N], r[N], x[N];

    lu_decomposition_n(Ad, Al, Au, n);
    std::cout << "n VECTOR:" << std::endl;
    print_vector(n);
    std::cout << std::endl;

    lu_decomposition_r(b, Al, n, r);
    std::cout << "r VECTOR:" << std::endl;
    print_vector(r);
    std::cout << std::endl;

    solve_equation(n, r, Au, x);
    std::cout << "SOLUTION:" << std::endl;
    for (int i = 0; i < N; i++) printf("x%d: %.2f\n", i + 1, x[i]);

    return 0;
}

void lu_decomposition_n(double d[N], double l[N], double u[N], double n[N]) {
    n[0] = d[0];
    for (int i = 1; i < N; i++) n[i] = d[i] - l[i - 1] * u[i - 1] / n[i - 1];
}

void lu_decomposition_r(double b[N], double l[N], double n[N], double r[N]) {
    r[0] = b[0];
    for (int i = 1; i < N; i++) r[i] = b[i] - l[i - 1] * r[i - 1] / n[i - 1];
}

void solve_equation(double n[N], double r[N], double u[N], double x[N]) {
    x[N - 1] = r[N - 1] / n[N - 1];
    for (int i = N - 2; i >= 0; i--) x[i] = (r[i] - u[i] * x[i + 1]) / n[i];
}

void print_vector(double v[N]) {
    for (int i = 0; i < N; i++) printf("%6.2f ", v[i]);
    std::cout << std::endl;
}
