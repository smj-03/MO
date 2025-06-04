/*
 * Napisz program w języku „C/C++”, demonstrujący zjawisko Rungego w interpolacji
 * wielomianowej Lagrange’a, na przykładzie interpolacji funkcji f(x) = x/(1 + x^4),
 * określonej na przedziale [−1, 1]. Zastosuj bazę Lagrange’a do konstrukcji wielomianów interpolacyjnych.
 * Porównaj wyniki interpolacji na węzłach równoodległych z wynikami interpolacji na węzłach
 * Czebyszewa. Wykonaj wykres interpolowanej funkcji oraz uzyskanych wielomianów interpolacyjnych.
 */

#include <cmath>
#include <iostream>

#define A -1
#define B 1
#define N 60
#define X_NUM 10000

double function(double);

void linear_space(double, double, int, double *);

void czebyszew_space(double, double, int, double *);

double polynomial_interpolation(double, int, double *, double *);

void print_array(double array[], int lenght) {
    for (int i = 0; i < lenght; i++)
        std::cout << array[i] << " ";
}

void plot_results();

int main() {
    double space[X_NUM];
    linear_space(A, B, X_NUM, space);

    double linear_nodes[N], linear_values[N], czebyszew_nodes[N], czebyszew_values[N];
    linear_space(A, B, N, linear_nodes);
    czebyszew_space(A, B, N, czebyszew_nodes);
    for (int i = 0; i < N; i++) {
        linear_values[i] = function(linear_nodes[i]);
        czebyszew_values[i] = function(czebyszew_nodes[i]);
    };

    FILE *file1 = fopen("../data/Lab12/punkty_funkcji.txt", "w");
    FILE *file2 = fopen("../data/Lab12/punkty_interpolacji_lagrange.txt", "w");
    FILE *file3 = fopen("../data/Lab12/punkty_interpolacji_czebyszew.txt", "w");
    for (int i = 0; i < X_NUM; i++) {
        fprintf(file1, "%.6fe\t%.6fe\n", space[i], function(space[i]));
        fprintf(file2, "%.6fe\t%.6fe\n", space[i],
                polynomial_interpolation(space[i], N, linear_nodes, linear_values));
        fprintf(file3, "%.6fe\t%.6fe\n", space[i],
                polynomial_interpolation(space[i], N, czebyszew_nodes, czebyszew_values));
    }

    return 0;
}

double function(const double x) {
    return x / (1.0 + x * x * x * x);
}

double polynomial_interpolation(double x, int n, double *nodes, double *values) {
    double p[n][n];
    for (int i = 0; i < n; i++) p[i][0] = values[i];
    for (int j = 1; j < n; j++)
        for (int i = 0; i < n - j; i++)
            p[i][j] = ((x - nodes[i + j]) * p[i][j - 1] +
                       (nodes[i] - x) * p[i + 1][j - 1]) /
                      (nodes[i] - nodes[i + j]);
    return p[0][n - 1];
}

void linear_space(double a, double b, int n, double *output) {
    if (n <= 0) return;
    if (n == 1) {
        output[0] = n;
        return;
    }

    double step = (b - a) / (static_cast<double>(n) - 1.0);
    for (int i = 0; i < n; i++) output[i] = a + i * step;
}

void czebyszew_space(double a, double b, int n, double *output) {
    if (n <= 0) return;
    if (n == 1) {
        output[0] = n;
        return;
    }

    for (int i = 0; i < n; i++)
        output[i] = (b + a) / 2.0 + (b - a) / 2.0 *
                    std::cos((2.0 * i + 1.0) / (2.0 * n + 2.0) * std::numbers::pi);
}
