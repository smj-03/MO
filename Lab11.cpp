#include <iostream>
#include <cmath>

#define D 1.0
#define A 0
#define B 1
#define X_NUM 1000

double analytical_solution(double, double);

void linear_space(double, double, int, double *);

int main() {
    std::cout << analytical_solution(0.5, 0.2);
    double space[X_NUM];
    linear_space(A, B, X_NUM, space);

    FILE *file = fopen("../data/Lab11/wykres_analityczny", "w");
    constexpr double t_values[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
    constexpr int t_count = std::size(t_values);
    for (int i = 0; i < X_NUM; i++) {
        fprintf(file, "%.6fe", space[i]);
        for (int j = 0; j < t_count; j++)
            fprintf(file, "\t%.6fe", analytical_solution(space[i], t_values[j]));
        fprintf(file, "\n");
    }

    return 0;
}

double analytical_solution(const double x, const double t) {
    return (1.0 - std::exp(-std::numbers::pi * std::numbers::pi * D * t)) * std::sin(std::numbers::pi * x);
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
