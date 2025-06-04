#include <iostream>
#include <cmath>

#define A 0.0L
#define B 1.0L
#define D 1.0L
#define H 0.05L
#define DELTA_T 0.001L
#define X_NUM 1000

double analytical_solution(const long double x, const long double t) {
    return (1.0 - std::expl(-std::numbers::pi * std::numbers::pi * D * t)) * std::sinl(std::numbers::pi * x);
}

void linear_space(long double, long double, int, long double *);

void finite_difference_method(long double, long double, long double);

int main() {
    std::cout << analytical_solution(0.5, 0.2);
    long double space[X_NUM];
    linear_space(A, B, X_NUM, space);

    FILE *file = fopen("../data/Lab11/wykres_analityczny", "w");
    constexpr long double t_values[] = {0.0L, 0.1L, 0.2L, 0.3L, 0.4L, 0.5L};
    constexpr int t_count = std::size(t_values);
    for (int i = 0; i < X_NUM; i++) {
        fprintf(file, "%.6fe", space[i]);
        for (int j = 0; j < t_count; j++)
            fprintf(file, "\t%.6fe", analytical_solution(space[i], t_values[j]));
        fprintf(file, "\n");
    }

    return 0;
}

double finite_difference_method(double u) {
}


void linear_space(const long double a, const long double b, const int n, long double *output) {
    if (n <= 0) return;
    if (n == 1) {
        output[0] = n;
        return;
    }

    const long double step = (b - a) / (static_cast<long double>(n) - 1.0);
    for (int i = 0; i < n; i++) output[i] = a + i * step;
}
