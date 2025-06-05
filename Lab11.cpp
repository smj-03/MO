#include <iostream>
#include <cmath>

#define T_MAX 1
#define A 0.0L
#define B 1.0L
#define D 1.0L
#define DT 0.001L

#define DX_EXPLICIT 0.05L
#define T_NUM_EXPLICIT static_cast<int>(1 / DT) + 1
#define X_NUM_EXPLICIT static_cast<int>((B - A) / DX_EXPLICIT) + 1

#define DX_IMPLICIT 0.05L

#define X_NUM 1000

long double analytical_solution(const long double x, const long double t) {
    return (1.0L - std::expl(-std::numbers::pi * std::numbers::pi * D * t)) * std::sinl(std::numbers::pi * x);
}

void explicit_method(long double (&)[T_NUM_EXPLICIT][X_NUM_EXPLICIT]);

void linear_space(long double, long double, int, long double *);

int main() {
    long double space[X_NUM];
    linear_space(A, B, X_NUM, space);

    FILE *file1 = fopen("../data/Lab11/wykres_analityczny", "w");
    constexpr long double t_values[] = {0.0L, 0.1L, 0.2L, 0.3L, 0.4L, 0.5L};
    constexpr int t_count = std::size(t_values);
    for (int i = 0; i < X_NUM; i++) {
        fprintf(file1, "%.12LE", space[i]);
        for (int j = 0; j < t_count; j++)
            fprintf(file1, "\t%.12LE", analytical_solution(space[i], t_values[j]));
        fprintf(file1, "\n");
    }

    long double u_explicit[T_NUM_EXPLICIT][X_NUM_EXPLICIT] = {0.0};
    explicit_method(u_explicit);

    FILE *file2 = fopen("../data/Lab11/wykres_numeryczny", "w");
    for (int i = 0; i < X_NUM_EXPLICIT; i++) {
        const long double xi = A + i * DX_EXPLICIT;
        fprintf(file2, "%.12LE", xi);
        for (int j = 0; j < t_count; j++) {
            const int t_index = static_cast<int>(t_values[j] / DT); // Map time to matrix index
            fprintf(file2, "\t%.12LE", u_explicit[t_index][i]);
        }
        fprintf(file2, "\n");
    }

    fclose(file1);
    fclose(file2);

    return 0;
}

void explicit_method(long double (&u)[T_NUM_EXPLICIT][X_NUM_EXPLICIT]) {
    constexpr long double lambda = D * DT / (DX_EXPLICIT * DX_EXPLICIT);
    for (int k = 0; k < T_NUM_EXPLICIT - 1; k++) {
        for (int i = 1; i < X_NUM_EXPLICIT - 1; i++) {
            const long double xi = A + i * DX_EXPLICIT;
            u[k + 1][i] = lambda * u[k][i - 1] + (1.0L - 2.0L * lambda) * u[k][i] + lambda * u[k][i + 1] +
                          DT * D * std::numbers::pi * std::numbers::pi * std::sinl(std::numbers::pi * xi);
        }
    }
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
