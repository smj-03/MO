#include <iostream>
#include <cmath>

#define T_MAX 1
#define A 0.0L
#define B 1.0L
#define D 1.0L

/**
 * PARAMETRY DLA METODY KMB
 */
#define DT_EXPLICIT 0.001L
#define DX_EXPLICIT 0.05L
#define T_NUM_EXPLICIT static_cast<int>(1 / DT_EXPLICIT) + 1
#define X_NUM_EXPLICIT static_cast<int>((B - A) / DX_EXPLICIT) + 1

/**
 * PARAMETRY DLA METODY LAASONEN
 */
#define DT_IMPLICIT 0.0025L
#define DX_IMPLICIT 0.05L
#define T_NUM_IMPLICIT static_cast<int>(1 / DT_IMPLICIT) + 1
#define X_NUM_IMPLICIT static_cast<int>((B - A) / DX_IMPLICIT) + 1

#define X_NUM 1000

long double analytical_solution(const long double x, const long double t) {
    return (1.0L - std::expl(-std::numbers::pi * std::numbers::pi * D * t)) * std::sinl(std::numbers::pi * x);
}

void explicit_method(long double (&u)[T_NUM_EXPLICIT][X_NUM_EXPLICIT]);

void implicit_method(long double (&u)[T_NUM_IMPLICIT][X_NUM_IMPLICIT]);

void thomas_algorithm(long double *d, const long double *l, const long double *u, long double *b, long double *x,
                      int n);

void linear_space(long double a, long double b, int n, long double *output);

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
    fclose(file1);

    long double u_explicit[T_NUM_EXPLICIT][X_NUM_EXPLICIT] = {0.0};
    explicit_method(u_explicit);

    FILE *file2 = fopen("../data/Lab11/wykres_numeryczny", "w");
    for (int i = 0; i < X_NUM_EXPLICIT; i++) {
        const long double xi = A + i * DX_EXPLICIT;
        fprintf(file2, "%.12LE", xi);
        for (int j = 0; j < t_count; j++) {
            const int t_index = static_cast<int>(t_values[j] / DT_EXPLICIT);
            fprintf(file2, "\t%.12LE", u_explicit[t_index][i]);
        }
        fprintf(file2, "\n");
    }
    fclose(file2);

    long double u_implicit[T_NUM_IMPLICIT][X_NUM_IMPLICIT] = {0.0};
    implicit_method(u_implicit);

    FILE *file3 = fopen("../data/Lab11/wykres_numeryczny_laasonen", "w");
    for (int i = 0; i < X_NUM_IMPLICIT; i++) {
        const long double xi = A + i * DX_IMPLICIT;
        fprintf(file3, "%.12LE", xi);
        for (int j = 0; j < t_count; j++) {
            const int t_index = static_cast<int>(t_values[j] / DT_IMPLICIT);
            fprintf(file3, "\t%.12LE", u_implicit[t_index][i]);
        }
        fprintf(file3, "\n");
    }
    fclose(file3);

    return 0;
}

void explicit_method(long double (&u)[T_NUM_EXPLICIT][X_NUM_EXPLICIT]) {
    constexpr long double lambda = D * DT_EXPLICIT / (DX_EXPLICIT * DX_EXPLICIT);
    for (int k = 0; k < T_NUM_EXPLICIT - 1; k++)
        for (int i = 1; i < X_NUM_EXPLICIT - 1; i++) {
            const long double x_i = A + i * DX_EXPLICIT;
            u[k + 1][i] = lambda * u[k][i - 1] + (1.0L - 2.0L * lambda) * u[k][i] + lambda * u[k][i + 1] +
                          DT_EXPLICIT * D * std::numbers::pi * std::numbers::pi * std::sinl(std::numbers::pi * x_i);
        }
}

void implicit_method(long double (&u)[T_NUM_IMPLICIT][X_NUM_IMPLICIT]) {
    constexpr long double lambda = D * DT_IMPLICIT / (DX_IMPLICIT * DX_IMPLICIT);
    constexpr int N = X_NUM_IMPLICIT;
    constexpr long double diagonal_value = -(1.0L + 2.0L * lambda);

    long double diagonal[N], lower[N - 1], upper[N - 1];
    diagonal[0] = 1.0L;
    upper[0] = 0.0L;
    diagonal[N - 1] = 1.0L;
    lower[N - 2] = 0.0L;
    for (int i = 1; i < N - 1; ++i) diagonal[i] = diagonal_value;
    for (int i = 0; i < N - 2; ++i) lower[i] = lambda;
    for (int i = 0; i < N - 1; ++i) upper[i] = lambda;

    for (int t = 1; t < T_NUM_IMPLICIT; t++) {
        long double d[N], b[N], x[N];
        for (int i = 0; i < N; ++i) d[i] = diagonal[i];

        b[0] = 0.0L;
        b[N - 1] = 0.0L;
        for (int i = 1; i < N - 1; ++i) {
            const long double x_i = A + i * DX_IMPLICIT;
            const long double source_term = D * DT_IMPLICIT * std::numbers::pi * std::numbers::pi * sinl(
                                                std::numbers::pi * x_i);
            b[i] = -u[t - 1][i] - source_term;
        };

        thomas_algorithm(d, lower, upper, b, x, N);
        for (int i = 0; i < N; ++i) u[t][i] = x[i];
    }
}

void thomas_algorithm(long double *d, const long double *l, const long double *u, long double *b, long double *x,
                      const int n) {
    for (int i = 1; i < n; i++) {
        const double m = l[i - 1] / d[i - 1];
        d[i] = d[i] - m * u[i - 1];
        b[i] = b[i] - m * b[i - 1];
    }

    x[n - 1] = b[n - 1] / d[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
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
