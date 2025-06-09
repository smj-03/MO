#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

using Vector = std::vector<long double>;
using Matrix = std::vector<Vector>;

constexpr long double PI = std::numbers::pi_v<long double>;

constexpr long double T_MAX = 1.0L;
constexpr long double D = 1.0L;
constexpr long double X_A = 0.0L;
constexpr long double X_B = 1.0L;

constexpr long double LAMBDA_EXPLICIT = 0.4L;
constexpr long double LAMBDA_IMPLICIT = 1.0L;

/**
 * PARAMETRY DLA METODY KMB
 */
#define DT_EXPLICIT 0.001L
#define DX_EXPLICIT 0.05L
#define T_NUM_EXPLICIT static_cast<int>(1 / DT_EXPLICIT) + 1
#define X_NUM_EXPLICIT static_cast<int>((X_B - X_A) / DX_EXPLICIT) + 1

/**
 * PARAMETRY DLA METODY LAASONEN
 */
#define DT_IMPLICIT 0.0025L
#define DX_IMPLICIT 0.05L
#define T_NUM_IMPLICIT static_cast<int>(1 / DT_IMPLICIT) + 1
#define X_NUM_IMPLICIT static_cast<int>((X_B - X_A) / DX_IMPLICIT) + 1

#define X_NUM 1000

void thomas_algorithm(long double *d, const long double *l, const long double *u, long double *b, long double *x,
                      int n);

void lu_decomposition(long double M[X_NUM_IMPLICIT][X_NUM_IMPLICIT], int n);

void solve_lu_equation(long double M[X_NUM_IMPLICIT][X_NUM_IMPLICIT], long double *b, int n);

long double calculate_max_error(const Matrix &u_matrix, long double h);

long double analytical_solution(const long double x, const long double t) {
    return (1.0L - std::exp(-PI * PI * D * t)) * std::sin(PI * x);
}

Matrix explicit_method(const long double h) {
    const long double dt = LAMBDA_EXPLICIT * h * h / D;
    const int x_num = static_cast<int>((X_B - X_A) / h) + 1;
    const int t_num = static_cast<int>(T_MAX / dt) + 1;

    Matrix u(t_num, Vector(x_num, 0.0L));

    for (int k = 0; k < t_num - 1; k++)
        for (int i = 1; i < x_num - 1; i++) {
            const long double x_i = X_A + i * h;
            const long double source_term = dt * D * PI * PI * std::sin(PI * x_i);
            u[k + 1][i] = LAMBDA_EXPLICIT * u[k][i - 1] + (1.0L - 2.0L * LAMBDA_EXPLICIT) * u[k][i] + LAMBDA_EXPLICIT *
                          u[k][i + 1] + source_term;
        }

    return u;
}

Matrix implicit_method_thomas(const long double h) {
    const long double dt = LAMBDA_IMPLICIT * h * h / D;
    const int x_num = static_cast<int>((X_B - X_A) / h) + 1;
    const int t_num = static_cast<int>(T_MAX / dt) + 1;

    Matrix u(t_num, Vector(x_num, 0.0L));

    constexpr long double diagonal_value = -(1.0L + 2.0L * LAMBDA_IMPLICIT);
    long double diagonal[x_num], lower[x_num - 1], upper[x_num - 1];
    diagonal[0] = 1.0L;
    upper[0] = 0.0L;
    diagonal[x_num - 1] = 1.0L;
    lower[x_num - 2] = 0.0L;
    for (int i = 1; i < x_num - 1; ++i) diagonal[i] = diagonal_value;
    for (int i = 0; i < x_num - 2; ++i) lower[i] = LAMBDA_IMPLICIT;
    for (int i = 1; i < x_num - 1; ++i) upper[i] = LAMBDA_IMPLICIT;

    for (int t = 1; t < t_num; t++) {
        long double d[x_num], b[x_num], x[x_num];
        for (int i = 0; i < x_num; ++i) d[i] = diagonal[i];

        b[0] = 0.0L;
        b[x_num - 1] = 0.0L;
        for (int i = 1; i < x_num - 1; ++i) {
            const long double x_i = X_A + i * h;
            const long double source_term = dt * D * PI * PI * std::sin(PI * x_i);
            b[i] = -u[t - 1][i] - source_term;
        };

        thomas_algorithm(d, lower, upper, b, x, x_num);
        for (int i = 0; i < x_num; ++i) u[t][i] = x[i];
    }

    return u;
}

void implicit_method_thomas(long double u[T_NUM_IMPLICIT][X_NUM_IMPLICIT]);

void implicit_method_lu(long double u[T_NUM_IMPLICIT][X_NUM_IMPLICIT]);

int main() {
    constexpr int precision = 16;

    std::ofstream f_explicit("../data/Lab11/errors_explicit.txt");
    f_explicit << std::scientific << std::setprecision(precision);

    for (long double h = 0.005; h < T_MAX - 0.1; h += 0.0001) {
        Matrix u_explicit = explicit_method(h);
        const long double err_explicit = calculate_max_error(u_explicit, h);
        f_explicit << h << "\t" << err_explicit << std::endl;
    }

    return 0;
}

void thomas_algorithm(long double *d, const long double *l, const long double *u, long double *b, long double *x,
                      const int n) {
    for (int i = 1; i < n; i++) {
        const long double m = l[i - 1] / d[i - 1];
        d[i] = d[i] - m * u[i - 1];
        b[i] = b[i] - m * b[i - 1];
    }

    x[n - 1] = b[n - 1] / d[n - 1];
    for (int i = n - 2; i >= 0; i--) x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
}

void implicit_method_lu(long double u[T_NUM_IMPLICIT][X_NUM_IMPLICIT]) {
    constexpr long double lambda = D * DT_IMPLICIT / (DX_IMPLICIT * DX_IMPLICIT);
    constexpr int N = X_NUM_IMPLICIT;
    constexpr long double diagonal_value = -(1.0L + 2.0L * lambda);

    long double A[N][N] = {0.0}, b[N] = {0.0};

    for (int i = 0; i < N; ++i) {
        A[i][i] = diagonal_value;
        if (i > 0) A[i][i - 1] = lambda;
        if (i < N - 1) A[i][i + 1] = lambda;
    }

    A[0][0] = 1.0L;
    A[0][1] = 0.0L;
    A[N - 1][N - 1] = 1.0L;
    A[N - 1][N - 2] = 0.0L;

    lu_decomposition(A, N);

    for (int t = 1; t <= T_NUM_IMPLICIT; t++) {
        for (int i = 1; i < N - 1; i++) {
            const long double x_i = X_A + i * DX_IMPLICIT;
            const long double source_term = D * DT_IMPLICIT * std::numbers::pi * std::numbers::pi * sinl(
                                                std::numbers::pi * x_i);
            b[i] = -u[t - 1][i] - source_term;
        }
        b[0] = 0.0L;
        b[N - 1] = 0.0L;

        solve_lu_equation(A, b, N);
        for (int i = 0; i < N; i++) u[t][i] = b[i];
    }
}

void lu_decomposition(long double M[X_NUM_IMPLICIT][X_NUM_IMPLICIT], const int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int k = i + 1; k < n; k++) {
            const long double quotient = M[k][i] / M[i][i];
            M[k][i] = quotient;
            for (int j = i + 1; j < n; j++)
                M[k][j] -= M[i][j] * quotient;
        }
    }
}

void solve_lu_equation(long double M[X_NUM_IMPLICIT][X_NUM_IMPLICIT], long double *b, const int n) {
    long double y[n];

    for (int i = 0; i < n; i++) {
        long double sum = 0.0;
        for (int j = 0; j < i; j++) sum += M[i][j] * y[j];
        y[i] = b[i] - sum;
    }

    for (int i = n - 1; i >= 0; i--) {
        long double sum = 0.0;
        for (int j = i + 1; j < n; j++) sum += M[i][j] * b[j];
        b[i] = (y[i] - sum) / M[i][i];
    }
}

long double calculate_max_error(const Matrix &u_matrix, const long double h) {
    const Vector &final_u = u_matrix.back();
    const int x_num = final_u.size();

    long double max_error = 0.0L;

    for (int i = 0; i < x_num; i++) {
        const long double x_i = X_A + i * h;
        const long double numerical_val = final_u[i];
        const long double analytical_val = analytical_solution(x_i, T_MAX);
        max_error = std::max(max_error, std::abs(numerical_val - analytical_val));
    }

    return max_error;
}
