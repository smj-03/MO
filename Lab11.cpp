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

void thomas_algorithm(Vector d, const Vector &l, const Vector &u, Vector b, Vector x);

void lu_decomposition(Matrix M);

void solve_lu_equation(const Matrix &M, Vector b);

long double calculate_max_error(const Matrix &u, long double h);

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
    Vector diagonal(x_num), lower(x_num - 1), upper(x_num - 1);
    diagonal[0] = 1.0L;
    upper[0] = 0.0L;
    diagonal[x_num - 1] = 1.0L;
    lower[x_num - 2] = 0.0L;
    for (int i = 1; i < x_num - 1; ++i) diagonal[i] = diagonal_value;
    for (int i = 0; i < x_num - 2; ++i) lower[i] = LAMBDA_IMPLICIT;
    for (int i = 1; i < x_num - 1; ++i) upper[i] = LAMBDA_IMPLICIT;

    for (int t = 1; t < t_num; t++) {
        Vector d(x_num), b(x_num), x(x_num);
        for (int i = 0; i < x_num; ++i) d[i] = diagonal[i];

        b[0] = 0.0L;
        b[x_num - 1] = 0.0L;
        for (int i = 1; i < x_num - 1; ++i) {
            const long double x_i = X_A + i * h;
            const long double source_term = dt * D * PI * PI * std::sin(PI * x_i);
            b[i] = -u[t - 1][i] - source_term;
        }

        thomas_algorithm(d, lower, upper, b, x);
        for (int i = 0; i < x_num; ++i) u[t][i] = x[i];
    }

    return u;
}

Matrix implicit_method_lu(const long double h) {
    const long double dt = LAMBDA_IMPLICIT * h * h / D;
    const int x_num = static_cast<int>((X_B - X_A) / h) + 1;
    const int t_num = static_cast<int>(T_MAX / dt) + 1;

    Matrix u(t_num, Vector(x_num, 0.0L));

    constexpr long double diagonal_value = -(1.0L + 2.0L * LAMBDA_IMPLICIT);

    Matrix A(x_num, Vector(x_num, 0.0L));
    Vector b(x_num);

    for (int i = 0; i < x_num; ++i) {
        A[i][i] = diagonal_value;
        if (i > 0) A[i][i - 1] = LAMBDA_IMPLICIT;
        if (i < x_num - 1) A[i][i + 1] = LAMBDA_IMPLICIT;
    }

    A[0][0] = 1.0L;
    A[0][1] = 0.0L;
    A[x_num - 1][x_num - 1] = 1.0L;
    A[x_num - 1][x_num - 2] = 0.0L;

    lu_decomposition(A);

    for (int t = 1; t <= t_num; t++) {
        for (int i = 1; i < x_num - 1; i++) {
            const long double x_i = X_A + i * h;
            const long double source_term = dt * D * PI * PI * std::sin(PI * x_i);
            b[i] = -u[t - 1][i] - source_term;
        }
        b[0] = 0.0L;
        b[x_num - 1] = 0.0L;

        solve_lu_equation(A, b);
        for (int i = 0; i < x_num; i++) u[t][i] = b[i];
    }

    return u;
}

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

void thomas_algorithm(Vector d, const Vector &l, const Vector &u, Vector b, Vector x) {
    const int n = d.size();

    for (int i = 1; i < n; i++) {
        const long double m = l[i - 1] / d[i - 1];
        d[i] = d[i] - m * u[i - 1];
        b[i] = b[i] - m * b[i - 1];
    }

    x[n - 1] = b[n - 1] / d[n - 1];
    for (int i = n - 2; i >= 0; i--) x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
}

void lu_decomposition(Matrix M) {
    const int n = M.size();
    for (int i = 0; i < n - 1; i++) {
        for (int k = i + 1; k < n; k++) {
            const long double quotient = M[k][i] / M[i][i];
            M[k][i] = quotient;
            for (int j = i + 1; j < n; j++)
                M[k][j] -= M[i][j] * quotient;
        }
    }
}

void solve_lu_equation(const Matrix &M, Vector b) {
    const int n = M.size();
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

long double calculate_max_error(const Matrix &u, const long double h) {
    const Vector &final_u = u.back();
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
