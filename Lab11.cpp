#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

using Real = long double;
using Vector = std::vector<Real>;
using Matrix = std::vector<Vector>;

constexpr Real PI = std::numbers::pi_v<Real>;

constexpr Real T_MAX = 1.0L;
constexpr Real D = 1.0L;
constexpr Real X_A = 0.0L;
constexpr Real X_B = 1.0L;

constexpr Real LAMBDA_EXPLICIT = 0.4L;
constexpr Real LAMBDA_IMPLICIT = 1.0L;

void thomas_algorithm(Vector d, const Vector &l, const Vector &u, Vector &b, Vector &x);

void lu_decomposition(Matrix &M);

void solve_lu_equation(const Matrix &M, Vector &b);

Real analytical_solution(const Real x, const Real t) {
    return (1.0L - std::exp(-PI * PI * D * t)) * std::sin(PI * x);
}

Real source_term(const Real dt, const Real h, const int i) {
    const Real x_i = X_A + i * h;
    return dt * D * PI * PI * std::sin(PI * x_i);
}

Matrix explicit_method(const Real h) {
    const Real dt = LAMBDA_EXPLICIT * h * h / D;
    const int x_num = static_cast<int>(std::round((X_B - X_A) / h)) + 1;
    const int t_num = static_cast<int>(std::round(T_MAX / dt)) + 1;

    Matrix u(t_num, Vector(x_num, 0.0L));

    for (int k = 0; k < t_num - 1; k++)
        for (int i = 1; i < x_num - 1; i++)
            u[k + 1][i] =
                    LAMBDA_EXPLICIT * u[k][i - 1] +
                    (1.0L - 2.0L * LAMBDA_EXPLICIT) * u[k][i] +
                    LAMBDA_EXPLICIT * u[k][i + 1] +
                    source_term(dt, h, i);

    return u;
}

Matrix implicit_method_thomas(const Real h) {
    const Real dt = LAMBDA_IMPLICIT * h * h / D;
    constexpr Real diagonal_value = -(1.0L + 2.0L * LAMBDA_IMPLICIT);
    const int x_num = static_cast<int>(std::round((X_B - X_A) / h)) + 1;
    const int t_num = static_cast<int>(std::round(T_MAX / dt)) + 1;

    Matrix u(t_num, Vector(x_num, 0.0L));

    Vector diagonal(x_num, diagonal_value);
    diagonal[0] = 1.0L;
    diagonal[x_num - 1] = 1.0L;

    Vector lower(x_num - 1, LAMBDA_IMPLICIT);
    lower[x_num - 2] = 0.0L;

    Vector upper(x_num - 1, LAMBDA_IMPLICIT);
    upper[0] = 0.0L;

    for (int t = 1; t < t_num; t++) {
        Vector b(x_num), x(x_num);

        b[0] = 0.0L;
        b[x_num - 1] = 0.0L;
        for (int i = 1; i < x_num - 1; ++i) b[i] = -u[t - 1][i] - source_term(dt, h, i);

        thomas_algorithm(diagonal, lower, upper, b, x);
        for (int i = 0; i < x_num; ++i) u[t][i] = x[i];
    }

    return u;
}

Matrix implicit_method_lu(const Real h) {
    const Real dt = LAMBDA_IMPLICIT * h * h / D;
    constexpr Real diagonal_value = -(1.0L + 2.0L * LAMBDA_IMPLICIT);
    const int x_num = static_cast<int>(std::round((X_B - X_A) / h)) + 1;
    const int t_num = static_cast<int>(std::round(T_MAX / dt)) + 1;

    Matrix u(t_num, Vector(x_num, 0.0L));
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

    for (int t = 1; t < t_num; t++) {
        b[0] = 0.0L;
        b[x_num - 1] = 0.0L;
        for (int i = 1; i < x_num - 1; i++) b[i] = -u[t - 1][i] - source_term(dt, h, i);

        solve_lu_equation(A, b);
        for (int i = 0; i < x_num; i++) u[t][i] = b[i];
    }

    return u;
}

int main() {
    constexpr int precision = 16;
    constexpr Real h = 0.05L;
    constexpr Real dt = LAMBDA_EXPLICIT * h * h / D;
    std::cout << dt << std::endl;
    Vector t_values = {0.001L, 0.01L, 0.1L, 1.0L};
    constexpr int x_num = static_cast<int>((X_B - X_A) / h) + 1;

    std::ofstream file("../data/Lab11/explicit_values.txt");
    file << std::scientific << std::setprecision(precision);

    // for (int i = 0; i <= 10000; i++) {
    //     const Real x_i = X_A + i * 0.0001;
    //     file << x_i;
    //     for (const Real t: t_values) file << "\t" << analytical_solution(x_i, t);
    //     file << std::endl;
    // }

    Matrix u = explicit_method(h);
    for (int i = 0; i < x_num; i++) {
        const Real x_i = X_A + i * h;
        file << x_i;
        for (const Real t: t_values) {
            const int t_index = static_cast<int>(std::round(t / dt));
            std::cout << t_index << std::endl;
            file << "\t" << u[t_index][i];
        }
        file << std::endl;
    }

    return 0;
}

void thomas_algorithm(Vector d, const Vector &l, const Vector &u, Vector &b, Vector &x) {
    const int n = static_cast<int>(d.size());
    if (n == 0) return;

    for (int i = 1; i < n; i++) {
        const Real m = l[i - 1] / d[i - 1];
        d[i] = d[i] - m * u[i - 1];
        b[i] = b[i] - m * b[i - 1];
    }

    x[n - 1] = b[n - 1] / d[n - 1];
    for (int i = n - 2; i >= 0; i--) x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
}

void lu_decomposition(Matrix &M) {
    const int n = static_cast<int>(M.size());
    if (n == 0) return;

    for (int i = 0; i < n - 1; i++) {
        for (int k = i + 1; k < n; k++) {
            const Real quotient = M[k][i] / M[i][i];
            M[k][i] = quotient;
            for (int j = i + 1; j < n; j++)
                M[k][j] -= M[i][j] * quotient;
        }
    }
}

void solve_lu_equation(const Matrix &M, Vector &b) {
    const int n = static_cast<int>(b.size());
    Real y[n];

    for (int i = 0; i < n; i++) {
        Real sum = 0.0;
        for (int j = 0; j < i; j++) sum += M[i][j] * y[j];
        y[i] = b[i] - sum;
    }

    for (int i = n - 1; i >= 0; i--) {
        Real sum = 0.0;
        for (int j = i + 1; j < n; j++) sum += M[i][j] * b[j];
        b[i] = (y[i] - sum) / M[i][i];
    }
}
