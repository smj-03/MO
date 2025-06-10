﻿#include <iostream>
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

Real calculate_max_error(const Matrix &u, Real h);

Real analytical_solution(const Real x, const Real t) {
    return (1.0L - std::exp(-PI * PI * D * t)) * std::sin(PI * x);
}

Real source_term(const Real dt, const Real h, const int i) {
    const Real x_i = X_A + i * h;
    return dt * D * PI * PI * std::sin(PI * x_i);
}

Matrix explicit_method(const Real h) {
    const Real dt = LAMBDA_EXPLICIT * h * h / D;
    const int x_num = static_cast<int>((X_B - X_A) / h) + 1;
    const int t_num = static_cast<int>(T_MAX / dt) + 1;

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
    const int x_num = static_cast<int>((X_B - X_A) / h) + 1;
    const int t_num = static_cast<int>(T_MAX / dt) + 1;

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
    const int x_num = static_cast<int>((X_B - X_A) / h) + 1;
    const int t_num = static_cast<int>(T_MAX / dt) + 1;

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
    constexpr Real DT_EXPLICIT = 0.001L;
    constexpr Real DX_EXPLICIT = 0.05L;
    constexpr Real X_NUM_EXPLICIT = static_cast<int>((X_B - X_A) / DX_EXPLICIT) + 1;

    constexpr Real DT_IMPLICIT = 0.0025L;
    constexpr Real DX_IMPLICIT = 0.05L;
    constexpr Real X_NUM_IMPLICIT = static_cast<int>((X_B - X_A) / DX_IMPLICIT) + 1;

    constexpr Real t_values[] = {0.0L, 0.1L, 0.2L, 0.3L, 0.4L, 0.5L};

    const Matrix u_explicit = explicit_method(DX_EXPLICIT);

    FILE *file2 = fopen("../data/Lab11/wykres_numeryczny", "w");
    for (int i = 0; i < X_NUM_EXPLICIT; i++) {
        const Real xi = X_A + i * DX_EXPLICIT;
        fprintf(file2, "%.12LE", xi);
        for (const Real t_value: t_values) {
            const int t_index = static_cast<int>(t_value / DT_EXPLICIT);
            fprintf(file2, "\t%.12LE", u_explicit[t_index][i]);
        }
        fprintf(file2, "\n");
    }
    fclose(file2);

    const Matrix u_implicit_thomas = implicit_method_thomas(DX_IMPLICIT);

    FILE *file3 = fopen("../data/Lab11/wykres_numeryczny_laasonen_thomas", "w");
    for (int i = 0; i < X_NUM_IMPLICIT; i++) {
        const Real xi = X_A + i * DX_IMPLICIT;
        fprintf(file3, "%.12LE", xi);
        for (const Real t_value: t_values) {
            const int t_index = static_cast<int>(t_value / DT_IMPLICIT);
            fprintf(file3, "\t%.12LE", u_implicit_thomas[t_index][i]);
        }
        fprintf(file3, "\n");
    }
    fclose(file3);

    const Matrix u_implicit_lu = implicit_method_lu(DX_IMPLICIT);

    FILE *file4 = fopen("../data/Lab11/wykres_numeryczny_laasonen_lu", "w");
    for (int i = 0; i < X_NUM_IMPLICIT; i++) {
        const Real xi = X_A + i * DX_IMPLICIT;
        fprintf(file4, "%.12LE", xi);
        for (const Real t_value: t_values) {
            const int t_index = static_cast<int>(t_value / DT_IMPLICIT);
            fprintf(file4, "\t%.12LE", u_implicit_lu[t_index][i]);
        }
        fprintf(file4, "\n");
    }
    fclose(file4);

    return 0;

    // constexpr int precision = 16;
    //
    // std::ofstream f_explicit("../data/Lab11/errors_explicit.txt");
    // f_explicit << std::scientific << std::setprecision(precision);
    //
    // for (Real h = 0.005; h < T_MAX - 0.1; h += 0.0001) {
    //     Matrix u_explicit = explicit_method(h);
    //     const Real err_explicit = calculate_max_error(u_explicit, h);
    //     f_explicit << h << "\t" << err_explicit << std::endl;
    // }
    //
    // return 0;
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

Real calculate_max_error(const Matrix &u, const Real h) {
    const Vector &final_u = u.back();
    const int x_num = static_cast<int>(final_u.size());

    Real max_error = 0.0L;

    for (int i = 0; i < x_num; i++) {
        const Real x_i = X_A + i * h;
        const Real numerical_val = final_u[i];
        const Real analytical_val = analytical_solution(x_i, T_MAX);
        max_error = std::max(max_error, std::abs(numerical_val - analytical_val));
    }

    return max_error;
}
