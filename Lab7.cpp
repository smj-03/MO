/*
* Napisz program w języku „C/C++”, rozwiązujący układ czterech równań liniowych metodami
* iteracyjnymi: (a) Jacobiego, (b) Gaussa-Seidela, (c) SOR z parametrem ω = 1/2, a następnie zastosuj
* ten program do rozwiązania układu równań liniowych Ax = b, gdzie
*
*     [ 50  5  4  3  2 ]      [ 140 ]                                        [ 6 ]
*     [  1 40  1  2  3 ]      [  67 ]                                        [ 6 ]
* A = [  4  5 30 -5 -4 ], b = [  62 ]. Przyjmij przybliżenie początkowe x0 = [ 6 ]
*     [ -3 -2 -1 20  0 ]      [  89 ]                                        [ 6 ]
*     [  1  2  3  4 30 ]      [ 153 ]                                        [ 6 ]
*
* Zastosuj trzy niezależne kryteria zakończenia iteracji. Zadbaj o to, aby wyprowadzać na konsolę
* wyniki pośrednie obliczeń dla każdej iteracji, tak aby możliwe było obserwowanie zbieżności
* kolejnych przybliżeń pierwiastków i porównanie liczby iteracji niezbędnych do uzyskania (za
* pomocą różnych metod) rozwiązania o zadanej dokładności bezwzględnej. W szczególności oblicz
* jak zmienia się estymator błędu rozwiązania oraz residuum układu w trakcie kolejnych iteracji.
*/

#include <iostream>
#include <cmath>

#define N 5

#define N_MAX 100
#define TOLX 1e-8
#define TOLF 1e-8

void jacobi(double [N][N], double [N], double [N]);

void gauss_seidel(double [N][N], double [N], double [N]);

void sor(double [N][N], double [N], double [N], double);

void calculate_residuum(double [N][N], double [N], double [N], double [N]);

double vector_norm(double [N]);

void print_vector(double [N]);

int main() {
    double M[N][N] = {
        {50.0, 5.0, 4.0, 3.0, 2.0},
        {1.0, 40.0, 1.0, 2.0, 3.0},
        {4.0, 5.0, 30.0, -5.0, -4.0},
        {-3.0, -2.0, -1.0, 20.0, 0.0},
        {1.0, 2.0, 3.0, 4.0, 30.0}
    };

    double b[N] = {140.0, 67.0, 62.0, 89.0, 153.0};
    double x[N] = {6.0, 6.0, 6.0, 6.0, 6.0};

    jacobi(M, b, x);
    gauss_seidel(M, b, x);
    sor(M, b, x, 0.5);

    return 0;
}

void jacobi(double M[N][N], double b[N], double x0[N]) {
    double x[N], u[N];
    for (int i = 0; i < N; i++) x[i] = x0[i];

    std::cout << "---------------------- Jacobi ---------------------" << std::endl;
    std::printf("%-4s%-36s%-24s%-17s\n", "i", "Wektor x", "Estymator", "Residuum");

    for (int iteration = 1; iteration < N_MAX; iteration++) {
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 0; j < N; j++)
                if (j != i) sum += M[i][j] * x[j];
            u[i] = (b[i] - sum) / M[i][i];
        }

        double estimator[N], residuum[N];
        for (int i = 0; i < N; i++) estimator[i] = u[i] - x[i];
        calculate_residuum(M, b, x, residuum);
        const double estimator_norm = vector_norm(estimator);
        const double residuum_norm = vector_norm(residuum);
        if (estimator_norm < TOLX && residuum_norm< TOLF) break;

        std::cout << iteration << ". ";
        print_vector(x);
        std::printf("%2.10f\t\t%2.10f", estimator_norm, residuum_norm);
        std::cout << std::endl;

        for (int i = 0; i < N; i++) x[i] = u[i];
    }
}

void gauss_seidel(double M[N][N], double b[N], double x0[N]) {
    double x[N], x_old[N];
    for (int i = 0; i < N; i++) x[i] = x0[i];

    std::cout << "------------------- Gauss-Seidel ------------------" << std::endl;
    std::printf("%-4s%-36s%-24s%-17s\n", "i", "Wektor x", "Estymator", "Residuum");

    for (int iteration = 1; iteration < N_MAX; iteration++) {
        for (int i = 0; i < N; i++) x_old[i] = x[i];

        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 0; j < N; j++)
                if (j != i) sum += M[i][j] * x[j];
            x[i] = (b[i] - sum) / M[i][i];
        }

        double estimator[N], residuum[N];
        for (int i = 0; i < N; i++) estimator[i] = x[i] - x_old[i];
        calculate_residuum(M, b, x, residuum);
        const double estimator_norm = vector_norm(estimator);
        const double residuum_norm = vector_norm(residuum);
        if (estimator_norm < TOLX && residuum_norm< TOLF) break;

        std::cout << iteration << ". ";
        print_vector(x);
        std::printf("%2.10f\t\t%2.10f", estimator_norm, residuum_norm);
        std::cout << std::endl;
    }
}

void sor(double M[N][N], double b[N], double x0[N], const double w) {
    double x[N], x_old[N];
    for (int i = 0; i < N; i++) x[i] = x0[i];

    std::cout << "--------------------- SOR ------------------------" << std::endl;
    std::printf("%-4s%-36s%-24s%-17s\n", "i", "Wektor x", "Estymator", "Residuum");

    for (int iteration = 1; iteration < N_MAX; iteration++) {
        for (int i = 0; i < N; i++) x_old[i] = x[i];

        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 0; j < N; j++)
                if (j != i) sum += M[i][j] * x[j];
            x[i] = (1.0 - w) * x[i] + w * (b[i] - sum) / M[i][i];
        }

        double estimator[N], residuum[N];
        for (int i = 0; i < N; i++) estimator[i] = x[i] - x_old[i];
        calculate_residuum(M, b, x, residuum);
        const double estimator_norm = vector_norm(estimator);
        const double residuum_norm = vector_norm(residuum);
        if (estimator_norm < TOLX && residuum_norm< TOLF) break;

        std::cout << iteration << ". ";
        print_vector(x);
        std::printf("%2.10f\t\t%2.10f", estimator_norm, residuum_norm);
        std::cout << std::endl;
    }
}

void calculate_residuum(double M[N][N], double b[N], double x[N], double residuum[N]) {
    for (int i = 0; i < N; ++i) {
        double Mx_i = 0.0;
        for (int j = 0; j < N; ++j) {
            Mx_i += M[i][j] * x[j];
        }
        residuum[i] = b[i] - Mx_i;
    }
}

double vector_norm(double x[N]) {
    double sum = 0.0;
    for (int i = 0; i < N; i++)
        sum += x[i] * x[i];
    return std::sqrt(sum);
}

void print_vector(double x[5]) {
    std::printf("%5.2f %5.2f %5.2f %5.2f %5.2f\t", x[0], x[1], x[2], x[3], x[4]);
}
