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

#define N 5

#define N_MAX 100
#define TOLX 1e-12
#define TOLF 1e-12

void jacobi(double [N][N], double [N][N], double [N][N], double [N], double [N]);

void ldu_decomposition(double [N][N], double [N][N], double [N][N], double [N][N]);

bool invert_matrix(double [5][5], double [5][5]);

void add_matrices(double [N][N], double [N][N], double [N][N]);

void print_matrix(double [N][N]);

int main() {
    double M[N][N] = {
        {50.0, 5.0, 4.0, 3.0, 2.0},
        {1.0, 40.0, 1.0, 2.0, 3.0},
        {4.0, 5.0, 30.0, -5.0, -4.0},
        {-3.0, -2.0, -1.0, 20.0, 0.0},
        {1.0, 2.0, 3.0, 4.0, 30.0}
    };

    double b[N] = {140.0, 67.0, 62.0, 89.0, 153.0};

    double L[N][N], D[N][N], U[N][N];

    ldu_decomposition(M, L, D, U);

    std::cout << "L MATRIX:" << std::endl;
    print_matrix(L);
    std::cout << std::endl;

    std::cout << "D MATRIX:" << std::endl;
    print_matrix(D);
    std::cout << std::endl;

    std::cout << "U MATRIX:" << std::endl;
    print_matrix(U);
    std::cout << std::endl;

    return 0;
}

void jacobi(double M[N][N], double b[N], double x[N]) {
    double L[N][N], D[N][N], U[N][N];
    ldu_decomposition(M, L, D, U);

    double inv_D[N][N];
    invert_matrix(D, inv_D);

    for (int i = 0; i <)

        // for (int i = 1; i < N_MAX; i++) {
        //
        // }



}

void ldu_decomposition(double M[N][N], double L[N][N], double D[N][N], double U[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (j < i) L[i][j] = M[i][j];
            else L[i][j] = 0.0;

            if (i == j) D[i][j] = M[i][j];
            else D[i][j] = 0.0;

            if (j > i) U[i][j] = M[i][j];
            else U[i][j] = 0.0;
        }
    }
}

bool invert_matrix(double M[5][5], double inv_M[5][5]) {
}

void add_matrices(double A[N][N], double B[N][N], double R[N][N]) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            R[i][j] = A[i][j] + B[i][j];
}


void print_matrix(double M[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            printf("%5.2f ", M[i][j]);
        std::cout << std::endl;
    }
}
