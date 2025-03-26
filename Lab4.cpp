/*
* Napisz program w języku „C/C++”, realizujący uogólnioną metodę Newtona rozwiązywania układu
* trzech algebraicznych równań nieliniowych, i zastosuj ten program do przykładu z zadania 1.
* Przyjmij takie przybliżenie początkowe, aby uzyskać zbieżność metody. Zastosuj trzy niezależne
* kryteria zakończenia iteracji. Zadbaj o to, aby wyprowadzać na konsolę wyniki pośrednie obliczeń
* dla każdej iteracji, tak aby możliwe było obserwowanie zbieżności kolejnych przybliżeń
* pierwiastków. W szczególności oblicz jak zmienia się estymator błędu rozwiązania oraz residuum
* układu w trakcie iteracji.
*/

#include <iostream>
#include <cmath>
#include <functional>
#include <cstdbool>

#define N_MAX 100
#define TOLX 1e-12
#define TOLF 1e-12

double f1(double x, double y, double z) { return x * x + y * y + z * z - 4.0; }
double f2(double x, double y, double z) { return x * x + y * y / 2.0 - 1.0; }
double f3(double x, double y, double z) { return x * y - 0.5; }
double df1_dx(double x, double y, double z) { return 2.0 * x; }
double df1_dy(double x, double y, double z) { return 2.0 * y; }
double df1_dz(double x, double y, double z) { return 2.0 * z; }
double df2_dx(double x, double y, double z) { return 2 * x; }
double df2_dy(double x, double y, double z) { return y; }
double df2_dz(double x, double y, double z) { return 0; }
double df3_dx(double x, double y, double z) { return y; }
double df3_dy(double x, double y, double z) { return x; }
double df3_dz(double x, double y, double z) { return 0; }

bool invert_matrix(double [3][3], double [3][3]);

void newton(double, double, double);

int main() {
    newton(5.0, 12.0, -3.0);
    return 0;
}

void newton(double x, double y, double z) {
    std::cout << "----------------------------------- Newton -----------------------------------" << std::endl;
    std::printf("%-4s%-15s%-15s%-15s%-15s%-15s\n", "i", "Wartosc x", "Wartosc y", "Wartosc z", "Estymator", "Residuum");

    for (int i = 1; i <= N_MAX; i++) {
        double J[3][3] = {
            {df1_dx(x, y, z), df1_dy(x, y, z), df1_dz(x, y, z)},
            {df2_dx(x, y, z), df2_dy(x, y, z), df2_dz(x, y, z)},
            {df3_dx(x, y, z), df3_dy(x, y, z), df3_dz(x, y, z)}
        };
        double inv_J[3][3];

        if (!invert_matrix(J, inv_J)) {
            std::cout << "MACIERZ JACOBIEGO JEST OSOBLIWA!" << std::endl;
            return;
        }

        double F[3] = {f1(x, y, z), f2(x, y, z), f3(x, y, z)};
        double delta[3] = {
            inv_J[0][0] * F[0] + inv_J[0][1] * F[1] + inv_J[0][2] * F[2],
            inv_J[1][0] * F[0] + inv_J[1][1] * F[1] + inv_J[1][2] * F[2],
            inv_J[2][0] * F[0] + inv_J[2][1] * F[1] + inv_J[2][2] * F[2]
        };

        x -= delta[0];
        y -= delta[1];
        z -= delta[2];

        double estimator = std::sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
        double residuum = std::sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);

        std::printf("%-3d %-10.12f %-10.12f %-10.12f %-10.12f %-10.12f\n", i, x, y, z, estimator, residuum);

        if (estimator < TOLX || residuum < TOLF) break;
    }
    std::cout << std::endl;
}

bool invert_matrix(double M[3][3], double inv_M[3][3]) {
    double det = M[0][0] * M[1][1] * M[2][2] + M[0][1] * M[1][2] * M[2][0] + M[0][2] * M[1][0] * M[2][1] -
                 M[0][2] * M[1][1] * M[2][0] - M[0][1] * M[1][0] * M[2][2] - M[0][0] * M[1][2] * M[2][1];

    if (det == 0) return false;
    double inv_det = 1.0 / det;

    inv_M[0][0] = inv_det * (M[1][1] * M[2][2] - M[1][2] * M[2][1]);
    inv_M[0][1] = inv_det * (M[0][2] * M[2][1] - M[0][1] * M[2][2]);
    inv_M[0][2] = inv_det * (M[0][1] * M[1][2] - M[0][2] * M[1][1]);
    inv_M[1][0] = inv_det * (M[1][2] * M[2][0] - M[1][0] * M[2][2]);
    inv_M[1][1] = inv_det * (M[0][0] * M[2][2] - M[0][2] * M[2][0]);
    inv_M[1][2] = inv_det * (M[0][2] * M[1][0] - M[0][0] * M[1][2]);
    inv_M[2][0] = inv_det * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    inv_M[2][1] = inv_det * (M[0][1] * M[2][0] - M[0][0] * M[2][1]);
    inv_M[2][2] = inv_det * (M[0][0] * M[1][1] - M[0][1] * M[1][0]);

    return true;
}
