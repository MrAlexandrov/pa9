#include "matrix.hpp"
#include "plotter.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <string>
// Параметры элементов схемы

using namespace NMatrix;
using namespace NPlotter;

constexpr long double EPS = 1e-13;

constexpr long double R2 = 1000.0;
constexpr long double R21 = 1e+5;
constexpr long double C1 = 1e-6;
constexpr long double C2 = 1e-6;
constexpr long double L2 = 1e-3;
constexpr long double Cb = 2e-12;
constexpr long double Ru = 1e+6;
constexpr long double Rb = 20;
constexpr long double MFt = 0.026;
constexpr long double It = 1e-12;
constexpr long double Re2 = 1e-6;
constexpr long double T = 1e-4;
constexpr long double A = 1e+9;

long double I2(long double t) {
    return (10 / Re2) * std::sin(2 * M_PI * t / T);
}

long double Id(long double p3, long double p5) {
    return It * (std::exp((p3 - p5) / MFt) - 1);
}

long double dId_dp3(long double p3, long double p5) {
    return It * std::exp((p3 - p5) / MFt) / MFt;
}

long double dId_dp5(long double p3, long double p5) {
    return -It * std::exp((p3 - p5) / MFt) / MFt;
}

template<typename T>
struct TMaximum {
    T value;
    int row;
    int col;
};

auto FindAbsMax(const TMatrix<>& matrix) {
    int resutlRow = 0;
    int resutlCol = 0;
    long double maximum = matrix[0][0];
    for (int row = 0, endRow = matrix.Rows(); row < endRow; ++row) {
        for (int col = 0, endCol = matrix.Cols(); col < endCol; ++col) {
            auto value = std::fabs(matrix[row][col]);
            if (maximum < value) {
                resutlRow = row;
                resutlCol = col;
                maximum = value;
            }
        }
    }
    return TMaximum{maximum, resutlRow, resutlCol};
};

// Метод Гаусса c выбором главного элемента
TMatrix<> Gauss(
    const TMatrix<>& coef,
    const TMatrix<>& rightPart
    ) 
{
    if (coef.Rows() != coef.Cols() || coef.Rows() != rightPart.Rows() || rightPart.Cols() != 1) {
        throw std::runtime_error("Error in function Gauss(): size of coef-matrix must be [n]x[n] and size of rightPart-matrix mustbe [n]x[1].");
    }
    TMatrix<> a = coef;
    TMatrix<> b = rightPart;
    TMatrix<> result(b.Rows(), 1);
    std::vector<int> p(a.Rows());  // вектор перестановок столбцов матрицы коэффициентов и
                                  // соответствующихим строк решения

    // Прямой ход метода Гаусса
    for (int i = 0; i < a.Cols() - 1; i++) {
        auto [max, rowMax, colMax] = 
            FindAbsMax(GetSubMatrix(a, i, a.Rows(), i,  a.Cols()));
        rowMax += i;
        colMax += i;

        if (i != rowMax) {  // Перестановка строк
            std::swap(a[i], a[rowMax]);
            std::swap(b[i], b[rowMax]);
        }
        if (i != colMax) {  // Перестановка столбцов
            for (int k = 0; k < a.Rows(); k++) {
                std::swap(a[k][i], a[k][colMax]);
                p[i] = colMax;
            }
        }
        if (std::abs(a[i][i]) < EPS) {
            throw std::runtime_error("Error in fuction Gauss(): coef-matrix must not be degenerate.");
        }
        for (int j = i + 1; j < a.Rows(); j++) {
            long double m = a[j][i] / a[i][i];
            for (int k = i; k < a.Cols(); k++) {
                a[j][k] -= m * a[i][k];
            }
            b[j][0] -= m * b[i][0];
        }
    }
    // Обратный ход метода Гаусса
    if (std::abs(a[a.Rows() - 1][a.Cols() - 1]) < EPS) {
        throw std::runtime_error("Error in fuction Gauss(): coef-matrix must not be degenerate.");
    }
    for (int i = b.Rows() - 1; i >= 0; --i) {
        long double sum = 0.0;
        for (int j = i + 1; j < a.Rows(); ++j) {
            sum += a[i][j] * result[j][0];
        }

        result[i][0] = (b[i][0] - sum) / a[i][i];
    }
    for (int i = p.size() - 1; i >= 0; --i) {  // обратная перестановка строк решения
        if (p[i] != 0) {
            std::swap(result[i][0], result[p[i]][0]);
        }
    }
    return result;
}

// Функция заполнения матрицы проводимости и вектора токов
void Initialize(
    const int timeIteration,
    const double t,
    const double dt,
    TMatrix<>& nodeAdmittance,
    TMatrix<>& residualVector,
    const TMatrix<>& basis,
    const std::vector<double>& uc1,
    const std::vector<double>& uc2,
    const std::vector<double>& ucb,
    const std::vector<double>& il2
    )
{
    residualVector[0][0] = basis[0][0] / R2 +
                    C1 * (basis[0][0] - basis[1][0] - uc1[timeIteration - 1]) / dt +
                    (basis[0][0] - basis[4][0]) / R21 + I2(t) - (basis[3][0] - basis[0][0]) / Re2;
    residualVector[1][0] = -C1 * (basis[0][0] - basis[1][0] - uc1[timeIteration - 1]) / dt +
                    (basis[1][0] - basis[2][0]) / Rb;
    residualVector[2][0] = -(basis[1][0] - basis[2][0]) / Rb +
                    Cb * (basis[2][0] - basis[4][0] - ucb[timeIteration - 1]) / dt +
                    (basis[2][0] - basis[4][0]) / Ru + Id(basis[2][0], basis[4][0]);
    residualVector[3][0] = -I2(t) + (basis[3][0] - basis[0][0]) / Re2 +
                    (il2[timeIteration - 1] + dt * (basis[3][0] - basis[4][0]) / L2);
    residualVector[4][0] = -Cb * (basis[2][0] - basis[4][0] - ucb[timeIteration - 1]) / dt -
                    (basis[2][0] - basis[4][0]) / Ru - Id(basis[2][0], basis[4][0]) -
                    (basis[0][0] - basis[4][0]) / R21 -
                    (il2[timeIteration - 1] + dt * (basis[3][0] - basis[4][0]) / L2) +
                    C2 * (basis[4][0] - uc2[timeIteration - 1]) / dt;
    nodeAdmittance[0][0] = 1 / R2 + C1 / dt + 1 / R21 + 1 / Re2;
    nodeAdmittance[0][1] = -C1 / dt;
    nodeAdmittance[0][3] = -1 / Re2;
    nodeAdmittance[0][4] = -1 / R21;
    nodeAdmittance[1][0] = -C1 / dt;
    nodeAdmittance[1][1] = C1 / dt + 1 / Rb;
    nodeAdmittance[1][2] = -1 / Rb;
    nodeAdmittance[2][1] = -1 / Rb;
    nodeAdmittance[2][2] = 1 / Rb + Cb / dt + 1 / Ru + dId_dp3(basis[2][0], basis[4][0]);
    nodeAdmittance[2][4] = -Cb / dt - 1 / Ru + dId_dp5(basis[2][0], basis[4][0]);
    nodeAdmittance[3][0] = -1 / Re2;
    nodeAdmittance[3][3] = 1 / Re2 + dt / L2;
    nodeAdmittance[3][4] = -dt / L2;
    nodeAdmittance[4][0] = -1 / R21;
    nodeAdmittance[4][2] = -Cb / dt - 1 / Ru - dId_dp3(basis[2][0], basis[4][0]);

    nodeAdmittance[4][3] = -dt / L2;
    nodeAdmittance[4][4] =
        Cb / dt + 1 / Ru - dId_dp5(basis[2][0], basis[4][0]) + C2 / dt + dt / L2 + 1 / R21;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& other) {
    for (const auto& item : other) {
        out << item << " ";
    }
    return out;
}

int main() {
    constexpr int MAX_STEPS = 10;                                       // максимальное число итераций метода Ньютона
    constexpr long double DT_MIN = 1e-12;                               // минимальный шаг интегрирования по времени
    constexpr long double eps = 1e-9;                                   // максимальное значение поправок
    constexpr long double eps_min = 1e-3;                               // нижняя граница для оценки локальной точности
    constexpr long double eps_max = 5e-2;                               // верхняя граница для оценки локальной точности
    constexpr long double MAX_TIME = 1e-3;                              // время расчёта
    long double t = 0;                                                  // время
    long double dt = DT_MIN;                                            // шаг интегрирования по времени
    long double dt_prev1 = dt;
    long double dt_prev2 = dt;                                          // предыдущие шаги интегрирования по времени
    TMatrix<> basis(5, 1);                                  // вектор базиса метода (узловые потенциалы)

    std::vector<TMatrix<>> previousBasis(3, TMatrix<>(basis.Rows(), 1)); // предыдущие значения базиса

    TMatrix<> nodeAdmittance(basis.Rows(), basis.Rows());      // матрица узловых проводимостей
    TMatrix<> residualVector(basis.Rows(), 1);                     // вектор невязок
    std::vector<long double> time;
    // Переменные состояния
    std::vector<double> uc1 = {0};
    std::vector<double> uc2 = {0};
    std::vector<double> ucb = {0};
    std::vector<double> il2 = {0};
    std::vector<std::vector<long double>> phi(5);
    int timeIteration = 1;
    int currentIteration = 1;
    while (dt >= DT_MIN && t <= MAX_TIME) {
        TMatrix<> delta(basis.Rows(), 1, eps_max + 1);  // вектор поправок
        int n = 0;
        basis = ((previousBasis[0] - previousBasis[1]) * (dt_prev1 + dt) / dt) + previousBasis[1]; // начальные приближения
        while (std::abs(FindAbsMax(delta).value) > eps && n < MAX_STEPS) {  // метод Ньютона
            Initialize(timeIteration, t, dt, nodeAdmittance, residualVector, basis, uc1, uc2, ucb, il2);
            delta = Gauss(nodeAdmittance, -residualVector);
            basis += delta;
            n++;
            currentIteration++;
        }
        if (n > MAX_STEPS) {
            dt /= 2;
            continue;
        }
        if (timeIteration > 2) {  // оценка локальной точности
            long double d = 0.5 * dt * dt * std::abs(FindAbsMax(((previousBasis[0] - previousBasis[1]) * (1 / (dt_prev1 * dt_prev1)) - (previousBasis[1] - previousBasis[2]) * (1 / (dt_prev1 * dt_prev2)))).value);
            if (d < eps_min) {
                t += dt;
                dt_prev2 = dt_prev1;
                dt_prev1 = dt;
                dt *= 2;
            } else if (d < eps_max) {
                t += dt;
                dt_prev2 = dt_prev1;
                dt_prev1 = dt;
            } else {
                dt /= 2;
                continue;
            }
        } else {
            t += dt;
            dt_prev2 = dt_prev1;
            dt_prev1 = dt;
        }
        previousBasis[2] = previousBasis[1];
        previousBasis[1] = previousBasis[0];
        previousBasis[0] = basis;
        time.push_back(t);
        for (int i = 0; i < 5; ++i) {
            phi[i].push_back(basis[i][0]);
        }
        uc1.push_back(basis[0][0] - basis[1][0]);
        uc2.push_back(basis[4][0]);
        ucb.push_back(basis[2][0] - basis[4][0]);
        il2.push_back(il2[timeIteration - 1] + dt * (basis[3][0] - basis[4][0]) / L2);
        timeIteration++;
    }
    if (dt < DT_MIN) {
        std::cerr << "dt < DT_MIN at " << timeIteration << " time iteration!\n";
    }
    std::cout << "Итераций по времени: " << timeIteration << std::endl;
    std::cout << "Всего итераций: " << currentIteration << std::endl;
    for (int i = 0; i < 5; ++i) {
        std::string title = "phi_" + std::to_string(i + 1);
        TPlotter graphic(title);
        graphic.SetXValues(time);
        graphic.AddGraphic(title, phi[i]);
    }
    return 0;
}
