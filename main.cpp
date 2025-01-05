#include "matrix.hpp"
#include "plotter.hpp"
#include "interfase.hpp"

#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <string>
// Параметры элементов схемы

using namespace NMatrix;

// constexpr long double EPS = 1e-13;

constexpr long double R2 = 1000.0;
// constexpr long double R21 = 1e+5;
constexpr long double C1 = 1e-6;
constexpr long double C2 = 1e-6;
constexpr long double L2 = 1e-3;
constexpr long double Cb = 2e-12;
constexpr long double Ru = 1e+6;
constexpr long double Rb = 20;
constexpr long double MFt = 0.026;
constexpr long double It = 1e-12;
constexpr long double Re1 = 1e-6;
constexpr long double T = 1e-4;
constexpr long double C_NEW = 1e-6;

long double I2(long double currentTime, long double p0, long double p4) {
    // return (1e1 / Re1) * std::sin(2 * M_PI * currentTime / T);
    return (p0 - p4 + 15) / Re1;
}

long double dI2_dp0() {
    return 1.0 / Re1;
}

long double dI2_dp4() {
    return -1.0 / Re1;
}

long double Id(long double p2, long double p4) {
    return It * (std::exp((p2 - p4) / MFt) - 1);
}

long double dId_dp2(long double p2, long double p4) {
    return It * std::exp((p2 - p4) / MFt) / MFt;
}

long double dId_dp4(long double p2, long double p4) {
    return -It * std::exp((p2 - p4) / MFt) / MFt;
}

// Функция заполнения матрицы проводимости и вектора токов
void Initialize(
    const int timeIteration,
    const long double currentTime,
    const long double dt,
    TMatrix<>& nodeAdmittance,                              // матрица узловых проводимостей
    TMatrix<>& residualVector,                              // вектор невязок
    const TMatrix<>& basis,
    const std::vector<long double>& uc1,
    const std::vector<long double>& uc2,
    const std::vector<long double>& ucb,
    const std::vector<long double>& il2,
    const std::vector<long double>& u_new
) {
    auto& p0 = basis[0][0];
    auto& p1 = basis[1][0];
    auto& p2 = basis[2][0];
    auto& p3 = basis[3][0];
    auto& p4 = basis[4][0];
    residualVector = {
        {
            basis[0][0] / R2
            + C1 * (basis[0][0] - basis[1][0] - uc1.back()) / dt
            + C_NEW * (basis[0][0] - basis[4][0] - u_new.back()) / dt
            + I2(currentTime, basis[0][0], basis[4][0])
            - (basis[3][0] - basis[0][0]) / Re1
        },
        {
            -C1 * (basis[0][0] - basis[1][0] - uc1.back()) / dt
            + (basis[1][0] - basis[2][0]) / Rb
        },
        {
            -(basis[1][0] - basis[2][0]) / Rb
            + Cb * (basis[2][0] - basis[4][0] - ucb.back()) / dt
            + (basis[2][0] - basis[4][0]) / Ru
            + Id(basis[2][0], basis[4][0])
        },
        {
            -I2(currentTime, basis[0][0], basis[4][0])
            + (basis[3][0] - basis[0][0]) / Re1
            + (il2.back() + dt * (basis[3][0] - basis[4][0]) / L2)
        },
        {
            -Cb * (basis[2][0] - basis[4][0] - ucb.back()) / dt
            - (basis[2][0] - basis[4][0]) / Ru
            - Id(basis[2][0], basis[4][0])
            - C_NEW * (basis[0][0] - basis[4][0] - u_new.back()) / dt
            - (il2.back() + dt * (basis[3][0] - basis[4][0]) / L2)
            + C2 * (basis[4][0] - uc2.back()) / dt
        }
    };

    nodeAdmittance = {
        {
            1 / R2 + C1 / dt + (C_NEW / dt) + 1 / Re1,
            -C1 / dt,
            0,
            - 1 / Re1,
            -(C_NEW / dt)
        },
        {
            -C1 / dt,
            C1 / dt + 1 / Rb,
            -1 / Rb,
            0,
            0
        },
        {
            0,
            -1.0 / Rb,
            1.0 / Rb + Cb / dt + 1 / Ru - dId_dp2(basis[2][0], basis[4][0]),
            0,
            -Cb / dt - 1.0 / Ru + dId_dp4(basis[2][0], basis[4][0])
        },
        {
            -1 / Re1,
            0,
            0,
            dt / L2,
            -dt / L2 + 1 / Re1
        },
        {
            -(C_NEW / dt),
            0,
            -Cb / dt - 1.0 / Ru - dId_dp2(basis[2][0], basis[4][0]),
            -dt / L2,
            Cb / dt + 1.0 / Ru - dId_dp4(basis[2][0], basis[4][0]) + C2 / dt + dt / L2 + (C_NEW / dt)
        }
    };
    // std::cout << Cb / dt << std::endl;
    // std::cout << 1 / Ru << std::endl;
    // std::cout << nodeAdmittance[2][2] << std::endl;
    // std::cout << nodeAdmittance[2][4] << std::endl;
    // std::cout << nodeAdmittance[4][2] << std::endl;
    // std::cout << nodeAdmittance[4][4] << std::endl;
    // std::cout << 1 / Rb + Cb / dt + 1 / Ru - dId_dp2(basis[2][0], basis[4][0]) << std::endl;
    // std::cout << -Cb / dt - 1 / Ru + dId_dp4(basis[2][0], basis[4][0]) << std::endl;
    // std::cout << -Cb / dt - 1 / Ru - dId_dp2(basis[2][0], basis[4][0]) << std::endl;
    // std::cout << Cb / dt + 1 / Ru - dId_dp4(basis[2][0], basis[4][0]) + C2 / dt + dt / L2 + (C_NEW / dt) << std::endl;
    
    // std::cout << 1.0 / Rb << std::endl;
    // std::cout << Cb / dt << std::endl;
    // std::cout << 1 / Ru << std::endl;
    // std::cout << dId_dp2(basis[2][0], basis[4][0]) << std::endl;
    
    // std::cout << It << std::endl;
    // std::cout << p2 << std::endl;
    // std::cout << p4 << std::endl;
    // std::cout << "MFt: " << MFt << std::endl;
    // std::cout << "(p2 - p4): " << (p2 - p4) << std::endl;
    // std::cout << "(p2 - p4) / MFt: " << (p2 - p4) / MFt << std::endl;
    // std::cout << "std::exp((p2 - p4) / MFt): " << std::exp((p2 - p4) / MFt) << std::endl;
    // std::cout << "It * std::exp((p2 - p4) / MFt) / MFt: " << It * std::exp((p2 - p4) / MFt) / MFt << std::endl;
    // std::cout << std::endl;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& other) {
    for (const auto& item : other) {
        out << item << " ";
    }
    return out;
}

bool PerformNewtonIteration(
    int timeIteration,
    long double currentTime,
    long double dt,
    NMatrix::TMatrix<>& basis,
    const std::vector<long double>& uc1,
    const std::vector<long double>& uc2,
    const std::vector<long double>& ucb,
    const std::vector<long double>& il2,
    const std::vector<long double>& u_new,
    int& currentIteration
) {
    constexpr int MAX_STEPS = 10;   // максимальное число итераций Ньютона
    constexpr long double eps = 1e-9;

    NMatrix::TMatrix<> nodeAdmittance(basis.Rows(), basis.Rows());
    NMatrix::TMatrix<> residualVector(basis.Rows(), 1);

    // Начнём с "большой" поправки, чтобы зайти в цикл
    NMatrix::TMatrix<> delta(basis.Rows(), 1, 1e-6);

    int n = 0; // счётчик итераций Ньютона

    // Выполняем итерационный процесс
    while (std::fabs(FindAbsMax(delta).value) > eps && n < MAX_STEPS) {
        Initialize(timeIteration, currentTime, dt, nodeAdmittance, residualVector, basis, uc1, uc2, ucb, il2, u_new);
        std::cout << std::setprecision(10) << nodeAdmittance << std::endl;

        delta = Gauss(nodeAdmittance, -residualVector);
        basis += delta;
        ++n;
        ++currentIteration;
    }

    // Проверяем результат
    if (n > MAX_STEPS) {
        // Не удалось достичь нужной точности за отведённое число итераций
        return false;
    }
    return true;
}

int main() {
    constexpr long double DT_MIN = 1e-12;                               // минимальный шаг интегрирования по времени
    constexpr long double EPS_MIN = 1e-8;                               // нижняя граница для оценки локальной точности
    constexpr long double EPS_MAX = 1e-2;                               // верхняя граница для оценки локальной точности
    constexpr long double TIME_MAX = 2e-5;                              // время расчёта
    long double currentTime = 0;                                        // время
    long double dt = 1e-8;                                            // шаг интегрирования по времени
    long double dt_prev1 = dt;
    long double dt_prev2 = dt;                                          // предыдущие шаги интегрирования по времени
    TMatrix<> basis(5, 1);                                   // вектор базиса метода (узловые потенциалы)

    std::vector<TMatrix<>> previousBasis(3, basis); // предыдущие значения базиса

    std::vector<long double> time;
    // Переменные состояния
    std::vector<long double> uc1 = {0};
    std::vector<long double> uc2 = {0};
    std::vector<long double> ucb = {0};
    std::vector<long double> il2 = {0};
    std::vector<long double> u_new = {0};
    std::vector<std::vector<long double>> phi(5);
    int timeIteration = 1;
    int currentIteration = 1;
    while (dt >= DT_MIN && currentTime <= TIME_MAX) {
        int n = 0;
        basis = ((previousBasis[0] - previousBasis[1]) * (dt_prev1 + dt) / dt) + previousBasis[1]; // начальные приближения
        bool success = PerformNewtonIteration(timeIteration, currentTime, dt, basis, uc1, uc2, ucb, il2, u_new, currentIteration);

        if (!success) {
            dt /= 2;
            continue;
        }

        if (timeIteration > 2) {  // оценка локальной точности
            long double d = 0.5 * dt * dt * std::fabs(FindAbsMax(((previousBasis[0] - previousBasis[1]) * (1 / (dt_prev1 * dt_prev1)) - (previousBasis[1] - previousBasis[2]) * (1 / (dt_prev1 * dt_prev2)))).value);
            if (d < EPS_MIN) {
                currentTime += dt;
                dt_prev2 = dt_prev1;
                dt_prev1 = dt;
                dt *= 2;
            } else if (d < EPS_MAX) {
                currentTime += dt;
                dt_prev2 = dt_prev1;
                dt_prev1 = dt;
            } else {
                dt /= 2;
                continue;
            }
        } else {
            currentTime += dt;
            dt_prev2 = dt_prev1;
            dt_prev1 = dt;
        }
        previousBasis[2] = previousBasis[1];
        previousBasis[1] = previousBasis[0];
        previousBasis[0] = basis;
        time.push_back(currentTime);
        for (int i = 0; i < 5; ++i) {
            phi[i].push_back(basis[i][0]);
        }
        uc1.push_back(basis[0][0] - basis[1][0]);
        uc2.push_back(basis[4][0]);
        ucb.push_back(basis[2][0] - basis[4][0]);
        il2.push_back(il2.back() + dt * (basis[3][0] - basis[4][0]) / L2);
        u_new.push_back(basis[0][0] - basis[4][0]);
        ++timeIteration;
    }
    if (dt < DT_MIN) {
        // std::cerr << "dt < DT_MIN at " << timeIteration << " time iteration!\n";
    }
    // std::cout << "Итераций по времени: " << timeIteration << std::endl;
    // std::cout << "Всего итераций: " << currentIteration << std::endl;

    std::vector<std::string> colors = {"red", "orange", "green", "blue", "violet"};

    {
        for (int i = 0; i < 5; ++i) {
            std::string title = "phi_" + std::to_string(i + 1);
            NPlotter::TPlotter graphic(title);
            graphic.SetXValues(time);
            graphic.AddGraphic(title, colors[i], phi[i]);
        }
    }
    // const int derivativeSize = time.size() - 1;
    // std::vector<long double> derivativeTime;
    // std::vector<std::vector<long double>> derivativePhi(5);

    // for (int i = 0; i < derivativeSize; ++i) {
    //     derivativeTime.emplace_back(time[i] + (time[i + 1] - time[i]) / 2);
    //     for (int j = 0; j < 5; ++j) {
    //         derivativePhi[j].emplace_back((phi[j][i + 1] - phi[j][i]) / (time[i + 1] - time[i]));
    //     }
    // }
    // {
    //     int i = 0;
    //     NPlotter::TPlotter graphic("derivative" + std::to_string(i + 1));
    //     for (int i = 0; i < 5; ++i) {
    //         std::string title = "derivative_phi_" + std::to_string(i + 1);
    //         graphic.SetXValues(derivativeTime);
    //         graphic.AddGraphic(title, colors[i], derivativePhi[i]);
    //     }
    // }
    return 0;
}
