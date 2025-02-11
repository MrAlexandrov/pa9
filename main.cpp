// #define YANDEX

#if defined(YANDEX)
#include <junk/mrralexandrov/test_directory/pa9/include/matrix.hpp>
#include <junk/mrralexandrov/test_directory/pa9/include/plotter.hpp>
// #include <junk/mrralexandrov/test_directory/pa9/include/interface.hpp>
#else
#include "matrix.hpp"
#include "plotter.hpp"
// #include "interface.hpp"
#endif // YANDEX

#include <iostream>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
// Параметры элементов схемы

#define CHANGE_RESISTOR_TO_CAPACITOR
#define CHANGE_E1

using namespace NMatrix;

constexpr long double EPS = 1e-13;

constexpr long double R2 = 1000.0;
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
#if defined(CHANGE_RESISTOR_TO_CAPACITOR)
constexpr long double C_NEW = 1e-6;
#else
constexpr long double R21 = 1e+5;
#endif // CHANGE_RESISTOR_TO_CAPACITOR

long double Id(long double p3, long double p5) {
    return It * (std::exp((p3 - p5) / MFt) - 1);
}

long double DiffIdDiffp3(long double p3, long double p5) {
    return It * std::exp((p3 - p5) / MFt) / MFt;
}

long double DiffIdDiffp5(long double p3, long double p5) {
    return -It * std::exp((p3 - p5) / MFt) / MFt;
}

#if defined(CHANGE_E1)
long double NewI2(long double phi0, long double phi4) {
    return (15 + (phi0 - phi4)) / Re1;
}

long double NewI2DiffPhi0() {
    return 1 / Re1;
}

long double NewI2DiffPhi4() {
    return - 1 / Re1;
}
#else
long double I2(long double currentTime) {
    return (1e1 / Re1) * std::sin(2 * M_PI * currentTime / T);
}
#endif // CHANGE_E1

template<typename T>
struct TMaximum {
    T Value;
    int Row;
    int Col;
};

auto FindAbsMax(const TMatrix<>& matrix) -> TMaximum<long double> {
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
    const TMatrix<>& coef
    , const TMatrix<>& rightPart
) {
    if (coef.Rows() != coef.Cols() || coef.Rows() != rightPart.Rows() || rightPart.Cols() != 1) {
        throw std::runtime_error("Error in function Gauss(): size of coef-matrix must be [n]x[n] and size of rightPart-matrix mustbe [n]x[1].");
    }
    TMatrix<> a = coef;
    TMatrix<> b = rightPart;
    TMatrix<> result(b.Rows(), 1);
    std::vector<int> p(a.Rows());                       // вектор перестановок столбцов матрицы коэффициентов и
                                                          // соответствующихим строк решения

    // Прямой ход метода Гаусса
    for (int i = 0; i < static_cast<int>(a.Cols()) - 1; ++i) {
        auto [max, rowMax, colMax] = 
            FindAbsMax(GetSubMatrix(a, i, a.Rows(), i,  a.Cols()));
        rowMax += i;
        colMax += i;

        if (i != rowMax) {  // Перестановка строк
            std::swap(a[i], a[rowMax]);
            std::swap(b[i], b[rowMax]);
        }
        if (i != colMax) {  // Перестановка столбцов
            for (int k = 0; k < static_cast<int>(a.Rows()); ++k) {
                std::swap(a[k][i], a[k][colMax]);
                p[i] = colMax;
            }
        }
        if (std::fabs(a[i][i]) < EPS) {
            throw std::runtime_error("Error in fuction Gauss(): coef-matrix must not be degenerate.");
        }
        for (int j = i + 1; j < static_cast<int>(a.Rows()); ++j) {
            long double m = a[j][i] / a[i][i];
            for (int k = i; k < static_cast<int>(a.Cols()); ++k) {
                a[j][k] -= m * a[i][k];
            }
            b[j][0] -= m * b[i][0];
        }
    }
    // Обратный ход метода Гаусса
    if (std::fabs(a[a.Rows() - 1][a.Cols() - 1]) < EPS) {
        throw std::runtime_error("Error in fuction Gauss(): coef-matrix must not be degenerate.");
    }
    for (int i = b.Rows() - 1; i >= 0; --i) {
        long double sum = 0.0;
        for (int j = i + 1; j < static_cast<int>(a.Rows()); ++j) {
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
    const int timeIteration
    , const long double currentTime
    , const long double dt
    , TMatrix<>& nodeAdmittance                 // матрица узловых проводимостей
    , TMatrix<>& residualVector                 // вектор невязок
    , const TMatrix<>& basis
    , const std::vector<long double>& uc1
    , const std::vector<long double>& uc2
    , const std::vector<long double>& ucb
    , const std::vector<long double>& il2
    #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
    , const std::vector<long double>& u_new
    #endif // CHANGE_RESISTOR_TO_CAPACITOR
) {
    residualVector = {
        {
            (basis[0][0] / R2) - 0                                  // resistor R2                  0, ground
            + C1 * (basis[0][0] - basis[1][0] - uc1.back()) / dt        // capasitor C1                 0, 1
            #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
            + C_NEW * (basis[0][0] - basis[4][0] - u_new.back()) / dt   // capasitor C_NEW              0, 4
            #else
            + (basis[0][0] - basis[4][0]) / R21
            #endif // CHANGE_RESISTOR_TO_CAPACITOR
            #if defined(CHANGE_E1)
            + NewI2(basis[0][0], basis[4][0])
            #else
            + I2(currentTime)                                           // inductor I2 from EDS         0, 3
            #endif // CHANGE_E1
            + (basis[0][0] - basis[3][0]) / Re1                         // resistor Re1 from EDS        0, 3
        },
        {
            - C1 * (basis[0][0] - basis[1][0] - uc1.back()) / dt    // capasitor C1                 1, 0   
            + (basis[1][0] - basis[2][0]) / Rb                          // resistor Rb from diode       1, 2
        },
        {
            - (basis[1][0] - basis[2][0]) / Rb                      // resistor Rb from diode       2, 1
            + Cb * (basis[2][0] - basis[4][0] - ucb.back()) / dt        // capasitor Cb from diode      2, 4
            + (basis[2][0] - basis[4][0]) / Ru                          // resistor Ru from diode       2, 4
            + Id(basis[2][0], basis[4][0])                       // inductor Id from diode       2, 4
        },
        {
            #if defined(CHANGE_E1)
            - NewI2(basis[0][0], basis[4][0])
            #else
            - I2(currentTime)                                           // inductor I2 from EDS         3, 0
            #endif // CHANGE_E1
            - (basis[0][0] - basis[3][0]) / Re1                         // resistor Re1 from EDS        3, 0
            + (il2.back() + dt * (basis[3][0] - basis[4][0]) / L2)      // inductor L2                  3, 4
        },
        {
            - Cb * (basis[2][0] - basis[4][0] - ucb.back()) / dt    // capasitor Cb from diode      4, 2
            - (basis[2][0] - basis[4][0]) / Ru                          // resistor Ru from diode       4, 2
            - Id(basis[2][0], basis[4][0])                       // inductor Id from diode       4, 2
            #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
            - C_NEW * (basis[0][0] - basis[4][0] - u_new.back()) / dt   // capasitor C_NEW              4, 0
            #else
            - (basis[0][0] - basis[4][0]) / R21
            #endif // CHANGE_RESISTOR_TO_CAPACITOR
            - (il2.back() + dt * (basis[3][0] - basis[4][0]) / L2)      // inductor L2                  4, 3
            + (C2 * (basis[4][0] - uc2.back()) / dt) - 0                // capasitor C2                 4, ground
        }
    };

    nodeAdmittance = {
        {
            #if defined(CHANGE_E1)
            + NewI2DiffPhi0()   // added
            #endif // CHANGE_E1
            + 1 / R2                                                    // resistor R2
            + C1 / dt                                                   // capasitor C1
            #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
            + (C_NEW / dt)                                              // capasitor C_NEW
            #else
            + 1 / R21                                                   // resistor R21
            #endif // CHANGE_RESISTOR_TO_CAPACITOR
            + 1 / Re1,                                                  // resistor Re1
            - C1 / dt,                                              // capasitor C1
            0,
            - 1 / Re1,                                              // resistor Re1
            #if defined(CHANGE_E1)
            + NewI2DiffPhi4()   // added
            #endif // CHANGE_E1
            #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
            - (C_NEW / dt)                                              // capasitor C_NEW
            #else
            - 1 / R21                                                   // resistor R21
            #endif // CHANGE_RESISTOR_TO_CAPACITOR
        },
        {
            - C1 / dt,                                              // capasitor C1
            + C1 / dt                                               // capasitor C1
            + 1 / Rb,                                                   // resistor Rb from diode
            - 1 / Rb,                                               // resistor Rb from diode
            0,
            0
        },
        {                                                           
            0,                                                      
            - 1 / Rb,                                               // resistor Rb from diode
            1 / Rb                                                  // resistor Rb from diode
            + Cb / dt                                                   // capasitor Cb from diode
            + 1 / Ru                                                    // resistor Ru from diode
            + DiffIdDiffp3(basis[2][0], basis[4][0]),            // inductor Id from diode
            0,
            - Cb / dt                                               // capasitor Cb from diode
            - 1 / Ru                                                    // resistor Ru from diode
            + DiffIdDiffp5(basis[2][0], basis[4][0])             // inductor Id from diode
        },
        {
            #if defined(CHANGE_E1)
            - NewI2DiffPhi0()   // added
            #endif // CHANGE_E1
            - 1 / Re1,                                                  // resistor Re1
            0,
            0,
            1 / Re1                                                 // resistor Re1
            + dt / L2,                                                  // inductor L2
            #if defined(CHANGE_E1)
            - NewI2DiffPhi4()   // added
            #endif // CHANGE_E1
            - dt / L2                                                   // inductor L2
        },
        {
            #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
            - (C_NEW / dt),                                         // capasitor C_NEW
            #else
            - 1 / R21,                                                  // resistor R21
            #endif // CHANGE_RESISTOR_TO_CAPACITOR
            0,
            - Cb / dt                                               // capasitor Cb from diode
            - 1 / Ru                                                    // resistor Ru from diode
            - DiffIdDiffp3(basis[2][0], basis[4][0]),            // inductor Id from diode
            dt / L2,                                                // inductor L2
            Cb / dt                                                 // capasitor Cb from diode
            + 1 / Ru                                                    // resistor Ru from diode
            - DiffIdDiffp5(basis[2][0], basis[4][0])             // inductor Id from diode
            + C2 / dt                                                   // capasitor C2
            + dt / L2                                                   // inductor L2
            #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
            + (C_NEW / dt)                                              // capasitor C2
            #else
            + (1 / R21)                                                 // capasitor C2
            #endif // CHANGE_RESISTOR_TO_CAPACITOR
        }
    };
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& other) {
    for (const auto& item : other) {
        out << item << " ";
    }
    return out;
}

bool PerformNewtonIteration(
    int timeIteration
    , long double currentTime
    , long double dt
    , NMatrix::TMatrix<>& basis
    , const std::vector<long double>& uc1
    , const std::vector<long double>& uc2
    , const std::vector<long double>& ucb
    , const std::vector<long double>& il2
    #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
    , const std::vector<long double>& u_new
    #endif // CHANGE_RESISTOR_TO_CAPACITOR
    , int& currentIteration
) {
    constexpr int MAX_STEPS = 10;   // максимальное число итераций Ньютона
    constexpr long double eps = 1e-9;

    NMatrix::TMatrix<> nodeAdmittance(basis.Rows(), basis.Rows());
    NMatrix::TMatrix<> residualVector(basis.Rows(), 1);

    // Начнём с "большой" поправки, чтобы зайти в цикл
    NMatrix::TMatrix<> delta(basis.Rows(), 1, 1e-6);

    int n = 0; // счётчик итераций Ньютона

    // Выполняем итерационный процесс
    while (std::fabs(FindAbsMax(delta).Value) > eps && n < MAX_STEPS) {
        Initialize(
            timeIteration
            , currentTime
            , dt
            , nodeAdmittance
            , residualVector
            , basis
            , uc1
            , uc2
            , ucb
            , il2
            #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
            , u_new
            #endif // CHANGE_RESISTOR_TO_CAPACITOR
        );
        // std::cout << nodeAdmittance << std::endl;
        delta = Gauss(nodeAdmittance, -residualVector);
        // std::cout << nodeAdmittance << std::endl;
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

int main(int argc, char* argv[]) {
    constexpr long double DT_MIN = 1e-12;                               // минимальный шаг интегрирования по времени
    constexpr long double EPS_MIN = 1e-6;                               // нижняя граница для оценки локальной точности
    constexpr long double EPS_MAX = 5e-2;                               // верхняя граница для оценки локальной точности
    constexpr long double TIME_MAX = 1e-3;                              // время расчёта
    long double currentTime = 0;                                        // время
    long double dt = 1e-8;                                              // шаг интегрирования по времени
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
    #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
    std::vector<long double> u_new = {0};
    #endif // CHANGE_RESISTOR_TO_CAPACITOR
    std::vector<std::vector<long double>> phi(5);
    int timeIteration = 1;
    int currentIteration = 1;
    while (dt >= DT_MIN && currentTime <= TIME_MAX) {
        // int n = 0;
        basis = ((previousBasis[0] - previousBasis[1]) * (dt_prev1 + dt) / dt) + previousBasis[1]; // начальные приближения
        bool success = PerformNewtonIteration(
            timeIteration
            , currentTime
            , dt
            , basis
            , uc1
            , uc2
            , ucb
            , il2
            #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
            , u_new
            #endif // CHANGE_RESISTOR_TO_CAPACITOR
            , currentIteration
        );

        if (!success) {
            dt /= 2;
            continue;
        }

        if (timeIteration > 2) {  // оценка локальной точности
            const auto diff1 = (previousBasis[0] - previousBasis[1]) * (1 / (dt_prev1 * dt_prev1));
            const auto diff2 = (previousBasis[1] - previousBasis[2]) * (1 / (dt_prev1 * dt_prev2));
            long double d = 0.5 * dt * dt * std::fabs(FindAbsMax(diff1 - diff2).Value);
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
        #if defined(CHANGE_RESISTOR_TO_CAPACITOR)
        u_new.push_back(basis[0][0] - basis[4][0]);
        #endif // CHANGE_RESISTOR_TO_CAPACITOR
        ++timeIteration;
    }
    if (dt < DT_MIN) {
        std::cerr << "dt < DT_MIN at " << timeIteration << " time iteration!\n";
    }
    std::cout << "Итераций по времени: " << timeIteration << std::endl;
    std::cout << "Всего итераций: " << currentIteration << std::endl;

    std::vector<std::string> colors = {"red", "orange", "green", "blue", "violet"};

    for (int i = 0; i < 5; ++i) {
        std::string title = "phi_" + std::to_string(i + 1);
        NPlotter::TPlotter<decltype(time)::value_type> graphic(title);
        graphic.SetXValues(time);
        graphic.EmplaceGraphic<std::decay_t<decltype(phi[i])>::value_type>(title, phi[i], colors[i]);
    }
    return 0;
}
