#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <string>
// Параметры элементов схемы

const double R2 = 1000.0;
const double R21 = 1e+5;
const double C1 = 1e-6;
const double C2 = 1e-6;
const double L2 = 1e-3;
const double Cb = 2e-12;
const double Ru = 1e+6;
const double Rb = 20;
const double MFt = 0.026;
const double It = 1e-12;
const double Re2 = 1e-6;
const double T = 1e-4;
const double A = 1e+9;

double I2(double t) {
    return (10 / Re2) * std::sin(2 * M_PI * t / T);
}
double Id(double p3, double p5) {
    double result = It * (std::exp((p3 - p5) / MFt) - 1);
    return result;
};
double dId_dp3(double p3, double p5) {
    double result = It * std::exp((p3 - p5) / MFt) / MFt;
    return result;
};
double dId_dp5(double p3, double p5) {
    double result = -It * std::exp((p3 - p5) / MFt) / MFt;
    return result;
};

// Класс матриц
class Matrix {
private:
    std::vector<std::vector<double>> m = {{0.}};
    int nrow = 1;
    int ncol = 1;
    struct Max {
        double value;
        int row;
        int col;
    };

public:
    Matrix() = default;
    Matrix(int row, int col, double d = 0) : nrow{row}, ncol{col} {
        m.resize(row);
        for (int i = 0; i < row; i++) {
            m[i].resize(col);
            for (int j = 0; j < col; j++) {
                m[i][j] = d;
            }
        }
    }
    Matrix(std::vector<std::vector<double>> a)
        : m{a}, nrow{int(a.size())}, ncol{int(a[0].size())} {}
    Matrix(const Matrix& a) : m{a.m}, nrow{a.nrow}, ncol{a.ncol} {}
    ~Matrix() = default;
    std::vector<double>& operator[](int i) { return m[i]; }
    std::vector<double> operator[](int i) const { return m[i]; }
    friend std::ostream& operator<<(std::ostream& out, const Matrix& a) {
        out << std::fixed;
        for (int i = 0; i < a.nrow; i++) {
            for (int j = 0; j < a.ncol; j++)
                out << '\t' << a.m[i][j];
            out << '\n';
        }
        return out << std::defaultfloat;
    }
    Matrix operator-() const {
        Matrix result(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                result[i][j] = -m[i][j];
            }
        }
        return result;
    }
    Matrix operator+(const Matrix& a) const {
        Matrix result(nrow, ncol);
        if (a.nrow != nrow || a.ncol != ncol) {
            std::cerr << "Error in operator +: matrices must be the same size!";
            std::exit(4);
        }
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                result[i][j] = m[i][j] + a[i][j];
            }
        }
        return result;
    }
    Matrix operator-(const Matrix& a) const {
        Matrix result(nrow, ncol);
        if (a.nrow != nrow || a.ncol != ncol) {
            std::cerr << "Error in operator -: matrices must be the same size!";
            std::exit(4);
        }
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                result[i][j] = m[i][j] - a[i][j];
            }
        }
        return result;
    }
    Matrix operator*(const Matrix& a) const {
        Matrix result(nrow, a.ncol);
        if (a.nrow != ncol || a.ncol != nrow) {
            std::cerr << "Error in operator *: matrices must be suitable for multiplication!";
            std::exit(4);
        }
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < a.ncol; j++) {
                for (int k = 0; k < ncol; k++) {
                    result[i][j] += m[i][k] * a[k][j];
                }
            }
        }
        return result;
    }
    Matrix operator*(double k) const {
        Matrix result(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                result[i][j] = k * m[i][j];
            }
        }
        return result;
    }
    friend Matrix operator*(double k, const Matrix& a) {
        Matrix result(a.nrow, a.ncol);
        for (int i = 0; i < a.nrow; i++) {
            for (int j = 0; j < a.ncol; j++) {
                result[i][j] = k * a[i][j];
            }
        }
        return result;
    }
    Matrix operator/(double k) {
        Matrix result(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                result[i][j] = m[i][j] / k;
            }
        }
        return result;
    }
    int row() const { return nrow; }
    int col() const { return ncol; }
    Max abs_max() const {
        double value = m[0][0];
        int row = 0;
        int col = 0;
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                if (std::abs(m[i][j]) > value) {
                    value = std::abs(m[i][j]);
                    row = i;
                    col = j;
                }
            }
        }
        return {value, row, col};
    }
    Matrix transpose() const {
        Matrix result(ncol, nrow);
        for (int i = 0; i < ncol; i++)
            for (int j = 0; j < nrow; j++)
                result.m[i][j] = m[j][i];
        return result;
    }

    Matrix submatrix(int row_begin, int col_begin, int row_end, int col_end) const {
        if (row_begin < 0 || row_begin > row_end || row_end > nrow - 1 || col_begin < 0 ||
            col_begin > col_end || col_end > ncol - 1) {
            std::cerr << "Error in function submatrix: wrong submatrix size!\n";
            std::exit(5);
        }
        Matrix result(row_end - row_begin + 1, col_end - col_begin + 1);
        for (int i = 0; i < row_end - row_begin + 1; i++) {
            for (int j = 0; j < col_end - col_begin + 1; j++) {
                result[i][j] = m[row_begin + i][col_begin + j];
            }
        }
        return result;
    }
};
// Метод Гаусса c выбором главного элемента
Matrix gauss(const Matrix& coef, const Matrix& right_part) {
    if (coef.row() != coef.col() || coef.row() != right_part.row() || right_part.col() != 1) {
        std::cerr << "Error in function gauss(): size of coef-matrix must be [n]x[n] and size of "
                     "right_part-matrix mustbe [n]x[1].\n";
        std::exit(-1);
    }
    Matrix a = coef;
    Matrix b = right_part;
    Matrix result(b.row(), 1);
    std::vector<int> p(a.row());  // вектор перестановок столбцов матрицы коэффициентов и
                                  // соответствующихим строк решения
    // Прямой ход метода Гаусса
    for (int i = 0; i < a.col() - 1; i++) {
        auto [max, row_max, col_max] = a.submatrix(i, i, a.row() - 1, a.col() - 1).abs_max();
        row_max += i;
        col_max += i;

        if (i != row_max) {  // Перестановка строк
            std::swap(a[i], a[row_max]);
            std::swap(b[i], b[row_max]);
        }
        if (i != col_max) {  // Перестановка столбцов
            for (int k = 0; k < a.row(); k++) {
                std::swap(a[k][i], a[k][col_max]);
                p[i] = col_max;
            }
        }
        if (std::abs(a[i][i]) < 1e-13) {
            std::cerr << "Error in fuction gauss(): coef-matrix must not be degenerate.\n";
            std::exit(-2);
        }
        for (int j = i + 1; j < a.row(); j++) {
            double m = a[j][i] / a[i][i];
            for (int k = i; k < a.col(); k++) {
                a[j][k] -= m * a[i][k];
            }
            b[j][0] -= m * b[i][0];
        }
    }
    // Обратный ход метода Гаусса
    if (std::abs(a[a.row() - 1][a.col() - 1]) < 1e-13) {
        std::cerr << "Error in fuction gauss(): coef-matrix must not be degenerate.\n";
        std::exit(-2);
    }
    for (int i = b.row() - 1; i > -1; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < a.row(); j++) {
            sum += a[i][j] * result[j][0];
        }

        result[i][0] = (b[i][0] - sum) / a[i][i];
    }
    for (int i = p.size() - 1; i >= 0; i--) {  // обратная перестановка строк решения
        if (p[i]) {
            std::swap(result[i][0], result[p[i]][0]);
        }
    }
    return result;
}
// Функция для вывода графиков с помощью gnuplot
void plot(std::vector<double> t, std::vector<double> p, std::string s, int j) {
    FILE* pipe = popen("gnuplot -persist", "w");
    if (!pipe) {
        std::cerr << "Gnuplot not found\n";
        std::exit(6);
    }
    fprintf(pipe, "$p << EOD\n");
    for (int i = 0; i < p.size(); i++) {
        fprintf(pipe, "%f\t%f\n", t[i], p[i]);
    }
    fprintf(pipe, "EOD\n");
    fprintf(pipe, "set xlabel 't, с'\nset ylabel '");
    fprintf(pipe, s.c_str());
    fprintf(pipe, "'\n");
    fprintf(pipe, "set title '");
    fprintf(pipe, s.c_str());
    fprintf(pipe, "'\n");
    fprintf(pipe, "set terminal png size 640, 480 \n");
    fprintf(pipe, "set output 'φ_%d.png'\n", j);

    fprintf(pipe, "plot '$p' using 1:2 with lines notitle\n");
    fflush(pipe);  //, 'sonya_all_phi_1_5.tb9' u 1:%d w lines title 'PA9'\n", j + 1);
    pclose(pipe);
}
// Функция заполнения матрицы проводимости и вектора токов
void init(const int time_iteration, const double t, const double dt, Matrix& node_admittance,
          Matrix& current, const Matrix& basis, const std::vector<double>& uc1,
          const std::vector<double>& uc2, const std::vector<double>& ucb,
          const std::vector<double>& il2) {
    current[0][0] = basis[0][0] / R2 +
                    C1 * (basis[0][0] - basis[1][0] - uc1[time_iteration - 1]) / dt +
                    (basis[0][0] - basis[4][0]) / R21 + I2(t) - (basis[3][0] - basis[0][0]) / Re2;
    current[1][0] = -C1 * (basis[0][0] - basis[1][0] - uc1[time_iteration - 1]) / dt +
                    (basis[1][0] - basis[2][0]) / Rb;
    current[2][0] = -(basis[1][0] - basis[2][0]) / Rb +
                    Cb * (basis[2][0] - basis[4][0] - ucb[time_iteration - 1]) / dt +
                    (basis[2][0] - basis[4][0]) / Ru + Id(basis[2][0], basis[4][0]);
    current[3][0] = -I2(t) + (basis[3][0] - basis[0][0]) / Re2 +
                    (il2[time_iteration - 1] + dt * (basis[3][0] - basis[4][0]) / L2);
    current[4][0] = -Cb * (basis[2][0] - basis[4][0] - ucb[time_iteration - 1]) / dt -
                    (basis[2][0] - basis[4][0]) / Ru - Id(basis[2][0], basis[4][0]) -
                    (basis[0][0] - basis[4][0]) / R21 -
                    (il2[time_iteration - 1] + dt * (basis[3][0] - basis[4][0]) / L2) +
                    C2 * (basis[4][0] - uc2[time_iteration - 1]) / dt;
    node_admittance[0][0] = 1 / R2 + C1 / dt + 1 / R21 + 1 / Re2;
    node_admittance[0][1] = -C1 / dt;
    node_admittance[0][3] = -1 / Re2;
    node_admittance[0][4] = -1 / R21;
    node_admittance[1][0] = -C1 / dt;
    node_admittance[1][1] = C1 / dt + 1 / Rb;
    node_admittance[1][2] = -1 / Rb;
    node_admittance[2][1] = -1 / Rb;
    node_admittance[2][2] = 1 / Rb + Cb / dt + 1 / Ru + dId_dp3(basis[2][0], basis[4][0]);
    node_admittance[2][4] = -Cb / dt - 1 / Ru + dId_dp5(basis[2][0], basis[4][0]);
    node_admittance[3][0] = -1 / Re2;
    node_admittance[3][3] = 1 / Re2 + dt / L2;
    node_admittance[3][4] = -dt / L2;
    node_admittance[4][0] = -1 / R21;
    node_admittance[4][2] = -Cb / dt - 1 / Ru - dId_dp3(basis[2][0], basis[4][0]);

    node_admittance[4][3] = -dt / L2;
    node_admittance[4][4] =
        Cb / dt + 1 / Ru - dId_dp5(basis[2][0], basis[4][0]) + C2 / dt + dt / L2 + 1 / R21;
}
int main() {
    const int n_max = 10;  // максимальное число итераций метода Ньютона
    const double dt_min = 1e-12;  // минимальный шаг интегрирования по времени
    const double eps = 1e-9;  // максимальное значение поправок
    const double eps_min = 1e-3;  // нижняя граница для оценки локальной точности
    const double eps_max = 5e-2;  // верхняя граница для оценки локальной точности
    const double t_max = 1e-3;  // время расчёта
    double t = 0;               // время
    double dt = dt_min;         // шаг интегрирования по времени
    double dt_prev1 = dt, dt_prev2 = dt;  // предыдущие шаги интегрирования по времени
    Matrix basis(5, 1);  // вектор базиса метода (узловые потенциалы)
    Matrix basis_prev1(basis.row(), 1), basis_prev2(basis.row(), 1),
        basis_prev3(basis.row(), 1);  // предыдущие значения базиса
    Matrix node_admittance(basis.row(), basis.row());  // матрица узловых проводимостей
    Matrix current(basis.row(), 1);                    // вектор невязок
    std::vector<double> time;
    time.push_back(t);
    // Переменные состояния
    std::vector<double> uc1;
    uc1.push_back(0);
    std::vector<double> uc2;
    uc2.push_back(0);
    std::vector<double> ucb;
    ucb.push_back(0);
    std::vector<double> il2;
    il2.push_back(0);
    std::vector<double> phi1, phi2, phi3, phi4, phi5;
    int time_iteration = 1, iter = 1;
    while (dt >= dt_min && t <= t_max) {
        Matrix delta(basis.row(), 1, eps_max + 1);  // вектор поправок
        int n = 0;
        basis = ((dt_prev1 + dt) / dt) * (basis_prev1 - basis_prev2) +
                basis_prev2;  // начальные приближения
        while (std::abs(delta.abs_max().value) > eps && n < n_max) {  // метод Ньютона
            init(time_iteration, t, dt, node_admittance, current, basis, uc1, uc2, ucb, il2);
            delta = gauss(node_admittance, -current);
            basis = basis + delta;
            n++;
            iter++;
        }
        if (n > n_max) {
            dt /= 2;
            continue;
        }
        if (time_iteration > 2) {  // оценка локальной точности
            double d = 0.5 * dt * dt *
                       std::abs(((basis_prev1 - basis_prev2) * (1 / (dt_prev1 * dt_prev1)) -
                                 (basis_prev2 - basis_prev3) * (1 / (dt_prev1 * dt_prev2)))
                                    .abs_max()
                                    .value);
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
        basis_prev3 = basis_prev2;
        basis_prev2 = basis_prev1;
        basis_prev1 = basis;
        time.push_back(t);
        phi1.push_back(basis[0][0]);
        phi2.push_back(basis[1][0]);
        phi3.push_back(basis[2][0]);
        phi4.push_back(basis[3][0]);
        phi5.push_back(basis[4][0]);
        uc1.push_back(basis[0][0] - basis[1][0]);
        uc2.push_back(basis[4][0]);
        ucb.push_back(basis[2][0] - basis[4][0]);
        il2.push_back(il2[time_iteration - 1] + dt * (basis[3][0] - basis[4][0]) / L2);
        time_iteration++;
    }
    if (dt < dt_min) {
        std::cerr << "dt < dt_min at " << time_iteration << " time iteration!\n";
    }
    std::cout << "Итераций по времени: " << time_iteration << std::endl;
    std::cout << "Всего итераций: " << iter << std::endl;
    plot(time, phi1, "φ_1", 1);
    plot(time, phi2, "φ_2", 2);
    plot(time, phi3, "φ_3", 3);
    plot(time, phi4, "φ_4", 4);
    plot(time, phi5, "φ_5", 5);
    return 0;
}
