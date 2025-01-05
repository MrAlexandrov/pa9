#pragma once

#include "matrix.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

constexpr long double EPS = 1e-13;

class IElement {
public:
    virtual ~IElement() = default;

    std::vector<int> GetNodes() const { return Nodes_; }

    virtual std::string GetType() const = 0;
    virtual void PrintInfo() const = 0;

protected:
    std::vector<int> Nodes_;
};


class TResistor : public IElement {
private:
    double resistance;

public:
    TResistor(double r, const std::vector<int>& n) 
        : resistance(r) 
    {
        Nodes_ = n;
    }

    double GetResistance() const { return resistance; }

    std::string GetType() const override { return "Resistor"; }

    void PrintInfo() const override {
        std::cout << "Resistor: R = " << resistance
                  << ", Nodes = (" << Nodes_[0] << ", " << Nodes_[1] << ")\n";
    }
};

class TCapacitor : public IElement {
private:
    double Capacitance_;

public:
    TCapacitor(double c, const std::vector<int>& n) 
        : Capacitance_(c) 
    {
        Nodes_ = n;
    }

    double GetCapacitance() const { return Capacitance_; }

    std::string GetType() const override { return "Capacitor"; }

    void PrintInfo() const override {
        std::cout << "Capacitor: C = " << Capacitance_
                  << ", Nodes = (" << Nodes_[0] << ", " << Nodes_[1] << ")\n";
    }
};

template<typename T>
struct TMaximum {
    T value;
    int row;
    int col;
};

auto FindAbsMax(const NMatrix::TMatrix<>& matrix) -> TMaximum<long double> {
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
NMatrix::TMatrix<> Gauss(
    const NMatrix::TMatrix<>& coef,
    const NMatrix::TMatrix<>& rightPart
) {
    if (coef.Rows() != coef.Cols() || coef.Rows() != rightPart.Rows() || rightPart.Cols() != 1) {
        throw std::runtime_error("Error in function Gauss(): size of coef-matrix must be [n]x[n] and size of rightPart-matrix mustbe [n]x[1].");
    }
    NMatrix::TMatrix<> a = coef;
    NMatrix::TMatrix<> b = rightPart;
    NMatrix::TMatrix<> result(b.Rows(), 1);
    std::vector<int> p(a.Rows());                       // вектор перестановок столбцов матрицы коэффициентов и
                                                          // соответствующихим строк решения

    // Прямой ход метода Гаусса
    for (int i = 0; i < a.Cols() - 1; ++i) {
        auto [max, rowMax, colMax] = 
            FindAbsMax(GetSubMatrix(a, i, a.Rows(), i,  a.Cols()));
        rowMax += i;
        colMax += i;

        if (i != rowMax) {  // Перестановка строк
            std::swap(a[i], a[rowMax]);
            std::swap(b[i], b[rowMax]);
        }
        if (i != colMax) {  // Перестановка столбцов
            for (int k = 0; k < a.Rows(); ++k) {
                std::swap(a[k][i], a[k][colMax]);
                p[i] = colMax;
            }
        }
        if (std::fabs(a[i][i]) < EPS) {
            throw std::runtime_error("Error in fuction Gauss(): coef-matrix must not be degenerate.");
        }
        for (int j = i + 1; j < a.Rows(); ++j) {
            long double m = a[j][i] / a[i][i];
            // if (a[i][i] < static_cast<long double>(1e-20)) {
            //     // m = std::numeric_limits<long double>::max();
            //     m = 1e10;
            //     // throw std::runtime_error("Errrorr");
            // }
            for (int k = i; k < a.Cols(); ++k) {
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

class TInductor : public IElement {
private:
    double inductance;

public:
    TInductor(double l, const std::vector<int>& n) : inductance(l) {
        Nodes_ = n;
    }

    double GetInductance() const { return inductance; }

    std::string GetType() const override { return "Inductor"; }

    void PrintInfo() const override {
        std::cout << "Inductor: L = " << inductance
                  << ", Nodes = (" << Nodes_[0] << ", " << Nodes_[1] << ")\n";
    }
};

class VoltageSource : public IElement {
private:
    double voltage;
public:
    VoltageSource(double v, const std::vector<int>& n) 
        : voltage(v) {
        Nodes_ = n;
    }

    std::string GetType() const override { return "VoltageSource"; }

    void PrintInfo() const override {
        std::cout << "VoltageSource: V = " << voltage
                  << ", Nodes = (" << Nodes_[0] << ", " << Nodes_[1] << ")\n";
    }

    double GetVoltage() const { return voltage; }
};

class CircuitSolver {
private:
    std::vector<std::unique_ptr<IElement>> elements; // Список элементов схемы
    NMatrix::TMatrix<> nodeAdmittance;               // Матрица проводимости
    NMatrix::TMatrix<> residualVector;               // Вектор невязки
    int nodeCount;                                   // Количество узлов в схеме

public:
    CircuitSolver(int nodeCount)
        : nodeAdmittance(nodeCount, nodeCount)
        , residualVector(nodeCount, 1)
        , nodeCount(nodeCount) {}

    void AddElement(std::unique_ptr<IElement> element) {
        elements.push_back(std::move(element));
    }

    void BuildSystem(double dt, const std::vector<double>& state) {
        // Очистка матрицы и вектора
        nodeAdmittance = {};
        residualVector = {};

        // Обработка каждого элемента
        for (const auto& element : elements) {
            const auto& nodes = element->GetNodes();
            if (element->GetType() == "Resistor") {
                const auto* resistor = static_cast<const TResistor*>(element.get());
                double g = 1.0 / resistor->GetResistance(); // Проводимость

                // Заполнение матрицы
                nodeAdmittance[nodes[0]][nodes[0]] += g;
                nodeAdmittance[nodes[0]][nodes[1]] -= g;
                nodeAdmittance[nodes[1]][nodes[0]] -= g;
                nodeAdmittance[nodes[1]][nodes[1]] += g;
            }
            else if (element->GetType() == "Capacitor") {
                const auto* capacitor = static_cast<const TCapacitor*>(element.get());
                double g = capacitor->GetCapacitance() / dt; // Эффективная проводимость

                // Заполнение матрицы
                nodeAdmittance[nodes[0]][nodes[0]] += g;
                nodeAdmittance[nodes[0]][nodes[1]] -= g;
                nodeAdmittance[nodes[1]][nodes[0]] -= g;
                nodeAdmittance[nodes[1]][nodes[1]] += g;

                // Заполнение вектора невязки
                residualVector[nodes[0]][0] += g * state[nodes[0]];
                residualVector[nodes[1]][0] -= g * state[nodes[1]];
            } else if (element->GetType() == "Inductor") {
                const auto* inductor = static_cast<const TInductor*>(element.get());
                double g = dt / inductor->GetInductance(); // Эффективная проводимость

                // Заполнение матрицы
                nodeAdmittance[nodes[0]][nodes[0]] += g;
                nodeAdmittance[nodes[0]][nodes[1]] -= g;
                nodeAdmittance[nodes[1]][nodes[0]] -= g;
                nodeAdmittance[nodes[1]][nodes[1]] += g;

                // Заполнение вектора невязки
                residualVector[nodes[0]][0] += g * state[nodes[0]];
                residualVector[nodes[1]][0] -= g * state[nodes[1]];
            }
        }
    }

    NMatrix::TMatrix<> SolveSystem() {
        // Решение системы методом Гаусса
        return Gauss(nodeAdmittance, residualVector);
    }

    void PrintElements() {
        for (const auto& element : elements) {
            element->PrintInfo();
        }
    }
};
