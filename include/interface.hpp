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
    double Resistance_;

public:
    TResistor(double r, const std::vector<int>& n) 
        : Resistance_(r) 
    {
        Nodes_ = n;
    }

    double GetResistance() const { return Resistance_; }

    std::string GetType() const override { return "Resistor"; }

    void PrintInfo() const override {
        std::cout << "Resistor: R = " << Resistance_
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

class TInductor : public IElement {
private:
    double Inductance_;

public:
    TInductor(double l, const std::vector<int>& n) : Inductance_(l) {
        Nodes_ = n;
    }

    double GetInductance() const { return Inductance_; }

    std::string GetType() const override { return "Inductor"; }

    void PrintInfo() const override {
        std::cout << "Inductor: L = " << Inductance_
                  << ", Nodes = (" << Nodes_[0] << ", " << Nodes_[1] << ")\n";
    }
};

class TVoltageSource : public IElement {
private:
    double Voltage_;
public:
    TVoltageSource(double v, const std::vector<int>& n) 
        : Voltage_(v) {
        Nodes_ = n;
    }

    std::string GetType() const override { return "VoltageSource"; }

    void PrintInfo() const override {
        std::cout << "VoltageSource: V = " << Voltage_
                  << ", Nodes = (" << Nodes_[0] << ", " << Nodes_[1] << ")\n";
    }

    double GetVoltage() const { return Voltage_; }
};

template<typename T>
struct TMaximum {
    T Value;
    int Row;
    int Col;
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

class CircuitSolver {
private:
    std::vector<std::unique_ptr<IElement>> Elements_; // Список элементов схемы
    NMatrix::TMatrix<> NodeAdmittance_;               // Матрица проводимости
    NMatrix::TMatrix<> ResidualVector_;               // Вектор невязки
    // int NodeCount_;                                   // Количество узлов в схеме

public:
    CircuitSolver(int nodeCount)
        : NodeAdmittance_(nodeCount, nodeCount)
        , ResidualVector_(nodeCount, 1)
        // , NodeCount_(nodeCount) 
        {
        }

    void AddElement(std::unique_ptr<IElement> element) {
        Elements_.push_back(std::move(element));
    }

    void BuildSystem(double dt, const std::vector<double>& state) {
        // Очистка матрицы и вектора
        NodeAdmittance_ = {};
        ResidualVector_ = {};

        // Обработка каждого элемента
        for (const auto& element : Elements_) {
            const auto& nodes = element->GetNodes();
            if (element->GetType() == "Resistor") {
                const auto* resistor = static_cast<const TResistor*>(element.get());
                double g = 1.0 / resistor->GetResistance(); // Проводимость

                // Заполнение матрицы
                NodeAdmittance_[nodes[0]][nodes[0]] += g;
                NodeAdmittance_[nodes[0]][nodes[1]] -= g;
                NodeAdmittance_[nodes[1]][nodes[0]] -= g;
                NodeAdmittance_[nodes[1]][nodes[1]] += g;
            }
            else if (element->GetType() == "Capacitor") {
                const auto* capacitor = static_cast<const TCapacitor*>(element.get());
                double g = capacitor->GetCapacitance() / dt; // Эффективная проводимость

                // Заполнение матрицы
                NodeAdmittance_[nodes[0]][nodes[0]] += g;
                NodeAdmittance_[nodes[0]][nodes[1]] -= g;
                NodeAdmittance_[nodes[1]][nodes[0]] -= g;
                NodeAdmittance_[nodes[1]][nodes[1]] += g;

                // Заполнение вектора невязки
                ResidualVector_[nodes[0]][0] += g * state[nodes[0]];
                ResidualVector_[nodes[1]][0] -= g * state[nodes[1]];
            } else if (element->GetType() == "Inductor") {
                const auto* inductor = static_cast<const TInductor*>(element.get());
                double g = dt / inductor->GetInductance(); // Эффективная проводимость

                // Заполнение матрицы
                NodeAdmittance_[nodes[0]][nodes[0]] += g;
                NodeAdmittance_[nodes[0]][nodes[1]] -= g;
                NodeAdmittance_[nodes[1]][nodes[0]] -= g;
                NodeAdmittance_[nodes[1]][nodes[1]] += g;

                // Заполнение вектора невязки
                ResidualVector_[nodes[0]][0] += g * state[nodes[0]];
                ResidualVector_[nodes[1]][0] -= g * state[nodes[1]];
            }
        }
    }

    NMatrix::TMatrix<> SolveSystem() {
        // Решение системы методом Гаусса
        return Gauss(NodeAdmittance_, ResidualVector_);
    }

    void PrintElements() {
        for (const auto& element : Elements_) {
            element->PrintInfo();
        }
    }
};
