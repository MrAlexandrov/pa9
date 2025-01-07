#pragma once
#include <cassert>
#include <cstdio>
#include <exception>
#include <iomanip>
#include <optional>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

namespace NPlotter {

namespace {

class TGraph final {
public:
    TGraph(const std::string& title, const std::vector<long double>& data)
        : Title_(title)
        , Data_(data)
        {}

    TGraph(const std::string& title, const std::string& color, const std::vector<long double>& data)
        : Title_(title)
        , Color_(color)
        , Data_(data)
        {}

    void SetTitle(std::string&& title) {
        Title_ = std::move(title);
    }

    void SetColor(const std::string& color) {
        Color_ = color;
    }

    void SetData(const std::vector<long double>& data) {
        Data_ = data;
    }

    std::string GetTitle() const {
        return Title_;
    }

    std::optional<std::string> GetColor() const {
        return Color_;
    }

    std::vector<long double> GetData() const {
        return Data_;
    }

    long double GetData(int index) const {
        assert(0 <= index && index < static_cast<int>(Data_.size()));
        return Data_[index];
    }

private:
    std::string Title_;
    std::optional<std::string> Color_;
    std::vector<long double> Data_;
};

} // namespace

class TPlotter final {
public:
    explicit TPlotter(const std::string& filename)
        : Width_(800)
        , Height_(600)
        , ImageName_(filename + ".png")
        , DataFile_(filename + ".csv")
        , Title_("Title")
        , XLabel_("Time, seconds")
        , YLabel_("Phi value")
        {}

    void SetXValues(const std::vector<long double>& values) {
        assert(!DataSize_.has_value() || DataSize_ == values.size());
        if (!DataSize_.has_value()) {
            DataSize_ = values.size();
        }
        XValues_ = values;
    }

    void AddGraphic(const std::string& name, const std::vector<long double>& data) {
        assert(!DataSize_.has_value() || DataSize_ == data.size());
        if (!DataSize_.has_value()) {
            DataSize_ = data.size();
        }
        Datas_.emplace_back(name, data);
    }

    void AddGraphic(const std::string& name, const std::string& color, const std::vector<long double>& data) {
        assert(!DataSize_.has_value() || DataSize_ == data.size());
        if (!DataSize_.has_value()) {
            DataSize_ = data.size();
        }
        Datas_.emplace_back(name, color, data);
    }

    ~TPlotter() {
        try {
            Plot();
        } catch (...) {
            std::cerr << "Unknown error in TPlotter destructor" << std::endl;
        }
    }

private:
    long double RoundIfApproxZero(long double value) const {
        constexpr long double EPS = 1e-300;
        return (-EPS < value && value < EPS ? 0 : value);
    }

    void SaveDataToFile() const {
        std::ofstream outFile(DataFile_);
        if (!outFile.is_open()) {
            throw std::ios_base::failure("Failed to open the output file: " + DataFile_);
        }
        for (int row = 0; row < DataSize_; ++row) {
            outFile << std::setw(12) << RoundIfApproxZero(XValues_[row]) << "\t";
            for (auto& data : Datas_) {
                outFile << std::setw(12) << RoundIfApproxZero(data.GetData(row)) << "\t";
            }
            outFile << "\n";
        }
        outFile.flush();
        outFile.close();
    }

    std::string GenerateGnuplotCommands() const {
        std::ostringstream commands;
        commands << "set term png size " << Width_ << "," << Height_ << "\n";
        commands << "set output '" << ImageName_ << "'\n";
        commands << "set title '" << Title_ << "'\n";
        commands << "set xlabel '" << XLabel_ << "'\n";
        commands << "set ylabel '" << YLabel_ << "'\n";
        commands << "set xrange [" << XValues_.front() << ":" << XValues_.back() << "]\n";
        commands << "plot ";

        int index = 2;
        for (const auto& data : Datas_) {
            commands << "'" << DataFile_ << "' using 1:" << index++ << " with lines ";
            auto color = data.GetColor();
            if (color.has_value()) {
                commands << "linecolor rgb \"" << color.value() << "\" ";
            }
            commands << "title '" << data.GetTitle() << "', ";
        }
        std::string result = commands.str();
        result.pop_back();
        result.pop_back();
        return result;
    }

    void Plot() const {
        SaveDataToFile();
        try {
            FILE* pipe{popen("gnuplot -persist", "w")};
            if (!pipe) {
                throw std::runtime_error("Failed to open pipe to gnuplot.");
            }
            std::string commands = GenerateGnuplotCommands();
            if (fwrite(commands.c_str(), sizeof(char), commands.size(), pipe) != commands.size()) {
                throw std::runtime_error("Failed to write to Gnuplot pipe.");
            }
            pclose(pipe);
        } catch (std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Unknow error, while plotting" << std::endl;
        }
    }

private:
    unsigned int Width_;
    unsigned int Height_;
    std::string ImageName_;
    std::string DataFile_;
    std::string Title_;
    std::string XLabel_;
    std::string YLabel_;
    std::vector<long double> XValues_;
    std::optional<int> DataSize_;
    std::vector<TGraph> Datas_;
};

} // NPlotter
