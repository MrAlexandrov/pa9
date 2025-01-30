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
#include <memory>
#include <typeindex>
#include <typeinfo>

namespace NPlotter {

template <typename YType>
class TGraph;

class IGraph {
public:
    virtual ~IGraph() = default;
    
    virtual std::string GetTitle() const = 0;
    virtual std::optional<std::string> GetColor() const = 0;
    
    virtual std::type_index GetDataType() const = 0;

    virtual std::vector<double> GetDataAsDouble() const = 0;

    template <typename T>
    const std::vector<T>& GetOriginalData() const {
        if (GetDataType() != typeid(T)) {
            throw std::bad_cast();
        }
        return dynamic_cast<const TGraph<T>*>(this)->GetData();
    }

    virtual size_t GetSize() const = 0;
};

template <typename YType>
class TGraph : public IGraph {
public:
    using value_type = YType;

    TGraph(const std::string& title, const std::vector<YType>& data)
        : Title_(title)
        , Data_(data)
        {
        }

    TGraph(const std::string& title, const std::vector<YType>& data, const std::string& color)
        : Title_(title)
        , Data_(data)
        , Color_(color)
        {
        }

    std::string GetTitle() const override { return Title_; }
    std::optional<std::string> GetColor() const override { return Color_; }
    
    size_t GetSize() const override { return Data_.size(); }

    std::type_index GetDataType() const override { return typeid(YType); }

    std::vector<double> GetDataAsDouble() const override {
        return std::vector<double>(Data_.begin(), Data_.end());
    }

    const std::vector<YType>& GetData() const {
        return Data_;
    }

private:
    std::string Title_;
    std::vector<YType> Data_;
    std::optional<std::string> Color_;
};

template <typename XType>
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
        {
        }

    void SetXValues(const std::vector<XType>& values) {
        assert(!DataSize_.has_value() || DataSize_ == values.size());
        if (!DataSize_.has_value()) {
            DataSize_ = values.size();
        }
        XValues_ = values;
    }

    template <typename YType, typename... Args>
    void EmplaceGraphic(Args&&... args) {
        Graphs_.emplace_back(std::make_unique<TGraph<YType>>(std::forward<Args>(args)...));
    }

    void Plot() const {
        this->SaveDataToFile();
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
    void SaveDataToFile() const {
        std::ofstream outFile(DataFile_);
        if (!outFile.is_open()) {
            throw std::ios_base::failure("Failed to open the output file: " + DataFile_);
        }

        for (int row = 0; row < DataSize_; ++row) {
            outFile << std::setw(12) << XValues_[row] << "\t";
            for (const auto& graph : Graphs_) {
                outFile << std::setw(12) << graph->GetDataAsDouble()[row] << "\t";
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
        for (const auto& data : Graphs_) {
            commands << "'" << DataFile_ << "' using 1:" << index++ << " with lines title '" << data->GetTitle() << "', ";
        }
        std::string result = commands.str();
        result.pop_back();
        result.pop_back();
        return result;
    }

private:
    unsigned int Width_;
    unsigned int Height_;
    std::string ImageName_;
    std::string DataFile_;
    std::string Title_;
    std::string XLabel_;
    std::string YLabel_;
    std::vector<XType> XValues_;
    std::optional<int> DataSize_;
    std::vector<std::unique_ptr<IGraph>> Graphs_;
};

}  // namespace NPlotter
