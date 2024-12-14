#include <gtest/gtest.h>
#include <stdexcept>
#include "../include/matrix.hpp"

using namespace NMatrix;

TEST(MatrixConstructorTest, DefaultConstructor) {
    NMatrix::TMatrix<> mat;
    EXPECT_EQ(mat.rows(), 0.0);
    EXPECT_EQ(mat.cols(), 0.0);
}

TEST(MatrixConstructorTest, ParameterizedConstructor) {
    TMatrix<> mat(1, 3);

    TMatrix<> expected{{0, 0, 0}};
    EXPECT_EQ(mat, expected);
}

TEST(MatrixConstructorTest, ParameterizedConstructorWithDefaultValue) {
    TMatrix<> mat(3, 3, 1.0);
    EXPECT_EQ(mat.rows(), 3);
    EXPECT_EQ(mat.cols(), 3);

    TMatrix<> expected{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
    EXPECT_EQ(mat, expected);
}

TEST(MatrixAccessTest, OutOfBoundsAccessRows) {
    TMatrix<> mat(3, 3, 1.0);
    EXPECT_THROW(mat[3][3], std::out_of_range);
}

TEST(MatrixAccessTest, OutOfBoundsAccessCols) {
    TMatrix<> mat(3, 3, 1.0);
    EXPECT_THROW(mat[0][3], std::out_of_range);
}

TEST(MatrixAccessTest, Rows) {
    TMatrix<> mat(1, 1);
    EXPECT_EQ(mat.rows(), 1);
}

TEST(MatrixAccessTest, Cols) {
    TMatrix<> mat(1, 1);
    EXPECT_EQ(mat.cols(), 1);
}

TEST(MatrixOperationsTest, AdditionCorrect) {
    TMatrix<> mat1(2, 2, 1.0);
    TMatrix<> mat2(2, 2, 2.0);
    TMatrix<> result = mat1 + mat2;

    TMatrix<> expected{{3, 3}, {3, 3}};
    EXPECT_EQ(result, expected);
}

TEST(MatrixOperationsTest, AdditionIncorrect) {
    TMatrix<> mat1(1, 1);
    TMatrix<> mat2(1, 2);

    EXPECT_THROW(mat1 + mat2, std::invalid_argument);
}


TEST(MatrixOperationsTest, SubtractionCorrect) {
    TMatrix<> mat1(2, 2, 3.0);
    TMatrix<> mat2(2, 2, 2.0);
    TMatrix<> result = mat1 - mat2;

    TMatrix<> expected{{1, 1}, {1, 1}};
    EXPECT_EQ(result, expected);
}

TEST(MatrixOperationsTest, SubtractionInorrect) {
    TMatrix<> mat1(1, 1, 3.0);
    TMatrix<> mat2(2, 2, 2.0);

    EXPECT_THROW(mat1 - mat2, std::invalid_argument);
}

TEST(MatrixOperationsTest, MultiplicationCorrect) {
    TMatrix<> mat1{{2, -3, 1}, {5, 4, -2}};
    TMatrix<> mat2{{-7, 5}, {2, -1}, {4, 3}};
    TMatrix<> result = mat1 * mat2;

    TMatrix<> expected{{-16, 16}, {-35, 15}};
    EXPECT_EQ(result, expected);
}


TEST(MatrixOperationsTest, MultiplicationIncorrect) {
    TMatrix<> mat1{{1}, {2}};
    TMatrix<> mat2{{2, 0}, {1, 2}};

    EXPECT_THROW(mat1 * mat2, std::invalid_argument);
    
    TMatrix<> result = mat2 * mat1;
    TMatrix<> expected{{2}, {5}};
    EXPECT_EQ(result, expected);
}


TEST(MatrixOperationsTest, AdittionEqualScalarCorrect) {
    TMatrix<> mat{{1, 2, 3}, {4, 5, 6}};
    mat += 1;

    TMatrix<> expected{{2, 3, 4}, {5, 6, 7}};
    EXPECT_EQ(mat, expected);
}

TEST(MatrixOperationsTest, SubtractionEqualScalarCorrect) {
    TMatrix<> mat{{1, 2, 3}, {4, 5, 6}};
    mat -= 1;

    TMatrix<> expected{{0, 1, 2}, {3, 4, 5}};
    EXPECT_EQ(mat, expected);
}

TEST(MatrixOperationsTest, MultiplicationEqualScalarCorrect) {
    TMatrix<> mat{{1, 2, 3}, {4, 5, 6}};
    mat *= 2;

    TMatrix<> expected{{2, 4, 6}, {8, 10, 12}};
    EXPECT_EQ(mat, expected);
}

TEST(MatrixOperationsTest, DivisionEqualScalarCorrect) {
    TMatrix<> mat{{2, 4, 6}, {8, 10, 12}};
    mat /= 2;

    TMatrix<> expected{{1, 2, 3}, {4, 5, 6}};
    EXPECT_EQ(mat, expected);
}

TEST(MatrixOperationsTest, DivisionEqualScalarIncorrect) {
    TMatrix<> mat{{1, 2, 3}, {4, 5, 6}};
    EXPECT_THROW(mat /= 0, std::runtime_error);
}

TEST(MatrixOperationsTest, AdittionScalarCorrect) {
    TMatrix<> mat{{1, 2, 3}, {4, 5, 6}};
    TMatrix<> result = mat + 1;

    TMatrix<> expected{{2, 3, 4}, {5, 6, 7}};
    EXPECT_EQ(result, expected);
}

TEST(MatrixOperationsTest, SubtractionScalarCorrect) {
    TMatrix<> mat{{1, 2, 3}, {4, 5, 6}};
    TMatrix<> result = mat - 1;

    TMatrix<> expected{{0, 1, 2}, {3, 4, 5}};
    EXPECT_EQ(result, expected);
}

TEST(MatrixOperationsTest, MultiplicationScalarCorrect) {
    TMatrix<> mat{{1, 2, 3}, {4, 5, 6}};
    TMatrix<> result = mat * 2;

    TMatrix<> expected{{2, 4, 6}, {8, 10, 12}};
    EXPECT_EQ(result, expected);
}

TEST(MatrixOperationsTest, DivisionScalarCorrect) {
    TMatrix<> mat{{2, 4, 6}, {8, 10, 12}};
    TMatrix<> result = mat / 2;

    TMatrix<> expected{{1, 2, 3}, {4, 5, 6}};
    EXPECT_EQ(result, expected);
}

TEST(MatrixOperationsTest, DivisionScalarIncorrect) {
    TMatrix<> mat{{1, 2, 3}, {4, 5, 6}};

    EXPECT_THROW(TMatrix<> result = mat / 0, std::runtime_error);
}

TEST(MatrixSpecialFunctionsTest, Transpose) {
    TMatrix<> mat{{1, 2, 3}, {4, 5, 6}};
    TMatrix<> result = mat.Transpose();

    TMatrix<> expected{{1, 4}, {2, 5}, {3, 6}};
    EXPECT_EQ(result, expected);
}

TEST(MatrixSpecialFunctionsTest, Equality) {
    TMatrix<> mat1{{1, 2}, {3, 4}};
    TMatrix<> mat2{{1, 2}, {3, 4}};
    TMatrix<> mat3{{4, 3}, {2, 1}};

    EXPECT_TRUE(mat1 == mat2);
    EXPECT_FALSE(mat1 == mat3);
}

// TEST(MatrixSpecialFunctionsTest, Printing) {
//     NMatrix::TMatrix<> matrix(2, 3, 1);

//     std::stringstream output;
//     output << matrix;

//     std::string expected_output = 
//         "1 1 1 \n"
//         "1 1 1 \n";

//     EXPECT_EQ(output.str(), expected_output);
// }

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
