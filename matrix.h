#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

class Matrix{
        private:
                std::vector <std::vector<double>> matrix;
                int rows, cols;
        
        public:
                // Constructor
                Matrix(int x, int y){
                        rows = x;
                        cols = y;
                        matrix.resize(rows, std::vector<double>(cols, 0));
                }
                
                // Matrix Functions

                void printMatrix(){ // Prints the matrix
                        for (const auto& row : matrix) {
                                for (const auto& matrixElement : row){
                                        std::cout << matrixElement << " ";
                                }
                                std::cout << std::endl;
                        }
                }

                void setElement(int x, int y, double value){ // Sets the element of the matrix to a specific value
                        if (x >= 0 && x < rows && y >=0 && y < cols){
                                matrix[x][y] = value;
                        }
                        else{
                                std::cout << "Out of bounds.";
                        }
                }

                double getElement(int x, int y){ // Prints the element of the matrix
                        if (x >= 0 && x < rows && y >=0 && y < cols){
                                return matrix[x][y];
                        }
                        else{
                                std::cout << "Out of bounds.";
                                return 0;
                        }
                }

                void unitMatrix(){ // Converts the matrix to a unit matrix
                        for (size_t i = 0; i < rows; i++){
                                for (size_t j = 0; j < cols; j++){
                                        if (i == j){
                                                matrix[i][j] = 1.0;
                                        }
                                        else{
                                                matrix[i][j] = 0.0;
                                        }
                                }
                        }
                }

                void ones(){ // Fills the matrix with ones
                        for (size_t i = 0; i < rows; i++){
                                for (size_t j = 0; j < cols; j++){
                                        matrix[i][j] = 1.0;
                                }
                        }
                }

                void zeros(){ // Fills the matrix with zeros
                        for (size_t i = 0; i < rows; i++){
                                for (size_t j = 0; j < cols; j++){
                                        matrix[i][j] = 0.0;
                                }
                        }
                }

                void T(){ //Transposes the matrix
                        Matrix temp(rows, cols);
                        for (size_t i = 0; i < rows; i++){
                                for (size_t j = 0; j < cols; j++){
                                        temp.matrix[i][j] = matrix[i][j];
                                }
                        }
                        for (size_t i = 0; i < rows; i++){
                                for (size_t j = 0; j < cols; j++){
                                        matrix[i][j] = temp.matrix[j][i];
                                }
                        }
                        std::swap(rows, cols);
                }

                double trace(){ // Calculates the trace of the matrix
                        double trace = 0.0;
                        for (size_t i = 0; i < rows; i++){
                                for (size_t j = 0; j < cols; j++){
                                        if (i == j){
                                                trace += matrix[i][j];
                                        }
                                }
                        }

                        return trace;
                }

                double det(){ // Calculate the determinant
                        if (rows != cols) {
                                std::cout << "Error: The matrix is not square.\n";
                                return 0.0;
                        }
                        if (rows == 1) { // For 1x1 matrix
                                return matrix[0][0];
                        }
                        
                        double determinant = 0.0;
                        int sign = 1;
                        Matrix submatrix(rows - 1, cols - 1);

                        for (int col = 0; col < cols; ++col) {
                                int subi = 0;
                                for (int i = 1; i < rows; ++i) {
                                        int subj = 0;
                                        for (int j = 0; j < cols; ++j) {
                                                if (j != col) {
                                                        submatrix.matrix[subi][subj] = matrix[i][j];
                                                        ++subj;
                                                }
                                        }
                                        ++subi;
                                }

                                determinant += sign * matrix[0][col] * submatrix.det();
                                sign = -sign;
                        }

                        return determinant;
                }

                void absoluteVal() { // Calculates the absolute value of each matrix element
                        for (size_t i = 0; i < rows; i++) {
                                for (size_t j = 0; j < cols; j++) {
                                        matrix[i][j] = abs(matrix[i][j]);
                                }
                        }
                }

                // Operations between two matrices

                Matrix operator+(Matrix const& secMatrix){
                        Matrix result(rows, cols);
                        for (size_t i = 0; i < rows; i++){
                                for (size_t j = 0; j < cols; j++){
                                        result.matrix[i][j] = matrix[i][j] + secMatrix.matrix[i][j];
                                }
                        }

                        return result;
                }

                Matrix operator-(Matrix const& secMatrix){
                        Matrix result(rows, cols);
                        for (size_t i = 0; i < rows; i++){
                                for (size_t j = 0; j < cols; j++){
                                        result.matrix[i][j] = matrix[i][j] - secMatrix.matrix[i][j];
                                }
                        }

                        return result;
                }

                Matrix operator*(Matrix const& secMatrix){
                        Matrix result(rows, cols);
                        for (size_t i = 0; i < rows; i++){
                                for (size_t j = 0; j < cols; j++){
                                        result.matrix[i][j] = matrix[i][j] * secMatrix.matrix[i][j];
                                }
                        }

                        return result;
                }

                Matrix operator/(Matrix const& secMatrix){
                        Matrix result(rows, cols);
                        for (size_t i = 0; i < rows; i++){
                                for (size_t j = 0; j < cols; j++){
                                        result.matrix[i][j] = matrix[i][j] / secMatrix.matrix[i][j];
                                }
                        }

                        return result;
                }

                // Matrix multiplication

                Matrix dot(Matrix const& secMatrix){
                        if (cols != secMatrix.rows) {
                                std::cout << "Invalid matrix dimensions for multiplication.";
                                return Matrix(0, 0);
                        }

                        Matrix result(rows, secMatrix.cols);
                        
                        for (size_t i = 0; i < rows; i++){
                                for (size_t j = 0; j < secMatrix.cols; j++){
                                        double sum = 0.0;
                                        for (size_t k = 0; k < cols; k++){
                                                sum += matrix[i][k] * secMatrix.matrix[k][j];
                                        }
                                        result.matrix[i][j] = sum;
                                }
                        }

                        return result;
                }

                // Operations between matrices and scalars

                Matrix operator^(double power){
                        Matrix result(rows, cols);
                        for (size_t i = 0; i < rows; i++) {
                                for (size_t j = 0; j < cols; j++) {
                                        result.matrix[i][j] = pow(matrix[i][j], power);
                                }
                        }

                        return result;
                }

                Matrix operator+(double scalar) {
                        Matrix result(rows, cols);
                        for (size_t i = 0; i < rows; i++) {
                                for (size_t j = 0; j < cols; j++) {
                                        result.matrix[i][j] = matrix[i][j] + scalar;
                                }
                        }

                        return result;
                }

                friend Matrix operator+(double scalar, Matrix const& matrix) {
                        Matrix result(matrix.rows, matrix.cols);
                        for (size_t i = 0; i < matrix.rows; i++) {
                                for (size_t j = 0; j < matrix.cols; j++) {
                                        result.matrix[i][j] = scalar + matrix.matrix[i][j];
                                }
                        }

                        return result;
                }

                Matrix operator-(double scalar) {
                        Matrix result(rows, cols);
                        for (size_t i = 0; i < rows; i++) {
                                for (size_t j = 0; j < cols; j++) {
                                        result.matrix[i][j] = matrix[i][j] - scalar;
                                }
                        }

                        return result;
                }

                friend Matrix operator-(double scalar, Matrix const& matrix) {
                        Matrix result(matrix.rows, matrix.cols);
                        for (size_t i = 0; i < matrix.rows; i++) {
                                for (size_t j = 0; j < matrix.cols; j++) {
                                        result.matrix[i][j] = scalar - matrix.matrix[i][j];
                                }
                        }

                        return result;
                }

                Matrix operator*(double scalar) {
                        Matrix result(rows, cols);
                        for (size_t i = 0; i < rows; i++) {
                                for (size_t j = 0; j < cols; j++) {
                                        result.matrix[i][j] = matrix[i][j] * scalar;
                                }
                        }

                        return result;
                }

                friend Matrix operator*(double scalar, Matrix const& matrix) {
                        Matrix result(matrix.rows, matrix.cols);
                        for (size_t i = 0; i < matrix.rows; i++) {
                                for (size_t j = 0; j < matrix.cols; j++) {
                                        result.matrix[i][j] = scalar * matrix.matrix[i][j];
                                }
                        }

                        return result;
                }

                Matrix operator/(double scalar) {
                        Matrix result(rows, cols);
                        for (size_t i = 0; i < rows; i++) {
                                for (size_t j = 0; j < cols; j++) {
                                        result.matrix[i][j] = matrix[i][j] / scalar;
                                }
                        }

                        return result;
                }

                friend Matrix operator/(double scalar, Matrix const& matrix) {
                        Matrix result(matrix.rows, matrix.cols);
                        for (size_t i = 0; i < matrix.rows; i++) {
                                for (size_t j = 0; j < matrix.cols; j++) {
                                        result.matrix[i][j] = scalar / matrix.matrix[i][j];
                                }
                        }

                        return result;
                }
};


#endif