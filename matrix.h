#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>


class Matrix{
        private:
                std::vector <std::vector<double>> matrix;
                int rows, cols;

                std::vector<double> projection(std::vector<double> const& vectorV, std::vector<double> const& vectorU){
                        std::vector<double> vectorProj;
                        double dotProduct = 0.0;
                        double norm = 0.0;
                        for (size_t i = 0; i < vectorV.size(); i++){
                                dotProduct += vectorU[i] * vectorV[i];
                                norm += vectorU[i] * vectorU[i];
                        }
                        for (size_t i = 0; i < vectorV.size(); i++){
                                vectorProj.push_back(dotProduct * vectorU[i]/norm);
                        }

                        return vectorProj;
                }
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

                void absoluteVal(){ // Calculates the absolute value of each matrix element
                        for (size_t i = 0; i < rows; i++) {
                                for (size_t j = 0; j < cols; j++) {
                                        matrix[i][j] = abs(matrix[i][j]);
                                }
                        }
                }

                void orthogonalize(){ // Modified Gram-Schmidt orthogonalization
                        Matrix matrixOld(rows, cols);
                        matrixOld.matrix = matrix;
                        std::vector<std::vector<double>> vectorsOld;
                        for (size_t k = 0; k < cols; k++){ // Store all vectors in a 2d vector
                                std::vector<double> vectorV; // The old vector
                                for (size_t i = 0; i < rows; i++){
                                        vectorV.push_back(matrixOld.matrix[i][k]);
                                }
                                vectorsOld.push_back(vectorV);
                        }
                        for (size_t k = 0; k < cols; k++){ // Calculate the orthogonalized vectors
                                std::vector<double> vectorU; // The new vector
                                for (size_t i = 0; i < vectorsOld[k].size(); i++){
                                        vectorU.push_back(vectorsOld[k][i]);
                                }
                                if (k > 0){
                                        for (size_t kPrime = 0; kPrime < k; kPrime++){
                                                std::vector<double> vectorPrev;
                                                for (size_t i = 0; i < vectorsOld[k -1].size(); i++){
                                                        vectorPrev.push_back(matrixOld.matrix[i][k - 1]);
                                                }
                                                std::vector<double> proj = projection(vectorU, vectorPrev);
                                                for (size_t i = 0; i < vectorU.size(); i++){
                                                        vectorU[i] -= proj[i];
                                                }
                                        }
                                }
                                for (size_t i = 0; i < rows; i++){
                                        matrix[i][k] = vectorU[i];
                                }
                        }
                }




                std::vector<double> eigenValues(double tol = pow(10, -10)){ // Calculates the eigenvalues of a matrix and returns them as vector components
                        Matrix eigenMatrix(rows, cols);
                        Matrix matrixQ(rows, cols), matrixR(rows, cols);
                        Matrix matrixA(rows, cols);
                        Matrix matrixAPrev(rows, cols);
                        Matrix matrixTemp(rows, cols);
                        std::vector<double> eigenValues;
                        double diff = 2.0 * tol; // Difference between the matrix elements of the matrix A after each iteration
                        
                        matrixA.matrix = matrix;
                        matrixTemp.matrix = matrix;
                        // RNG setup
                        std::random_device rd;
                        std::mt19937 gen(rd());
                        std::uniform_real_distribution<double> dist(0.0, 1.0);

                        matrixTemp.orthogonalize();
                        matrixTemp.printMatrix();
                        matrixQ.matrix = matrixTemp.matrix;
                        matrixTemp.T();
                        matrixR = matrixTemp.dot(matrixA);
                        
                        while (diff > tol){
                                matrixAPrev.matrix = matrixA.matrix;
                                matrixA = matrixR.dot(matrixQ);
                                
                                double sum = 0.0;
                                for (size_t i = 0; i < rows; i++){
                                        for (size_t j = 0; j < cols; j++){
                                                sum += abs(matrixA.matrix[i][j] - matrixAPrev.matrix[i][j]);
                                        }
                                }
                                // std::cout << sum << std::endl;
                                diff = sum;
                                matrixTemp.matrix = matrixA.matrix;
                                matrixTemp.orthogonalize();
                                matrixQ.matrix = matrixTemp.matrix;
                                matrixTemp.T();
                                matrixR = matrixTemp.dot(matrixA);
                        }
                


                        
                        for (size_t i = 0; i < rows; i++){
                                eigenValues.push_back(matrixA.matrix[i][i]);
                        }
                        return eigenValues;
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