#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "matrix.h"

int main(){
        int n = 3;
        Matrix matrixOne(n, n);
        Matrix matrixTwo(n, n);
        Matrix matrixOp(n, n);
        // matrix.printMatrix();

        for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                        matrixOne.setElement(i, j, i + j);
                        matrixTwo.setElement(i, j, i - j);
                }
        }
        
        std::cout << "Initial Matrix A: " << std::endl;
        matrixOne.printMatrix();
        std::cout << "Orthogonalized Matrix Q: " << std::endl;
        matrixOne.orthogonalize();
        matrixOne.printMatrix();
        // std::vector<double> testVector;
        // testVector = matrixOne.eigenValues();
        // for (size_t i = 0; i < n; i++){
        //         std::cout << testVector[i] << " " << std::endl;
        // }
        return 0;
}