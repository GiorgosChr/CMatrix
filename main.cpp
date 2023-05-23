#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "matrix.h"

int main(){
        Matrix matrixOne(4, 4);
        Matrix matrixTwo(4, 4);
        Matrix matrixOp(4, 4);
        // matrix.printMatrix();

        for (int i = 0; i < 4; i++){
                for (int j = 0; j < 4; j++){
                        matrixOne.setElement(i, j, i + j);
                        matrixTwo.setElement(i, j, i - j);
                }
        }
        matrixOne.printMatrix();
        std::cout << std::endl;

        matrixTwo.printMatrix();
        std::cout << std::endl;
        
        matrixOp = matrixTwo.dot(matrixOne);
        matrixOp.printMatrix();
        std::cout << std::endl;

        matrixOp.T();
        matrixOp.printMatrix();
        std::cout << std::endl;
        return 0;
}