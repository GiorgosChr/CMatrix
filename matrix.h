#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix{
        private:
                std::vector <std::vector<double>> matrix;
                int rows, cols;
        
        public:
                Matrix(int x, int y){
                        rows = x;
                        cols = y;
                }             
};

#endif