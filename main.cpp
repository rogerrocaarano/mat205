#include <iostream>
#include "lib/matrix.h"

int main() {
    Matrix A = Matrix("matrix.txt");
    A.solveGaussSeidel(A, 100);
    Matrix B = Matrix("matrix.txt");
    B.gaussElim();
    B.print();
// Ejemplo descomposición LU
    Matrix AB = Matrix("matrix.txt");
    AB.solveLU();

    return 0;
}
