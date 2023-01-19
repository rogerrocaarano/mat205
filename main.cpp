#include <iostream>
#include "lib/matrix.h"

int main() {
//    Matrix gs = Matrix("matrix.txt");
//    gs.solveGaussSeidel(4);
//    Matrix jb = Matrix("matrix.txt");
//    jb.solveJacobi(4);
// Ejemplo descomposición LU
//    Matrix AB = Matrix("matrix.txt");
//    AB.solveLU();
    Matrix AB = Matrix("matrix.txt");
    AB.solveLU();

    return 0;
}
