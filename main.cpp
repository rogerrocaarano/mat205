#include <iostream>
#include "lib/matrix.h"

int main() {
    Matrix A = Matrix("matrix.txt");
    A.gaussSeidel(A, 10);


    return 0;
}
