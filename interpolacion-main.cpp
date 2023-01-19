#include "lib/interpolacion.h"

int main() {
    std::vector<double> x;
    std::vector<double> y;
    x = {1, -3, 5, 7};
    y = {-2, 1, 2, -3};

    polinomio R;
    R = lagrange::interpolar(x, y);
    R.print();
//    std::cout << R.eval(-1);
    return 0;
}
