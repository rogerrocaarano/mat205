#include <iostream>
#include <cmath>

/******************************************************************/
/*************   DEFINIR LOS VALORES DE ENTRADA   *****************/

double y(double x) {
    // y=f(x):
    return (2) / (pow(x, 2) - 4);
}

double a = 0; // Límite inferior de la integral
double b = 0.35; // Límite superior de la integral
int n = 4;    // Cantidad de partes en la que se subdivide el área

/******************************************************************/

bool par(int x) {
    if (x == 0) return false;
    return x % 2 == 0;
}

int main() {
    double h = (b - a) / n;

    // Construir tabla de valores:
    double xi = a;
    double sumImpares = 0;
    double sumPares = 0;
    std::cout << "TABLA DE VALORES" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "i\t" << "x\t\t" << "y" << std::endl;
    for (int i = 0; i <= n; i++) {
        double yi = y(xi);
        std::cout << i << "\t" << xi << "\t\t" << yi << std::endl;
        if (i != 0 && i != n) {
            if (par(i)) {
                sumPares += yi;
            } else {
                sumImpares += yi;
            }
        }
        xi += h;
    }

    // Calcular el área de la integral:
    double y0 = y(a);
    double yn = y(b);
    double A = h / 3 * (y0 + yn + (4 * sumImpares) + (2 * sumPares));
    // Devolver los valores calculados:
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Suma Pares" << ": " << sumPares << std::endl;
    std::cout << "Suma Impares" << ": " << sumImpares << std::endl;
    std::cout << "h: " << h << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "AREA" << ": " << A << std::endl;
    std::cout << "-------------------------------------------" << std::endl;


    return 0;
}
