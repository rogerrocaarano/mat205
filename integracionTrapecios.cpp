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

int main() {
    double h = (b - a) / n;
    // Construir tabla de valores:
    double xi = a;
    std::cout << "TABLA DE VALORES" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "i\t" << "x\t\t" << "y" << std::endl;
    for (int i = 0; i <= n; i++) {
        double yi = y(xi);
        std::cout << i << "\t" << xi << "\t\t" << yi << std::endl;
        xi += h;
    }
    // Evaluar y para xa y xb;
    double ya = y(a);
    double yb = y(b);
    double extremos = ya + yb;
    // calcular los valores medios:
    double medios = 0;
    xi = a;
    for (int i = 1; i < n; i++) {
        xi += h;
        medios += y(xi);
    }

    // calcular el área de la integral:
    double A = (h / 2) * (extremos + 2 * medios);
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Extremos: " << extremos << std::endl;
    std::cout << "Medios: " << medios << std::endl;
    std::cout << "h: " << h << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "INTEGRAL APROXIMADA: " << A << std::endl;
    std::cout << "-------------------------------------------" << std::endl;

    return 0;
}
