//
// Created by rogerroca on 27/12/2022.
//

#include "interpolacion.h"

polinomio lagrange::interpolar(std::vector<double> &x, std::vector<double> &y) {
    polinomio Pn;
    if (x.size() != y.size()) {
        std::cout << "Revisar datos." << std::endl;
        return Pn;
    }
    int n = x.size();
    for (int i = 0; i < n; i++) {
        polinomio yTerm;
        yTerm.poner_termino(y[i], 0);
        polinomio LiTerm;
        LiTerm = Li(i, x);
        std::cout << "y" << i << ": ";
        yTerm.print();
        std::cout << "L" << i << ": ";
        LiTerm.print();
        polinomio PnTerm;
        PnTerm.multiplicar(yTerm, LiTerm);
        polinomio PnSum;
        PnSum.sumar(Pn, PnTerm);
        Pn = PnSum;
    }
    return Pn;
}

polinomio lagrange::Li(int k, std::vector<double> x) {
    int n = x.size();
    // Cálculo del denominador de la productiva (x-xi)/(xi-xj)
    double denominador = 1;
    for (int i = 0; i < n; i++) {
        if (i != k) {
            denominador = denominador * (x[k] - x[i]);
        }
    }
    // Crear un polinomio y calcular la productiva de (x-xi)/(xi-xj)
    polinomio Li;
    Li.poner_termino(1, 0);
    for (int i = 0; i < n; i++) {
        if (i != k) {
            polinomio term;
            term.poner_termino(1, 1);
            term.poner_termino(-x[i], 0);
            polinomio aux;
            aux.multiplicar(Li, term);
            Li = aux;
        }
    }
    int grado = Li.Grado();
    for (int i = 0; i <= grado; i++) {
        double coef = Li.coeficiente(i);
        coef = coef / denominador;
        Li.AsignarCoeficiente(coef, i);
    }
    return Li;
}

polinomio newton::interpolar(std::vector<double> &x, std::vector<double> &y) {
    polinomio Pn;
    if (x.size() != y.size()) {
        std::cout << "Revisar datos." << std::endl;
        return Pn;
    }
    // Crear una matriz para guardar los datos:
    matriz M_difDivididas;
    M_difDivididas.dimensionar(x.size(), x.size());
    // Cargar los valores de y en la primera fila
    for (int i = 1; i <= y.size(); i++) {
        M_difDivididas.poner(1, i, y[i - 1]);
    }
    // Calcular las diferencias divididas:
    for (int i = 2; i <= M_difDivididas.dimension_fila(); i++) {
        cargarFila_difDividida(M_difDivididas, x, i);
    }
    M_difDivididas.mostrar();
    // Sumar los términos de Pn para obtener el polinomio:
    for (int i = 0; i < x.size(); i++) {
        polinomio Term;
        Term = PnTerm(M_difDivididas, x, i);
        polinomio sum;
        sum.sumar(Pn, Term);
        Pn = sum;
    }
    return Pn;
}

double newton::difDividida(double x0, double x1, double y0, double y1) {
    return (y1 - y0) / (x1 - x0);
}

void newton::cargarFila_difDividida(matriz &m, std::vector<double> x, int fila) {
    int i_x0 = 0;
    int i_x1 = fila - 1;
    for (int i = fila; i <= m.dimension_columna(); i++) {
        double x0 = x[i_x0];
        double x1 = x[i_x1];
        double y0 = m.elemento(fila - 1, i - 1);
        double y1 = m.elemento(fila - 1, i);
        m.poner(fila, i, difDividida(x0, x1, y0, y1));
        i_x0++;
        i_x1++;
    }
}

polinomio newton::PnTerm(matriz &m, std::vector<double> x, int grado) {
    double a = m.elemento(grado + 1, grado + 1); // Normalizando el índice con +1 para la estructura de matriz.
    polinomio PnTerm;
    PnTerm.poner_termino(a, 0);
    if (grado == 0) {
        return PnTerm;
    }
    for (int i = 0; i < grado; i++) {
        polinomio aux;
        aux.poner_termino(1, 1);
        aux.poner_termino(-x[i], 0);
        polinomio mult;
        mult.multiplicar(PnTerm, aux);
        PnTerm = mult;
    }
    return PnTerm;
}
