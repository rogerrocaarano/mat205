// Librería de interpolación polinómica
// para la asignatura MAT205 - Métodos numéricos.
// Incluir este archivo y llamar a las funciones según su namespace.
//
// Created by rogerroca on 27/12/2022.
//

#ifndef MAT205_INTERPOLACION_INTERPOLACION_H
#define MAT205_INTERPOLACION_INTERPOLACION_H


#include "polinomio.h"
#include "matriz.h"
#include "vector"
#include "iostream"

namespace lagrange {
    polinomio interpolar(std::vector<double> &x, std::vector<double> &y);

    polinomio Li(int k, std::vector<double> x);
}
namespace newton {
    polinomio interpolar(std::vector<double> &x, std::vector<double> &y);

    double difDividida(double x0, double x1, double y0, double y1);

    void cargarFila_difDividida(matriz &m, std::vector<double> x, int fila);

    polinomio PnTerm(matriz &m, std::vector<double> x, int grado);

}
#endif //MAT205_INTERPOLACION_INTERPOLACION_H
