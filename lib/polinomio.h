//
// Created by rogerroca on 20/10/2022.
//

#ifndef DATASTRUCTURE_POLINOMIOVECTORDOUBLE_H
#define DATASTRUCTURE_POLINOMIOVECTORDOUBLE_H

#include "vector" // Para soportar el uso de vectores como datos.

const int MAX_SIZE_POL_VECT_DOUBLE = 100;

class polinomio {
private:
    double vCoef[MAX_SIZE_POL_VECT_DOUBLE]{};
    int vExp[MAX_SIZE_POL_VECT_DOUBLE]{};
    int length;

    int BuscarExponente(int exp);

    void rmTerm(int pos);

public:
    polinomio();

    bool EsCero();

    int Grado();

    double coeficiente(int exp);

    void AsignarCoeficiente(double coef, int exp);

    void poner_termino(double coef, int exp);

    int numero_terminos();

    int exponente(int term);

    void sumar(polinomio p1, polinomio p2);

    void restar(polinomio p1, polinomio p2);

    void multiplicar(polinomio p1, polinomio p2);

    void opuesto(polinomio p1, polinomio p2);

    void print();

    void derivar(polinomio p1);

    double integrar(double a, double b, double dx);

    double eval(double x);

};


#endif //DATASTRUCTURE_POLINOMIOVECTORDOUBLE_H
