//
// Created by rogerroca on 20/10/2022.
//

#include "polinomio.h"
#include "iostream"
#include "cmath"

using namespace std;

polinomio::polinomio() {
    length = 0;
    for (int i = 0; i < MAX_SIZE_POL_VECT_DOUBLE; i++) {
        vExp[i] = 0;
        vCoef[i] = 0;
    }
}

bool polinomio::EsCero() {
    return length == 0;
}

int polinomio::Grado() {
    if (this->length > 0) {
        int grade = vExp[0];
        for (int i = 0; i <= length - 1; i++) {
            if (vExp[i] > grade) grade = vExp[i];
        }
        return grade;
    } else cout << "Polynomial is empty" << endl;
}

double polinomio::coeficiente(int exp) {
    int pos = BuscarExponente(exp);
    if (pos != -1) {
        return vCoef[pos];
    } else cout << "Exponent not found." << endl;
}

void polinomio::AsignarCoeficiente(double coef, int exp) {
    int pos = BuscarExponente(exp);
    if (pos != -1) {
        if (coef == 0) {
            rmTerm(pos);
        } else vCoef[pos] = coef;
    }
}

void polinomio::poner_termino(double coef, int exp) {
    int pos = BuscarExponente(exp);
    if (pos != -1) {
        vCoef[pos] = vCoef[pos] + coef;
        if (vCoef[pos] == 0) rmTerm(pos);
    } else {
        if (length < MAX_SIZE_POL_VECT_DOUBLE) {
            vExp[length] = exp;
            vCoef[length] = coef;
            length++;
        } else cout << "Not enough space for new term." << endl;
    }
}

int polinomio::numero_terminos() {
    return length;
}

int polinomio::exponente(int term) {
    if (term > 0 && term <= length) {
        return vExp[term - 1];
    } else cout << "Invalid position." << endl;
}

void polinomio::sumar(polinomio p1, polinomio p2) {
    for (int i = 1; i <= p1.numero_terminos(); i++) {
        int exp = p1.exponente(i);
        double coef = p1.coeficiente(exp);
        poner_termino(coef, exp);
    }
    for (int i = 1; i <= p2.numero_terminos(); i++) {
        int exp = p2.exponente(i);
        double coef = p2.coeficiente(exp);
        poner_termino(coef, exp);
    }
}

void polinomio::restar(polinomio p1, polinomio p2) {
    for (int i = 1; i <= p1.numero_terminos(); i++) {
        int exp = p1.exponente(i);
        double coef = p1.coeficiente(exp);
        poner_termino(coef, exp);
    }
    for (int i = 1; i <= p1.numero_terminos(); i++) {
        int exp = p2.exponente(i);
        double coef = p2.coeficiente(exp) * (-1);
        poner_termino(coef, exp);
    }
}

void polinomio::multiplicar(polinomio p1, polinomio p2) {
    for (int i = 1; i <= p1.numero_terminos(); i++) {
        for (int j = 1; j <= p2.numero_terminos(); j++) {
            int exp = p1.exponente(i) + p2.exponente(j);
            double coef = p1.coeficiente(p1.exponente(i))
                          * p2.coeficiente(p2.exponente(j));
            poner_termino(coef, exp);
        }
    }
}

void polinomio::opuesto(polinomio p1, polinomio p2) {
    sumar(p1, p2);
    if (EsCero()) {
        cout << "Son polinomios opuestos." << endl;
    } else {
        cout << "No son polinomios opuestos." << endl;
    }
}

void polinomio::print() {
    if (EsCero()) {
        cout << "0" << endl;
        return;
    }
    for (int i = 0; i < length; i++) {
        if (i == 0) {
            cout << vCoef[i];
            if (vExp[i] != 0) {
                cout << "*x^" << vExp[i];
            }
        } else {
            if (vCoef[i] > 0) {
                cout << "+" << vCoef[i];
            } else {
                cout << vCoef[i];
            }
            if (vExp[i] != 0) {
                cout << "*x^" << vExp[i];
            }
        }
    }
    cout << endl;
}

void polinomio::derivar(polinomio p1) {
    for (int i = 1; i <= p1.numero_terminos(); i++) {
        if (p1.exponente(i) != 0) {
            int exp = p1.exponente(i) - 1;
            int coef = p1.coeficiente(p1.exponente(i))
                       * p1.exponente(i);
            poner_termino(coef, exp);
        }
    }
}

int polinomio::BuscarExponente(int exp) {
    int pos = -1;
    int i = 0;
    while (pos == -1 && i <= length - 1) {
        if (vExp[i] == exp) {
            pos = i;
        }
        i++;
    }
    return pos;
}

void polinomio::rmTerm(int pos) {
    if (pos < length) {
        for (int i = pos; i < length - 1; i++) {
            vExp[i] = vExp[i + 1];
            vCoef[i] = vCoef[i + 1];
        }
        vExp[length - 1] = 0;
        vCoef[length - 1] = 0;
        length--;
    } else cout << "Invalid direction" << endl;
}

double polinomio::integrar(double a, double b, double dx) {
    if (b > a) {
        double area = 0;
        double y;
        while ((a + dx) < b) {
            y = eval((a + dx) / 2);
            y = abs(y);
            area = area + y * dx;
            a = a + dx;
        }
        y = eval((b - a) / 2);
        y = abs(y);
        area = area + y * ((b - a) / 2);
        return area;
    } else cout << "Invalid range.";
}

double polinomio::eval(double x) {
    double y = 0;
    for (int i = 1; i <= length; i++) {
        int exp = exponente(i);
        double coef = coeficiente(exp);
        double term = coef * pow(x, exp);
        y = y + term;
    }
    return y;
}
