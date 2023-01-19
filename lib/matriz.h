//
// Created by rogerroca on 14/12/2022.
//

#ifndef INF220_MATRIZDISPERSAVECTOR_H
#define INF220_MATRIZDISPERSAVECTOR_H

const int MAX_VECTOR_SIZE = 256;

class matriz {
    typedef double DATA_TYPE;

private:
    int Vf[MAX_VECTOR_SIZE]; // Filas
    int Vc[MAX_VECTOR_SIZE]; // Columnas
    DATA_TYPE Vd[MAX_VECTOR_SIZE];  // elementos
    int df; // Dimensión filas
    int dc; // Dimensión columnas
    int nt; // Términos
    DATA_TYPE repe; // es el elemento que se repetirá en la matriz

    int obtenerPosicion(int f, int c);

    void eliminar(int pos);

public:
    matriz();

    void dimensionar(int m, int n);

    int dimension_fila();

    int dimension_columna();

    void poner(int f, int c, DATA_TYPE valor);

    DATA_TYPE elemento(int f, int c);

    void definir_valor_repetido(DATA_TYPE valor);

    void mostrar();
};


#endif //INF220_MATRIZDISPERSAVECTOR_H
