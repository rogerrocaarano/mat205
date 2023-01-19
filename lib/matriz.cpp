//
// Created by rogerroca on 14/12/2022.
//

#include <iostream>
#include "matriz.h"

matriz::matriz() {
    df = 0;
    dc = 0;
    repe = 0;
    nt = 0;
}

void matriz::dimensionar(int n, int m) {
    df = n;
    dc = m;
}

int matriz::dimension_fila() {
    return df;
}

int matriz::dimension_columna() {
    return dc;
}

void matriz::poner(int f, int c, matriz::DATA_TYPE valor) {
    int pos = obtenerPosicion(f, c);
    if (pos >= 0) {
        // Asignamos el valor y si el valor es igual a "repe" hay que eliminarlo:
        Vd[pos] = valor;
        if (valor == repe) {
            eliminar(pos);
        }
        return;
    }
    // No existe el valor en la matriz.
    // Crear un nuevo elemento:
    if (nt < MAX_VECTOR_SIZE) {
        Vd[nt] = valor;
        Vf[nt] = f;
        Vc[nt] = c;
        nt++;
    } else std::cout << "Error: No existe suficiente espacio." << std::endl;
}

matriz::DATA_TYPE matriz::elemento(int f, int c) {
    if (f < 0 || f > df) {
        std::cout << "Error: Fila fuera de rango." << std::endl;
        return repe;
    }
    if (c < 0 || c > dc) {
        std::cout << "Error: Columna fuera de rango." << std::endl;
        return repe;
    }
    int pos = obtenerPosicion(f, c);
    return pos == -1 ? repe : Vd[pos];
}

void matriz::definir_valor_repetido(matriz::DATA_TYPE valor) {
    // En tiempo de ejecución se eliminaría de los vectores de la matriz
    // los valores coincidentes con el nuevo valor para "repe"
    int i = 0;
    while (i < nt) {
        if (Vd[i] == valor) {
            eliminar(i);
            break;
        }
        i++;
    }
    // definir el nuevo valor para "repe"
    repe = valor;
}

int matriz::obtenerPosicion(int f, int c) {
    for (int i = 0; i < nt; i++) {
        if (Vf[i] == f && Vc[i] == c) {
            return i;
        }
    }
    return -1;
}

void matriz::eliminar(int pos) {
    if (pos >= nt || pos < 0) {
        std::cout << "Error: Posicion fuera de rango" << std::endl;
        return;
    }
    // Desplazar los valores siguientes del vector.
    for (int i = pos + 1; i < nt - 1; i++) {
        Vf[i - 1] = Vf[i];
        Vc[i - 1] = Vc[i];
        Vd[i - 1] = Vd[i];
    }
    // limpiar el último valor:
    Vf[nt - 1] = repe;
    Vc[nt - 1] = repe;
    Vd[nt - 1] = repe;
    // disminuir "nt".
    nt--;
}

void matriz::mostrar() {
    for (int i = 1; i <= df; i++) {
        for (int j = 1; j <= dc; j++) {
            std::cout << "[" << elemento(i, j) << "]";
        }
        std::cout << std::endl;
    }

}