//
//  matrix.cpp
//
//  Created by Furkanicus on 12/04/15.
//  Copyright (c) 2015 Furkan. All rights reserved.
//

#include "matrix.h"

using namespace std;

// Constructor for Any Matrix
Matrix::Matrix(unsigned rowSize, unsigned colSize, double initial) {
    m_rowSize = rowSize;
    m_colSize = colSize;
    m_matrix.resize(rowSize);
    for (unsigned i = 0; i < m_matrix.size(); i++) {
        m_matrix[i].resize(colSize, initial);
    }
}

// Constructor for Given Matrix
Matrix::Matrix(const char *fileName) {
    ifstream file_A(fileName); // input file stream to open the file A.txt

    // Task 1
    // Keeps track of the Column and row sizes
    int colSize = 0;
    int rowSize = 0;

    // read it as a vector
    string line_A;
    int idx = 0;
    double element_A;
    double *vector_A = nullptr;

    if (file_A.is_open() && file_A.good()) {
        // cout << "File A.txt is open. \n";
        while (getline(file_A, line_A)) {
            rowSize += 1;
            stringstream stream_A(line_A);
            colSize = 0;
            while (1) {
                stream_A >> element_A;
                if (!stream_A)
                    break;
                colSize += 1;
                double *tempArr = new double[idx + 1];
                copy(vector_A, vector_A + idx, tempArr);
                tempArr[idx] = element_A;
                vector_A = tempArr;
                idx += 1;
            }
        }
    } else {
        cout << " WTF! failed to open. \n";
    }

    int j;
    idx = 0;
    m_matrix.resize(rowSize);
    for (unsigned i = 0; i < m_matrix.size(); i++) {
        m_matrix[i].resize(colSize);
    }
    for (int i = 0; i < rowSize; i++) {
        for (j = 0; j < colSize; j++) {
            this->m_matrix[i][j] = vector_A[idx];
            idx++;
        }
    }
    m_colSize = colSize;
    m_rowSize = rowSize;
    delete[] vector_A; // Tying up loose ends


}

// Copy Constructor
Matrix::Matrix(const Matrix &B) {
    this->m_colSize = B.getCols();
    this->m_rowSize = B.getRows();
    this->m_matrix = B.m_matrix;

}

Matrix::~Matrix() {

}

// Addition of Two Matrices
Matrix Matrix::operator+(Matrix &B) {
    Matrix sum(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            sum(i, j) = this->m_matrix[i][j] + B(i, j);
        }
    }
    return sum;
}

// Subtraction of Two Matrices
Matrix Matrix::operator-(Matrix &B) {
    Matrix diff(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            diff(i, j) = this->m_matrix[i][j] - B(i, j);
        }
    }

    return diff;
}

// Multiplication of Two Matrices
Matrix Matrix::operator*(Matrix &B) {
    Matrix multip(m_rowSize, B.getCols(), 0.0);
    if (m_colSize == B.getRows()) {
        unsigned i, j, k;
        double temp = 0.0;
        for (i = 0; i < m_rowSize; i++) {
            for (j = 0; j < B.getCols(); j++) {
                temp = 0.0;
                for (k = 0; k < m_colSize; k++) {
                    temp += m_matrix[i][k] * B(k, j);
                }
                multip(i, j) = temp;
                //cout << multip(i,j) << " ";
            }
            //cout << endl;
        }
        return multip;
    } else {
        return "Error";
    }
}

// Scalar Addition
Matrix Matrix::operator+(double scalar) {
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            result(i, j) = this->m_matrix[i][j] + scalar;
        }
    }
    return result;
}

// Scalar Subraction
Matrix Matrix::operator-(double scalar) {
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            result(i, j) = this->m_matrix[i][j] - scalar;
        }
    }
    return result;
}

// Scalar Multiplication
Matrix Matrix::operator*(double scalar) {
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            result(i, j) = this->m_matrix[i][j] * scalar;
        }
    }
    return result;
}

// Scalar Division
Matrix Matrix::operator/(double scalar) {
    Matrix result(m_rowSize, m_colSize, 0.0);
    unsigned i, j;
    for (i = 0; i < m_rowSize; i++) {
        for (j = 0; j < m_colSize; j++) {
            result(i, j) = this->m_matrix[i][j] / scalar;
        }
    }
    return result;
}


// Returns value of given location when asked in the form A(x,y)
double &Matrix::operator()(const unsigned &rowNo, const unsigned &colNo) {
    return this->m_matrix[rowNo][colNo];
}

// No brainer - returns row #
unsigned Matrix::getRows() const {
    return this->m_rowSize;
}

// returns col #
unsigned Matrix::getCols() const {
    return this->m_colSize;
}

// Take any given matrices transpose and returns another matrix
Matrix Matrix::transpose() {
    Matrix Transpose(m_colSize, m_rowSize, 0.0);
    for (unsigned i = 0; i < m_colSize; i++) {
        for (unsigned j = 0; j < m_rowSize; j++) {
            Transpose(i, j) = this->m_matrix[j][i];
        }
    }
    return Transpose;
}

// Prints the matrix beautifully
void Matrix::print() const {
    for (unsigned i = 0; i < m_rowSize; i++) {
        for (unsigned j = 0; j < m_colSize; j++) {
            cout << "[" << m_matrix[i][j] << "] ";
        }
        cout << endl;
    }
}

// Returns 3 values
//First: Eigen Vector
//Second: Eigen Value
//Third: Flag
tuple<Matrix, double, int> Matrix::powerIter(unsigned rowNum, double tolerance) {
    // Picks a classic X vector
    Matrix X(rowNum, 1, 1.0);
    // Initiates X vector with values 1,2,3,4
    for (unsigned i = 1; i <= rowNum; i++) {
        X(i - 1, 0) = i;
    }
    int errorCode = 0;
    double difference = 1.0; // Initiall value greater than tolerance
    unsigned j = 0;
    unsigned location;
    // Defined to find the value between last two eigen values
    vector<double> eigen;
    double eigenvalue = 0.0;
    eigen.push_back(0.0);

    while (abs(difference) > tolerance) // breaks out when reached tolerance
    {
        j++;
        // Normalize X vector with infinite norm
        for (int i = 0; i < rowNum; ++i) {
            eigenvalue = X(0, 0);
            if (abs(X(i, 0)) >= abs(eigenvalue)) {
                // Take the value of the infinite norm as your eigenvalue
                eigenvalue = X(i, 0);
                location = i;
            }
        }
        if (j >= 5e5) {
            cout << "Oops, that was a nasty complex number wasn't it?" << endl;
            cout << "ERROR! Returning code black, code black!";
            errorCode = -1;
            return make_tuple(X, 0.0, errorCode);
        }
        eigen.push_back(eigenvalue);
        difference = eigen[j] - eigen[j - 1];
        // Normalize X vector with its infinite norm
        X = X / eigenvalue;

        // Multiply The matrix with X vector
        X = (*this) * X;
    }

    // Take the X vector and what you've found is an eigenvector!
    X = X / eigenvalue;
    return make_tuple(X, eigenvalue, errorCode);
}

Matrix Matrix::deflation(Matrix &X, double &eigenvalue) {
    // Deflation formula exactly applied
    double denominator = eigenvalue / (X.transpose() * X)(0, 0);
    Matrix Xtrans = X.transpose();
    Matrix RHS = (X * Xtrans);
    Matrix RHS2 = RHS * denominator;
    Matrix A2 = *this - RHS2;
    return A2;
}

// Ampliado por rogerroca para la asignatura de métodos numéricos
// 12/12/2022

void Matrix::multRow(unsigned int row, double mult) {
    for (int i = 0; i < m_colSize; i++) {
        this->m_matrix[row][i] = this->m_matrix[row][i] * mult;
    }
}

void Matrix::divRow(unsigned int row, double div) {
    for (int i = 0; i < m_colSize; i++) {
        this->m_matrix[row][i] = this->m_matrix[row][i] / div;
    }
}

void Matrix::moveRows(unsigned int row1, unsigned int row2) {
    if (row1 <= m_rowSize && row2 <= m_colSize && row1 != row2) {
        for (int i = 0; i < m_colSize; i++) {
            float aux = this->m_matrix[row1][i];
            this->m_matrix[row1][i] = this->m_matrix[row2][i];
            this->m_matrix[row2][i] = aux;
        }
    }
}

void Matrix::zeroCol(unsigned int pivotRow, unsigned int pivotCol, unsigned int row, bool rowDown) {
    if (row >= m_rowSize || row < 0)
        return;

    if (row != pivotRow) {
        double pivot = this->m_matrix[pivotRow][pivotCol];
        double div = pivot != 1 ? pivot : 1;
        double x = this->m_matrix[row][pivotCol];
        double mult = x / div;
        for (int i = 0; i < m_colSize; i++) {
            double c1 = this->m_matrix[pivotRow][i];
            double c2 = this->m_matrix[row][i];
            this->m_matrix[row][i] = c2 - mult * c1;
        }
    }

    if (rowDown)
        Matrix::zeroCol(pivotRow, pivotCol, row + 1, rowDown);
    else
        Matrix::zeroCol(pivotRow, pivotCol, row - 1, rowDown);
}

void Matrix::gaussElim() {
    std::cout << "ELIMINACION GAUSSIANA" << std::endl;
    if (m_colSize >= m_rowSize + 1) {
        for (int i = 0; i < m_colSize - 1; i++) {
            sort(i, i);
            if (m_matrix[i][i] != 1)
                divRow(i, m_matrix[i][i]);
            zeroCol(i, i, i, true);
            print();
        }
    }
}

void Matrix::sort(unsigned int row, unsigned int col) {
    if (m_matrix[row][col] != 0)
        return;

    unsigned int testRow = row + 1;
    while (testRow < m_rowSize) {
        if (m_matrix[testRow][col] != 0) {
            moveRows(row, testRow);
        }
        testRow++;
    }
}

Matrix Matrix::gaussSeidelIter() {
    // Inicializar una matriz con la iteración inicial de resultados.
    Matrix result(1, m_rowSize, 0);
    // Preparar la matriz, asumiendo que es diagonal dominante
    for (int i = 0; i < m_colSize - 1; i++) {
        // obtener Ai,i:
        double div = m_matrix[i][i];
        // dividir toda la fila entre div:
        divRow(i, div);
    }
    int var = 0;
    while (var < m_rowSize) {
        for (int col = 0; col < m_colSize; col++) {
            if (col != var) {
                if (col == m_colSize - 1) {
                    result(0, var) = result(0, var) + m_matrix[var][col];
                } else {
                    result(0, var) = result(0, var) - m_matrix[var][col] * result(0, col);
                }
            }
        }
        var++;
    }
    return result;
}

Matrix Matrix::gaussSeidelIter(Matrix &iterResult) {
    Matrix result(1, m_rowSize, 0);
    int var = 0;
    while (var < m_rowSize) {
        for (int col = 0; col < m_colSize; col++) {
            if (col != var) {
                if (col == m_colSize - 1) {
                    result(0, var) = result(0, var) + m_matrix[var][col];
                } else {
                    double mult = result(0, col) == 0 ? iterResult(0, col) : result(0, col);
                    result(0, var) = result(0, var) - m_matrix[var][col] * mult;
                }
            }
        }
        var++;
    }
    return result;
}

void Matrix::solveGaussSeidel(unsigned int iterations) {
    cout << "RESOLUCION DEL SISTEMA POR EL METODO DE GAUSS-SEIDEL" << endl;
    print();
    int i = 1;
    cout << "Iter " << i << endl;
    Matrix previousIter = gaussSeidelIter();
    previousIter.print();
    i++;
    while (i <= iterations) {
        cout << "Iter " << i << endl;
        Matrix newIter = gaussSeidelIter(previousIter);
        newIter.print();
        Er(previousIter, newIter);
        previousIter = newIter;
        i++;
    }
}

void Matrix::Er(Matrix initialIter, Matrix newIter) {
    if (initialIter.m_colSize != newIter.m_colSize || initialIter.m_rowSize != newIter.m_rowSize) {
        cout << "Error: Las matrices no son iguales." << endl;
        return;
    }
    if (initialIter.m_rowSize > 1) {
        cout << "Error: Matrices de tamaño incorrecto." << endl;
        return;
    }
    cout << "Er:";
    for (int i = 0; i < initialIter.m_colSize; i++) {
        double Er = (newIter(0, i) - initialIter(0, i)) / newIter(0, i);
        cout << " [" << Er << "]";
    }
    cout << endl;
}

void Matrix::factLU(Matrix &L, Matrix &U) {
    Matrix A = getAfromAB();
    U = A;
    // La diagonal principal de la matriz L debe ser 1:
    for (int i = 0; i < L.m_rowSize; i++) {
        L(i, i) = 1;
    }
    cout << "La diagonal principal de la matriz L debe ser 1:" << endl;
    L.print();
    for (int i = 0; i < m_colSize - 1; i++) {
        double divisor = U(i, i);
        for (int i_L = i + 1; i_L < L.m_rowSize; i_L++) {
            double dividendo = U(i_L, i);
            L(i_L, i) = dividendo / divisor;
            cout << "---------------------------------------" << endl;
            L.print();
        }
        U.zeroCol(i, i, i, true);
        cout << "---------------------------------------" << endl;
        U.print();
    }
    cout << "FACTORIZACION LU" << endl;
    cout << "-----MATRIZ L-----" << endl;
    L.print();
    cout << "-----MATRIZ U-----" << endl;
    U.print();
}

Matrix Matrix::getAfromAB() {
    Matrix A(m_rowSize, m_colSize - 1, 0);
    for (int i = 0; i < m_rowSize; i++) {
        for (int j = 0; j < m_colSize - 1; j++) {
            A(i, j) = m_matrix[i][j];
        }
    }
    return A;
}

Matrix Matrix::getBfromAB() {
    Matrix B(m_rowSize, 1, 0);
    for (int i = 0; i < m_rowSize; i++) {
        B(i, 0) = m_matrix[i][m_colSize - 1];
    }
    return B;
}

void Matrix::solveLU() {
    Matrix B = getBfromAB();
    Matrix U = Matrix(getRows(), getRows(), 0);
    Matrix L = Matrix(getRows(), getRows(), 0);
    factLU(L, U);
// Eliminación hacia adelante:
    cout << "Eliminacion hacia adelante" << endl;
    for (int i = 1; i < B.m_rowSize; i++) {
        for (int j = 0; j < i; j++) {
            B(i, 0) = B(i, 0) - (L(i, j) * B(j, 0));
        }
    }
    B.print();
    // Eliminación hacia atras:
    cout << "Eliminacion hacia atras" << endl;
    B(B.m_rowSize - 1, 0) = B(B.m_rowSize - 1, 0) / U(U.m_rowSize - 1, U.m_colSize - 1);
    for (int i = B.m_rowSize - 2; i >= 0; i--) {
        for (int j = U.m_colSize - 1; j >= i; j--) {
            if (j == i) {
                B(i, 0) = B(i, 0) / U(i, j);
            } else {
                B(i, 0) = B(i, 0) - (U(i, j) * B(j, 0));
            }
        }
    }
    B.print();
}

Matrix Matrix::jacobiIter() {
    // Inicializar una matriz con la iteración inicial de resultados.
    Matrix seed(1, m_rowSize, 0);
    Matrix result(1, m_rowSize, 0);
    // Preparar la matriz, asumiendo que es diagonal dominante
    for (int i = 0; i < m_colSize - 1; i++) {
        // obtener Ai,i:
        double div = m_matrix[i][i];
        // dividir toda la fila entre div:
        divRow(i, div);
    }
    int var = 0;
    while (var < m_rowSize) {
        for (int col = 0; col < m_colSize; col++) {
            if (col != var) {
                if (col == m_colSize - 1) {
                    result(0, var) = result(0, var) + m_matrix[var][col];
                } else {
                    result(0, var) = result(0, var) - m_matrix[var][col] * seed(0, col);
                }
            }
        }
        var++;
    }
    return result;
}

Matrix Matrix::jacobiIter(Matrix &iterResult) {
    Matrix result(1, m_rowSize, 0);
    int var = 0;
    while (var < m_rowSize) {
        for (int col = 0; col < m_colSize; col++) {
            if (col != var) {
                if (col == m_colSize - 1) {
                    result(0, var) = result(0, var) + m_matrix[var][col];
                } else {
                    result(0, var) = result(0, var) - m_matrix[var][col] * iterResult(0, col);
                }
            }
        }
        var++;
    }
    return result;
}

void Matrix::solveJacobi(unsigned int iterations) {
    cout << "RESOLUCION DEL SISTEMA POR EL METODO DE JACOBI" << endl;
    print();
    int i = 1;
    cout << "Iter " << i << endl;
    Matrix previousIter = jacobiIter();
    previousIter.print();
    i++;
    while (i <= iterations) {
        cout << "Iter " << i << endl;
        Matrix newIter = jacobiIter(previousIter);
        newIter.print();
        Er(previousIter, newIter);
        previousIter = newIter;
        i++;
    }
}











