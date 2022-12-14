//
//  matrix.h
//
//  Created by Furkanicus on 12/04/15.
//  Copyright (c) 2015 Furkan. All rights reserved.
// Learnrope.com

#ifndef __EE_242_Project_2__matrix__
#define __EE_242_Project_2__matrix__

#include <stdio.h>
#include <fstream> // for file access
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

using std::vector;
using std::tuple;

class Matrix {
private:
    unsigned m_rowSize;
    unsigned m_colSize;
    vector<vector<double> > m_matrix;
public:
    Matrix(unsigned, unsigned, double);

    Matrix(const char *);

    Matrix(const Matrix &);

    ~Matrix();

    // Matrix Operations
    Matrix operator+(Matrix &);

    Matrix operator-(Matrix &);

    Matrix operator*(Matrix &);

    Matrix transpose();

    // Scalar Operations
    Matrix operator+(double);

    Matrix operator-(double);

    Matrix operator*(double);

    Matrix operator/(double);

    // Aesthetic Methods
    double &operator()(const unsigned &, const unsigned &);

    void print() const;

    unsigned getRows() const;

    unsigned getCols() const;

    // Power Iteration
    tuple<Matrix, double, int> powerIter(unsigned, double);

    // Deflation
    Matrix deflation(Matrix &, double &);

    // Ampliado por rogerroca para la asignatura de métodos numéricos
    // 12/12/2022

    void multRow(unsigned int row, double mult);

    void divRow(unsigned int row, double div);

    void moveRows(unsigned int row1, unsigned int row2);

    void zeroCol(unsigned int pivotRow, unsigned int pivotCol, unsigned int row, bool rowDown);

    void gaussElim();

    void sort(unsigned int row, unsigned int col);

    Matrix getMatrixU(Matrix &A);

    Matrix getMatrixL(Matrix &A, Matrix &U);

    Matrix gaussSeidelIter(Matrix &A);

    Matrix gaussSeidelIter(Matrix &A, Matrix &iterResult);

    void gaussSeidel(Matrix &A, unsigned int iterations);

    void Er(Matrix initialIter, Matrix newIter);


};

#endif /* defined(__EE_242_Project_2__matrix__) */
