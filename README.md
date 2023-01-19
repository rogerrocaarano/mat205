# Algoritmos de resolución para métodos numéricos

## Ecuaciones Lineales

Incluir un archivo matrix.txt con el contenido de la matriz que se desea utilizar junto con el ejecutable compilado.

```
4 -2 -1 9  
5 1 -1 7  
1 2 -1 12
```

Modificar el archivo `ecuaciones-lineales.cpp` con los datos del problema a analizar.

```
int main() {
	Matrix M = Matrix("matrix.txt");
	AB.solveLU();
}
```

- Eliminación de Gauss `M.gaussElim()`
- Resolución por método de Gauss-Seidel `M.solveGaussSeidel(iterations);`
- Resolución por método de Factorización LU `M.solveLU();`
- Resolución por método de Jacobi `M.solveJacobi(iterations);`

## Interpolación

Modificar el archivo `interpolacion-main.cpp`, los métodos para interpolar disponibles están divididos en namespaces

- Interpolación de LaGrange `lagrange::interpolar(x,y)`
- Interpolación de Newton `newton::interpolar(x,y)`

## Integración numérica

Se incluyen los siguientes métodos:

- Integración por Trapecios
- Integración por método de Simpson
  Modificar los datos en los archivos `integracionTrapecios.cpp` e `integracionSimpson.cpp` respectivamente.