#include <stdlib.h>
#include <stdio.h>

#include "matriz_essencial.h"


// Aluno: Luiz Fernando Giongo dos Santos
// GRR: 20203965




/***********************
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * i,j: coordenadas do elemento a ser calculado (0<=i,j<n)
 * k: numero de diagonais da matriz A
 ***********************/
inline double generateRandomA( unsigned int i, unsigned int j, unsigned int k ) {
  double invRandMax = 1.0 / (double)RAND_MAX;

  return ( (i==j)?(double)(k<<1) : 1.0 )  * (double)rand() * invRandMax;
}



/***********************
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * k: numero de diagonais da matriz A
 ***********************/
inline double generateRandomB( unsigned int k ) {
  double invRandMax = 1.0 / (double)RAND_MAX;

  return (double)(k<<2) * (double)rand() * invRandMax;
}



//--------------------------------------------------------------------------------



// aloca o espaco na memoria para uma matriz n por n
matriz_quad* alloca_matriz_quadrada(int n, int k) {
  matriz_quad* A = malloc(sizeof(matriz_quad));
  if (A == NULL) {
    return NULL;
  }

  A->matriz = malloc(n * sizeof(double*));
  if (A->matriz == NULL) {
    free(A);
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    A->matriz[i] = malloc(n * sizeof(double));

    if (A->matriz[i] == NULL) {
      free(A);
      free(A->matriz);
      for (int x = i-1; x > 0; x--) {
        free(A->matriz[x]);
      }

      return NULL;
    }
  }

  A->n = n;
  A->k = k;
  return A;
}

// libera o espaco alocado pela matriz A
void libera_matriz_quadrada(matriz_quad* A) {
  for (int i = 0; i < A->n; i++) {
    free(A->matriz[i]);
  }
  free(A->matriz);

  free(A);
}


// usando generateRandomA preencho a matriz quadrada A
void gera_valores_matriz_quadrada(matriz_quad* A) {
  for (int i = 0; i < A->n; i++) {
    for (int j = 0; j < A->n; j++) {
      A->matriz[i][j] = generateRandomA(i, j, A->k);
    }
  }
}



//--------------------------------------------------------------------------------



// aloca o espaco usado por um vetor de tamanho n
vetor* alloca_vetor(int n) {
  vetor* b = malloc(sizeof(vetor));
  if (b == NULL) {
    return NULL;
  }

  b->vetor = malloc(n * sizeof(double));
  if (b->vetor == NULL) {
    free(b);
    return NULL;
  }

  b->n = n;
  return b;
}


// libera o espaco alocado pelo vetor b
void libera_vetor(vetor* b) {
  free(b->vetor);
  free(b);
}


// usando generateRandomB preencho o vetor com os valores independentes
void gera_valores_independentes_vetor(vetor* b, int k) {
  for (int i = 0; i < b->n; i++) {
    b->vetor[i] = generateRandomB(k);
  }
}


// copia um vetor origem em outro destino
void copia_vetores(vetor* origem, vetor* destino) {
  for (int i = 0; i < origem->n; i++) {
    destino->vetor[i] = origem->vetor[i];
  }
}


// coleta a diagonal principal de A
// salva em um vetor
vetor* obtem_diag_principal(matriz_quad* A) {
  vetor* Diag = alloca_vetor(A->n);
  if (Diag == NULL) {
    return NULL;
  }

  for (int i = 0; i < A->n; i++) {
    Diag->vetor[i] = A->matriz[i][i];
  }

  Diag->n = A->n;
  return Diag;
}



//--------------------------------------------------------------------------------



// aloca o espaco usado por um vetor de tamanho n por i iteracoes
mult_vetor* alloca_multiplo_vetor(int n, int iter) {
  mult_vetor* V = malloc(sizeof(mult_vetor));
  if (V == NULL) {
    return NULL;
  }

  V->vetores = malloc(iter * sizeof(vetor*));
  if (V->vetores == NULL) {
    free(V);
    return NULL;
  }

  for (int i = 0; i < iter; i++) {
    V->vetores[i] = alloca_vetor(n);
    if (V->vetores[i] == NULL) {
      for (int x = i-1; x > 0; x--) {
        libera_vetor(V->vetores[x]);
      }
      free(V->vetores);
      free(V);
      return NULL;
    }
  }

  V->i = iter;
  return V;
}


// libera o espaco alocado pelo mult vetor
void libera_mult_vetor(mult_vetor* V) {
  for (int i = 0; i < V->i; i++) {
    libera_vetor(V->vetores[i]);
  }
  free(V->vetores);

  free(V);
}



//--------------------------------------------------------------------------------