#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "matriz_essencial.h"
#include "matriz_math.h"


// Aluno: Luiz Fernando Giongo dos Santos
// GRR: 20203965




// calcula o escalar alpha
//             rT * r
//          ---------------  = a
//            pT * A * p
// P = vetor direcional
// R = vetor residuo
double escalar_alpha(vetor* R, vetor* P, matriz_quad* A) {
    // obtem o valor que vai ser dividido, salvo em a
    double a = 0;
    for (int i = 0; i < R->n; i++) {
        a += R->vetor[i] * R->vetor[i];
    }

    // obtem o valor que vai dividir, salvo em b
    double b = 0;
    double temp2;
    for (int i = 0; i < P->n; i++) {
        temp2 = 0;

        for (int k = 0; k < A->n; k++) {
            temp2 += A->matriz[i][k] * P->vetor[k];
        }

        b += P->vetor[i] * temp2;
    }

    // se b for 0 ele resultara em divisao por 0
    if (b == 0) {
        fprintf(stderr, "ERRO, divisao por 0 ao calcular escalar alpha");
        return 0.0;
    }

    // retorna a divisao
    return a / b;
}


// calcula o escalar beta
//             r1T * r1
//          ---------------  = b
//             r0T * r0
double escalar_beta(vetor* R0, vetor* R1) {
    // obtem o valor que vai ser divido, salvo em a
    double a = 0;
    for (int i = 0; i < R1->n; i++) {
        a += R1->vetor[i] * R1->vetor[i];
    }

    // obtem o valor que vai dividir, salvo em b
    double b = 0;
    for (int i = 0; i < R0->n; i++) {
        b += R0->vetor[i] * R0->vetor[i];
    }

    // se b for igual a 0 dara divisao por 0
    if (b == 0) {
        fprintf(stderr, "ERRO, divisao por 0 ao calcular escalar alpha");
        return 0.0;
    }

    // retorna a divisao
    return a / b;
}



//-------------------------------------------------------------------



// calcula o residuo inicial sem escalar alpha e salva no vetor residuo
// r = b - A*X0
// R = vetor residuo
void calcula_residuo_inicial(vetor* R, vetor* X0, vetor* B, matriz_quad* A) {
    for (int i = 0; i < R->n; i++) {
        R->vetor[i] = B->vetor[i];

        for (int k = 0; k < R->n; k++) {
            R->vetor[i] -= A->matriz[i][k] * X0->vetor[k];
        }
    }
}


// com o escalar alpha podemos calcular um residuo melhor, do residuo anterior
// r1 = r0 - alpha * A * P
void calcula_residuo(vetor* R1, vetor* R0, double alpha, matriz_quad* A, vetor* P) {
    double temp;
    for (int i = 0; i < R1->n; i++) {
        temp = 0;

        for (int k = 0; k < A->n; k++) {
            temp += A->matriz[i][k] * P->vetor[k];
        }

        R1->vetor[i] = R0->vetor[i] - alpha * temp;
    }
}



//-------------------------------------------------------------------



// calcula o proximo vetor direcional
// p1 = r1 + beta * p0
void calcula_prox_direcional(vetor* P1, vetor* R1, double beta, vetor* P0) {
    for (int i =0; i < P1->n; i++) {
        P1->vetor[i] = R1->vetor[i] + beta * P0->vetor[i];
    }
}



//-------------------------------------------------------------------



// com o uso do direcional e o escalar alpha calculamos o proximo X
// x1 = x0 + alpha * P
void calcula_prox_x(vetor* X1, vetor* X0, double alpha, vetor* P) {
    for (int i = 0; i < X1->n; i++) {
        X1->vetor[i] = X0->vetor[i] + alpha * P->vetor[i];
    }
}



//-------------------------------------------------------------------



// calcula o erro aproximado absoluto maximo
// max( |xi - xi-1| / |xi|) = e
double calculo_de_erro(mult_vetor* MV, int i) {
    double maxima = 0.0;
    double temp;

    // caso nao haja um vetor anterior a esse para calcular o erro
    if (i <= 0) {
        return 0.0;
    }

    for (int k = 0; k < MV->vetores[i]->n; k++) {
        // caso o valor do vetor i seja 0 ele dara divisao por zero
        if (MV->vetores[i]->vetor[k] == 0) {
            fprintf(stderr, "ERRO, divisao por 0 ao calcular erro");
            return 0.0;
        }

        temp = fabs(MV->vetores[i]->vetor[k] - MV->vetores[i-1]->vetor[k]) / fabs(MV->vetores[i]->vetor[k]);

        if (maxima < temp) {
            maxima = temp;
        }
    }


    return maxima;
}


// calcula o iter k desejado
// iter k = max(| xi - xi-1 |)
double calculo_iter_k(mult_vetor* X) {
    double maxima = 0.0;
    double temp;

    for (int i = 1; i < X->i; i++) {
        temp = 0.0;

        for (int n = 0; n < X->vetores[i]->n; n++) {
            temp += fabs(X->vetores[i]->vetor[n] - X->vetores[i-1]->vetor[n]);
        }

        if (maxima < temp) {
            maxima = temp;
        }
    }

    return maxima;
}


// calcula a norma maxima de multiplo vetor
double max_norma_mult_vetor(mult_vetor* MV) {
    double maximo = 0.0;
    double temp;

    for (int i = 0; i < MV->i; i++) {
        temp = 0.0;

        for (int k = 0; k < MV->vetores[i]->n; k ++) {
            temp += fabs(MV->vetores[i]->vetor[k]);
        }

        if (maximo < temp) {
            maximo = temp;
        }
    }

    return maximo;
}


// calcula a norma euclidiana de um vetor isolado
double norma_euclidiana_vetor(vetor* V) {
    double val = 0.0;

    for (int i = 0; i < V->n; i ++) {
        val += fabs(V->vetor[i]);
    }

    return val;
}



//--------------------------------------------------------



// de um vetor V, calculamos sua inversa
// assumindo que o vetor V representa uma diagonal principal
void calcula_inversa_diagonal_principal(vetor* V) {
    for (int i = 0; i < V->n; i++) {
        if (V->vetor[i] == 0) {
            continue;
        }

        V->vetor[i] = 1/V->vetor[i];
    }
}


// multiplica uma matriz A pela diagonal principal diag
void multiplica_matriz_por_diagonal(matriz_quad* A, vetor* diag) {
    for (int i = 0; i < A->n; i++) {
        for (int k = 0; k < A->n; k++) {
            A->matriz[i][k] *= diag->vetor[i];
        }
    }
}


// multiplica um vetor B pela diagonal principal diag
void multiplica_vetor_por_diagonal(vetor* B, vetor* diag) {
    for (int i = 0; i < B->n; i++) {
        B->vetor[i] *= diag->vetor[i];
    }
}