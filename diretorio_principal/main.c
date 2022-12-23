#include <stdlib.h>
#include <stdio.h>

#include "matriz_essencial.h"
#include "matriz_math.h"
#include "time_manip.h"
#include "io_handler.h"


// Aluno: Luiz Fernando Giongo dos Santos
// GRR: 20203965




// funcao que roda todas as funcoes criadas
int main(int argc, char *argv[]) {
    srand(20222);

    //------------------------------------------------------------------------
    // inicializacao de variaveis
    //------------------------------------------------------------------------

    double tempTEMPO, tempoIter, tempoPC, tempoResid, iterk, residuo;

    //------------------------------------------------------------------------
    // obtemos os parametros primeiros
    //------------------------------------------------------------------------

    int n, k, p, i;
    double e;
    char* o;
    FILE* output_file;


    n = get_int_parameter(argc, argv, "-n");
    if (n <= 10) {
        fprintf(stderr, "ERRO CRITICO, n deve ser maior que 10\n");
        return 1;
    }


    k = get_int_parameter(argc, argv, "-k");
    if (k <= 1 || k % 2 == 0) {
        fprintf(stderr, "ERRO CRITICO, k deve ser maior que 1 e impar\n");
        return 1;
    }


    p = get_int_parameter(argc, argv, "-p");
    if (p == -1) {
        fprintf(stderr, "ERRO CRITICO, parametro -p deve ser igual ou maior que 0\n");
        return 1;
    }


    i = get_int_parameter(argc, argv, "-i");
    if (i == -1) {
        fprintf(stderr, "ERRO CRITICO, parametro -i obrigatorio nao inserido\n");
        return 1;
    }


    e = get_double_parameter(argc, argv, "-e");
    // e == -1, nao inserido


    o = get_string_parameter(argc, argv, "-o");
    if (o == NULL) {
        fprintf(stderr, "ERRO CRITICO, arquivo de saida nao especificado\n");
        return 1;
    }


    output_file = fopen(o, "w");
    if (output_file == NULL) {
        fprintf(stderr, "ERRO CRITICO, ao abrir arquivo de saida");
        return 3;
    }

    //------------------------------------------------------------------------
    // apos obter os parametros criamos as matrizes
    //------------------------------------------------------------------------

    matriz_quad* A = alloca_matriz_quadrada(n, k);
    if (A == NULL) {
        fprintf(stderr, "ERRO CRITICO, matriz A nao foi alocada\n");
        return 2;
    }
    gera_valores_matriz_quadrada(A);


    vetor* B = alloca_vetor(n);
    if (B == NULL) {
        fprintf(stderr, "ERRO CRITICO, vetor B nao foi alocado\n");
        return 2;
    }
    gera_valores_independentes_vetor(B, A->k);


    mult_vetor* X = alloca_multiplo_vetor(n, i);
    if (X == NULL) {
        fprintf(stderr, "ERRO CRITICO, vetor X nao foi alocado\n");
        return 2;
    }


    mult_vetor* R = alloca_multiplo_vetor(n, i);
    if (R == NULL) {
        fprintf(stderr, "ERRO CRITICO, R nao foi alocado\n");
        return 2;
    }


    mult_vetor* P = alloca_multiplo_vetor(n, i);
    if (P == NULL) {
        fprintf(stderr, "ERRO CRITICO, P nao foi alocado\n");
        return 2;
    }


    vetor* erros = alloca_vetor(n);
    if (erros == NULL) {
        fprintf(stderr, "ERRO CRITICO, vetor erros nao foi alocado\n");
        return 2;
    }

    //------------------------------------------------------------------------
    // caso o parametro tenha sido definido, prepara precondicionador
    //------------------------------------------------------------------------

    tempTEMPO = timestamp();

    if (p > 0) {    // pre-condicionador de Jacobi
        vetor* inversa = obtem_diag_principal(A);
        if (inversa == NULL) {
            fprintf(stderr, "ERRO CRITICO, vetor inversa nao foi alocado\n");
            return 2;
        }

        calcula_inversa_diagonal_principal(inversa);

        multiplica_matriz_por_diagonal(A, inversa);
        multiplica_vetor_por_diagonal(B, inversa);

        libera_vetor(inversa);
    }

    tempoPC = timestamp() - tempTEMPO;

    //------------------------------------------------------------------------
    // tudo preparado, agora podemos fazer os calculos
    //------------------------------------------------------------------------

    double alpha, beta;
    int itera;


    tempTEMPO = timestamp();


    // a primeira iteracao do codigo deve ser especial para contar com a ausencia
    // do direcional P
    // X0 deve ser inicialmente um vetor de 0, per enunciado
    calcula_residuo_inicial(R->vetores[0], X->vetores[0], B, A);
    copia_vetores(R->vetores[0], P->vetores[0]);
    alpha = escalar_alpha(R->vetores[0], P->vetores[0], A);

    tempoIter = timestamp() - tempTEMPO;


    for (itera = 1; itera < i; itera++) {
        tempTEMPO = timestamp();


        calcula_prox_x(X->vetores[itera], X->vetores[itera-1], alpha, P->vetores[itera-1]);
        erros->vetor[itera] = calculo_de_erro(X, itera);

        if (e > 0.0 && erros->vetor[itera] < e) {
            itera++;
            break;
        }

        calcula_residuo(R->vetores[itera], R->vetores[itera-1], alpha, A, P->vetores[itera-1]);
        beta = escalar_beta(R->vetores[itera-1], R->vetores[itera]);
        calcula_prox_direcional(P->vetores[itera], R->vetores[itera], beta, P->vetores[itera-1]);
        alpha = escalar_alpha(R->vetores[itera], P->vetores[itera], A);


        tempoIter += timestamp() - tempTEMPO;
    }

    // para pegar a media do tempo que cada iteracao tomou dividimos por i
    tempoIter /= i;

    //------------------------------------------------------------------------
    // desejamos a norma do residuo agora
    //------------------------------------------------------------------------

    tempTEMPO = timestamp();

    residuo = norma_euclidiana_vetor(R->vetores[i-1]);

    tempoResid = timestamp() - tempTEMPO;

    //------------------------------------------------------------------------
    // calculo do iter k, que o enunciado pede
    //------------------------------------------------------------------------

    iterk = calculo_iter_k(X);

    //------------------------------------------------------------------------
    // calculos finalizamos, agora enviamos a resposta para o arquivo de saida
    //------------------------------------------------------------------------

    saida_do_programa(output_file, tempoPC, tempoIter, tempoResid, erros, iterk, residuo, X, itera);

    //------------------------------------------------------------------------
    // sanitizacao do codigo, liberacao de memorias e etc
    //------------------------------------------------------------------------

    libera_matriz_quadrada(A);
    libera_vetor(B);

    libera_mult_vetor(X);

    libera_mult_vetor(R);

    libera_mult_vetor(P);

    libera_vetor(erros);

    return 0;   // sem erros
}