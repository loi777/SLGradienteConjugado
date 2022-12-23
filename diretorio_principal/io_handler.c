#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "matriz_essencial.h"


// Aluno: Luiz Fernando Giongo dos Santos
// GRR: 20203965




// obtem o parametro determinado por parameter_indicator
// apenas para valores inteiros
int get_int_parameter(int argc, char *argv[], char* parameter_indicator) {
    
    // percorre todos os parametros procurando o indicador
    for (int i = 0; i < argc; i ++) {
        if (!strcmp(argv[i], parameter_indicator)) {

            // achado o indicador retorna ele
            if (i+1 < argc) {
                return (int) strtol(argv[i+1], NULL, 10);
            } else {
                // se o indicador for a ultima coisa ele retorna erro
                return -1;  // valor de erro
            }
        }
    }
    
    return -1;  // valor de erro
}


// obtem o parametro determinado por parameter_indicator
// apenas para valores double
double get_double_parameter(int argc, char *argv[], char* parameter_indicator) {
    
    // percorre todos os parametros procurando o indicador
    for (int i = 0; i < argc; i ++) {
        if (!strcmp(argv[i], parameter_indicator)) {

            // achado o indicador retorna ele
            if (i+1 < argc) {
                return strtod(argv[i+1], NULL);
            } else {
                // se o indicador for a ultima coisa ele retorna erro
                return -1;  // valor de erro
            }
        }
    }
    
    return -1;  // valor de erro
}


// obtem o parametro determinado por parameter_indicator
// apenas para textos
char* get_string_parameter(int argc, char *argv[], char* parameter_indicator) {

    // percorre todos os parametros procurando o indicador
    for (int i = 0; i < argc; i ++) {

        // achado o indicador retorna ele
        if (!strcmp(argv[i], parameter_indicator)) {
            if (i+1 < argc) {
                return argv[i+1];
            } else {
                // se o indicador for a ultima coisa ele retorna erro
                return NULL;    // valor de erro
            }
        }
    }

    return NULL;  // valor de erro
}


//-------------------------------------------------------------------


/*
iter k: norma máxima do erro aproximado em x após a k-ésima iteração (max|xi - xi-1|);
*/

// coleta as informacoes calculadas e imprime no padrao do enunciado
void saida_do_programa(FILE* out_file, double tempoPC, double tempoIter, double tempoResid, vetor* iter, double iterk, double residuo, mult_vetor* respostas, int itera) {

    fprintf(out_file, "# lfgs20 Luiz Fernando Giongo dos Santos\n#\n");

    for (int i = 0; i < itera; i++) {
        fprintf(out_file, "# iter %d: %.15g\n", i, iter->vetor[i]);
    }
    fprintf(out_file, "# iter k: %.15g\n", iterk);
    
    fprintf(out_file, "# residuo: %.15g\n", residuo);

    fprintf(out_file, "# Tempo PC: %.15g\n", tempoPC);
    fprintf(out_file, "# Tempo iter: %.15g\n", tempoIter);
    fprintf(out_file, "# Tempo residuo: %.15g\n#\n", tempoResid);

    // imprime o vetor do ultimo chute, da resposta
    fprintf(out_file, "%d\n", respostas->vetores[itera-1]->n);
    for (int i = 0; i < respostas->vetores[itera-1]->n; i++) {
        fprintf(out_file, "%.15g ", respostas->vetores[itera-1]->vetor[i]);
    }
}