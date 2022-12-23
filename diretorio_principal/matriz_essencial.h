typedef struct matriz_quad_struct {
    int n;
    int k;
    double** matriz;
} matriz_quad;

typedef struct vetorial_struct {
    int n;
    double* vetor;
} vetor;

typedef struct mult_vetorial_struct {
    int i;
    vetor** vetores;
} mult_vetor;






// aloca o espaco na memoria para uma matriz n por n
matriz_quad* alloca_matriz_quadrada(int n, int k);

// libera o espaco alocado pela matriz A
void libera_matriz_quadrada(matriz_quad* A);

// usando generateRandomA preencho a matriz quadrada A
void gera_valores_matriz_quadrada(matriz_quad* A);




// aloca o espaco usado por um vetor de tamanho n
vetor* alloca_vetor(int n);

// libera o espaco alocado pelo vetor b
void libera_vetor(vetor* b);

// usando generateRandomB preencho o vetor com os valores independentes
void gera_valores_independentes_vetor(vetor* b, int k);

// copia um vetor em outro
void copia_vetores(vetor* origem, vetor* destino);

// coleta a diagonal principal de A
// salva em um vetor
vetor* obtem_diag_principal(matriz_quad* A);




// aloca o espaco usado por um vetor de tamanho n por i iteracoes
mult_vetor* alloca_multiplo_vetor(int n, int iter);

// libera o espaco alocado pelo mult vetor
void libera_mult_vetor(mult_vetor* V);




// gera os valores para usar na matriz A
double generateRandomA( unsigned int i, unsigned int j, unsigned int k );

// gera os valores para usar na matriz B
double generateRandomB( unsigned int k );