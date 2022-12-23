
// calcula o escalar alpha
//             rT * r
//          ---------------  = a
//            pT * A * p
// P = vetor direcional
// R = vetor residuo
double escalar_alpha(vetor* R, vetor* P, matriz_quad* A);

// calcula o escalar beta
//             r1T * r1
//          ---------------  = b
//             r0T * r0
double escalar_beta(vetor* R0, vetor* R1);



// calcula o residuo inicial sem escalar alpha e salva no vetor residuo
// r = b - A*X0
// R = vetor residuo
void calcula_residuo_inicial(vetor* R, vetor* X0, vetor* B, matriz_quad* A);

// com o escalar alpha podemos calcular um residuo melhor, do residuo anterior
// r1 = r0 - alpha * A * P
void calcula_residuo(vetor* R1, vetor* R0, double alpha, matriz_quad* A, vetor* P);



// calcula o proximo vetor direcional
// p1 = r1 + beta * p0
void calcula_prox_direcional(vetor* P1, vetor* R1, double beta, vetor* P0);

// com o uso do direcional e o escalar alpha calculamos o proximo X
void calcula_prox_x(vetor* X1, vetor* X0, double alpha, vetor* P);



// calcula o erro aproximado absoluto maximo
// max( |xi - xi-1| / |xi|) = e
double calculo_de_erro(mult_vetor* MV, int i);

// calcula o iter k desejado
double calculo_iter_k(mult_vetor* X);

// calcula a norma maxima de multiplo vetor
double max_norma_mult_vetor(mult_vetor* MV);

// calcula a norma euclidiana de um vetor isolado
double norma_euclidiana_vetor(vetor* V);



// de um vetor V, calculamos sua inversa
// assumindo que o vetor V representa uma diagonal principal
void calcula_inversa_diagonal_principal(vetor* V);

// multiplica uma matriz A pela diagonal principal diag
void multiplica_matriz_por_diagonal(matriz_quad* A, vetor* diag);

// multiplica um vetor B pela diagonal principal diag
void multiplica_vetor_por_diagonal(vetor* B, vetor* diag);