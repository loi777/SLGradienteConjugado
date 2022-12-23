// obtem o parametro determinado por parameter_indicator
// apenas para valores inteiros
int get_int_parameter(int argc, char *argv[], char* parameter_indicator);

// obtem o parametro determinado por parameter_indicator
// apenas para valores double
double get_double_parameter(int argc, char *argv[], char* parameter_indicator);

// obtem o parametro determinado por parameter_indicator
// apenas para textos
char* get_string_parameter(int argc, char *argv[], char* parameter_indicator);


// coleta as informacoes calculadas e imprime no padrao do enunciado
void saida_do_programa(FILE* out_file, double tempoPC, double tempoIter, double tempoResid, vetor* iter, double iterk, double residuo, mult_vetor* respostas, int itera);