#include <time.h>


// Aluno: Luiz Fernando Giongo dos Santos
// GRR: 20203965




// obtem o tempo no momento
double timestamp(void) {
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);

  return ((double)(tp.tv_sec + tp.tv_nsec*1.0e-9));
}