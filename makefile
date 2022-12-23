DIRECTORY = diretorio_principal/
OBJECTS = main.o matriz_essential.o time_manip.o io_handler.o matriz_math.o

FLAGS = -Wall -lm

all: clean cgSolver
	@echo "\n\n\nCompilacao finalizada\nExecute cgSolver com os parametros desejados\n-n | -k | -p | -i | -e | -o\nParametro -e opcional\n\n"

#-----------------


clean :	purge
	@rm -f $(addprefix $(DIRECTORY), $(OBJECTS))

purge :
	@rm -f cgSolver

#-----------------

cgSolver : $(addprefix $(DIRECTORY), $(OBJECTS))
	gcc -o cgSolver $(FLAGS) $(addprefix $(DIRECTORY), $(OBJECTS))

#-----------------

$@.o : $@.c
	gcc -c $(FLAGS) $@