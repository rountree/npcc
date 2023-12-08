# Remember to upgrade your version of gcc.  On LLNL machines, use:
# 	module load gcc/12.1.1-magic

debug:
	gcc -std=c2x -Wall -Wextra -Og -fsanitize=address -o np.debug nanopond.c

fast:
	gcc -std=c2x -Wall -Wextra -Ofast -g -o np.fast nanopond.c

clean:
	rm -f *.o nanopond *.dSYM




