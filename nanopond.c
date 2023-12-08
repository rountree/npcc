/* ----------------------------------------------------------------------- */
/* Tunable parameters                                                      */
/* ----------------------------------------------------------------------- */

/* Frequency of comprehensive reports-- lower values will provide more
 * info while slowing down the simulation. Higher values will give less
 * frequent updates. */
/* This is also the frequency of screen refreshes if SDL is enabled. */
#define REPORT_FREQUENCY 200000

/* Mutation rate -- range is from 0 (none) to 0xffffffff (all mutations!) */
/* To get it from a float probability from 0.0 to 1.0, multiply it by
 * 4294967295 (0xffffffff) and round. */
#define MUTATION_RATE 5000

/* How frequently should random cells / energy be introduced?
 * Making this too high makes things very chaotic. Making it too low
 * might not introduce enough energy. */
#define INFLOW_FREQUENCY 100

/* Base amount of energy to introduce per INFLOW_FREQUENCY ticks */
#define INFLOW_RATE_BASE 600

/* A random amount of energy between 0 and this is added to
 * INFLOW_RATE_BASE when energy is introduced. Comment this out for
 * no variation in inflow rate. */
#define INFLOW_RATE_VARIATION 1000

/* Size of pond in X and Y dimensions. */
#define POND_SIZE_X 800
#define POND_SIZE_Y 600

/* Depth of pond in four-bit codons -- this is the maximum
 * genome size. This *must* be a multiple of 16! */
#define POND_DEPTH 1024

/* This is the divisor that determines how much energy is taken
 * from cells when they try to KILL a viable cell neighbor and
 * fail. Higher numbers mean lower penalties. */
#define FAILED_KILL_PENALTY 3

/* Define this to use SDL. To use SDL, you must have SDL headers
 * available and you must link with the SDL library when you compile. */
/* Comment this out to compile without SDL visualization support. */
#define USE_SDL 1

/* ----------------------------------------------------------------------- */

#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

volatile uint64_t prngState[2];
static inline uintptr_t getRandom()
{
	// https://en.wikipedia.org/wiki/Xorshift#xorshift.2B
	uint64_t x = prngState[0];
	const uint64_t y = prngState[1];
	prngState[0] = y;
	x ^= x << 23;
	const uint64_t z = x ^ y ^ (x >> 17) ^ (y >> 26);
	prngState[1] = z;
	return (uintptr_t)(z + y);
}

/* Pond depth in machine-size words.  This is calculated from
 * POND_DEPTH and the size of the machine word. (The multiplication
 * by two is due to the fact that there are two four-bit values in
 * each eight-bit byte.) */
#define POND_DEPTH_SYSWORDS (POND_DEPTH / (sizeof(uintptr_t) * 2))

/* Number of bits in a machine-size word */
#define SYSWORD_BITS (sizeof(uintptr_t) * 8)

/* Constants representing neighbors in the 2D grid. */
#define N_LEFT 0
#define N_RIGHT 1
#define N_UP 2
#define N_DOWN 3

/* Word and bit at which to start execution */
/* This is after the "logo" */
#define EXEC_START_WORD 0
#define EXEC_START_BIT 4

/* Number of bits set in binary numbers 0000 through 1111 */
static const uintptr_t BITS_IN_FOURBIT_WORD[16] = { 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4 };

/**
 * Structure for a cell in the pond
 */
struct Cell
{
	/* Globally unique cell ID */
	uint64_t ID;

	/* ID of the cell's parent */
	uint64_t parentID;

	/* Counter for original lineages -- equal to the cell ID of
	 * the first cell in the line. */
	uint64_t lineage;

	/* Generations start at 0 and are incremented from there. */
	uintptr_t generation;

	/* Energy level of this cell */
	uintptr_t energy;

	/* Memory space for cell genome (genome is stored as four
	 * bit instructions packed into machine size words) */
	uintptr_t genome[POND_DEPTH_SYSWORDS];
};

/* The pond is a 2D array of cells */
static struct Cell pond[POND_SIZE_X][POND_SIZE_Y];

/* This is used to generate unique cell IDs */
static volatile uint64_t cellIdCounter = 0;

volatile struct {
	/* Counts for the number of times each instruction was
	 * executed since the last report. */
	double instructionExecutions[16];

	/* Number of cells executed since last report */
	double cellExecutions;

	/* Number of viable cells replaced by other cells' offspring */
	uintptr_t viableCellsReplaced;

	/* Number of viable cells KILLed */
	uintptr_t viableCellsKilled;

	/* Number of successful SHARE operations */
	uintptr_t viableCellShares;
} statCounters;

static void doReport(const uint64_t clock)
{
	static uint64_t lastTotalViableReplicators = 0;

	uintptr_t x,y;

	uint64_t totalActiveCells = 0;
	uint64_t totalEnergy = 0;
	uint64_t totalViableReplicators = 0;
	uintptr_t maxGeneration = 0;

	for(x=0;x<POND_SIZE_X;++x) {
		for(y=0;y<POND_SIZE_Y;++y) {
			struct Cell *const c = &pond[x][y];
			if (c->energy) {
				++totalActiveCells;
				totalEnergy += (uint64_t)c->energy;
				if (c->generation > 2)
					++totalViableReplicators;
				if (c->generation > maxGeneration)
					maxGeneration = c->generation;
			}
		}
	}

	/* Look here to get the columns in the CSV output */

	/* The first five are here and are self-explanatory */
	printf("%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",",
		(uint64_t)clock,
		(uint64_t)totalEnergy,
		(uint64_t)totalActiveCells,
		(uint64_t)totalViableReplicators,
		(uint64_t)maxGeneration,
		(uint64_t)statCounters.viableCellsReplaced,
		(uint64_t)statCounters.viableCellsKilled,
		(uint64_t)statCounters.viableCellShares
		);

	/* The next 16 are the average frequencies of execution for each
	 * instruction per cell execution. */
	double totalMetabolism = 0.0;
	for(x=0;x<16;++x) {
		totalMetabolism += statCounters.instructionExecutions[x];
		printf(",%.4f",(statCounters.cellExecutions > 0.0) ? (statCounters.instructionExecutions[x] / statCounters.cellExecutions) : 0.0);
	}

	/* The last column is the average metabolism per cell execution */
	printf(",%.4f\n",(statCounters.cellExecutions > 0.0) ? (totalMetabolism / statCounters.cellExecutions) : 0.0);
	fflush(stdout);

	if ((lastTotalViableReplicators > 0)&&(totalViableReplicators == 0))
		fprintf(stderr,"[EVENT] Viable replicators have gone extinct. Please reserve a moment of silence.\n");
	else if ((lastTotalViableReplicators == 0)&&(totalViableReplicators > 0))
		fprintf(stderr,"[EVENT] Viable replicators have appeared!\n");

	lastTotalViableReplicators = totalViableReplicators;

	/* Reset per-report stat counters */
	for(x=0;x<sizeof(statCounters);++x)
		((uint8_t *)&statCounters)[x] = (uint8_t)0;
}

static inline struct Cell *getNeighbor(const uintptr_t x,const uintptr_t y,const uintptr_t dir)
{
	/* Space is toroidal; it wraps at edges */
	switch(dir) {
		case N_LEFT:
			return (x) ? &pond[x-1][y] : &pond[POND_SIZE_X-1][y];
		case N_RIGHT:
			return (x < (POND_SIZE_X-1)) ? &pond[x+1][y] : &pond[0][y];
		case N_UP:
			return (y) ? &pond[x][y-1] : &pond[x][POND_SIZE_Y-1];
		case N_DOWN:
			return (y < (POND_SIZE_Y-1)) ? &pond[x][y+1] : &pond[x][0];
	}
	return &pond[x][y]; /* This should never be reached */
}

static inline int accessAllowed(struct Cell *const c2,const uintptr_t c1guess,int sense)
{
	/* Access permission is more probable if they are more similar in sense 0,
	 * and more probable if they are different in sense 1. Sense 0 is used for
	 * "negative" interactions and sense 1 for "positive" ones. */
	return sense ? (((getRandom() & 0xf) >= BITS_IN_FOURBIT_WORD[(c2->genome[0] & 0xf) ^ (c1guess & 0xf)])||(!c2->parentID)) : (((getRandom() & 0xf) <= BITS_IN_FOURBIT_WORD[(c2->genome[0] & 0xf) ^ (c1guess & 0xf)])||(!c2->parentID));
}

volatile int exitNow = 0;

static void *run(void *targ)
{
	const uintptr_t threadNo = (uintptr_t)targ;
	uintptr_t x,y,i;
	uintptr_t clock = 0;

	/* Buffer used for execution output of candidate offspring */
	uintptr_t outputBuf[POND_DEPTH_SYSWORDS];

	/* Miscellaneous variables used in the loop */
	uintptr_t currentWord,wordPtr,shiftPtr,inst,tmp;
	struct Cell *pptr,*tmpptr;

	/* Virtual machine memory pointer register (which
	 * exists in two parts... read the code below...) */
	uintptr_t ptr_wordPtr;
	uintptr_t ptr_shiftPtr;

	/* The main "register" */
	uintptr_t reg;

	/* Which way is the cell facing? */
	uintptr_t facing;

	/* Virtual machine loop/rep stack */
	uintptr_t loopStack_wordPtr[POND_DEPTH];
	uintptr_t loopStack_shiftPtr[POND_DEPTH];
	uintptr_t loopStackPtr;

	/* If this is nonzero, we're skipping to matching REP */
	/* It is incremented to track the depth of a nested set
	 * of LOOP/REP pairs in false state. */
	uintptr_t falseLoopDepth;

	/* If this is nonzero, cell execution stops. This allows us
	 * to avoid the ugly use of a goto to exit the loop. :) */
	int stop;

	/* Main loop */
	while (!exitNow) {
		/* Increment clock and run reports periodically */
		/* Clock is incremented at the start, so it starts at 1 */
		++clock;
		if ((threadNo == 0)&&(!(clock % REPORT_FREQUENCY))) {
			doReport(clock);
		}

		/* Introduce a random cell somewhere with a given energy level */
		/* This is called seeding, and introduces both energy and
		 * entropy into the substrate. This happens every INFLOW_FREQUENCY
		 * clock ticks. */
		if (!(clock % INFLOW_FREQUENCY)) {
			x = getRandom() % POND_SIZE_X;
			y = getRandom() % POND_SIZE_Y;
			pptr = &pond[x][y];

			pptr->ID = cellIdCounter;
			pptr->parentID = 0;
			pptr->lineage = cellIdCounter;
			pptr->generation = 0;
#ifdef INFLOW_RATE_VARIATION
			pptr->energy += INFLOW_RATE_BASE + (getRandom() % INFLOW_RATE_VARIATION);
#else
			pptr->energy += INFLOW_RATE_BASE;
#endif /* INFLOW_RATE_VARIATION */
			for(i=0;i<POND_DEPTH_SYSWORDS;++i)
				pptr->genome[i] = getRandom();
			++cellIdCounter;

			/* Update the random cell on SDL screen if viz is enabled */
		}

		/* Pick a random cell to execute */
		i = getRandom();
		x = i % POND_SIZE_X;
		y = ((i / POND_SIZE_X) >> 1) % POND_SIZE_Y;
		pptr = &pond[x][y];

		/* Reset the state of the VM prior to execution */
		for(i=0;i<POND_DEPTH_SYSWORDS;++i)
			outputBuf[i] = ~((uintptr_t)0); /* ~0 == 0xfffff... */
		ptr_wordPtr = 0;
		ptr_shiftPtr = 0;
		reg = 0;
		loopStackPtr = 0;
		wordPtr = EXEC_START_WORD;
		shiftPtr = EXEC_START_BIT;
		facing = 0;
		falseLoopDepth = 0;
		stop = 0;

		/* We use a currentWord buffer to hold the word we're
		 * currently working on.  This speeds things up a bit
		 * since it eliminates a pointer dereference in the
		 * inner loop. We have to be careful to refresh this
		 * whenever it might have changed... take a look at
		 * the code. :) */
		currentWord = pptr->genome[0];

		/* Keep track of how many cells have been executed */
		statCounters.cellExecutions += 1.0;

		/* Core execution loop */
		while ((pptr->energy)&&(!stop)) {
			/* Get the next instruction */
			inst = (currentWord >> shiftPtr) & 0xf;

			/* Randomly frob either the instruction or the register with a
			 * probability defined by MUTATION_RATE. This introduces variation,
			 * and since the variation is introduced into the state of the VM
			 * it can have all manner of different effects on the end result of
			 * replication: insertions, deletions, duplications of entire
			 * ranges of the genome, etc. */
			if ((getRandom() & 0xffffffff) < MUTATION_RATE) {
				tmp = getRandom(); /* Call getRandom() only once for speed */
				if (tmp & 0x80) /* Check for the 8th bit to get random boolean */
					inst = tmp & 0xf; /* Only the first four bits are used here */
				else reg = tmp & 0xf;
			}

			/* Each instruction processed costs one unit of energy */
			--pptr->energy;

			/* Execute the instruction */
			if (falseLoopDepth) {
				/* Skip forward to matching REP if we're in a false loop. */
				if (inst == 0x9) /* Increment false LOOP depth */
					++falseLoopDepth;
				else if (inst == 0xa) /* Decrement on REP */
					--falseLoopDepth;
			} else {
				/* If we're not in a false LOOP/REP, execute normally */

				/* Keep track of execution frequencies for each instruction */
				statCounters.instructionExecutions[inst] += 1.0;

				switch(inst) {
					case 0x0: /* ZERO: Zero VM state registers */
						reg = 0;
						ptr_wordPtr = 0;
						ptr_shiftPtr = 0;
						facing = 0;
						break;
					case 0x1: /* FWD: Increment the pointer (wrap at end) */
						if ((ptr_shiftPtr += 4) >= SYSWORD_BITS) {
							if (++ptr_wordPtr >= POND_DEPTH_SYSWORDS)
								ptr_wordPtr = 0;
							ptr_shiftPtr = 0;
						}
						break;
					case 0x2: /* BACK: Decrement the pointer (wrap at beginning) */
						if (ptr_shiftPtr)
							ptr_shiftPtr -= 4;
						else {
							if (ptr_wordPtr)
								--ptr_wordPtr;
							else ptr_wordPtr = POND_DEPTH_SYSWORDS - 1;
							ptr_shiftPtr = SYSWORD_BITS - 4;
						}
						break;
					case 0x3: /* INC: Increment the register */
						reg = (reg + 1) & 0xf;
						break;
					case 0x4: /* DEC: Decrement the register */
						reg = (reg - 1) & 0xf;
						break;
					case 0x5: /* READG: Read into the register from genome */
						reg = (pptr->genome[ptr_wordPtr] >> ptr_shiftPtr) & 0xf;
						break;
					case 0x6: /* WRITEG: Write out from the register to genome */
						pptr->genome[ptr_wordPtr] &= ~(((uintptr_t)0xf) << ptr_shiftPtr);
						pptr->genome[ptr_wordPtr] |= reg << ptr_shiftPtr;
						currentWord = pptr->genome[wordPtr]; /* Must refresh in case this changed! */
						break;
					case 0x7: /* READB: Read into the register from buffer */
						reg = (outputBuf[ptr_wordPtr] >> ptr_shiftPtr) & 0xf;
						break;
					case 0x8: /* WRITEB: Write out from the register to buffer */
						outputBuf[ptr_wordPtr] &= ~(((uintptr_t)0xf) << ptr_shiftPtr);
						outputBuf[ptr_wordPtr] |= reg << ptr_shiftPtr;
						break;
					case 0x9: /* LOOP: Jump forward to matching REP if register is zero */
						if (reg) {
							if (loopStackPtr >= POND_DEPTH)
								stop = 1; /* Stack overflow ends execution */
							else {
								loopStack_wordPtr[loopStackPtr] = wordPtr;
								loopStack_shiftPtr[loopStackPtr] = shiftPtr;
								++loopStackPtr;
							}
						} else falseLoopDepth = 1;
						break;
					case 0xa: /* REP: Jump back to matching LOOP if register is nonzero */
						if (loopStackPtr) {
							--loopStackPtr;
							if (reg) {
								wordPtr = loopStack_wordPtr[loopStackPtr];
								shiftPtr = loopStack_shiftPtr[loopStackPtr];
								currentWord = pptr->genome[wordPtr];
								/* This ensures that the LOOP is rerun */
								continue;
							}
						}
						break;
					case 0xb: /* TURN: Turn in the direction specified by register */
						facing = reg & 3;
						break;
					case 0xc: /* XCHG: Skip next instruction and exchange value of register with it */
						if ((shiftPtr += 4) >= SYSWORD_BITS) {
							if (++wordPtr >= POND_DEPTH_SYSWORDS) {
								wordPtr = EXEC_START_WORD;
								shiftPtr = EXEC_START_BIT;
							} else shiftPtr = 0;
						}
						tmp = reg;
						reg = (pptr->genome[wordPtr] >> shiftPtr) & 0xf;
						pptr->genome[wordPtr] &= ~(((uintptr_t)0xf) << shiftPtr);
						pptr->genome[wordPtr] |= tmp << shiftPtr;
						currentWord = pptr->genome[wordPtr];
						break;
					case 0xd: /* KILL: Blow away neighboring cell if allowed with penalty on failure */
						tmpptr = getNeighbor(x,y,facing);
						if (accessAllowed(tmpptr,reg,0)) {
							if (tmpptr->generation > 2)
								++statCounters.viableCellsKilled;

							/* Filling first two words with 0xfffff... is enough */
							tmpptr->genome[0] = ~((uintptr_t)0);
							tmpptr->genome[1] = ~((uintptr_t)0);
							tmpptr->ID = cellIdCounter;
							tmpptr->parentID = 0;
							tmpptr->lineage = cellIdCounter;
							tmpptr->generation = 0;
							++cellIdCounter;
						} else if (tmpptr->generation > 2) {
							tmp = pptr->energy / FAILED_KILL_PENALTY;
							if (pptr->energy > tmp)
								pptr->energy -= tmp;
							else pptr->energy = 0;
						}
						break;
					case 0xe: /* SHARE: Equalize energy between self and neighbor if allowed */
						tmpptr = getNeighbor(x,y,facing);
						if (accessAllowed(tmpptr,reg,1)) {
							if (tmpptr->generation > 2)
								++statCounters.viableCellShares;
							tmp = pptr->energy + tmpptr->energy;
							tmpptr->energy = tmp / 2;
							pptr->energy = tmp - tmpptr->energy;
						}
						break;
					case 0xf: /* STOP: End execution */
						stop = 1;
						break;
				}
			}

			/* Advance the shift and word pointers, and loop around
			 * to the beginning at the end of the genome. */
			if ((shiftPtr += 4) >= SYSWORD_BITS) {
				if (++wordPtr >= POND_DEPTH_SYSWORDS) {
					wordPtr = EXEC_START_WORD;
					shiftPtr = EXEC_START_BIT;
				} else shiftPtr = 0;
				currentWord = pptr->genome[wordPtr];
			}
		}

		/* Copy outputBuf into neighbor if access is permitted and there
		 * is energy there to make something happen. There is no need
		 * to copy to a cell with no energy, since anything copied there
		 * would never be executed and then would be replaced with random
		 * junk eventually. See the seeding code in the main loop above. */
		if ((outputBuf[0] & 0xff) != 0xff) {
			tmpptr = getNeighbor(x,y,facing);
			if ((tmpptr->energy)&&accessAllowed(tmpptr,reg,0)) {
				/* Log it if we're replacing a viable cell */
				if (tmpptr->generation > 2)
					++statCounters.viableCellsReplaced;

				tmpptr->ID = ++cellIdCounter;
				tmpptr->parentID = pptr->ID;
				tmpptr->lineage = pptr->lineage; /* Lineage is copied in offspring */
				tmpptr->generation = pptr->generation + 1;

				for(i=0;i<POND_DEPTH_SYSWORDS;++i)
					tmpptr->genome[i] = outputBuf[i];
			}
		}

	}

	return (void *)0;
}

/**
 * Main method
 *
 * @param argc Number of args
 * @param argv Argument array
 */
int main(int ,char **)
{
	uintptr_t i,x,y;

	/* Seed and init the random number generator */
	prngState[0] = 13;
	srand(13);
	prngState[1] = (uint64_t)rand();

	/* Reset per-report stat counters */
	for(x=0;x<sizeof(statCounters);++x)
		((uint8_t *)&statCounters)[x] = (uint8_t)0;

	/* Clear the pond and initialize all genomes
	 * to 0xffff... */
	for(x=0;x<POND_SIZE_X;++x) {
		for(y=0;y<POND_SIZE_Y;++y) {
			pond[x][y].ID = 0;
			pond[x][y].parentID = 0;
			pond[x][y].lineage = 0;
			pond[x][y].generation = 0;
			pond[x][y].energy = 0;
			for(i=0;i<POND_DEPTH_SYSWORDS;++i)
				pond[x][y].genome[i] = ~((uintptr_t)0);
		}
	}
	run((void *)0);
	return 0;
}
