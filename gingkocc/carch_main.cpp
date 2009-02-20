#include <iostream>
#include "carch_defs.h"
#include "greet.h"


class Population;
class Cell;

typedef std::vector<float> Genotype;


class Individual 
	{
		//static list<std::vector<Individual> > realIndividual;
		//static std::stack<Individual *> freeIndividuals;
	
	public:
		/*
		static Individual * GetNewIndividual() {
			if (freeIndividuals.empty()) {
				realIndividual.push_back(std::vector<Individual>());
				list<std::vector<Individual> >::iterator rvIt = realIndividual.rbegin();
				rvIt->resize(10000);
				...
				
			}
		}
		static ReleaseNewIndividual(Individual * used)  {
			freeIndividuals.push(used);
		}
		*/
		const Genotype & getGenotype() const {
			return this->nonNeutralGenotype;
		}
	
	private:
		Individual(const Genotype & g)
			:nonNeutralGenotype(g) {
		}
		
		Genotype nonNeutralGenotype;	
		
	};
	
class Species 
	{
	public:
	
	private:
		std::list<Populations *> extantPopulations;
	};

class Population 
	{
	public:
	
	private:
		Species & species;
		Cell * cell;
		std::vector<Individual> individuals;		
	};

class Cell
	{
	public:
		Cell(unsigned nSpecies)
			:population(nSpecies)
			{}
	private:
		std::vector<Population> population; // each population is of a different species
		
		Cell * northNeighbor;
		Cell * southNeighbor;
		Cell * westNeighbor;
		Cell * eastNeighbor;
	};

class World
	{
	public:
		World(unsigned nCells, unsigned nSpecies)
			{
			Cell dummy;
			this->cells.assign(nCells, dummy);
			}
	private:
		std::vector<Cell> cells;
	};
	
	
int main(int, char * []) {

	greet("World");
	
	World world();
	for (unsigned i = 0; i < nIterations; ++i) {
	}
	return 0;
}
	
