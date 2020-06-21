#include "Evolution.h"

int evolution(int initPopulation, int generations)
{
	// Correct population to multiple of cores
	int cores = 4;
	int population = initPopulation - (initPopulation % cores);

	// Print timestamp
	cout << "Evolving neural networks with genetic algorithm\n"
		<< population << " individuals, " << generations << " generations\n"
		<< currentDateTime() << "\n\n";

	// Declare variables
	vector<shared_ptr<Network>> generation;
	generation.reserve(population);
	vector<double> fitnesses;
	fitnesses.reserve(generations);
	shared_ptr<Network> winner;
	vector<std::thread> threads(cores);
	double currentFitness = 0;
	double gap = 0.2;
	vector<int> durations(generations);
	int groups = population / cores;

	// Initialise start population with random parameters
	for (int i = 0; i < population; i++)
		generation.push_back(make_shared<Network>());

	// Define time variables
	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

	// Loop which iterates through generations
	for (int g = 0; g < generations; g++)
	{
		// Start clock
		start = std::chrono::high_resolution_clock::now();

		// Print status
		std::cout << "Generation " << g << "\n";

		// Single thread
		/*for (int i = 0; i < population; i++)
			generation[i]->operator()();*/

		// Informed multithread
		for (int group = 0; group < groups; group++)
		{
			std::cout << "Group " << group << "\n";

			for (int core = 0; core < cores; core++)
			{
				shared_ptr<Network> individual = generation[core + cores * group];
				threads[core] = std::thread(std::ref(*individual));
			}

			for (int core = 0; core < cores; core++)
				threads[core].join();
		}

		// Print status
		currentFitness = findFittest(generation)->getKS();
		fitnesses.push_back(currentFitness);
		std::cout << "\tBest fitness (KS): " << currentFitness;

		// Stop clock
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
		std::cout << "\t" << duration.count() << " ms";
		durations[g] = (int)duration.count();

		// End evolution if last generation
		if (g == (generations - 1))
		{
			std::cout << "\t-Evolution complete-\n\n";
			break;
		}

		std::cout << "\t-selection commenced-\n\n";

		// Initialise new population
		generation = selectDescendants(generation, gap);

		// Clear performance
		for (const auto individual : generation)
			individual->clear();
	}

	// Select highest performing Network from generation
	winner = findFittest(generation);

	// Assess strategy by visualization
	winner->visualize();
	winner->exportSynapses();
	winner->exportParameters();
	std::cout << "\n\nResulting fitness (KS): " << winner->getKS() << "\n";

	// Produce fitness statistics
	saveVector("fitness-vector.dat", fitnesses);
	std::cout << "Winner alpha: " << winner->getAlpha() << "\n";

	return 0;
}

vector<shared_ptr<Network>> selectDescendants(vector<shared_ptr<Network>> &generation, double gap)
{
	// Newest implementation
	vector<shared_ptr<Network>> elite;
	vector<shared_ptr<Network>> newGeneration;
	newGeneration.reserve(generation.size());
	vector<shared_ptr<Network>> parents;
	parents.reserve(2);

	int numElite = int(round(generation.size() * gap));
	elite.reserve(numElite);

	//std::sort(generation.begin(), generation.end(), Network::sortByKappa);
	std::sort(generation.begin(), generation.end(), Network::sortByKS);

	for (int i = 0; i < numElite; i++)
	{
		elite.push_back(make_shared<Network>(*generation[i]));
		newGeneration.push_back(make_shared<Network>(*generation[i]));
	}
	for (unsigned i = numElite; i < generation.size(); i++)
	{
		parents = selectParents(elite);
		newGeneration.push_back(make_shared<Network>(*parents[0], *parents[1]));
		parents.clear();
	}
	return newGeneration;
}

vector<shared_ptr<Network>> selectParents(vector<shared_ptr<Network>> &elite)
{
	vector<int> indeces = randomIntVector(elite.size());

	vector<shared_ptr<Network>> parents;
	parents.reserve(2);

	parents.push_back(make_shared<Network>(*elite[indeces[0]]));
	parents.push_back(make_shared<Network>(*elite[indeces[1]]));
	return parents;
}

shared_ptr<Network> findFittest(vector<shared_ptr<Network>> &generation)
{
	//std::sort(generation.begin(), generation.end(), Network::sortByKappa);
	std::sort(generation.begin(), generation.end(), Network::sortByKS);
	return generation[0];
}

void testNetworkConstructor()
{
	Network Adam;
	Adam();
	std::cout << "Network kappa: " << Adam.getKappa() << std::endl;

	Network Eve(Adam);
	Eve();
	std::cout << "Network kappa: " << Eve.getKappa() << std::endl;

	Network Datsa = Adam;
	Datsa();
	std::cout << "Network kappa: " << Datsa.getKappa() << std::endl;
}

void testKappaComparison()
{
	Network Adam;
	Adam();
	std::cout << "Network kappa: " << Adam.getKappa() << std::endl;

	Network Eve(Adam);
	Eve();
	std::cout << "Network kappa: " << Eve.getKappa() << std::endl;

	Network Datsa = Adam;
	Datsa();
	std::cout << "Network kappa: " << Datsa.getKappa() << std::endl;
}

void testSpeed()
{
	int pop = 3;
	std::vector<std::shared_ptr<Network>> Eden;
	std::vector<std::thread> threads(pop);

	// Define time variables
	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

	// Initialise start population with random parameters
	for (int seed = 0; seed < pop; seed++)
		Eden.push_back(make_shared<Network>(seed));

	// Loop which iterates through generations
	for (int g = 0; g < 1; g++)
	{
		// Print status
		std::cout << "Generation " << g;

		// Start clock
		start = std::chrono::high_resolution_clock::now();

		/*Multi-thread approach*/
		// Iterate through current population
		/*for (int i = 0; i < pop; i++)
		{
			shared_ptr<Network> individual = Eden[i];
			threads[i] = std::thread(std::ref(*individual));
		}

		for (int i = 0; i < pop; i++)
			threads[i].join();*/

		/*Single-thread approach*/
		for (int i = 0; i < pop; i++)
			Eden[i]->operator()();

		// Stop clock
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
		std::cout << "\t" << duration.count() << " ms";
	}
	std::cout << "\n\n";
}