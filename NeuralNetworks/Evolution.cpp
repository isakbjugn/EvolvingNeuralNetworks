#include "Evolution.h"

int evolution(int population, int generations)
{
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
	vector<std::thread> threads(population);
	double currentFitness = 0;
	double gap = 0.2;
	vector<int> durations(generations);

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
		std::cout << "Generation " << g;

		/*Multi-thread approach*/
		// Iterate through current population
		for (int i = 0; i < population; i++)
		{
			shared_ptr<Network> individual = generation[i];
			threads[i] = std::thread(std::ref(*individual));
		}

		for (int i = 0; i < population; i++)
			threads[i].join();

		/*Single-thread approach*/
		//for (int i = 0; i < population; i++)
		//	generation[i]->operator()();

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
			break;

		std::cout << "\t-selection commenced-\n";

		// Initialise new population
		generation = selectDescendants(generation, gap);

		// Clear performance
		for (const auto individual : generation)
			individual->clear();
	}

	// Select highest performing Network from generation
	winner = findFittest(generation);

	// Save duration vector
	/*
	string fileStr = "durations";
	fileStr += "-p" + std::to_string(population)
		+ "-g" + std::to_string(generations)
		+ "-n" + std::to_string(winner->getNeurons())
		+ "-d" + std::to_string(winner->getDuration())
		+ ".dat";
	const char* filename = fileStr.c_str();
	saveVector(filename, durations);
	*/

	// Assess strategy by visualization
	winner->visualize();
	winner->exportSynapses();
	winner->exportParameters();
	std::cout << "\n\nResulting fitness (KS): " << winner->getKS() << "\n";

	// Produce fitness statistics
	std::cout << "Historic fitness: ";
	for (const double& k :fitnesses)
		std::cout << k << ", ";
	std::cout << "\n";

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