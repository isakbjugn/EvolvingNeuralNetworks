#include "Connectivity.h"

void makeConnectivityMatrix()
{
	// Initilize parameters
	int N = 1024; // Number of neurons
	int H = 3; // Number of hierarchical levels
	int M = (int)pow(2, H - 1); // Number of modules;
	double balance = 0.3; // Percentage of inhibitory neurons
	int n = N / M; // Number of neurons per module

	// Indicate indeces of inhibitory neurons
	std::vector<int> inhibNeurons(N, 0);
	for (int m = 0; m < M; m++)
		for (int i = 0; i < floor(n*balance); i++)
			inhibNeurons[(m + 1) * n - i - 1] = 1; 

	// Define connectivity variables
	int degree;
	int minDegree = 20;
	int maxDegree = n;
	double alpha = 1.0;
	double gamma = 2.0;

	// Define data structures
	int preModule;
	int postModule;
	std::vector<int> degrees = powerLawVector(N, alpha, maxDegree, minDegree);
	std::vector<int> preModules;
	std::vector<int> preNeurons;
	preNeurons.reserve(n);
	std::map<int, int> moduleCount;

	// Create connectivity matrix
	Matrix synapses = Matrix(N, N, 0.0);

	// Add connections weighted by hierarchy
	for (int postNeuron = 0; postNeuron < N; postNeuron++)
	{
		// Find degree
		degree = degrees.back();
		degrees.pop_back();

		// Find presynaptic modeles
		postModule = getPostModule(postNeuron, N, M);
		preModules = getPreModuleVector(degree, postModule, H, gamma);
		moduleCount = countOccurences(preModules);

		for (auto it = moduleCount.begin(); it != moduleCount.end(); it++)
		{
			preModule = it->first;

			// Treat excitatory and inhibitory neurons equally
			preNeurons = randomIntVector(n, preModule*n);

			for (int i = 0; i < it->second; i++)
			{
				if (preNeurons.empty())
					continue;

				if (synapses(preNeurons.back(), postNeuron) != 0)
					cout << "Existing neuron!\n";
				
				synapses(preNeurons.back(), postNeuron) = 1;
				preNeurons.pop_back();
			}
		}

		// Alternative treatment of inhibitory neurons
		// - Initialize inhibitory synapses independently
		// - Skip these when looking for *presynaptic* neurons
		// - Specifically: Alter input to "randomIntVector()"
	}

	// Remove auto connections
	//for (int i = 0; i < N; i++)
	//	synapses(i, i) = 0;

	synapses.save("synapses.dat"); // This takes *A LOT OF TIME*
}


int getPostModule(int postNeuron, int neurons, int modules)
{
	// 0-indexed modeles

	if (neurons < modules || (neurons % modules))
		cout << "getPostModule():\tMismatch between number of neurons and modules\n";

	return (int)floor(postNeuron * modules / neurons);
}

int getPreNeuron(int preModule, int neurons, int modules)
{
	// 0-indexed modules

	if (neurons < modules || !(neurons % modules))
		cout << "getPreNeuron():\tMismatch between number of neurons and modules\n";

	int neuron = randomIntWithLimits(0, (int)(neurons / modules));
	return neuron + preModule * neurons / modules;
}

int getPreModule(int postModule, int hierarchyDepth)
{
	int hierarchyLevel = randomPowerDiscrete(2, 1, hierarchyDepth);
	//int hierarchyLevel = randomIntWithLimits(1, hierarchyDepth);
	//int hierarchyLevel = hierarchyDepth;
	//if (hierarchyLevel == 5)
	//	cout << "Kvekk deluxe!";
	//cout << "getPreModule() calls getPreModule(): H = " << hierarchyDepth << ", h = " << hierarchyLevel << "\n";
	return getPreModule(postModule, hierarchyDepth, hierarchyLevel);
}

int getPreModule(int postModule, int hierarchyDepth, int hierarchyLevel)
{
	if (hierarchyLevel > hierarchyDepth)
		cout << "getPreModule():\tHierarchy level exceeds legal bounds\n";

	// Should this be called the same as upper getPreModule()?

	// 0-indexed modules

	using std::vector;

	int rowIdx; // Vertical Index
	int colIdx; // Horisontal index

	switch (hierarchyLevel)
	{
	case 1:
		return postModule;
	case 2:
		return postModule + (int)pow(-1, postModule % 2);
	case 3:
	{
		static vector<vector<int>> m3{ {2, 0, 6, 4, 10, 8, 14, 12},
									   {3, 1, 7, 5, 11, 9, 15, 13} };
		rowIdx = randomIntWithLimits(0, 1);
		colIdx = (int)floor(postModule / 2);
		return m3[rowIdx][colIdx];
	}
	case 4:
	{
		static vector<vector<int>> m4{ {4, 0, 12,  8},
									   {5, 1, 13,  9},
									   {6, 2, 14, 10},
									   {7, 3, 15, 11} };
		rowIdx = randomIntWithLimits(0, 3);
		colIdx = (int)floor(postModule / 4);
		return m4[rowIdx][colIdx];
	}
	case 5:
	{
		static vector<vector<int>> m5{ { 8, 0},
									   { 9, 1},
									   {10, 2},
									   {11, 3},
									   {12, 4},
									   {13, 5},
									   {14, 6},
									   {15, 7} };
		rowIdx = randomIntWithLimits(0, 7);
		colIdx = (int)floor(postModule / 8);
		return m5[rowIdx][colIdx];
	}
	default:
		return 0;
	};
}

std::vector<int> getPreModuleVector(int connections, int postModule, int hierarchyDepth, double gamma)
{
	std::vector<int> preModules;
	preModules.reserve(connections);

	std::vector<int> levels = powerLawVector(connections, gamma, hierarchyDepth);
	for (auto it = levels.begin(); it != levels.end(); it++)
		preModules.push_back(getPreModule(postModule, hierarchyDepth, *it));

	return preModules;

	// Alternative: Draw independent h's
	//for (int i = 0; i < connections; i++)
	//	preModules.push_back(getPreModule(postModule, hierarchyDepth));
	//return preModules;
}



Matrix readMatrixFromFile(const char* filename)
{
	int N = 1024;
	Matrix result(1024, 0);

	std::string line;
	std::ifstream myfile(filename);
	if (myfile.is_open())
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; i < N; i++)
				myfile >> result(i, j);
		myfile.close();
	}

	else cout << "Unable to open file";
	return result;
}

void oldMakeConnectivityMatrix()
{
	// Initilize parameters
	int N = 1024; // Number of neurons
	int H = 4; // Number of hierarchical levels
	int M = (int)pow(2, H - 1); // Number of modules;
	double balance = 0.3; // Percentage of inhibitory neurons
	int n = N / M; // Number of neurons per module

	// Indicate indeces of inhibitory neurons
	std::vector<int> inhibNeurons(N, 0);
	for (int m = 0; m < M; m++)
		for (int i = 0; i < floor(n*balance); i++)
			inhibNeurons[(m + 1) * n - i - 1] = 1;

	// Define connectivity variables
	int degree;
	int minDegree = 10;
	int maxDegree = n / 8;
	double alpha = 2;

	// Define data structures
	int preModule;
	int postModule;
	std::vector<int> degrees = powerLawVector(N, alpha, maxDegree, minDegree);
	std::vector<int> synCounts;
	synCounts.reserve(N);
	int synCount;
	saveVector("super-degrees.dat", degrees);
	std::vector<int> preModules;
	std::vector<int> preNeurons;
	preNeurons.reserve(n);
	std::map<int, int> moduleCount;

	// Create connectivity matrix
	Matrix synapses = Matrix(N, N, 0.0);

	// Test: See if degrees are conserved. Yes, with certain precautions they are preserved
	/*
	for (int postNeuron = 0; postNeuron < N; postNeuron++)
	{
		// Get degree
		degree = degrees.back();
		degrees.pop_back();

		// Uniform outcoming degree
		for (int j = 0; j < degree; j++)
		{
			if (preNeurons.empty())
				preNeurons = randomIntVector(N);
			while (preNeurons.back() == postNeuron || synapses(preNeurons.back(), postNeuron) == 1)
				preNeurons.pop_back();

			synapses(preNeurons.back(), postNeuron) = 1;
			preNeurons.pop_back();
		}

		// Apparently Gaussian degree
		//if (i < N / 2)
		//	for (int j = 0; j < degree; j++)
		//		synapses(i + j + 1, i) = 1;
		//else
		//	for (int j = 0; j < degree; j++)
		//		synapses(i - j - 1, i) = 1;
	}
	*/

	// Add intramodular connections
	/*
	for (int m = 0; m < M; m++)
	{
		for (int pre = 0; pre < n; pre++)
		{
			for (int post = 0; post < n; post++)
			{
				// Treat all neurons equally
				synapses(pre + m*n, post + m*n) = 1;

				// Differential treatment of inhibitory neurons?
			}
		}
	}
	*/

	// Add connections weighted by hierarchy
	for (int postNeuron = 0; postNeuron < N; postNeuron++)
	{
		synCount = 0;

		// Find degree
		degree = degrees.back();
		degrees.pop_back();

		// Find presynaptic modeles
		postModule = getPostModule(postNeuron, N, M);
		preModules = getPreModuleVector(degree, postModule, H);
		moduleCount = countOccurences(preModules);
		//cout << "Difference: " << degree - preModules.size() << "\t" << degree - summarizeMap(moduleCount) << "\n";

		for (auto it = moduleCount.begin(); it != moduleCount.end(); it++)
		{
			preModule = it->first;

			// Treat excitatory and inhibitory neurons equally
			preNeurons = randomIntVector(n, preModule*n);

			for (int i = 0; i < it->second; i++)
			{
				if (synapses(preNeurons.back(), postNeuron) != 0)
					cout << "Existing neuron!\n";
				if ((postModule % 2) ^ (preModule % 2))
					synapses(preNeurons.back(), postNeuron) = 1;
				else
					synapses(preNeurons.back(), postNeuron) = -1;
				preNeurons.pop_back();
				synCount++;
			}
		}

		synCounts.push_back(synCount);

		// Alternative treatment of inhibitory neurons
		// - Initialize inhibitory synapses independently
		// - Skip these when looking for *presynaptic* neurons
		// - Specifically: Alter input to "randomIntVector()"
	}

	// Remove auto connections
	//for (int i = 0; i < N; i++)
	//	synapses(i, i) = 0;

	synapses.save("synapses.dat"); // This takes *A LOT OF TIME*

	std::sort(synCounts.begin(), synCounts.end());

	cout << "Kvekk!\n";
}


/*Tests*/

void testGetPreModule()
{
	cout << "Test getPreModule()\n";
	int M = 16;
	int H = 5;

	for (int postModule = 0; postModule < M; postModule++)
	{
		cout << postModule << ":\t";
		for (int h = 1; h <= 5; h++)
		{
			cout << getPreModule(postModule, H, h) << "\t";
		}
		cout << endl;
	}
}

void testNeuronModuleConversion()
{
	int N = 512;
	int M = 16;

	cout << "Neuron-to-module conversion:\n";
	std::vector<int> postNeurons{ 0, 1, 10, 20, 100, 257, 470, 511 };
	for (auto it = postNeurons.begin(); it != postNeurons.end(); it++)
	{
		cout << *it << ": " << getPostModule(*it, N, M) << "\t";
	}
	
	cout << "\n\nModule-to-neuron conversion:\n";
	for (int i = 0; i < M; i++)
	{
		cout << getPreNeuron(i, N, M) << "\t";
	}
	cout << endl;
}

void testPowerLawDiscrete()
{
	std::ofstream myfile;
	myfile.open("PowerDist.dat", std::ios::trunc);
	int number;

	for (int i = 0; i < 10000; i++)
	{
		number = randomPowerDiscrete(2, 1, 20);
		//cout << number << "\t";
		myfile << number << " ";
	}
	//cout << endl;
	myfile << "\n";
	myfile.close();
}

void testGetDegrees()
{
	int neurons = 1024;
	int minDegree = 10;
	int maxDegree = 128;
	double alpha = 2;
	std::vector<int> degrees;

	// Test sample-and-insert method
	degrees = powerLawVectorSampleAndInsert(neurons, alpha, maxDegree, minDegree);
	saveVector("sample-and-insert-degrees.dat", degrees);

	// Test fill-and-shuffle method
	degrees = powerLawVector(neurons, alpha, maxDegree, minDegree);
	saveVector("fill-and-shuffle-degrees.dat", degrees);
}

void testCountOccurences()
{
	std::vector<int> example { 1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 6, 7, 8, 9, 9, 9 };
	std::map<int, int> exampleCount = countOccurences(example);

	for (auto it = exampleCount.begin(); it != exampleCount.end(); it++)
	{
		cout << it->first << ":\t" << it->second << "\n";
	}
}

void testSynapseDegrees()
{
	Matrix synapses = readMatrixFromFile("synapses.dat");

}

void testBalance(double balance)
{
	// Initilize parameters
	int N = 1024; // Number of neurons
	int H = 6; // Number of hierarchical levels
	int M = (int)pow(2, H - 1); // Number of modules;
	int n = N / M; // Number of neurons per module
	int ni = (int)floor(n * balance);
	int Ni = M * ni;
	int Ne = N - Ni;

	// Indicate indeces of inhibitory neurons
	std::vector<int> inhibNeurons(N, 0);
	for (int m = 0; m < M; m++)
		for (int i = 0; i < ni; i++)
			inhibNeurons[(m + 1) * n - i - 1] = 1;

	int numExcitatory = 0;
	int numInhibitory = 0;

	for (unsigned int i = 0; i < inhibNeurons.size(); i++)
	{
		if (inhibNeurons[i])
			numInhibitory++;
		else
			numExcitatory++;
	}
	cout << "Total:" << N << "\tExc: " << numExcitatory << "\tInh: " << numInhibitory << "\tSum: " << numExcitatory + numInhibitory
		<< "\tNe: " << Ne << "\tNi: " << Ni << "\tSum: " << Ne + Ni << "\n";
}