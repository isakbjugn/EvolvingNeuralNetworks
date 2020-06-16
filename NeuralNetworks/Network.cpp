#include "Network.h"

/*Constructors*/

Network::Network(int seed)
{
	// Initialize random engine
	engine = startEngine(seed);

	// Make new genome
	importGenome(newGenome());

	// Initialize network structure
	initialize();
}

Network::Network(const Network &other)
{
	// Question: Should the current random engine be copied or remade? Currently remade
	
	// Copy members
	importGenome(other.exportGenome());
	defineNeurons();
	S = other.S;
	a = other.a;
	b = other.b;
	c = other.c;
	d = other.d;

	// Copy mutable member
	alpha = other.alpha;
	kappa = other.kappa;
	KS = other.KS;
	avalanches = other.avalanches;
	activeNeurons = other.activeNeurons;

	if (activeNeurons.size())
		cout << "Error in ctor Network(const Network &other): vector activeNeurons not empty\n";

	// Initialize engine
	engine = startEngine();
}

Network::Network(const Network& mom, const Network& dad, int seed)
{
	// Initialize engine
	engine = startEngine(seed);

	// Combine genomes
	combineGenomes(mom.exportGenome(), dad.exportGenome());		// Draw approachh
	//spliceGenomes(mom.exportGenome(), dad.exportGenome());	// Splice approach

	// Mutate genome	
	mutate();

	// Initialize network structure
	initialize();
}

/*Operators*/

Network Network::operator= (Network rhs)
{
	return rhs;
}

void Network::operator()()
{
	run();

	if (alive)
	{
		analyze();
		std::cout << "Thread " << std::this_thread::get_id() << " finished\n";

	}
}

/*Initialization methods*/

void Network::initialize()
{
	// Initilize network structure
	defineNeurons();
	initializeIzhikevich();
	createSynapses();
	addSynapseWeights();

	// Reset performance variables
	clear();
}

std::vector<double> Network::newGenome() const
{
	// Remember to check the limits of n after calling this

	return std::vector<double>
	{
		0.01 * randInt(20, 30),
		(double)randInt(1, 5),
		(double)randInt(1, 30),
		(double)randInt(1, N),
		randReal(1.0, 2.0),
		randReal(1.5, 3.0),
		randReal(4.0, 6.0),
		randReal(1.0, 3.0),
		randReal(0.1, 2.0),
		randReal(-2.0, -0.1)
	};
	// Izhikevich parameters a, b, c, d (mean, variance)
	// Incoming degree mean (rather than exponent, min and max)*/
}

void Network::importGenome(const std::vector<double> &genome)
{
	balance = genome[0];
	H = (int)genome[1];
	minDegree = (int)genome[2];
	maxDegree = (int)genome[3];
	beta = genome[4];
	gamma = genome[5];
	exNoise = genome[6];
	inNoise = genome[7];
	exWeight = genome[8];
	inWeight = genome[9];

	// Check limits of n
	n = N / (int)pow(2, H - 1);
	if (maxDegree > n || maxDegree < minDegree)
		maxDegree = randInt(minDegree, n);
}

void Network::combineGenomes(const std::vector<double> &momGenome, const std::vector<double> &dadGenome)
{
	std::vector<double> genome = momGenome;
	for (size_t i = 0; i < genome.size(); i++)
	{
		if (randReal() < 0.5)
			genome[i] = dadGenome[i];
	}

	importGenome(genome);
}

void Network::spliceGenomes(const std::vector<double> &momGenome, const std::vector<double> &dadGenome)
{
	std::vector<double> genome = momGenome;
	int genes = genome.size();
	int splice = randInt(0, genes - 1);
	for (int i = splice; i < genes; i++)
		genome[i] = dadGenome[i];

	importGenome(genome);
}

void Network::defineNeurons()
{
	// Calculate network parameters
	M = (int)pow(2, H - 1);
	n = N / M;
	ni = (int)floor(n * balance);
	ne = n - ni;
	Ni = M * ni;
	Ne = N - Ni;

	// Indicate indeces of inhibitory neurons
	inhibNeurons = std::vector<int>(N, 0);
	for (int m = 0; m < M; m++)
		for (int i = 0; i < ni; i++)
			inhibNeurons[(m + 1) * n - i - 1] = 1;

	// Make vector of the inhibitory neuron indeces?
}

void Network::mutate()
{
	// Incorporates one mutation into genome

	std::vector<double> mew = this->exportGenome();
	std::vector<double> mewtwo = this->newGenome();
	int locus = randInt(0, mew.size() - 1);
	mew[locus] = mewtwo[locus];
	importGenome(mew);
}

void Network::initializeIzhikevich()
{
	// Random vectors
	Vector re = Vector(Ne);
	for (int i = 0; i < Ne; i++)
		re(i) = randReal();

	Vector ri = Vector(Ni);
	for (int i = 0; i < Ni; i++)
		ri(i) = randReal();

	int exIdx = 0;
	int inIdx = 0;

	// Initialize Izhikevich parameter vectors
	a = Vector(N);
	b = Vector(N);
	c = Vector(N);
	d = Vector(N);

	for (int neuron = 0; neuron < N; neuron++)
	{
		// Excitatory parameters
		if (!inhibNeurons[neuron])
		{
			a(neuron) = 0.02;
			b(neuron) = 0.2;
			c(neuron) = -65 + 15 * pow(re(exIdx), 2);
			d(neuron) = 8 - 6 * pow(re(exIdx), 2);
			exIdx++;
		}

		// Inhibitory parameters
		else
		{
			a(neuron) = 0.02 + 0.08 * ri(inIdx);
			b(neuron) = 0.25 - 0.05 * ri(inIdx);
			c(neuron) = -65;
			d(neuron) = 2;
			inIdx++;
		}
	}
}

void Network::createSynapses()
{
	// Define data structures
	int degree;
	int preModule;
	int postModule;
	std::vector<int> degrees = randPowerVector(N, beta, minDegree, maxDegree);
	std::vector<int> preModules;
	std::vector<int> preNeurons;
	preNeurons.reserve(n);
	std::map<int, int> moduleCount;

	// Create connectivity matrix
	S = Matrix(N, N);

	// Add connections weighted by hierarchy
	for (int postNeuron = 0; postNeuron < N; postNeuron++)
	{
		// Find degree
		degree = degrees.back();
		degrees.pop_back();

		// Find presynaptic modeles
		postModule = getPostModule(postNeuron, N, M);
		preModules = getPreModuleVector(degree, postModule);
		moduleCount = countOccurences(preModules);

		for (auto it = moduleCount.begin(); it != moduleCount.end(); it++)
		{
			preModule = it->first;

			if (it->first == postModule && it->second >= n)
				cout << "Error in createSynapses(): More synapses than neurons in module\n";

			// Treat excitatory and inhibitory neurons equally
			preNeurons = randIntVector(n, preModule*n);

			for (int i = 0; i < it->second; i++)
			{
				// Avoid auto connections
				if (preNeurons.back() == postNeuron)
					preNeurons.pop_back();

				if (S(preNeurons.back(), postNeuron))
					cout << "Existing neuron!\n";
				S(preNeurons.back(), postNeuron) = 1;
				preNeurons.pop_back();
			}
		}

		/*Alternative treatment of inhibitory neurons*/
		// Initialize inhibitory synapses independently
		// Skip these when looking for *presynaptic* neurons
		// Specifically: Alter input to "randomIntVector()"
	}
}

void Network::addSynapseWeights()
{
	std::uniform_real_distribution<double> exDist(0, exWeight);
	std::uniform_real_distribution<double> inDist(inWeight, 0);

	for (int preNeuron = 0; preNeuron < N; preNeuron++)
	{
		for (int postNeuron = 0; postNeuron < N; postNeuron++)
		{
			if (S(preNeuron, postNeuron))
			{
				if (inhibNeurons[preNeuron])
					S(preNeuron, postNeuron) = inDist(engine);
				else
					S(preNeuron, postNeuron) = exDist(engine);
			}
		}
	}
}

/*Function methods*/

void Network::run() const
{
	// Activity variables
	Vector v = Vector(N, -65);
	Vector u = b.hadamard(v);

	// Record firings for a given time step
	std::vector<int> fired = std::vector<int>();
	fired.reserve(N);
	int numFired = 0;

	// Define input current vectors
	Vector I = Vector(N);

	// Start clock
	auto start = std::chrono::high_resolution_clock::now();
	int checkpoint = 100;
	if (checkpoint >= duration)
		std::cout << "Error in Network::run(): Checkpoint exceeds duration\n";

	for (int t = 0; t < checkpoint; t++)
	{
		// Stochastic thalamic input
		makeNoise(I);

		// Find spiking events
		fired = v.find(30);
		numFired += fired.size();

		// Update activity variables
		updateAvalanches(fired, t);
		v.update(c, fired);
		u.update(u + d, fired);
		synapticCurrent(I, fired);

		for (int i = 0; i < 2; i++)
			v += 0.5 * ((0.04 * (v ^ 2)) + (5 * v) + 140 - u + I);
		u += a.hadamard(b.hadamard(v) - u);
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	//t100 = time.count();
	//f100 = numFired;

	if (time.count() > checkpoint * survival)
	{
		kill();
		return;
	}

	for (int t = checkpoint; t < duration; t++)
	{

		// Stochastic thalamic input
		makeNoise(I);

		// Find spiking events
		fired = v.find(30);
		numFired += fired.size();

		// Update activity variables
		updateAvalanches(fired, t);
		v.update(c, fired);
		u.update(u + d, fired);
		synapticCurrent(I, fired);

		for (int i = 0; i < 2; i++)
			v += 0.5 * ((0.04 * (v ^ 2)) + (5 * v) + 140 - u + I);
		u += a.hadamard(b.hadamard(v) - u);
	}

	// Record last, ongoing avalanche
	fired.clear();
	updateAvalanches(fired, duration);
}

void Network::updateAvalanches(const std::vector<int> &fired, int t) const
{
	// Add most recent firing events to avalanche vector

	// If no firing: Add current number of active neurons to avalanche vector
	if (fired.empty())
	{
		if (!activeNeurons.size())
			return;

		avalanches.push_back(activeNeurons.size());
		activeNeurons.clear();
		return;
	}

	for (const int &neuron : fired)
		activeNeurons.insert(neuron);
}

void Network::visualize() const
{
	// Activity variables
	Vector v = Vector(N, -65);
	Vector u = b.hadamard(v);

	// Record all firings
	std::vector<int> fireTiming = std::vector<int>();
	fireTiming.reserve(10 * N);
	std::vector<int> fireNeuron = std::vector<int>();
	fireNeuron.reserve(10 * N);
	std::vector<int> fired = std::vector<int>();
	fired.reserve(N);

	// Define input current vectors
	Vector I = Vector(N);

	for (int t = 0; t < duration; t++)
	{
		// Stochastic thalamic input
		makeNoise(I);

		// Find spiking events
		fired = v.find(30);

		// Update activity variables
		updateFirings(fireTiming, fireNeuron, fired, t);
		v.update(c, fired);
		u.update(u + d, fired);
		synapticCurrent(I, fired);

		for (int i = 0; i < 2; i++)
			v += 0.5 * ((0.04 * (v ^ 2)) + (5 * v) + 140 - u + I);
		u += a.hadamard(b.hadamard(v) - u);
	}

	// Write firing events to file
	writeToFile("network-activity.dat", fireTiming, fireNeuron);
}

void Network::visualize2() const
{
	// Activity variables
	Vector v = Vector(N, -65);
	Vector u = b.hadamard(v);

	// Record firings for a given time step
	std::vector<int> fired = std::vector<int>();
	fired.reserve(N);

	// Define input current vectors
	Vector I = Vector(N);

	// Open file
	std::ofstream myfile;
	myfile.open("visualize2.dat", std::ios::trunc);

	for (int t = 0; t < duration; t++)
	{
		// Stochastic thalamic input
		makeNoise(I);

		// Find spiking events
		fired = v.find(30);

		// Write spikes to file
		for (const int& n : fired)
			myfile << t << " " << n << "\n";

		// Update activity variables
		v.update(c, fired);
		u.update(u + d, fired);
		synapticCurrent(I, fired);

		for (int i = 0; i < 2; i++)
			v += 0.5 * ((0.04 * (v ^ 2)) + (5 * v) + 140 - u + I);
		u += a.hadamard(b.hadamard(v) - u);
	}

	// Close file
	myfile.close();
}

double Network::analyze() const
{
	// Assumes we have avalanche size matrix
	if (avalanches.empty())
		cout << "Error in Network::analyze(): Avalanche size vector is empty\n";

	//std::pair<double, double> alphaKappa = searchAlpha(avalanches, 1);
	//std::pair<double, double> alphaFitness = clauset(avalanches);
	std::pair<double, double> alphaFitness = alphaFromKS(avalanches);
	alpha = alphaFitness.first;
	//kappa = alphaFitness.second;
	KS = alphaFitness.second;
	return KS;
}

void Network::makeNoise(Vector& I) const
{

	std::normal_distribution<double> exDist(0, exNoise);
	std::normal_distribution<double> inDist(0, exNoise);

	for (int neuron = 0; neuron < N; neuron++)
	{
		// Noise to excitatory neurons
		if (!inhibNeurons[neuron])
			I(neuron) = exDist(engine);

		// Noise to inhibitory neurons
		else
			I(neuron) = inDist(engine);
	}
}

void Network::synapticCurrent(Vector &I, const std::vector<int> &fired) const
{
	for (unsigned int i = 0; i < S.getRows(); i++)
	{
		for (auto it = fired.begin(); it != fired.end(); it++)
		{
			I(i) += S(i, *it);
		}
	}
}

std::vector<int> Network::getPreModuleVector(int connections, int postModule) const
{
	std::vector<int> preModules;
	preModules.reserve(connections);

	std::vector<int> levels = randPowerVector(connections, gamma, 1, H);
	for (auto it = levels.begin(); it != levels.end(); it++)
		preModules.push_back(getPreModule(postModule, H, *it));

	return preModules;

	// Alternative: Draw independent h's
	//for (int i = 0; i < connections; i++)
	//	preModules.push_back(getPreModule(postModule, hierarchyDepth));
	//return preModules;
}

void Network::clear() const
{
	avalanches.clear();
	activeNeurons.clear();

	// Also clear performance variables?
	alpha = -1;
	kappa = -1;
	KS = 1;
}

void Network::kill() const
{
	alive = false;
	std::cout << "Thread " << std::this_thread::get_id() << " killed\n";
	clear();
}

/*Export methods*/

std::vector<double> Network::exportGenome() const
{
	// Genome: Balance, H, minDegree, maxDegree, beta, gamma, exNoise, inNoise, exWeight, inWeight
	std::vector<double> genome = {
		balance,
		(double)H,
		(double)minDegree,
		(double)maxDegree,
		beta,
		gamma,
		exNoise,
		inNoise,
		exWeight,
		inWeight
	};
	return genome;
}

void Network::exportSynapses() const
{
	S.save("network-synapses.dat");
}

void Network::exportParameters() const
{
	std::ofstream myfile;
	myfile.open("network-parameters.txt", std::ios::trunc);
	myfile << currentDateTime();
	myfile << "\nNeurons: " << this->N
		<< "\nDuration: " << this->duration
		<< "\nInhibition ratio: " << this->balance
		<< "\nHierarchical levels: " << this->H
		<< "\nMinimum degree:" << minDegree
		<< "\nMaximum degree:" << maxDegree
		<< "\nBeta:" << beta
		<< "\nGamma:" << gamma
		<< "\nStd of noise to excitatory neurons: " << exNoise
		<< "\nStd of noise to inhibitory neurons: " << inNoise
		<< "\nMaximum weight of excitatory neurons: " << exWeight
		<< "\nMaximum weight of inhibitory neurons: " << inWeight
		<< "\n";
	myfile.close();
}

/*Static methods*/

bool Network::sortByKappa(const std::shared_ptr<Network> &a1, const std::shared_ptr<Network> &a2)
{
	// Question: is a1 better than a2?

	return abs(a1->kappa - 1) < abs(a2->kappa - 1);
}

bool Network::sortByKS(const std::shared_ptr<Network> &a1, const std::shared_ptr<Network> &a2)
{
	// Question: is a1 better than a2?

	return abs(a1->KS) < abs(a2->KS);
}

/*Random methods*/

std::mt19937 Network::startEngine(int seed)
{
	if (seed)
		return std::mt19937(seed);

	std::random_device rd;
	return std::mt19937(rd());
}

int Network::randInt(int a, int b) const
{
	std::uniform_int_distribution<> dist(a, b);
	return dist(engine);
}

double Network::randReal(double a , double b) const
{
	std::uniform_real_distribution<> dist(a, b);
	return dist(engine);
}

double Network::randNormal(double mu, double sigma) const
{
	std::normal_distribution<> dist(mu, sigma);
	return dist(engine);
}

std::vector<int> Network::randIntVector(int length, int min) const
{
	std::vector<int> result;
	result.reserve(length);
	for (int i = 0; i < length; i++)
	{
		result.push_back(i + min);
	}

	// Shuffle
	std::shuffle(result.begin(), result.end(), engine);

	return result;
}

std::vector<int> Network::randPowerVector(int length, double exponent, int min, int max) const
{
	std::vector<int> result;
	result.reserve(length + max);
	double normalization = 0;

	// Calculate normalization
	for (int x = min; x <= max; x++)
		normalization += pow(x, -exponent);
	normalization = 1 / normalization;

	double density;
	int number;

	// Fill vector
	for (int x = max; x >= min; x--)
	{
		density = pow(x, (-exponent)) * normalization;
		number = (int)floor(density * length); // Floor vs. round vs. ceil affects total length

		// Secure a tail
		if (number < 1)
			number = 1;

		for (int i = 0; i < number; i++)
			result.push_back(x);
	}

	// Correct number of elements (subject to rounding variations)
	while (result.size() > (size_t)length)
		result.pop_back();
	while (result.size() < (size_t)length)
		result.push_back(min);

	// Shuffle
	std::shuffle(result.begin(), result.end(), engine);

	return result;
}

/*External functions*/

std::vector<double> testGenome()
{
	return std::vector<double>
	{
		0.2,
		1,
		10,
		128,
		2,
		2,
		5,
		2,
		0.5,
		-1
	};
}

std::vector<double> goodGenome()
{
	return std::vector<double>
	{
		0.29,
		1,
		24,
		315,
		1.66352,
		1.7877,
		4.61092,
		2.37324,
		1.56954,
		-1
	};
}

std::vector<double> veryGoodGenome()
{
	return std::vector<double>
	{
		0.29,
		2,
		8,
		104,
		1.69871,
		1.5153,
		4.11085,
		2.55698,
		0.239477,
		-1.56041
	};
}

/*Tests*/

void testRandNormal()
{
	int n = 100000;
	std::vector<double> samples;
	samples.reserve(n);

	for (int i = 0; i < n; i++)
		samples.push_back(randomNormal(0, 5));
	
	/*static std::random_device device;
	static std::mt19937 engine(device());
	std::normal_distribution<> normal(0, 5); // This can't be declared static

	for (int i = 0; i < n; i++)
		samples.push_back(normal(engine));*/

	saveVector("testRandNormal.dat", samples);
}

void testVisualize()
{
	Network Moses;
	Moses.visualize2();
}

void testSeededCtor()
{
	Network Adam(1);
	//Adam.visualize();

	Network Eve(2);
	Network Cain;
	//Network Cain(Eve, Adam, 14);
	
	for (int seed = 1; seed < 1000; seed++)
	{
		std::cout << seed << "\n";
		Cain = Network(Eve, Adam, seed);
	}
}

void testRandPowerVector()
{
	int xMin = 1;
	int xMax = 20;
	double alpha = 2;
	int N = 100;
}