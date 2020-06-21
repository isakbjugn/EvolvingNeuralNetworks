#include "Network.h"

/*Constructors*/

Network::Network()
{
	// Make genome
	newGenome();
	//testGenome();

	// Initilize network structure
	initialize();
	initializeIzhikevich();
	createSynapses();
	addSynapseWeights();

	// Reset performance variables
	clear();
}

Network::Network(const Network &other)
{
	// Copy members
	importGenome(other.exportGenome());
	initialize();
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
}

Network::Network(int neurons)
{
	newGenome();
	
	// Initialize
	H = 1;
	n = N;
	Ne = n;
	Ni = 0;
	ne = n;
	ni = 0;
	inhibNeurons = std::vector<int>(neurons, 0);
	a = Vector(neurons, 1);
	b = Vector(neurons, 2);
	c = Vector(neurons, 3);
	d = Vector(neurons, 4);

	// Initialize connectivity matrix
	S = Matrix(neurons);
	S.populate();

	// Initialize performance variables
	clear();
	avalanches = std::vector<int>{ 1, 1, 1, 2, 2, 3 };
	activeNeurons = std::unordered_set<int>();
}

Network::Network(const Network& mom, const Network& dad)
{
	// Combine genomes
	std::vector<double> momGenome = mom.exportGenome();
	std::vector<double> dadGenome = dad.exportGenome();
	int genes = momGenome.size();
	int splice = randomIntWithLimits(0, genes - 1);

	std::vector<double> newGenome(genes);
	for (int i = 0; i < splice; i++)
		newGenome[i] = momGenome[i];
	for (int i = splice; i < genes; i++)
		newGenome[i] = dadGenome[i];

	// Incorporate genome, with mutations
	importGenome(newGenome);
	mutate();

	// Avoid more synapses than neurons in module (with high gamma)
	if (maxDegree > n)
		maxDegree = n;

	// Initialize network structure
	initialize();
	initializeIzhikevich();
	createSynapses();
	addSynapseWeights();

	// Reset performance variables
	clear();
}

/*Operators*/

Network Network::operator= (Network rhs)
{
	std::swap(this->a, rhs.a);
	std::swap(this->b, rhs.b);
	std::swap(this->c, rhs.c);
	std::swap(this->d, rhs.d);
	std::swap(this->S, rhs.S);
	std::swap(this->S, rhs.S);
	return rhs;
}

void Network::operator()()
{
	run(); // So far so good
	analyze();
}

/*Initialization methods*/

void Network::newGenome()
{
	// Seems to be error when maxDegree is equal to number of neurons total

	H = randomIntWithLimits(1, 5);
	M = (int)pow(2, H - 1);
	n = N / M;
	balance = 0.01 * randomIntWithLimits(20, 30);
	minDegree = randomIntWithLimits(1, 30);
	maxDegree = randomIntWithLimits(minDegree, n);

	//if (H == 1)
	//	maxDegree = randomIntWithLimits(minDegree, N/2);

	beta = randomRealWithLimits(1.0, 2.0); // Massobrio et al. (2015)
	gamma = randomRealWithLimits(1.5, 3.0); // Rubinov et al. (2011)
	exNoise = randomRealWithLimits(4.0, 6.0);
	inNoise = randomRealWithLimits(1.0, 3.0);
	exWeight = randomRealWithLimits(0.1, 2.0);
	inWeight = randomRealWithLimits(-2.0, -0.1);

	/*Other parameters discussed*/
	// Izhikevich parameters a, b, c, d (mean, variance)
	// Incoming degree mean (rather than exponent, min and max)	
}

void Network::testGenome()
{
	H = 1;
	M = (int)pow(2, H - 1);
	n = N / M;
	balance = 0.2;
	minDegree = 10;
	maxDegree = 128;
	beta = 0.0;
	gamma = 2;
	exNoise = 5;
	inNoise = 2;
	exWeight = 0.5;
	inWeight = -1;
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
}

void Network::initialize()
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
	double mutateChance = 0.02;

	if (randUniform() < mutateChance)
		H = randomIntWithLimits(1, 5);

	M = (int)pow(2, H - 1);
	n = N / M;

	if (randUniform() < mutateChance)
		balance = 0.01 * randomIntWithLimits(20, 30);
	if (randUniform() < mutateChance)
		minDegree = randomIntWithLimits(1, 30);
	if (randUniform() < mutateChance)
		maxDegree = randomIntWithLimits(minDegree, n);
	if (randUniform() < mutateChance)
		beta = randomRealWithLimits(1.0, 2.0); // Massobrio et al. (2015)
	if (randUniform() < mutateChance)
		gamma = randomRealWithLimits(1.5, 3.0); // Rubinov et al. (2011)
	if (randUniform() < mutateChance)
		exNoise = randomRealWithLimits(4.0, 6.0);
	if (randUniform() < mutateChance)
		inNoise = randomRealWithLimits(1.0, 3.0);
	if (randUniform() < mutateChance)
		exWeight = 0.5;
	if (randUniform() < mutateChance)
		inWeight = -1.0;
}

void Network::initializeIzhikevich()
{
	// Random vectors
	Vector re = Vector(Ne, 0, 1);
	Vector ri = Vector(Ni, 0, 1);

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
	std::vector<int> degrees = powerLawVector(N, beta, maxDegree, minDegree);
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
		preModules = getPreModuleVector(degree, postModule, H);
		moduleCount = countOccurences(preModules);

		for (auto it = moduleCount.begin(); it != moduleCount.end(); it++)
		{
			preModule = it->first;

			if (it->first == postModule && it->second >= n)
				cout << "Error in createSynapses(): More synapses than neurons in module\n";

			// Treat excitatory and inhibitory neurons equally
			preNeurons = randomIntVector(n, preModule*n);

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
	for (auto preNeuron = inhibNeurons.begin(); preNeuron != inhibNeurons.end(); preNeuron++)
	{
		for (int postNeuron = 0; postNeuron < N; postNeuron++)
		{
			if (S(*preNeuron, postNeuron))
			{
				if (*preNeuron)
					S(*preNeuron, postNeuron) = randomRealWithLimits(inWeight, 0);
				else
					S(*preNeuron, postNeuron) = randomRealWithLimits(0, exWeight);
			}
		}
	}		
}

/*Function methods*/

void Network::run() const
{
	// Reset activity variables?
	//clear();

	// Activity variables
	Vector v = Vector(N, -65);
	Vector u = b.hadamard(v);

	// Record firings for a given time step
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
		updateAvalanches(fired, t);
		v.update(c, fired);
		u.update(u + d, fired);
		I += synapticCurrent(S, fired);

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
		I += synapticCurrent(S, fired);

		for (int i = 0; i < 2; i++)
			v += 0.5 * ((0.04 * (v ^ 2)) + (5 * v) + 140 - u + I);
		u += a.hadamard(b.hadamard(v) - u);
	}

	// Write firing events to file
	writeToFile("visualize1.dat", fireTiming, fireNeuron);
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
		I += synapticCurrent(S, fired);

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
	std::pair<double, double> alphaFitness = clauset(avalanches);
	alpha = alphaFitness.first;
	//kappa = alphaFitness.second;
	KS = alphaFitness.second;
	return KS;
}

void Network::makeNoise(Vector& I) const
{
	for (int neuron = 0; neuron < N; neuron++)
	{
		// Noise to excitatory neurons
		if (!inhibNeurons[neuron])
			I(neuron) = randNormal(0, exNoise);

		// Noise to inhibitory neurons
		else
			I(neuron) = randNormal(0, inNoise);
	}
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
	myfile << &std::chrono::system_clock::now();
	myfile << "\nBalance: " << this->balance
		<< "\nH: " << this->H
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

	return a1->KS < a2->KS;
}


/*Tests*/

void testRandNormal()
{
	int n = 100000;
	std::vector<double> samples;
	samples.reserve(n);

	for (int i = 0; i < n; i++)
		samples.push_back(randNormal(0, 5));
	
	/*static std::random_device device;
	static std::mt19937 generator(device());
	std::normal_distribution<> normal(0, 5); // This can't be declared static

	for (int i = 0; i < n; i++)
		samples.push_back(normal(generator));*/

	saveVector("testRandNormal.dat", samples);
}

void testVisualize()
{
	Network Moses;
	Moses.visualize2();
}

void testNetworkCopyCtor()
{
	Network Adam(4);
	Network Abel = Network(Adam);
	std::vector<double> AdamGenome = Adam.exportGenome();
	std::vector<double> AbelGenome = Abel.exportGenome();

	cout << "Genome differences: ";
	for (size_t i = 0; i < AdamGenome.size(); i++)
		cout << AdamGenome[i] - AbelGenome[i] << " ";

	std::vector<int> AdamAvalanches = Adam.exportAvalanches();
	std::vector<int> AbelAvalanches = Abel.exportAvalanches();
	
	cout << "\nNum Avalanche entries: " << AdamAvalanches.size() << " vs. " << AbelAvalanches.size();
	cout << "\nAvalanche differences: ";
	for (size_t i = 0; i < AdamAvalanches.size(); i++)
		cout << AdamAvalanches[i] - AbelAvalanches[i] << " ";

	cout << "\nKappa difference:" << Adam.getKappa() - Abel.getKappa();
	cout << "\n\n";
}