#pragma once
#include "Vector.h"
#include "Connectivity.h"
#include "Izhikevich.h"
#include "Criticality.h"
#include <chrono>
#include <thread>

class Network
{
	private:
		/*Const members*/
		const int N = 512;				// Number of neurons
		const int duration = 2000;		// Duration of simulation [ms]
		const int survival = 150;		// Threshold for which models live
		/*Genome members*/
		double balance;					// Ratio of inhibitory neurons
		int H;							// Hierarchical levels
		int minDegree, maxDegree;		// Degree parameters
		double beta;					// Degree distribution exponent
		double gamma;					// Intermodule connection density
		double exNoise, inNoise;		// Noise variance parameters
		double exWeight, inWeight;		// Maximum synaptic weights
		/*Network parameter members*/
		int M, n;						// Modules and neurons per module
		int Ne, Ni, ne, ni;				// Number of excitatory and inhibitory neurons
		std::vector<int> inhibNeurons;	// Vector indicating inhibitory neurons		// Consider Vector vs std::vector<int>
		/*Izhikevich members*/
		Vector a, b, c, d;				// Vectors of izhikevich parameters
		Matrix S;						// Matrix of synaptic weights

		/*Members that change during simulation*/
		mutable std::mt19937 engine;					// Random number generator
		mutable double alpha;							// Observed power-law exponent
		mutable double kappa;							// Observed kappa value
		mutable double KS;								// Observed Kolmogorov-Smirnov distance
		mutable std::vector<int> avalanches;			// Avalanche sizes
		mutable std::unordered_set<int> activeNeurons;	// Active neurons in current avalanche
		//mutable int t100;								// Elapsed time for 100ms simulation
		//mutable int f100;
		mutable bool alive = true;						// Indicate whether killed during simulation
	public:
		/*Constructors*/
		explicit Network(int seed = 0);
		Network(const Network &other);
		Network(const Network& mom, const Network& dad, int seed = 0);
		/*Operators*/
		Network operator= (Network rhs);
		void operator()();
		/*Initialization methods*/
		void initialize();
		void importGenome(const std::vector<double> &genome);
		void combineGenomes(const std::vector<double> &momGenome, const std::vector<double> &dadGenome);
		void spliceGenomes(const std::vector<double> &momGenome, const std::vector<double> &dadGenome);
		std::vector<double> newGenome() const;
		void defineNeurons();
		void mutate();
		void initializeIzhikevich();
		void createSynapses();
		void addSynapseWeights();
		/*Function methods*/
		void run() const;
		void updateAvalanches(const std::vector<int> &fired, int t) const;
		void visualize() const;
		void visualize2() const;
		double analyze() const;
		void makeNoise(Vector& I) const;
		void synapticCurrent(Vector &I, const std::vector<int> &fired) const;
		std::vector<int> getPreModuleVector(int connections, int postModule) const;
		void clear() const;
		void kill() const;
		/*Get methods*/
		double getKappa() const { return kappa; }
		double getKS() const { return KS; }
		int getNeurons() const { return N; }
		int getDuration() const { return duration; }
		//int getT100() const { return t100; }
		//int getF100() const { return f100; }
		bool isAlive() const { return alive; }
		/*Export methods*/
		std::vector<double> exportGenome() const;
		void exportSynapses() const;
		void exportParameters() const;
		std::vector<int> exportAvalanches() const { return avalanches; }
		/*Static methods*/
		static bool sortByKappa(const std::shared_ptr<Network> &a1, const std::shared_ptr<Network> &a2);
		static bool sortByKS(const std::shared_ptr<Network> &a1, const std::shared_ptr<Network> &a2);
		/*Random methods*/
		std::mt19937 startEngine(int seed = 0);
		int randInt(int a = 0, int b = 1) const;
		double randReal(double a = 0, double b = 1) const;
		std::vector<int> randIntVector(int length, int min = 0) const;
		double randNormal(double mu = 0, double sigma = 1) const;
		std::vector<int> randPowerVector(int length, double exponent, int min, int max) const;
};

/*Tests*/
void testRandNormal();
void testVisualize();
void testSeededCtor();
std::vector<double> testGenome();
std::vector<double> goodGenome();
std::vector<double> veryGoodGenome();