#pragma once
#include "Vector.h"
#include "Connectivity.h"
#include "Izhikevich.h"
#include "Criticality.h"
#include <chrono>

class Network
{
	private:
		const int N = 512;				// Number of neurons
		double balance;					// Ratio of inhibitory neurons
		int Ne, Ni, ne, ni;				// Number of excitatory and inhibitory neurons
		int H, M, n;					// Hierarchy parameters
		int minDegree, maxDegree;		// Degree parameters
		double beta;					// Degree distribution exponent
		double gamma;					// Intermodule connection density
		Matrix S;						// Matrix of synaptic weights
		std::vector<int> inhibNeurons;	// Vector indicating inhibitory neurons		// Consider Vector vs std::vector<int>
		Vector a, b, c, d;				// Vectors of izhikevich parameters
		int duration = 400;			// Duration of simulation [ms]
		double exNoise, inNoise;		// Noise variance parameters
		double exWeight, inWeight;		// Maximum synaptic weights

		/*Members that change during simulation*/
		mutable double alpha;							// Observed power-law exponent
		mutable double kappa;							// Observed kappa value
		mutable double KS;								// Observed Kolmogorov-Smirnov distance
		mutable std::vector<int> avalanches;			// Avalanche sizes
		mutable std::unordered_set<int> activeNeurons;	// Active neurons in current avalanche
	public:
		/*Constructors*/
		Network();
		Network(const Network &other);
		Network(const Network& mom, const Network& dad);
		explicit Network(int neurons);
		/*Operators*/
		Network operator= (Network rhs);
		void operator()();
		/*Initialization methods*/
		void importGenome(const std::vector<double> &genome);
		void newGenome();
		void testGenome();
		void initialize();
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
		void clear() const;
		/*Get methods*/
		double getKappa() const { return kappa; }
		double getKS() const { return KS; }
		int getNeurons() const { return N; }
		int getDuration() const { return duration; }
		/*Export methods*/
		std::vector<double> exportGenome() const;
		void exportSynapses() const;
		void exportParameters() const;
		std::vector<int> exportAvalanches() const { return avalanches; }
		/*Static methods*/
		static bool sortByKappa(const std::shared_ptr<Network> &a1, const std::shared_ptr<Network> &a2);
		static bool sortByKS(const std::shared_ptr<Network> &a1, const std::shared_ptr<Network> &a2);
};

/*Tests*/
void testRandNormal();
void testVisualize();
void testNetworkCopyCtor();
