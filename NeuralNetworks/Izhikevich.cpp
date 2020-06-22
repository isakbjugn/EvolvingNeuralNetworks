#include "Izhikevich.h"

void izhikevich()
{
	int N = 100, Ne = 80, Ni = 20;
	int duration = 1000; // Number of milliseconds

	// Random vectors
	Vector re = Vector(Ne, 0, 1);
	Vector ri = Vector(Ni, 0, 1);

	// Excitatory parameters			Inhibitory parameters
	Vector a(Vector(Ne, 0.02),			Vector(0.02 + 0.08*ri));
	Vector b(Vector(Ne, 0.2),			Vector(0.25 - 0.05*ri));
	Vector c(Vector(-65 + 15 * (re^2)), Vector(Ni, -65));
	Vector d(Vector(8 - 6 * (re^2)),	Vector(Ni, 2));

	// Synaptic weights
	Matrix S = Matrix(N, Ne, 0, 0.5).expand(Matrix(N, Ni, -1, 0));
	std::cout << "S:\n" << S << "\n";

	// Activity variables
	Vector v = Vector(N, -65);
	Vector u = b.hadamard(v);

	// Record all firings
	std::vector<int> fireTiming = std::vector<int>();
	fireTiming.reserve(10*N);
	std::vector<int> fireNeuron = std::vector<int>();
	fireNeuron.reserve(10*N);
	std::vector<int> fired = std::vector<int>();
	fireTiming.reserve(N);

	// Define input current vectors
	Vector I = Vector(N);
	Vector I1 = Vector(Ne, 0.0);
	Vector I2 = Vector(Ni, 0.0);

	// Record single neuron
	std::vector<double> neuron = std::vector<double>();
	neuron.reserve(duration);

	// Simulation of 1000 ms
	for (int i = 0; i < duration; i++)
	{
		// Record from single neuron
		neuron.push_back(v(1));

		// Stochastic thalamic input
		I1.randNormal(0, 4);
		I2.randNormal(0, 1.5);
		I.copy(I1, I2);

		// Find spiking events
		fired = v.find(30);
		
		// Update activity variables
		updateFirings(fireTiming, fireNeuron, fired, i);
		v.update(c, fired);
		u.update(u + d, fired);
		I += synapticCurrent(S, fired);
		
		for (int i = 0; i < 2; i++)
			v += 0.5 * ((0.04 * (v ^ 2)) + (5 * v) + 140 - u + I);
		u += a.hadamard(b.hadamard(v) - u);
	}
	writeToFile("spikes.dat", fireTiming, fireNeuron);
	writeNeuronToFile("neuron.dat", neuron);
}

void updateFirings(Matrix & firings, const Vector & fired, int time)
{
	Matrix timestamp = Matrix(fired.getRows(), 1, time);
	timestamp.expand(fired);
	firings.extend(timestamp);
}

void updateFirings(std::vector<int> & fireTiming, std::vector<int>& fireNeuron, const std::vector<int> fired, int time)
{
	for (unsigned int i = 0; i < fired.size(); i++)
	{
		fireTiming.push_back(time);
		fireNeuron.push_back(fired[i]);
	}
}

Vector synapticCurrent(const Matrix & S, const std::vector<int> & fired)
{
	Vector I = Vector(S.getRows(), 0.0);

	// This is incorrect: Every neuron gets synaptic input, but from wrong weight
	/*for (unsigned int i = 0; i < S.getRows(); i++)
	{
		for (auto it = fired.begin(); it != fired.end(); it++)
		{
			I(i) += S(i, *it);
		}
	}*/

	// This is incorrect: Only firing neurons get synaptic input, from all others
	/*int postNeuron;
	for (unsigned int i = 0; i < fired.size(); i++)
	{
		postNeuron = int(fired[i]);
		for (unsigned int j = 0; j < S.getColumns(); j++)
		{
			I(postNeuron) += S(i, j);
		}
	}*/

	// Correct: Postsynaptic neurons get input from firing presynaptic neurons
	for (unsigned int postNeuron = 0; postNeuron < S.getRows(); postNeuron++)
	{
		for (auto preNeuron = fired.begin(); preNeuron != fired.end(); preNeuron++)
		{
			//I(postNeuron) += S(*preNeuron, postNeuron);
			I(postNeuron) += S(postNeuron, *preNeuron);
		}
	}

	return I;
}

void writeToFile(const char* filename, std::vector<int> & fireTiming, std::vector<int>& fireNeuron)
{
	std::ofstream myfile;
	myfile.open(filename, std::ios::trunc);

	for (unsigned int i = 0; i < fireTiming.size(); i++)
	{
		if (fireTiming[i] >= 0)
		{
			myfile << fireTiming[i] << " " << fireNeuron[i] << "\n";
		}
	}
	myfile.close();
}

void writeNeuronToFile(const char* filename, std::vector<double> & neuron)
{
	std::ofstream myfile;
	myfile.open(filename, std::ios::trunc);

	for (unsigned int i = 0; i < neuron.size(); i++)
	{
		myfile << neuron[i] << endl;
	}
	myfile.close();
}

void printVector(const std::vector<int> vec)
{
	for (unsigned int i = 0; i < vec.size(); i++)
	{
		cout << vec[i] << " ";
	}
	cout << endl;
}


// Tests

void testOperations()
{
	Vector c = Vector(4, -65);
	Vector d = Vector(4, 2);

	Vector v = Vector(4, 20);
	v(1) = 30;
	v(2) = 40;
	Vector u = 0.2 * v;
	cout << "v:\n" << v << endl << "u:\n" << u << endl;

	std::vector<int> fired = v.find(30);
	v.update(c, fired);
	u.update(u + d, fired);
	cout << "v:\n" << v << endl << "u:\n" << u << endl;

	Matrix S = Matrix(4, 4);
	S.randomize(0, 0.5);
	Vector I = Vector(4, 0, 5);
	I += synapticCurrent(S, fired);
	cout << "S:\n" << S << endl << "I:\n" << I << endl;

	for (int i = 0; i < 2; i++)
		v += 0.5 * ((0.04 * (v ^ 2)) + (5 * v) + 140 - u + I);
	u += 0.02*(0.2*v - u);
	cout << "v:\n" << v << endl << "u:\n" << u << endl;

}

void testSynapticCurrent()
{
	std::vector<int> fired = std::vector<int>();
	fired.push_back(1);
	fired.push_back(2);
	Matrix S = Matrix(4, 4);
	S.populate();
	S.populate();
	Vector I = Vector(4);
	I.populate();
	I.populate();
	cout << "I_thal:\n" << I << endl
		<< "S:\n" << S << endl
		<< "I_syn:\n" << synapticCurrent(S, fired) << endl;
	I += synapticCurrent(S, fired);
	cout << "I_total:\n" << I;
}