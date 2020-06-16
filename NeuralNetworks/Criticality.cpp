#include "Criticality.h"

std::vector<int> avalancheDetection(const std::vector<int> &times, const std::vector<int> &neurons)
{
	// Current implementation stricly uses time bin width of one time step (1 ms)

	// Check that data is correct
	if (times.size() != neurons.size())
		cout << "Error in avalancheDetection: Differently-sized arrays\n";

	// Initialize variables
	int t = -1;
	std::vector<int> sizes;
	sizes.reserve(1024);
	std::unordered_set<int> activeElectrodes;
	activeElectrodes.reserve(1024);

	// Iterate through all events
	for (size_t i = 0; i < times.size(); i++)
	{
		// Add electrodes for a given time steps
		if (times[i] == t)
			activeElectrodes.insert(neurons[i]);

		// Ongoing avalanche, new time step
		else if (times[i] == t + 1)
			activeElectrodes.insert(neurons[i]);

		// Jump in time step: Store current avalanche size
		else
		{
			sizes.push_back(activeElectrodes.size());
			activeElectrodes.clear();
		}
		t = times[i];
	}

	// Add last avalanche
	if (!activeElectrodes.empty())
		sizes.push_back(activeElectrodes.size());

	return sizes;
}

std::vector<double> readVector(const char* filename)
{
	std::vector<double> result;
	double input;

	std::string line;
	std::ifstream myfile(filename);
	if (myfile.is_open())
		while (myfile >> input)
			result.push_back(input);
	myfile.close();

	return result;
}

double calculateKappa(const std::vector<double> &pdf_theory, const std::vector<double> &pdf_emp, int x0, int x1)
{
	int m = 10;
	std::vector<int> beta = logSpacing(x0, x1, m);
	int n = x1 - x0 + 1;

	// Make CDF vectors
	std::vector<double> cdf_theory(n);
	std::vector<double> cdf_emp(n);

	for (int i = 0; i < n; i++)
	{
		cdf_theory[i] = partSum(pdf_theory, x0 - 1, x0 + i - 1);
		cdf_emp[i] = partSum(pdf_emp, x0 - 1, x0 + i - 1);
	}

	double kappa = 0;
	for (auto b = beta.begin(); b != beta.end(); b++)
		kappa += cdf_theory[*b] - cdf_emp[*b];
	
	return 1 + 1 / (double)m * kappa;
}

std::vector<int> logSpacing(int x0, int x1, int m)
{
	std::vector<int> result(m, 0);
	for (int i = 0; i < m; i++)
		result[i] = (int)round(pow(10, log10(x0) + i / (double)(m - 1) * log10(x1))) - 1; // -1 for 0-indexing
	
	return result;
}

double calculateKS(const std::vector<double> &pdf_theory, const std::vector<double> &pdf_emp, int x0, int x1)
{
	int n = x1 - x0 + 1;
	double max = 0;

	// Make CDF vectors
	std::vector<double> cdf_theory(n);
	std::vector<double> cdf_emp(n);

	for (int i = 0; i < n; i++)
	{
		cdf_theory[i] = partSum(pdf_theory, x0 - 1, x0 + i - 1);
		cdf_emp[i] = partSum(pdf_emp, x0 - 1, x0 + i - 1);

		if (abs(cdf_theory[i] - cdf_emp[i]) > max)
			max = abs(cdf_theory[i] - cdf_emp[i]);
	}

	return max;
}

double calculateKS(const std::vector<double> &pdf_theory, const std::vector<double> &pdf_emp)
{
	int n = pdf_emp.size();
	if (pdf_theory.size() != n)
		cout << "Error in calculateKS: PDF vectors not of equal length\n";
	double max = 0;

	// Make CDF vectors
	std::vector<double> cdf_theory(n);
	std::vector<double> cdf_emp(n);

	for (int i = 0; i < n; i++)
	{
		cdf_theory[i] = partSum(pdf_theory, 0, i);
		cdf_emp[i] = partSum(pdf_emp, 0, i);

		if (abs(cdf_theory[i] - cdf_emp[i]) > max)
			max = abs(cdf_theory[i] - cdf_emp[i]);
	}

	return max;
}

double mle(const std::vector<double> &pdf, int x0)
{
	// Uses pdf, converts into samples
	
	// Only alphas within a given range should be accepted, e.g. [1, 3]

	int n = pdf.size();
	double alpha = 0;

	// Scale up from pdf to histogram of observations
	std::vector<int> observations = histogramToSamples(pdf);

	for (auto it = observations.begin(); it != observations.end(); it++)
	{
		//alpha += log(*it / (double)(x0 - 0.5));	// Discrete version
		alpha += log(*it / (double)x0);				// Continuous version
	}

	return 1 + n * (1 / double(alpha));
}

double mle(const std::vector<int> &sizes, int x0)
{
	// Uses sizes, a list of all observations (not a histogram)

	// Only alphas within a given range should be accepted, e.g. [1, 3]

	int n = sizes.size();
	double alpha = 0;

	for (auto it = sizes.begin(); it != sizes.end(); it++)
	{
		//alpha += log(*it / (double)(x0 - 0.5));	// Discrete version
		alpha += log(*it / (double)x0);				// Continuous version
	}

	return 1 + n * (1 / double(alpha));
}

double mleSearch(const std::vector<int> &sizes, int x0)
{
	// Initialize parameters
	int x1 = -1;
	int n = sizes.size();
	double g = 0;

	// Calculate geometric mean
	for (const int& x : sizes)
	{
		g += log(x);
		if (x > x1)
			x1 = x;
	}
	g /= (double)n;

	return maxLikelihood(g, x0, x1);
}

double mleSearch(const std::vector<double> &sizes, int x0)
{
	// Initialize parameters
	int x1 = -1;
	int n = sizes.size();
	double g = 0;

	// Calculate geometric mean
	for (const double& x : sizes)
	{
		g += log(x);
		if (x > x1)
			x1 = (int)ceil(x);
	}
	g /= (double)n;

	return maxLikelihood(g, x0, x1);
}

double maxLikelihood(double g, int x0, int x1)
{
	// Initialize parameters
	double alphaHat = -1.0;
	double l = -1;
	double maxVal = -999999999;
	double start = 1.1;

	// Perform course scan
	double stepSize = 0.1;
	for (double alpha = 1.1; alpha < start + 2; alpha += stepSize)
	{
		l = log(alpha - 1) - log(pow(x0, 1 - alpha) - pow(x1, 1 - alpha)) - alpha * g;
		if (l > maxVal)
		{
			maxVal = l;
			alphaHat = alpha;
		}
	}

	// Perform fine scan
	stepSize = 0.01;
	start = alphaHat - 0.1;
	for (double alpha = start; alpha < start + 0.2; alpha += stepSize)
	{
		l = log(alpha - 1) - log(pow(x0, 1 - alpha) - pow(x1, 1 - alpha)) - alpha * g;
		if (l > maxVal)
		{
			maxVal = l;
			alphaHat = alpha;
		}
	}
	return alphaHat;
}

std::pair<double, double> searchAlpha(const std::vector<double> &pdf_emp, int x0)
{
	int n = pdf_emp.size();
	int xShift = 1;
	double alpha;
	double kappa;
	std::pair<double, double> alphaKappa( 0, 100.0 );
	std::vector<double> pdf_theory;

	for (int xmin = x0; xmin < x0 + xShift; xmin++) // Does not work for xmin > 1 !!!! 
	{
		alpha = mle(pdf_emp, xmin);
		pdf_theory = makePowerLawDistribution(alpha, n); // Currently assumes x0 = 1
		kappa = calculateKappa(pdf_theory, pdf_emp, xmin, n);
		if (abs(kappa - 1) < abs(alphaKappa.second - 1))
			alphaKappa = std::make_pair(alpha, kappa);
	}
	return alphaKappa;
}

std::pair<double, double> searchAlpha(const std::vector<int> &sizes, int x0)
{
	// Calculate optimal alpha via MLE, and corresponding Kolmogorov-Smirnov distance
	// Currently, only x0 = 1 is tested. Can expand to several x_min
	// Challenge: compare only parts of empirical PDF, not entire

	std::vector<double> pdf_emp = samplesToPdf(sizes);

	int n = pdf_emp.size();
	int xShift = 1;
	double alpha;
	double kappa;
	std::pair<double, double> alphaKappa(0, 100.0);
	std::vector<double> pdf_theory;

	for (int xmin = x0; xmin < x0 + xShift; xmin++) // Does not work for xmin > 1 !!!! 
	{
		alpha = mle(sizes, xmin);
		pdf_theory = makePowerLawDistribution(alpha, n); // Currently assumes x0 = 1
		kappa = calculateKappa(pdf_theory, pdf_emp, xmin, n);
		if (abs(kappa - 1) < abs(alphaKappa.second - 1))
			alphaKappa = std::make_pair(alpha, kappa);
	}
	return alphaKappa;
}

std::pair<double, double> clauset(const std::vector<int> &sizes)
{
	// Calculate optimal alpha via MLE, and corresponding Kolmogorov-Smirnov distance
	// Possible expansion: Scan several x_min
	// Challenge: compare only parts of empirical PDF, not entire

	std::vector<double> pdf_emp = samplesToPdf(sizes);

	int n = pdf_emp.size();
	double alpha = mle(sizes);
	std::pair<double, double> alphaKS(0, 100.0);
	std::vector<double> pdf_theory = makePowerLawDistribution(alpha, n);
	double KS = calculateKS(pdf_theory, pdf_emp);
	
	return std::make_pair(alpha, KS);
}

std::vector<double> makePowerLawDistribution(double alpha, int length)
{
	std::vector<double> pdf(length);
	for (int i = 0; i < length; i++)
		pdf[i] = pow(i + 1, -alpha);

	double vecSum = sumVector(pdf);
	for (int i = 0; i < length; i++)
		pdf[i] /= vecSum;

	return pdf;
}

std::pair<double, double> alphaFromKS(const std::vector<int>& sizes, int x0)
{
	// Initialize parameters
	std::pair<double, double> alphaKS(-1, 100.0);
	double KS;
	double start = 1.1;

	// Initialize variables
	std::vector<double> pdf_emp = samplesToPdf(sizes);
	int n = pdf_emp.size();
	std::vector<double> pdf_theory;

	// Perform course scan
	double stepSize = 0.1;
	for (double alpha = 1.1; alpha < start + 2; alpha += stepSize)
	{
		pdf_theory = makePowerLawDistribution(alpha, n); // Currently assumes x0 = 1
		KS = calculateKS(pdf_theory, pdf_emp, x0, n);
		if (abs(KS) < abs(alphaKS.second))
			alphaKS = std::make_pair(alpha, KS);
	}

	// Perform fine scan
	stepSize = 0.01;
	start = alphaKS.first - 0.1;
	for (double alpha = start; alpha < start + 0.2; alpha += stepSize)
	{
		pdf_theory = makePowerLawDistribution(alpha, n); // Currently assumes x0 = 1
		KS = calculateKS(pdf_theory, pdf_emp, x0, n);
		if (abs(KS) < abs(alphaKS.second))
			alphaKS = std::make_pair(alpha, KS);
	}
	return alphaKS;
}

/*Tests*/
void testReadVector()
{
	std::vector<double> cdf_theory = readVector("cdf_theory.dat");
	for (auto it = cdf_theory.begin(); it != cdf_theory.end(); it++)
		cout << *it << " ";
	cout << "\n";
}

void testKappa()
{
	std::vector<double> pdf_theory = readVector("pdf_theory.dat");
	std::vector<double> pdf_emp = readVector("pdf_emp.dat");
	double kappa = calculateKappa(pdf_theory, pdf_emp, 1, pdf_emp.size());
	cout << std::setprecision(5) << std::fixed;
	cout << "Kappa is " << kappa << "\n";
}

void testKS()
{
	std::vector<double> pdf_theory = readVector("pdf_theory.dat");
	std::vector<double> pdf_emp = readVector("pdf_emp.dat");
	double KS = calculateKS(pdf_theory, pdf_emp);
	cout << std::setprecision(5) << std::fixed;
	cout << "D_KS is " << KS << "\n";
}

void testSearchAlpha()
{
	std::vector<double> pdf_emp = readVector("pdf_emp.dat");
	std::pair<double, double> alphaKappa = searchAlpha(pdf_emp, 1);
	cout << "Alpha = " << alphaKappa.first << " gives kappa = " << alphaKappa.second << "\n";
}

void testAvalancheDetection()
{
	std::vector<int> fireTiming =  { 0, 1, 1, 2, 2, 2, 2, 3, 4, 4, 6, 7, 9, 10, 10, 10, 11 };
	std::vector<int> fireNeurons = { 0, 1, 2, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3,  0,  1,  2,  3 };
	std::vector<int> sizes = avalancheDetection(fireTiming, fireNeurons);

	cout << "Avalanche sizes: ";
	for (const int& s : sizes)
		cout << s << " ";
	cout << "\n";
}

void testMle()
{
	//std::vector<int> samples = powerLawVector(100, 2, 10, 1);
	std::vector<int> samples = powerLawVectorSampleAndInsert(1000, 2, 30, 1);
	saveVector("mle-vector.dat", samples);
	//std::vector<double> samples = readVector("samples.dat");
	//std::vector<double> samples = readVector("mle-vector.dat");
	double alpha = mleSearch(samples);
	cout << "Alpha: " << alpha << "\n";
}