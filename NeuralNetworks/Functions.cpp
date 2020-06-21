#include "Functions.h"

int randomIntWithLimits(int lower, int upper)
{
	static std::random_device device;
	static std::mt19937 generator(device());
	std::uniform_int_distribution<int> uniform(lower, upper);
	return uniform(generator);
}

double randomRealWithLimits(double lower, double upper)
{
	static std::random_device device;
	static std::mt19937 generator(device());
	std::uniform_real_distribution<> uniform(lower, upper);
	return uniform(generator);
}

double randUniform()
{
	return randomRealWithLimits(0, 1);
}

double randNormal(double mean, double stddev)
{
	static std::random_device device;
	static std::mt19937 generator(device());
	std::normal_distribution<> normal(mean, stddev); // This can't be declared static

	return normal(generator);
}

int randPowerDiscrete(double alpha, int xMin, int xMax)
{
	return (int)round(randPower(alpha, xMin - 0.5, xMax - 0.5) + 0.5);
	// Remark: Using "1/2" yields error from implicit rounding of ints
}

double randPower(double alpha, double xMin, double xMax)
{
	double u = randUniform();
	return pow((pow(xMax, 1 - alpha) - pow(xMin, 1 - alpha)) * u + pow(xMin, 1 - alpha), 1 / (double)(1 - alpha));
}

std::vector<int> powerLawVector(int n, double alpha, int x1, int x0)
{
	// Fill-and-shuffle method

	std::vector<int> result;
	result.reserve(n + x1);
	double normalization = 0;

	// Calculate normalization
	for (int x = x0; x <= x1; x++)
		normalization += pow(x, -alpha);
	normalization = 1 / normalization;

	double density;
	int number;

	// Fill vector
	for (int x = x0; x <= x1; x++)
	{
		density = pow(x, (-alpha)) * normalization;
		number = (int)floor(density * n); // Floor vs. round vs. ceil affects total length

		// Secure a tail
		if (number < 1)
			number = 1;

		for (int i = 0; i < number; i++)
			result.push_back(x);
	}

	// Correct number of elements (subject to rounding variations)
	while (result.size() > n)
		result.erase(result.begin());
	while (result.size() < n)
		result.push_back(x0);

	// Shuffle
	static std::random_device device;
	static std::default_random_engine generator(device());
	std::shuffle(result.begin(), result.end(), generator);

	return result;
}

std::vector<int> histogramToSamples(const std::vector<double>& pdf)
{
	// Converts histogram into possible observations

	// Find scaling factor by finding smallest number
	double smallest = 1.0;
	double epsilon = pow(10, -7);
	int n = pdf.size();
	for (auto it = pdf.begin(); it != pdf.end(); it++)
	{
		if (*it < epsilon)
			n--;
		else if (*it < smallest)
			smallest = *it;
	}

	// Correct number of elements for scaling factor
	n *= (int)(1.0 / smallest);

	// Fill-and-shuffle method
	std::vector<int> samples;
	samples.reserve(n);

	// Fill vector
	for (size_t x = 1; x <= pdf.size(); x++)
	{
		for (int i = 0; i < (int)(pdf[x-1] * n); i++)
			samples.push_back(x);
	}

	// Shuffle
	static std::random_device device;
	static std::default_random_engine generator(device());
	std::shuffle(samples.begin(), samples.end(), generator);

	return samples;
}

std::vector<double> samplesToPdf(const std::vector<int>& samples)
{
	// Make histogram from samples
	std::map<int, int> hist = countOccurences(samples);
	int norm = summarizeMap(hist);

	// Normalize histogram into PDF
	std::vector<double> pdf;
	pdf.reserve(hist.size());
	int idx = 1; // Assumes smallest possible sample is 1
	for (auto it = hist.begin(); it != hist.end(); it++)
	{
		while (it->first != idx)
		{
			pdf.push_back(0);
			idx++;
		}

		pdf.push_back(it->second / (double)norm);
		idx++;
	}

	return pdf;
}


std::vector<int> powerLawVectorSampleAndInsert(int n, double alpha, int x1, int x0)
{
	// Sample-and-insert method

	std::vector<int> degrees;
	degrees.reserve(n);

	for (int i = 0; i < n; i++)
		degrees.push_back(randPowerDiscrete(alpha, x0, x1));

	return degrees;
}

std::vector<int> unique(const std::vector<int>& source)
{
	// Takes in sorted vector
	std::vector<int> target;
	int last = -1;
	for (auto it = target.begin(); it < target.end(); it++)
	{
		if (*it != last)
		{
			target.push_back(*it);
			last = *it;
		}
	}
	return target;
}

std::vector<int> randomIntVector(int len, int shift)
{
	std::vector<int> result;
	result.reserve(len);
	for (int i = 0; i < len; i++)
	{
		result.push_back(i + shift);
	}

	// Shuffle
	static std::random_device device;
	static std::mt19937 generator(device());
	std::shuffle(result.begin(), result.end(), generator);

	return result;
}

std::map<int, int> countOccurences(const std::vector<int> &vec)
{
	std::map<int, int> freq;
	for (auto it = vec.begin(); it != vec.end(); it++)
		freq[*it]++;
	//for (auto const & x : vec)
	//	++freq[x];
	return freq;
}

int summarizeMap(std::map<int, int> mymap)
{
	int result = 0;
	for (auto it = mymap.begin(); it != mymap.end(); it++)
		result += it->second;
	return result;
}

std::string currentDateTime()
{
	std::time_t rawtime;
	struct tm timeinfo;
	char buffer[80];

	std::time(&rawtime);
	localtime_s(&timeinfo, &rawtime);

	std::strftime(buffer, 80, "%Y-%m-%d-%H-%M-%S", &timeinfo);
	return std::string(buffer);
}

/*Tests*/

void testRandPowerDiscrete()
{
	int n = 100000;
	std::vector<int> vec(n);
	for (int i = 0; i < n; i++)
	{
		vec[i] = randPowerDiscrete(2, 1, 100);
	}
	saveVector("testRandPowerDiscrete.dat", vec);
}