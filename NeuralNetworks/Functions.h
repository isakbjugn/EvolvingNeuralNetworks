#pragma once
#include <random>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <iterator>
#include <ctime>
#include <string>

int randomIntWithLimits(int lower, int upper);
double randomRealWithLimits(double lower, double upper);
double randUniform();
double randNormal(double mean, double stddev);
int randPowerDiscrete(double alpha, int xMin, int xMax);
double randPower(double alpha, double xMin, double xMax);
std::vector<int> powerLawVector(int n, double alpha, int x1, int x0 = 1);
std::vector<int> histogramToSamples(const std::vector<double>& pdf);
std::vector<double> samplesToPdf(const std::vector<int>& samples);
std::vector<int> powerLawVectorSampleAndInsert(int n, double alpha, int x1, int x0 = 1);
std::vector<int> unique(const std::vector<int> &source);
std::vector<int> randomIntVector(int len, int shift = 0);
std::map<int, int> countOccurences(const std::vector<int> &vec);
int summarizeMap(std::map<int, int> mymap);
std::string currentDateTime();

/*Template functions*/
template<typename T>
void saveVector(const char* filename, const std::vector<T>& vec);
template<typename T>
double partSum(const std::vector<T> &vec, int k0, int k1); // [k0, k1]
template<typename T>
T sumVector(const std::vector<T>& vec);

/*Implementation of template functions*/

template<typename T>
void saveVector(const char* filename, const std::vector<T> &vec)
{
	std::ofstream myfile;
	myfile.open(filename, std::ios::trunc);

	for (auto it = vec.begin(); it != vec.end(); it++)
	{
		myfile << *it << " ";
	}
	myfile << "\n";
	myfile.close();
}
template<typename T>
T sumVector(const std::vector<T>& vec)
{
	return partSum(vec, 0, vec.size() - 1);
}
template<typename T>
double partSum(const std::vector<T> &vec, int k0, int k1)
{
	// k1 is an index, not a length
	// From-and-including index k0 to-and-including k1

	T partialSum = 0;
	for (int i = k0; i < k1 + 1; i++)
		partialSum += vec[i];
	return partialSum;
}


/*Tests*/
void testRandPowerDiscrete();
