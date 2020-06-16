#include <chrono>
#include <iostream>
#include "Izhikevich.h"
#include "Connectivity.h"
#include "Network.h"
#include "Criticality.h"
#include "Evolution.h"


using namespace std::chrono;
using std::cout;

int main()
{
	// Start timer
	auto start = high_resolution_clock::now();

	// Run function
	evolution(50, 20);

	// End timer
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << "Time: " << duration.count() << " ms" << std::endl;

	return  0;
}

