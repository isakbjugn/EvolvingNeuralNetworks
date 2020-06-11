#pragma once
#include "Vector.h"
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

void izhikevich();
void updateFirings(Matrix & firings, const Vector & fired, int time);
void updateFirings(std::vector<int> & fireTiming, std::vector<int>& fireNeuron, const std::vector<int> fired, int time);
Vector synapticCurrent(const Matrix & S, const std::vector<int> & fired);
void writeToFile(const char* filename, std::vector<int> & fireTiming, std::vector<int>& fireNeuron);
void writeNeuronToFile(const char* filename, std::vector<double> & neuron);
//void printState(Vector v, Vector u, std::vector<int> fired);
void printVector(const std::vector<int> vec);

// Tests
void testOperations();
void testSynapticCurrent();