#pragma once
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include "Vector.h"
#include "Functions.h"

using std::cout;
using std::endl;

void makeConnectivityMatrix();
int getPostModule(int postNeuron, int neurons, int modules);
int getPreModule(int postModule, int hierarchyDepth);
int getPreModule(int postModule, int hierarchyDepth, int hierarchyLevel);
int getPreNeuron(int preModule, int neurons, int modules);
std::vector<int> getPreModuleVector(int connections, int postModule, int hierarchyDepth, double gamma = 2.0);
Matrix readMatrixFromFile(const char* filename);
void oldMakeConnectivityMatrix();

// Tests
void testGetPreModule();
void testNeuronModuleConversion();
void testPowerLawDiscrete();
void testGetDegrees();
void testCountOccurences();
void testBalance(double balance);
