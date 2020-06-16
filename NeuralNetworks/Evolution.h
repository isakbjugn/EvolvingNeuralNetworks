#pragma once
#include <queue>
#include <thread>
#include "Network.h"

using std::shared_ptr;
using std::make_shared;
using std::vector;
using std::move;
//using namespace std::chrono;
using std::string;

int evolution(int population, int generations);

vector<shared_ptr<Network>> selectDescendants(vector<shared_ptr<Network>> &generation, double gap);

vector<shared_ptr<Network>> selectParents(vector<shared_ptr<Network>> &generation);

shared_ptr<Network> findFittest(vector<shared_ptr<Network>> &generation);

//void recordRuntime(const vector<shared_ptr<Network>> &generation);

void testNetworkConstructor();

void testSpeed();
