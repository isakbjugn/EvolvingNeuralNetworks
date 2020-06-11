#pragma once
#include <iomanip>
#include <unordered_set>
#include <utility>
#include "Functions.h"
using std::cout;

// Note, for all functions:
// x0 and x1 are values, not indeces (cf. 0-indexing)

std::vector<int> avalancheDetection(const std::vector<int> &times, const std::vector<int> &neurons);
std::vector<double> readVector(const char* filename);
double calculateKappa(const std::vector<double> &pdf_theory, const std::vector<double> &pdf_emp, int x0, int x1);
double calculateKS(const std::vector<double> &pdf_theory, const std::vector<double> &pdf_emp, int x0, int x1);
double calculateKS(const std::vector<double> &pdf_theory, const std::vector<double> &pdf_emp);
std::vector<int> logSpacing(int x0, int x1, int m = 10);
std::vector<double> makePowerLawDistribution(double alpha, int length);
double mle(const std::vector<double> &pdf, int x0 = 1);
double mle(const std::vector<int> &sizes, int x0 = 1);
double mleSearch(const std::vector<int> &sizes, int x0 = 1);
double mleSearch(const std::vector<double> &sizes, int x0 = 1);
double maxLikelihood(double g, int x0, int x1);
std::pair<double, double> searchAlpha(const std::vector<double> &pdf_emp, int x0 = 1);
std::pair<double, double> searchAlpha(const std::vector<int> &sizes, int x0 = 1);
std::pair<double, double> clauset(const std::vector<int> &sizes);

/*Tests*/
void testReadVector();
void testKappa();
void testKS();
void testSearchAlpha();
void testAvalancheDetection();
void testMle();