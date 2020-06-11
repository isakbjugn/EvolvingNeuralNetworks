#pragma once
#include "Matrix.h"
#include <vector>

class Vector : public Matrix {
public:
	Vector()
		: Matrix() { 	}

	explicit Vector(unsigned int N, double fill = 0.0)
		: Matrix(N, 1, fill)
	{ }

	Vector(const Matrix & other); // To construct a Vector from a matrix, as for instance in the result matrix of a M*V operation.

	double & at(int index) {
		return Matrix::at(index, 0);
	}

	const double & at(int index) const {
		return Matrix::at(index, 0);
	}

	double get(int index) const {
		return at(index);
	}

	void set(int index, double value) {
		at(index) = value;
	}

	double dot(const Vector & rhs) const;
	double lengthSquared() const;
	double length() const;

	// Additional functionality for neural networks
	Vector(int rows, double lower, double upper);
	Vector(const Vector& top, const Vector& bottom);

	double & operator()(int row)
	{
		return Matrix::at(row, 0);
	}
	const double & operator()(int row) const
	{
		return Matrix::at(row, 0);
	}
	Vector & expand(Vector & rhs) = delete;
	void update(const Vector& source, const std::vector<int>& indeces);
	std::vector<int> find(double threshold) const;
	void copy(const Vector& top, const Vector& bottom);
};

/*Tests*/
void testVectorEfficiency();