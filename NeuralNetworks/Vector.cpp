#include "Vector.h"

Vector::Vector(const Matrix & other)
	: Matrix() // Start with this as invalid Vector.
{
	// If and only if the resultant matrix is valid, assign it.
	if (other.isValid() && other.getColumns() == 1)
		Matrix::operator=(other); // reuse matrix' operator=
}

double Vector::dot(const Vector & rhs) const {
	if (!this->isValid() || this->getRows() != rhs.getRows())
		return nan("");

	double sum = 0.0;
	for (unsigned int i = 0; i < this->getRows(); ++i)
		sum += this->get(i)*rhs.get(i);

	return sum;
}

double Vector::lengthSquared() const {
	return this->dot(*this);
}

double Vector::length() const {
	return sqrt(lengthSquared());
}

// Additional functionality for neural networks
Vector::Vector(int rows, double lower, double upper)
	: Matrix(rows, 1, lower, upper) {};

void Vector::update(const Vector & source, const std::vector<int> & indeces)
{
	if (this->getRows() != source.getRows())
	{
		return;
	}

	for (unsigned int i = 0; i < indeces.size(); i++)
	{
		set(indeces[i], source(indeces[i]));
	}
}

Vector::Vector(const Vector& top, const Vector& bottom)
	: Matrix(top.getRows() + bottom.getRows(), 1)
{
	copy(top, bottom);
}

std::vector<int> Vector::find(double threshold) const
{
	std::vector<int> trues = std::vector<int>();

	for (unsigned int i = 0; i < getRows(); i++)
	{
		if (get(i) >= threshold)
		{
			trues.push_back(i);
		}
	}
	return trues;
}

void Vector::copy(const Vector& top, const Vector& bottom)
{
	if (top.getRows() + bottom.getRows() != this->getRows())
	{
		invalidate();
	}

	for (unsigned int i = 0; i < this->getRows(); i++)
	{
		if (i < top.getRows())
		{
			set(i, top.get(i));
		}
		else
			set(i, bottom.get(i - top.getRows()));
	}
}


/*Tests*/

void testVectorEfficiency()
{
	int n = 4000000;
	
	// Fastest
	Vector a(n);
	for (int i = 0; i < n; i++)
	{
		a(i) = i;
	}

	// Second-fastest
	std::vector<int> b(n,0);
	for (int i = 0; i < n; i++)
	{
		b[i];
	}

	// Second-slowest
	std::vector<int> c;
	c.reserve(n);
	for (int i = 0; i < n; i++)
	{
		c.push_back(i);
	}

	// Slowest
	std::vector<int> d;
	for (int i = 0; i < n; i++)
	{
		d.push_back(i);
	}
}