#include "Matrix.h"
#include <ostream>

Matrix::Matrix() : rows(0), columns(0), data(nullptr) {
	// Ingenting her
}

Matrix::Matrix(unsigned int N) : Matrix(N, N, 0.0)
{
	for (unsigned int i = 0; i < N; ++i)
		set(i, i, 1.0); // Set equal to identity
}

Matrix::Matrix(unsigned int rows, unsigned int columns, const double fill)
	: rows(rows)
	, columns(columns)
	, data(nullptr)
{
	data = new double[rows*columns];
	for (unsigned int i = 0; i < rows*columns; ++i)
		data[i] = fill;
}

Matrix::Matrix(const Matrix & rhs)
	: rows(rhs.rows)
	, columns(rhs.columns)
	, data(0)
{
	// If rhs was not valid, we know that rhs.row and columns == 0, therefore we have an invalid matrix and can stop.
	if (!rhs.isValid())
		return;

	data = new double[rows*columns];
	for (unsigned int i = 0; i < rows*columns; ++i)
		data[i] = rhs.data[i];
}

Matrix::~Matrix() {
	invalidate();
}

void Matrix::invalidate() {
	delete[] data;
	data = nullptr;
	rows = 0;
	columns = 0;

}

bool Matrix::isValid() const {
	return !!data; // Forcing pointer to bool
}

Matrix& Matrix::operator=(Matrix rhs)
{
	std::swap(data, rhs.data);
	std::swap(rows, rhs.rows);
	std::swap(columns, rhs.columns);
	return *this;
}

Matrix & Matrix::operator +=(const Matrix & rhs) {
	if (this->rows == rhs.rows && this->columns == rhs.columns)
	{
		for (unsigned int row = 0; row < rows; ++row)
			for (unsigned int column = 0; column < columns; ++column)
				this->at(row, column) += rhs.at(row, column);
	}
	else
		invalidate();

	return *this;
}

Matrix Matrix::operator+(const Matrix & rhs) const {
	if (this->rows == rhs.rows && this->columns == rhs.columns)
		return Matrix(*this) += rhs;
	else
		return Matrix();
}

std::ostream & operator<<(std::ostream & out, const Matrix & elem)
{
	if (!elem.isValid())
		out << "The matrix is not valid.";
	else
	{
		for (unsigned int row = 0; row < elem.getRows(); ++row)
		{
			out << "| ";
			for (unsigned int column = 0; column < elem.getColumns(); ++column)
				out << elem.get(row, column) << ' ';
			out << '|' << std::endl;
		}
	}
	return out;
}

// The code below contains additional functionality, not part of the exercise.

Matrix & Matrix::operator-=(const Matrix & rhs) {
	if (this->rows == rhs.rows && this->columns == rhs.columns)
	{
		for (unsigned int row = 0; row < rows; ++row)
			for (unsigned int column = 0; column < columns; ++column)
				this->at(row, column) -= rhs.at(row, column);
	}
	else
		invalidate();

	return *this;
}

Matrix Matrix::operator-(const Matrix & rhs) const {
	if (this->rows == rhs.rows && this->columns == rhs.columns)
		return Matrix(*this) -= rhs;
	else
		return Matrix();
}

Matrix Matrix::operator *(const Matrix & rhs) const {
	if (this->isValid() && this->columns == rhs.rows)
	{
		Matrix temp(this->rows, rhs.columns);

		for (unsigned int row = 0; row < this->rows; ++row)
			for (unsigned int column = 0; column < rhs.columns; ++column)
			{
				temp.at(row, column) = 0.0;
				for (unsigned int i = 0; i < this->columns; ++i)
					temp.at(row, column) += this->at(row, i) * rhs.at(i, column);
			}
		return temp;
	}
	else
		return Matrix();
}

Matrix & Matrix::operator *=(const Matrix & rhs) {
	return *this = *this * rhs;
}

Matrix Matrix::operator*(double rhs) const {
	return Matrix(*this) *= rhs;
}

Matrix operator*(double lhs, const Matrix & rhs) {
	return rhs * lhs;
}

Matrix & Matrix::operator*=(double rhs) {
	for (unsigned int i = 0; i < this->rows * this->columns; ++i)
		this->data[i] *= rhs;
	return *this;
}

Matrix Matrix::operator-() const {
	Matrix temp(*this);

	for (unsigned int i = 0; i < this->rows * this->columns; ++i)
		temp.data[i] = -temp.data[i];

	return temp;
}

// Additional functionality for neural networks
Matrix::Matrix(int rows, int columns, double lower, double upper)
	: Matrix(rows, columns)
{
	randomize(lower, upper);
}

void Matrix::randomize(double lower, double upper)
{
	//static thread_local std::random_device device;
	static thread_local std::mt19937 generator;
	std::uniform_real_distribution<> uniform(lower, upper);

	for (unsigned int i = 0; i < rows*columns; i++)
		data[i] = uniform(generator);
}

void Matrix::randNormal(double mean, double stddev)
{
	//static thread_local std::random_device device;
	static thread_local std::mt19937 generator;
	std::normal_distribution<> normal(mean, stddev); // This can't be declared static

	for (unsigned int i = 0; i < rows*columns; i++)
		data[i] = normal(generator);
}

void Matrix::randNormal(std::mt19937 & generator, std::normal_distribution<> & dist)
{
	for (unsigned int i = 0; i < rows*columns; i++)
		data[i] = dist(generator);
}


Matrix::Matrix(const Matrix & top, const Matrix & bottom)
	: Matrix(top.rows + bottom.rows, top.columns)
{
	if (top.columns != bottom.columns)
		this->invalidate();

	for (unsigned int i = 0; i < top.rows + bottom.rows; i++)
	{
		for (unsigned int j = 0; j < this->columns; j++)
		{
			if (i < top.rows)
				this->at(i, j) = top.at(i, j);

			else
				this->at(i, j) = bottom.at(i - top.rows, j);
		}
	}
}

void Matrix::populate()
{
	std::vector<int> numbers = std::vector<int>(rows*columns);
	for (unsigned int i = 0; i < rows*columns; i++)
		numbers[i] = i;

	//std::random_device rd;
	static thread_local std::mt19937 generator;
	std::shuffle(numbers.begin(), numbers.end(), generator);

	for (unsigned int i = 0; i < rows*columns; i++)
	{
		data[i] = numbers.back();
		numbers.pop_back();
	}
}

Matrix & Matrix::operator+=(double rhs)
{
	for (unsigned int i = 0; i < rows*columns; i++)
		data[i] += rhs;

	return *this;
}

Matrix & Matrix::operator-=(double rhs)
{
	*this += (-rhs);
	return *this;
}

Matrix Matrix::operator+(double rhs) const
{
	return Matrix(*this) += rhs;
}

Matrix Matrix::operator-(double rhs) const
{
	return Matrix(*this) -= rhs;
}

Matrix & Matrix::extend(const Matrix & rhs)
{
	if (this->columns != rhs.columns)
	{
		this->invalidate();
		return *this;
	}

	Matrix temp(this->rows + rhs.rows, this->columns);
	for (unsigned int i = 0; i < this->rows + rhs.rows; i++)
	{
		for (unsigned int j = 0; j < this->columns; j++)
		{
			if (i < this->rows)
				temp.at(i, j) = this->at(i, j);
			else
				temp.at(i, j) = rhs.at(i - this->rows, j);
		}
	}
	return *this = temp;
}

Matrix & Matrix::expand(const Matrix & rhs)
{
	if (this->rows != rhs.rows)
	{
		this->invalidate();
		return *this;
	}

	Matrix temp(this->rows, this->columns + rhs.columns);
	for (unsigned int i = 0; i < this->rows; i++)
	{
		for (unsigned int j = 0; j < this->columns + rhs.columns; j++)
		{
			if (j < this->columns)
				temp.at(i, j) = this->at(i, j);

			else
				temp.at(i, j) = rhs.at(i, j - this->columns);
		}
	}
	return *this = temp;
}

Matrix Matrix::operator^(double rhs) const
{
	Matrix temp(rows, columns);
	for (unsigned int i = 0; i < rows*columns; i++)
		temp.data[i] = pow(this->data[i], rhs);

	return temp;
}

Matrix operator+(double lhs, const Matrix & rhs)
{
	return rhs + lhs;
}

Matrix operator-(double lhs, const Matrix & rhs)
{
	return (-rhs) + lhs;
}

Matrix Matrix::operator>=(double rhs) const
{
	Matrix result = Matrix(rows, columns, 0.0);
	for (unsigned int i = 0; i < rows*columns; i++)
	{
		if (data[i] >= rhs)
			result.data[i] = 1;
	}
	return result;
}

Matrix Matrix::hadamard(const Matrix & rhs) const
{
	if (this->rows != rhs.rows)
		return Matrix();

	Matrix temp(rows, columns);
	for (unsigned int i = 0; i < this->rows*this->columns; i++)
		temp.data[i] = this->data[i] * rhs.data[i];

	return temp;
}

void Matrix::save(const char* filename) const
{
	std::ofstream myfile;
	myfile.open(filename, std::ios::trunc);

	for (unsigned i = 0; i < this->rows*this->columns; i++)
	{
		myfile << data[i] << " ";
		if (!((i + 1) % this->columns))
			myfile << "\n";
	}

	// These methods use similar time

	/*
	for (unsigned int row = 0; row < this->getRows(); row++)
	{
		for (unsigned int column = 0; column < this->getColumns(); column++)
		{
			myfile << at(row, column) << " ";
		}
		myfile << std::endl;
	}
	*/

	myfile.close();
}
