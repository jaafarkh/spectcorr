#pragma once

#include <iostream>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
class polar {
	friend polar operator +(polar const& c1, polar const& c2);
	friend polar operator -(polar const& c1, polar const& c2);
	friend polar operator *(polar const& c1, polar const& c2);
	friend polar operator /(polar const& c1, polar const& c2);
public:
	polar(double r = 0.0, double ang = 0.0) :amp(r), angle(ang) {}
	double amp;
	double angle;
};
//define our own complex class to improve speed, out-of-the-box c++ complex class is slow
class complex {
	friend complex operator +(complex const& c1, complex const& c2);
	friend complex operator -(complex const& c1, complex const& c2);
	friend complex operator *(complex const& c1, complex const& c2);
	friend complex operator *(float const& c1, complex const& c2);
	friend complex operator *(double const& c1, complex const& c2);
	friend complex operator /(complex const& c1, complex const& c2);
public:
	complex(double x = 0.0, double y = 0.0) :real(x), imag(y) {}
	polar toPolar() const { return polar(sqrt(imag*imag + real*real), atan2(imag, real)); }
	double real;
	double imag;
	double amp() const { return sqrt(imag*imag + real*real); }
	double angle() const { return atan2(imag, real); }
	complex conj() const { return complex(real, -imag); }
};