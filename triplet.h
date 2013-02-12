// ----------------------------------------------------------------//
// Filename : triplet.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// header file of Triplet class, this is used for solving sparse 
// matrix
// ----------------------------------------------------------------//
// - Zigang Xiao - Fri Oct 29 16:44:11 CDT 2010
//   * created file
// - Zigang Xiao - Tue Dec  7 21:57:13 CST 2010
//   * modified to adapt complex number

#ifndef __TRIPLET_H__
#define __TRIPLET_H__

#include <vector>
#include <iostream>
#include "vec.h"
#include "algebra.h"
using namespace std;

class Vec;
class Triplet{
public:
	Triplet();
	~Triplet();
	void merge();
	void clear();
	void push_back(long i,long j,double x);
	size_t size() const;
	long get_row() const;
	void set_row(long row);
	void to_arrays(long * Ti, long * Tj, double * Tx) const;
	void diagonal_split(Triplet & L, Triplet & D, Triplet & U) const;

	friend ostream & operator <<(ostream & os, const Triplet & t);

	// vector matrix multiplication
	friend Vec operator *(const Vec & x, const Triplet & A);
	friend Vec operator *(const Triplet & A, const Vec & x);

	// scale
	Triplet & operator *= (double scale);
	friend Triplet operator *(const Triplet & A, double scale);
	friend Triplet operator *(double scale, const Triplet & A);

	// addition
	Triplet & operator += (const Triplet & B);
	Triplet operator + (const Triplet & B) const;

	friend class Algebra;
	friend class Circuit;
private:
	vector<long> Ti;
	vector<long> Tj;
	vector<double> Tx;
	long row;
};

inline long Triplet::get_row() const{
	return this->row ;
}
inline void Triplet::set_row(long row){
	this->row = row;
}

#endif
