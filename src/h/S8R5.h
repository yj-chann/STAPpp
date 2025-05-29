/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! 8-node quadratic degenerated shell element class
class CS8R5 : public CElement
{
public:

//!	Constructor
	CS8R5();

//!	Desconstructor
	~CS8R5();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress, for S8R5 element
	virtual void ElementStress(double* stress, double* Displacement);

//!	Calculate element non-homogeneous essential boundary conditions
	virtual void ElementNonHomo(double* Matrix, double* NonForce) override;

//!	Calculate element shape function matrix N at parent coordinate xi, eta, zeta
	void ElementShapeFunction(double (&N)[3][40], double xi, double eta, double zeta);

//!	Calculate derivative of element shape function matrix B(global frame) at parent coordinate xi, eta, zeta
	void ElementStrainFunction(double (&B)[5][40], double* det, double xi, double eta, double zeta);

//!	Calcualte detJ at parent coordinate xi, eta, zeta = 1
	void ElementDetTop(double* det, double xi, double eta);
};
