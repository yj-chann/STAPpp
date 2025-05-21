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

//! Beam31(Euler-Bernoulli) element class
class CB31EB : public CElement
{
public:

//!	Constructor
	CB31EB();

//!	Desconstructor
	~CB31EB();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//!	Calculate Q4 element shape function matrix N
	void ElementShapeFunction(double (&N)[2][8], double xi, double eta);

//!	Calculate derivative of Q4 element shape function matrix B and Jacobian determination
	void ElementStrainFunction(double (&B)[3][8], double* det, double xi, double eta);

//! Calcualte derivative of element shape function matrix DN at coordinate xt
	void ElementDerivativeShapeFunction(double (&DN)[2][4], double xi, double eta);
};
