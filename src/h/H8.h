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

//! H8 element class
class CH8 : public CElement
{
public:

//!	Constructor
	CH8();

//!	Desconstructor
	~CH8();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//!	Calculate element non-homogeneous essential boundary conditions
	virtual void ElementNonHomo(double* Matrix, double* NonForce) override;

//!	Calculate H8 element shape function matrix N
	void ElementShapeFunction(double(&N)[3][24], double xi, double eta , double zeta);

//!	Calculate derivative of H8 element shape function matrix B and Jacobian determination
	void ElementStrainFunction(double(&B)[6][24], double* det, double xi, double eta, double zeta);
};

