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

//! Beam21(Euler-Bernoulli) element class
class CB21EB : public CElement
{
public:

//!	Constructor
	CB21EB();

//!	Desconstructor
	~CB21EB();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress, for Beam element, they are moment and shear force in Gausspoints
	virtual void ElementStress(double* stress, double* Displacement);

//!	Calculate element shape function matrix N at parent coordinate xi
	void ElementShapeFunction(double* N, double xi);

//!	Calculate derivative of element shape function matrix B at parent coordinate xi(m = 2, actually second derivative)
	void ElementStrainFunction(double* B, double xi);

//!	Calculate second derivative of element shape function matrix S at parent coordinate xi
	void ElementDerivativeStrainFunction(double *S, double xi);
};
