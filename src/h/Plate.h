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

//! Mindlin-Reissner Plate element class
class CPlate : public CElement
{
public:

//!	Constructor
	CPlate();

//!	Desconstructor
	~CPlate();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//! Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
    virtual void GenerateLocationMatrix() override;

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//!	Calculate element non-homogeneous essential boundary conditions
	virtual void ElementNonHomo(double* Matrix, double* NonForce) override;

//!	Calculate Mindlin-Reissner Plate element shape function matrix N
	void ElementShapeFunction(double (&N)[1][12], double xi, double eta);

//!	Calculate Mindlin-Reissner Plate element bending kinematic matrix Bb, shear kinematic matrix Bs and Jacobian determination
	void ElementStrainFunction(double (&B)[3][12], double* det, double xi, double eta);
};
