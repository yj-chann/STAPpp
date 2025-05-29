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

#include "Outputter.h"

using namespace std;

//!	Node class
class CNode
{
public:

//!	Maximum number of degrees of freedom per node
/*!	For 3D bar and solid elements, NDF = 3. For 3D beam or shell elements, NDF = 5 or 6 */
	const static unsigned int NDF = 6;

//!	Node number
	unsigned int NodeNumber;

//!	x, y and z coordinates of the node
	double XYZ[6];

//!	Boundary code of each degree of freedom of the node
/*!		0: The corresponding degree of freedom is active (defined in the global system) */
/*!		1: The corresponding degree of freedom in nonactive (not defined) */
/*!	    2: The corresponding degree of freedom has Non-homogeneous essential boundary conditions */
/*!	    3: The corresponding degree of freedom is rotation freedom and is useless in solid element like H8 */
/*!	After call Domain::CalculateEquationNumber(), bcode stores the global equation number */
/*!	corresponding to each degree of freedom of the node */
	unsigned int bcode[NDF];

//! Node type, 0: Solid Element Node; 1: Beam\Plate\Shell Node
	unsigned int NodeType = 0;

//! Calculate non-homogeneous essential boundary conditions flag
/*!		0: Closed */
/*!		1: Open */
	unsigned int NonHomo = 0;

//! Non-homogeneous essential boundary conditions
	double BC[6] = {0};

//!	Constructor
	CNode(double X = 0, double Y = 0, double Z = 0);

//!	Read nodal point data from stream Input
	bool Read(ifstream& Input);

//!	Output nodal point data to stream
	void Write(COutputter& output);

//!	Output equation numbers of nodal point to stream OutputFile
	void WriteEquationNo(COutputter& OutputFile);

//!	Write nodal displacement
	void WriteNodalDisplacement(COutputter& OutputFile, double* Displacement);
};
