/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>

#include "Node.h"

CNode::CNode(double X, double Y, double Z)
{
    XYZ[0] = X;		// Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;

	XYZ[3] = 0;		// Rotation Angles of the node
    XYZ[4] = 0;
    XYZ[5] = 0;
    
    bcode[0] = 0;	// Boundary codes
    bcode[1] = 0;
    bcode[2] = 0;

	bcode[3] = -1;	// Boundary codes
    bcode[4] = -1;
    bcode[5] = -1;
};

//	Read element data from stream Input
bool CNode::Read(ifstream& Input)
{
	Input >> NodeNumber;	// node number
	Input >> bcode[0] >> bcode[1] >> bcode[2]
		  >> XYZ[0] >> XYZ[1] >> XYZ[2];

	if (Input.peek() == '\n') {
	} 
	else {
    Input >> bcode[3] >> bcode[4] >> bcode[5];
	NodeType = 1;
	}

	return true;
}

//	Output nodal point data to stream
void CNode::Write(COutputter& output)
{	
	if (NodeType) {
		output << setw(9) << NodeNumber << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) <<  bcode[2]
		<< setw(5) << bcode[3] << setw(5) << bcode[4] << setw(5) <<  bcode[5]
		<< setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
	}
	else {
		output << setw(9) << NodeNumber << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) <<  bcode[2]
		<< setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
	}
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output)
{
	output << setw(9) << NodeNumber << "       ";

	if (NodeType) {
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np, Beam/Plate/Shell node
		{
			output << setw(5) << bcode[dof];
		}
	}
	else {
		for (unsigned int dof = 0; dof < 3; dof++)	// Loop over for DOFs of node np, Solid Element node
		{
			output << setw(5) << bcode[dof];
		}
		for (unsigned int i = 0; i < 3; i++)	// Loop over for DOFs of node np, Solid Element node
		{
			output << setw(5) << "N/A";
		}
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, double* Displacement)
{
	output << setw(5) << NodeNumber << "        ";

	for (unsigned int j = 0; j < 3; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	if (NodeType) {
		for (unsigned int j = 3; j < 6; j++)
		{
			if (bcode[j] == 0)
			{
				output << setw(15) << 0.0;
			}
			else
			{
				output << setw(15) << Displacement[bcode[j] - 1];
			}
		}

	}
	else {
		output << setw(12) << "N/A"<< setw(12) << "N/A"<< setw(15) << "N/A";
	}

	output << endl;
}
