/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "CST.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CCST::CCST()
{
	NEN_ = 3;	// Each element has 3 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 6;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CCST::~CCST()
{
}

//	Read element data from stream Input
bool CCST::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3;	// first node number, second node number and third node number, Counterclockwise

	Input >> N1 >> N2 >> N3 >> MSet;
    ElementMaterial_ = dynamic_cast<CCSTMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];

	return true;
}

//	Write element data to stream
void CCST::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
           << setw(11) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CCST::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate CST element area
	double X[3];
    double Y[3];		//	x and y coordinates
	for (unsigned int i = 0; i < 3; i++)
	{
		X[i] = nodes_[i]->XYZ[0];
        Y[i] = nodes_[i]->XYZ[1];
    }    

    double A = 0.5*(X[0]*Y[1] - X[1]*Y[0] + X[2]*Y[0] - X[0]*Y[2] + X[1]*Y[2] - X[2]*Y[1]);

	double DX[3];	//	x3-x2, x1-x3, x2-x1
    double DY[3];	//	y2-y3, y3-y1, y1-y2
	for (unsigned int i = 0; i < 3; i++)
	{
        unsigned int j = (i + 1) % 3; 
        unsigned int k = (i + 2) % 3; 
    
        DX[i] = nodes_[k]->XYZ[0] - nodes_[j]->XYZ[0]; 
        DY[i] = nodes_[j]->XYZ[1] - nodes_[k]->XYZ[1]; 
    }

//	Calculate element strain matrix
	double B[3][6] = {0};
	B[0][0] = DY[0]/(2.0*A);
	B[0][2] = DY[1]/(2.0*A);
	B[0][4] = DY[2]/(2.0*A);
	B[1][1] = DX[0]/(2.0*A);
	B[1][3] = DX[1]/(2.0*A);
	B[1][5] = DX[2]/(2.0*A);
	B[2][0] = DX[0]/(2.0*A);
	B[2][1] = DY[0]/(2.0*A);
	B[2][2] = DX[1]/(2.0*A);
	B[2][3] = DY[1]/(2.0*A);
	B[2][4] = DX[2]/(2.0*A);
	B[2][5] = DY[2]/(2.0*A);

//	Calculate element stiffness matrix

	CCSTMaterial* material_ = dynamic_cast<CCSTMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E;
    double nu = material_->nu;
	double t = material_->t;
	int plane_strain = material_->plane_strain;	// Check if plane strain is active

	if (plane_strain)
	{
		E = E / (1.0 - nu * nu);
		nu = nu / (1.0 - nu);
	}

    double k = E / (1.0 - nu * nu);
	
	double D[3][3] = {0};	// Build the constitutive matrix D (3x3) for plane stress    
	D[0][0] = k;
    D[0][1] = k * nu;
    D[1][0] = k * nu;
    D[1][1] = k;
    D[2][2] = k * (1.0 - nu) / 2.0;

	double Ke[6][6] = {0};
    for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			for (int k = 0; k < 3; k++)
				for (int l = 0; l < 3; l++)
					Ke[i][j] += A * t * B[k][i] * D[k][l] * B[l][j];

	int index = 0;
    for (int j = 0; j < 6; j++)	// Column-major order
	{         
        for (int i = j; i >= 0; i--) {    // Upper triangle (i ≤ j)
            Matrix[index++] = Ke[i][j];   // Pack into one dimensional element stiffness matrix
        }
    }
}

//	Calculate element stress, for CST element, the strain and stress is constant 
void CCST::ElementStress(double* stress, double* Displacement)
{	
//	Calculate CST element area
	double X[3];
    double Y[3];		//	x and y coordinates
	for (unsigned int i = 0; i < 3; i++)
	{
		X[i] = nodes_[i]->XYZ[0];
        Y[i] = nodes_[i]->XYZ[1];
    }    

    double A = 0.5*(X[0]*Y[1] - X[1]*Y[0] + X[2]*Y[0] - X[0]*Y[2] + X[1]*Y[2] - X[2]*Y[1]);

	double DX[3];	//	x3-x2, x1-x3, x2-x1
    double DY[3];	//	y2-y3, y3-y1, y1-y2
	for (unsigned int i = 0; i < 3; i++)
	{
        unsigned int j = (i + 1) % 3; 
        unsigned int k = (i + 2) % 3; 
        
        DX[i] = nodes_[k]->XYZ[0] - nodes_[j]->XYZ[0]; 
        DY[i] = nodes_[j]->XYZ[1] - nodes_[k]->XYZ[1]; 
    }

//	Calculate element strain matrix
	double B[3][6] = {0};
	B[0][0] = DY[0]/(2.0*A);
	B[0][2] = DY[1]/(2.0*A);
	B[0][4] = DY[2]/(2.0*A);
	B[1][1] = DX[0]/(2.0*A);
	B[1][3] = DX[1]/(2.0*A);
	B[1][5] = DX[2]/(2.0*A);
	B[2][0] = DX[0]/(2.0*A);
	B[2][1] = DY[0]/(2.0*A);
	B[2][2] = DX[1]/(2.0*A);
	B[2][3] = DY[1]/(2.0*A);
	B[2][4] = DX[2]/(2.0*A);
	B[2][5] = DY[2]/(2.0*A);

	CCSTMaterial* material_ = dynamic_cast<CCSTMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E;
    double nu = material_->nu;
	int plane_strain = material_->plane_strain;	// Check if plane strain is active

	if (plane_strain)
	{
		E = E / (1.0 - nu * nu);
		nu = nu / (1.0 - nu);
	}

    double k = E / (1.0 - nu * nu);
	
	double D[3][3] = {0};	// Build the constitutive matrix D (3x3) for plane stress    
	D[0][0] = k;
    D[0][1] = k * nu;
    D[1][0] = k * nu;
    D[1][1] = k;
    D[2][2] = k * (1.0 - nu) / 2.0;

//  Get element displacement vector
	double d[6] = {0};
	for (unsigned int i = 0; i < 6; i++)
		if (LocationMatrix_[i])
			d[i] = Displacement[LocationMatrix_[i]-1];

	unsigned int Node_NDF = ND_ / NEN_;
// Resolving non-homogeneous essential boundary conditions
    if (NonHomo_)
        for (unsigned int i = 0; i < ND_; i++)
            if (!LocationMatrix_[i] && nodes_[i/Node_NDF]->BC[i%Node_NDF] != 0)
                d[i] = nodes_[i/Node_NDF]->BC[i%Node_NDF];


// Calculate element stress (constant), sigma = D * B * d
    double Bd[3] = {0};
    for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 6; j++)
            Bd[i] += B[i][j] * d[j];

	*stress = *(stress + 1) = *(stress + 2) = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            stress[i] += D[i][j] * Bd[j];
}

//	Calculate element non-homogeneous essential boundary conditions
void CCST::ElementNonHomo(double* Matrix, double* NonForce)
{   
    double Ke[6][6] = {0};
    unsigned int index = 0;
    for (int j = 0; j < 6; j++) 
        for (int i = j; i >= 0; i--) 
            Ke[i][j] = Matrix[index++];

    for (unsigned int i = 0; i < 6; i++) 
        for (unsigned int j = 0; j < i; j++) 
            Ke[i][j] = Ke[j][i];        

    double d[6] = {0};
    for (unsigned int i = 0; i < 3; i++)
    {
        d[2*i] = nodes_[i]->BC[0];
        d[2*i+1] = nodes_[i]->BC[1];
    }

    for (unsigned int i = 0; i < 6; i++)
    {   
        if (LocationMatrix_[i] == 0)
            continue;
        for (unsigned int j = 0; j < 6; j++)
        {
            if (LocationMatrix_[j] != 0)
                continue;
            NonForce[i] += Ke[i][j] * d[j];
        }
    }
}