/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "B21EB.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CB21EB::CB21EB()
{
	NEN_ = 2;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 6;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CB21EB::~CB21EB()
{
}

//	Read element data from stream Input
bool CB21EB::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// first node number, second node number

	Input >> N1 >> N2 >> MSet;
    ElementMaterial_ = dynamic_cast<CB21EBMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CB21EB::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CB21EB::ElementStiffness(double* Matrix)  // Local stiffness matrix can be derived analytically, then use DCM matrix
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate element stiffness matrix

	CB21EBMaterial* material_ = dynamic_cast<CB21EBMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E;
    double Area = material_->Area;
	double I = material_->I;

//	Calculate beam length
	double DX[2];		//	dx = x2-x1, dy = y2-y1
	for (unsigned int i = 0; i < 2; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[3];	//  Quadratic polynomial (dx^2, dy^2, dx*dy)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[0] * DX[1];

	double L2 = DX2[0] + DX2[1];
	double L = sqrt(L2);

    double Ke[6][6] = {0};

    Ke[0][0] = E * Area / L;
    Ke[0][3] = -E * Area / L;

    Ke[1][1] = 12.0 * E * I / (L * L2);
    Ke[1][2] = 6.0 * E * I / L2;
    Ke[1][4] = -12.0 * E * I / (L * L2);
    Ke[1][5] = 6.0 * E * I / L2;

    Ke[2][1] = 6.0 * E * I / L2;
    Ke[2][2] = 4.0 * E * I / L;
    Ke[2][4] = -6.0 * E * I / L2;
    Ke[2][5] = 2.0 * E * I / L;

    Ke[3][0] = -E * Area / L;
    Ke[3][3] = E * Area / L;

    Ke[4][1] = -12.0 * E * I / (L * L2);
    Ke[4][2] = 6.0 * E * I / L2;
    Ke[4][4] = 12.0 * E * I / (L * L2);
    Ke[4][5] = -6.0 * E * I / L2;

    Ke[5][1] = 6.0 * E * I / L2;
    Ke[5][2] = 2.0 * E * I / L;
    Ke[5][4] = -6.0 * E * I / L2;
    Ke[5][5] = 4.0 * E * I / L;

// Calculate direction cosine matrix
    double cos = DX[0] / L; // Cosine of angle between beam and x-axis
    double sin = DX[1] / L;

// Calculate transformation Matrix
    double Te[6][6] = {0};
    Te[0][0] = cos;
    Te[0][1] = sin;
    Te[1][0] = -sin;
    Te[1][1] = cos;
    Te[2][2] = 1.0;
    Te[3][3] = cos;
    Te[3][1] = sin;
    Te[4][3] = -sin;
    Te[4][4] = cos;
    Te[5][5] = 1.0;

    double K[6][6] = {0};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            for (int k = 0; k < 6; k++)
                for (int l = 0; l < 6; l++)
                    K[i][j] += Te[k][i] * Ke[k][l] * Te[l][j];

	int index = 0;
    for (int j = 0; j < 6; j++)	// Column-major order
	{         
        for (int i = j; i >= 0; i--) {    // Upper triangle (i ≤ j)
            Matrix[index++] = K[i][j];   // Pack into one dimensional element stiffness matrix
        }
    }
}

//	Calculate element stress, for CST element, the strain and stress is constant 
void CB21EB::ElementStress(double* stress, double* Displacement)
{	    
	CB21EBMaterial* material_ = dynamic_cast<CB21EBMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E;
    double Area = material_->Area;
	double I = material_->I;

// x and y coordinates, and rotation angle
    double C[6];
	for (unsigned int i = 0; i < 2; i++)
	{
		C[3*i] = nodes_[i]->XYZ[0];
        C[3*i+1] = nodes_[i]->XYZ[1];
        C[3*i+2] = nodes_[i]->XYZ[2];
    };   

//  Get element displacement vector
	double d[6] = {0};
	for (unsigned int i = 0; i < 6; i++)
		if (LocationMatrix_[i])
			d[i] = Displacement[LocationMatrix_[i]-1];

// Calculate element stress , sigma = D * B * d
    double GaussPoints[2] = {
            -sqrt(3) / 3.0,
            sqrt(3) / 3.0,
    };
    int index = 0;
    for (unsigned int gp = 0; gp < 2; gp++)
    {
        double B[4] = {0};
        double S[4] = {0};
        ElementStrainFunction(B, GaussPoints[gp]); 
        ElementDerivativeStrainFunction(S, GaussPoints[gp]);
        stress[4*index] += (C[0] + C[3]) / 2.0 + (C[3] - C[0]) / 2.0 * GaussPoints[gp];
        stress[4*index+1] += (C[1] + C[4]) / 2.0 + (C[4] - C[1]) / 2.0 * GaussPoints[gp];
        for (int i = 0; i < 4; i++)
            stress[4*index+2] += E * I * B[i] * d[i];
        for (int i = 0; i < 4; i++)
            stress[4*index+3] += -E * I * S[i] * d[i];
    }
}

//	Calculate element shape function matrix N at parent coordinate xi
void CB21EB::ElementShapeFunction(double* N, double xi)
{
    //	Calculate beam length
	double DX[2];		//	dx = x2-x1, dy = y2-y1
	for (unsigned int i = 0; i < 2; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[3];	//  Quadratic polynomial (dx^2, dy^2, dx*dy)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[0] * DX[1];

	double L2 = DX2[0] + DX2[1];
	double L = sqrt(L2);
    
    double xi3 = xi * xi * xi;
    double xi2 = xi * xi;
    N[0] = (xi3 - 3.0 * xi + 2.0) / 4.0;
    N[1] = L * (xi3 - xi2 - xi + 1.0) / 8.0;
    N[2] = (-xi3 + 3.0 * xi + 2.0) / 4.0;
    N[3] = L * (xi3 + xi2 - xi - 1.0) / 8.0;
}

//	Calculate derivative of element shape function matrix B at parent coordinate xi(m = 2, actually second derivative)
void CB21EB::ElementStrainFunction(double* B, double xi)
{   
    //	Calculate beam length
	double DX[2];		//	dx = x2-x1, dy = y2-y1
	for (unsigned int i = 0; i < 2; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[3];	//  Quadratic polynomial (dx^2, dy^2, dx*dy)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[0] * DX[1];

	double L2 = DX2[0] + DX2[1];
	double L = sqrt(L2);

    B[0] = 6.0 * xi / L2;
    B[1] = (3.0 * xi - 1) / L;
    B[2] = -6.0 * xi / L2;
    B[3] = (3.0 * xi + 1) / L;
}

//	Calculate second derivative of element shape function matrix S at parent coordinate xi
void CB21EB::ElementDerivativeStrainFunction(double *S, double xi)
{
    //	Calculate beam length
	double DX[2];		//	dx = x2-x1, dy = y2-y1
	for (unsigned int i = 0; i < 2; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[3];	//  Quadratic polynomial (dx^2, dy^2, dx*dy)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[0] * DX[1];

	double L2 = DX2[0] + DX2[1];
	double L = sqrt(L2);

    S[0] = 12.0 / (L2 * L);
    S[1] = 6 / L2;
    S[2] = -12.0 / (L2 * L);
    S[3] = 6 / L2;
}
