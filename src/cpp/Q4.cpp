/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Q4.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CQ4::CQ4()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 8;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CQ4::~CQ4()
{
}

//	Read element data from stream Input
bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// first node number, second node number, third node number and fourth node number, Counterclockwise

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void CQ4::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
           << setw(11) << nodes_[1]->NodeNumber
           << setw(11) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CQ4::ElementStiffness(double* Matrix)  // use full integration here, then is reduced integration with hourglass control
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate element stiffness matrix

	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);	// Pointer to material of the element

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

    double GaussPoints[2] = {
        -sqrt(3) / 3.0,
        sqrt(3) / 3.0,
    };
    double Ke[8][8] = {0};
    for (unsigned int gp1 = 0; gp1 < 2; gp1++)
    {
        for (unsigned int gp2 = 0; gp2 < 2; gp2++)
        {
            double B[3][8] = {0};
            double det;
            ElementStrainFunction(B, &det, GaussPoints[gp1], GaussPoints[gp2]);
            for (int i = 0; i < 8; i++)
                for (int j = 0; j < 8; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            Ke[i][j] += t * B[k][i] * D[k][l] * B[l][j] * det;  // Weight is 1
        }
    }

	int index = 0;
    for (int j = 0; j < 8; j++)	// Column-major order
	{         
        for (int i = j; i >= 0; i--) {    // Upper triangle (i ≤ j)
            Matrix[index++] = Ke[i][j];   // Pack into one dimensional element stiffness matrix
        }
    }
}

//	Calculate element stress, for CST element, the strain and stress is constant 
void CQ4::ElementStress(double* stress, double* Displacement)
{	    
	CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial_);	// Pointer to material of the element

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

//  Get element displacement vector
	double d[8] = {0};
	for (unsigned int i = 0; i < 8; i++)
		if (LocationMatrix_[i])
			d[i] = Displacement[LocationMatrix_[i]-1];
    
    unsigned int Node_NDF = ND_ / NEN_;
// Resolving non-homogeneous essential boundary conditions
    if (NonHomo_)
        for (unsigned int i = 0; i < ND_; i++)
            if (!LocationMatrix_[i] && nodes_[i/Node_NDF]->BC[i%Node_NDF] != 0)
                d[i] = nodes_[i/Node_NDF]->BC[i%Node_NDF];

// x and y coordinates
    double C[8];
	for (unsigned int i = 0; i < 4; i++)
	{
		C[2*i] = nodes_[i]->XYZ[0];
        C[2*i+1] = nodes_[i]->XYZ[1];
    };   

// Calculate element stress , sigma = D * B * d
    double GaussPoints[2] = {
            -sqrt(3) / 3.0,
            sqrt(3) / 3.0,
    };
    
    int index = 0;
    for (unsigned int gp1 = 0; gp1 < 2; gp1++)
    {
        for (unsigned int gp2 = 0; gp2 < 2; gp2++)
        {
            double B[3][8] = {0};
            double det;
            ElementStrainFunction(B, &det, GaussPoints[gp1], GaussPoints[gp2]);            
            for (unsigned int i = 0; i < 3; i++)
                for (unsigned int j = 0; j < 3; j++)
                    for (unsigned int k = 0; k < 8; k++)
                        stress[5*index+2+i] += D[i][j] * B[j][k] * d[k];    
            double N[2][8] = {0};
            ElementShapeFunction(N, GaussPoints[gp1], GaussPoints[gp2]);
            for (unsigned int i = 0; i < 2; i++)
                for (unsigned int j = 0; j < 8; j++)
                    stress[5*index+i] += N[i][j] * C[j];
            index++;
        }
    }
}

//	Calculate element non-homogeneous essential boundary conditions
void CQ4::ElementNonHomo(double* Matrix, double* NonForce)
{   
    double Ke[8][8] = {0};
    unsigned int index = 0;
    for (int j = 0; j < 8; j++) 
        for (int i = j; i >= 0; i--) 
            Ke[i][j] = Matrix[index++];

    for (unsigned int i = 0; i < 8; i++) 
        for (unsigned int j = 0; j < i; j++) 
            Ke[i][j] = Ke[j][i];        

    double d[8] = {0};
    for (unsigned int i = 0; i < 4; i++)
    {
        d[2*i] = nodes_[i]->BC[0];
        d[2*i+1] = nodes_[i]->BC[1];
    }

    for (unsigned int i = 0; i < 8; i++)
    {   
        if (LocationMatrix_[i] == 0)
            continue;
        for (unsigned int j = 0; j < 8; j++)
        {
            if (LocationMatrix_[j] != 0)
                continue;
            NonForce[i] += Ke[i][j] * d[j];
        }
    }
}

// Calculate Q4 element shape function matrix N
void CQ4::ElementShapeFunction(double (&N)[2][8], double xi, double eta)
{   
    double n[4] = {
        0.25 * (1.0 - xi) * (1.0 - eta),
	    0.25 * (1.0 + xi) * (1.0 - eta),
	    0.25 * (1.0 + xi) * (1.0 + eta),
	    0.25 * (1.0 - xi) * (1.0 + eta),
    };
    
    for (unsigned i = 0; i < 2; i++)
        for (unsigned j = 0; j < 4; j++)
            N[i][(2*j+i)] = n[j];
}

//!	Calculate derivative of Q4 element shape function matrix B and Jacobian determination
void CQ4::ElementStrainFunction(double (&B)[3][8], double* det, double xi, double eta)
{
    // Calculate the Grad(N) matrix, natural coordinate
    double GN[2][4] = {
        {0.25 * (eta - 1.0), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), 0.25 * (- eta - 1.0)},
        {0.25 * (xi - 1.0), 0.25 * (- xi - 1.0), 0.25 * (1.0 + xi), 0.25 * (1.0 - xi)},
    };
    // x and y coordinates
    double C[4][2];
	for (unsigned int i = 0; i < 4; i++)
	{
		C[i][0] = nodes_[i]->XYZ[0];
        C[i][1] = nodes_[i]->XYZ[1];
    };   

    // Calculate Jacobian matrix
    double J[2][2] = {0};
    for (unsigned int i = 0; i < 2; i++) 
        for (unsigned int j = 0; j < 2; j++) 
            for (unsigned int k = 0; k < 4; k++) 
                J[i][j] += GN[i][k] * C[k][j];
	// Calculate determinant of Jacobian
    *det = J[0][0] * J[1][1] - J[0][1] * J[1][0];

    // Calculate inverse matrix of Jacobian
    double invJ[2][2];
    double invDet = 1.0 / (*det);
    invJ[0][0] =  J[1][1] * invDet;
    invJ[0][1] = -J[0][1] * invDet;
    invJ[1][0] = -J[1][0] * invDet;
    invJ[1][1] =  J[0][0] * invDet;

    // Calculate the Grad matrix, physical coordinate
    double G[2][4] = {0};
    for (unsigned int i = 0; i < 2; i++)
        for (unsigned int j = 0; j < 4; j++)
            for (unsigned int k = 0; k < 2; k++) 
                G[i][j] += invJ[i][k] * GN[k][j];
    
    // Derivative of Q4 element shape function matrix B
    B[0][0] = G[0][0];
    B[0][2] = G[0][1];
    B[0][4] = G[0][2];
    B[0][6] = G[0][3];

    B[1][1] = G[1][0];
    B[1][3] = G[1][1];
    B[1][5] = G[1][2];
    B[1][7] = G[1][3];

    B[2][0] = G[1][0];
    B[2][1] = G[0][0];
    B[2][2] = G[1][1];
    B[2][3] = G[0][1];
    B[2][4] = G[1][2];
    B[2][5] = G[0][2];
    B[2][6] = G[1][3];    
    B[2][7] = G[0][3];
}

//! Calcualte derivative of element shape function matrix DN at coordinate xt
void CQ4::ElementDerivativeShapeFunction(double (&DN)[2][4], double xi, double eta)
{
    // Calculate the Grad(N) matrix, natural coordinate
    double GN[2][4] = {
        {0.25 * (eta - 1.0), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), 0.25 * (- eta - 1.0)},
        {0.25 * (xi - 1.0), 0.25 * (- xi - 1.0), 0.25 * (1.0 + xi), 0.25 * (1.0 - xi)},
    };
    // x and y coordinates
    double C[4][2];
	for (unsigned int i = 0; i < 4; i++)
	{
		C[i][0] = nodes_[i]->XYZ[0];
        C[i][1] = nodes_[i]->XYZ[1];
    };   

    // Calculate Jacobian matrix
    double J[2][2] = {0};
    for (unsigned int i = 0; i < 2; i++) 
        for (unsigned int j = 0; j < 2; j++) 
            for (unsigned int k = 0; k < 4; k++) 
                J[i][j] += GN[i][k] * C[k][j];
	// Calculate determinant of Jacobian
    double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];

    // Calculate inverse matrix of Jacobian
    double invJ[2][2];
    double invDet = 1.0 / det;
    invJ[0][0] =  J[1][1] * invDet;
    invJ[0][1] = -J[0][1] * invDet;
    invJ[1][0] = -J[1][0] * invDet;
    invJ[1][1] =  J[0][0] * invDet;

    // Calculate the Grad matrix, physical coordinate
    for (unsigned int i = 0; i < 2; i++)
        for (unsigned int j = 0; j < 4; j++)
            for (unsigned int k = 0; k < 2; k++) 
                DN[i][j] += invJ[i][k] * GN[k][j];
}

