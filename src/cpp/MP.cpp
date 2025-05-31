/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "MP.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CMP::CMP()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 12;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CMP::~CMP()
{
}

//	Read element data from stream Input
bool CMP::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// first node number, second node number, third node number and fourth node number, Counterclockwise

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CMPMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void CMP::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
           << setw(11) << nodes_[1]->NodeNumber
           << setw(11) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//! Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CMP::GenerateLocationMatrix()
{
    
    LocationMatrix_[0] = nodes_[0]->bcode[3];
    LocationMatrix_[1] = nodes_[0]->bcode[4];
    LocationMatrix_[2] = nodes_[0]->bcode[2];

    LocationMatrix_[3] = nodes_[1]->bcode[3];
    LocationMatrix_[4] = nodes_[1]->bcode[4];
    LocationMatrix_[5] = nodes_[1]->bcode[2];

    LocationMatrix_[6] = nodes_[2]->bcode[3];
    LocationMatrix_[7] = nodes_[2]->bcode[4];
    LocationMatrix_[8] = nodes_[2]->bcode[2];

    LocationMatrix_[9] = nodes_[3]->bcode[3];
    LocationMatrix_[10] = nodes_[3]->bcode[4];
    LocationMatrix_[11] = nodes_[3]->bcode[2];
    
    int break_all = 0;
    for (unsigned int N = 0; N < NEN_ && !break_all; N++) {
        for (unsigned int D = 2; D < 5 && !break_all; D++) {
            if (nodes_[N]->bcode[D] == 0 && nodes_[N]->BC[D] != 0) {
                NonHomo_ = 1;
                break_all = 1;
                break;
            }
        }                
    }      
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CMP::ElementStiffness(double* Matrix)  // Use selective reduced integration
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate element stiffness matrix

	CMPMaterial* material_ = dynamic_cast<CMPMaterial*>(ElementMaterial_);	// Pointer to material of the element

    double E = material_->E;
    double G = material_->G;
    double nu = material_->nu;
	double t = material_->t;
	double shcof = material_->k;  // Correction factor
    double D_0 = material_->D_0;

	double Db[3][3] = {0};	// Build the constitutive matrix D (3x3) for plane stress    
	Db[0][0] = D_0;
    Db[0][1] = D_0 * nu;
    Db[1][0] = D_0 * nu;
    Db[1][1] = D_0;
    Db[2][2] = D_0 * (1.0 - nu) / 2.0;

    double GaussPoints[2] = {
        -sqrt(3) / 3.0,
        sqrt(3) / 3.0,
    };

    double Ke[12][12] = {0};
    for (unsigned int gp1 = 0; gp1 < 2; gp1++)
    {
        for (unsigned int gp2 = 0; gp2 < 2; gp2++)
        {
            double Bb[3][12] = {0};
            double Bs[2][12] = {0};
            double det;
            ElementStrainFunction(Bb, Bs, &det, GaussPoints[gp1], GaussPoints[gp2]);
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            Ke[i][j] += Bb[k][i] * Db[k][l] * Bb[l][j] * det;  // Weight is 1
        }
    }

    double Ds[2][2] = {0};
    Ds[0][0] = shcof * G * t;
    Ds[1][1] = shcof * G * t;

    double Bb[3][12] = {0};
    double Bs[2][12] = {0};
    double det;
    ElementStrainFunction(Bb, Bs, &det, 0.0, 0.0);
    for (int i = 0; i < 12; i++)
        for (int j = 0; j < 12; j++)
            for (int k = 0; k < 2; k++)
                for (int l = 0; l < 2; l++)
                    Ke[i][j] += 2.0 * 2.0 * Bs[k][i] * Ds[k][l] * Bs[l][j] * det;  // Weight is 2

	int index = 0;
    for (int j = 0; j < 12; j++)	// Column-major order
	{         
        for (int i = j; i >= 0; i--) {    // Upper triangle (i ≤ j)
            Matrix[index++] = Ke[i][j];   // Pack into one dimensional element stiffness matrix
        }
    }		
}

//	Calculate element stress, for CST element, the strain and stress is constant 
void CMP::ElementStress(double* stress, double* Displacement)
{	    

}

//	Calculate element non-homogeneous essential boundary conditions
void CMP::ElementNonHomo(double* Matrix, double* NonForce)
{   
    double Ke[12][12] = {0};
    unsigned int index = 0;
    for (int j = 0; j < 12; j++) 
        for (int i = j; i >= 0; i--) 
            Ke[i][j] = Matrix[index++];

    for (unsigned int i = 0; i < 12; i++) 
        for (unsigned int j = 0; j < i; j++) 
            Ke[i][j] = Ke[j][i];        

    double d[12] = {0};
    for (unsigned int i = 0; i < 4; i++)
    {
        d[3*i] = nodes_[i]->BC[3];
        d[3*i+1] = nodes_[i]->BC[4];
        d[3*i+2] = nodes_[i]->BC[2];
    }

    for (unsigned int i = 0; i < 12; i++)
    {   
        if (LocationMatrix_[i] == 0)
            continue;
        for (unsigned int j = 0; j < 12; j++)
        {
            if (LocationMatrix_[j] != 0)
                continue;
            NonForce[i] += Ke[i][j] * d[j];
        }
    }
}

// Calculate Mindlin-Reissner Plate element shape function matrix N
void CMP::ElementShapeFunction(double (&N)[3][12], double xi, double eta)
{   
    double n[4] = {
        0.25 * (1.0 - xi) * (1.0 - eta),
	    0.25 * (1.0 + xi) * (1.0 - eta),
	    0.25 * (1.0 + xi) * (1.0 + eta),
	    0.25 * (1.0 - xi) * (1.0 + eta),
    };
    
    for (unsigned i = 0; i < 3; i++)
        for (unsigned j = 0; j < 4; j++)
            N[i][(3*j+i)] = n[j];
}

//!	Calculate Mindlin-Reissner Plate element bending kinematic matrix Bb, shear kinematic matrix Bs and Jacobian determination
void CMP::ElementStrainFunction(double (&Bb)[3][12], double (&Bs)[2][12], double* det, double xi, double eta)
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

    // Shape function
    double n[4] = {
        0.25 * (1.0 - xi) * (1.0 - eta),
	    0.25 * (1.0 + xi) * (1.0 - eta),
	    0.25 * (1.0 + xi) * (1.0 + eta),
	    0.25 * (1.0 - xi) * (1.0 + eta),
    };
    
    // Bending kinematic matrix Bb
    Bb[0][0] = G[0][0];
    Bb[0][3] = G[0][1];
    Bb[0][6] = G[0][2];
    Bb[0][9] = G[0][3];

    Bb[1][1] = G[1][0];
    Bb[1][4] = G[1][1];
    Bb[1][7] = G[1][2];
    Bb[1][10] = G[1][3];

    Bb[2][0] = G[1][0];
    Bb[2][1] = G[0][0];
    Bb[2][3] = G[1][1];
    Bb[2][4] = G[0][1];
    Bb[2][6] = G[1][2];
    Bb[2][7] = G[0][2];
    Bb[2][9] = G[1][3];    
    Bb[2][10] = G[0][3];

    // Shear kinematic matrix Bs
    Bs[0][0] = -n[0];
    Bs[1][1] = -n[0];
    Bs[0][2] = G[0][0];
    Bs[1][2] = G[1][0];

    Bs[0][3] = -n[1];
    Bs[1][4] = -n[1];
    Bs[0][5] = G[0][1];
    Bs[1][5] = G[1][1];

    Bs[0][6] = -n[2];
    Bs[1][7] = -n[2];
    Bs[0][8] = G[0][2];
    Bs[1][8] = G[1][2];

    Bs[0][9] = -n[3];
    Bs[1][10] = -n[3];
    Bs[0][11] = G[0][3];
    Bs[1][11] = G[1][3];
}
