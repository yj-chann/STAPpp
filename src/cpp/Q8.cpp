/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Q8.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CQ8::CQ8()
{
    NEN_ = 8;	// Each element has 8 nodes
    nodes_ = new CNode * [NEN_]; // pointer array to CNode objs 

    ND_ = 16; //    LM Height
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//	Desconstructor
CQ8::~CQ8()
{
}

//	Read element data from stream Input
bool CQ8::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;	// Material property set number
    unsigned int N1, N2, N3, N4, N5, N6, N7, N8;	// 1~8 node number, Counterclockwise

    Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;
    ElementMaterial_ = dynamic_cast<CQ8Material*>(MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];
    nodes_[4] = &NodeList[N5 - 1];
    nodes_[5] = &NodeList[N6 - 1];
    nodes_[6] = &NodeList[N7 - 1];
    nodes_[7] = &NodeList[N8 - 1];

    return true;
}

//	Write element data to stream
void CQ8::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
        << setw(11) << nodes_[1]->NodeNumber
        << setw(11) << nodes_[2]->NodeNumber
        << setw(11) << nodes_[3]->NodeNumber 
        << setw(11) << nodes_[4]->NodeNumber
        << setw(11) << nodes_[5]->NodeNumber
        << setw(11) << nodes_[6]->NodeNumber
        << setw(11) << nodes_[7]->NodeNumber
        <<setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CQ8::ElementStiffness(double* Matrix)  // use reduced integration (2*2)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    //	Calculate element stiffness matrix

    CQ8Material* material_ = dynamic_cast<CQ8Material*>(ElementMaterial_);	// Pointer to material of the element

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

    double D[3][3] = { 0 };	// Build the constitutive matrix D (3x3) for plane stress    
    D[0][0] = k;
    D[0][1] = k * nu;
    D[1][0] = k * nu;
    D[1][1] = k;
    D[2][2] = k * (1.0 - nu) / 2.0;

    double GaussPoints[2] = {
        -sqrt(3) / 3.0,
        sqrt(3) / 3.0,
    };
    double Ke[16][16] = { 0 };
    for (unsigned int i = 0; i < 2; i++)
    {
        for (unsigned int j = 0; j < 2; j++)
        {
            double B[3][16] = { 0 };
            double det;
            ElementStrainFunction(B, &det, GaussPoints[i], GaussPoints[j]);
            for (int loop1 = 0; loop1 < 16; loop1++)
                for (int loop2 = 0; loop2 < 16; loop2++)
                    for (int loop3 = 0; loop3 < 3; loop3++)
                        for (int loop4 = 0; loop4 < 3; loop4++)
                            Ke[loop1][loop2] += t * B[loop3][loop1] * D[loop3][loop4] * B[loop4][loop2] * det;  // Weight is 1
        }
    }

    int index = 0;
    for (int j = 0; j < 16; j++)	// Column-major order
    {
        for (int i = j; i >= 0; i--) {    // Upper triangle (i �� j)
            Matrix[index++] = Ke[i][j];   // Pack into one dimensional element stiffness matrix
        }
    }
}

//	Calculate element stress
void CQ8::ElementStress(double* stress, double* Displacement)
{
    CQ8Material* material_ = dynamic_cast<CQ8Material*>(ElementMaterial_);	// Pointer to material of the element

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

    double D[3][3] = { 0 };	// Build the constitutive matrix D (3x3) for plane stress    
    D[0][0] = k;
    D[0][1] = k * nu;
    D[1][0] = k * nu;
    D[1][1] = k;
    D[2][2] = k * (1.0 - nu) / 2.0;

    //  Get element displacement vector
    double d[16] = { 0 };
    for (unsigned int i = 0; i < 16; i++)
        if (LocationMatrix_[i])
            d[i] = Displacement[LocationMatrix_[i] - 1];

    // x and y coordinates
    double C[16];
    for (unsigned int i = 0; i < 8; i++)
    {
        C[2 * i] = nodes_[i]->XYZ[0];
        C[2 * i + 1] = nodes_[i]->XYZ[1];
    };

    // Calculate element stress , sigma = D * B * d
    double GaussPoints[2] = {
            -sqrt(3) / 3.0,
            sqrt(3) / 3.0,
    };
    int index = 0;
    for (unsigned int i = 0; i < 2; i++)
    {
        for (unsigned int j = 0; j < 2; j++)
        {
            double B[3][16] = { 0 };
            double det;
            ElementStrainFunction(B, &det, GaussPoints[i], GaussPoints[j]);
            for (unsigned int k = 0; k < 3; k++)
                for (unsigned int l = 0; l < 3; l++)
                    for (unsigned int m = 0; m < 16; m++)
                        stress[5 * index + 2 + k] += D[k][l] * B[l][m] * d[m];
            double N[2][16] = { 0 };
            ElementShapeFunction(N, GaussPoints[i], GaussPoints[j]);
            for (unsigned int k = 0; k < 2; k++)
                for (unsigned int l = 0; l < 16; l++)
                    stress[5 * index + k] += N[k][l] * C[l];
            index++;
        }
    }
}

//	Calculate element non-homogeneous essential boundary conditions
void CQ8::ElementNonHomo(double* Matrix, double* NonForce)
{   
    double Ke[16][16] = {0};
    unsigned int index = 0;
    for (int j = 0; j < 16; j++) 
        for (int i = j; i >= 0; i--) 
            Ke[i][j] = Matrix[index++];

    for (unsigned int i = 0; i < 16; i++) 
        for (unsigned int j = 0; j < i; j++) 
            Ke[i][j] = Ke[j][i];        

    double d[16] = {0};
    for (unsigned int i = 0; i < 8; i++)
    {
        d[2*i] = nodes_[i]->BC[0];
        d[2*i+1] = nodes_[i]->BC[1];
    }

    for (unsigned int i = 0; i < 16; i++)
    {   
        if (LocationMatrix_[i] == 0)
            continue;
        for (unsigned int j = 0; j < 16; j++)
        {
            if (LocationMatrix_[j] != 0)
                continue;
            NonForce[i] += Ke[i][j] * d[j];
        }
    }
}

// Calculate Q8 element shape function matrix N
void CQ8::ElementShapeFunction(double(&N)[2][16], double xi, double eta)
{
    double n[8] = {
        -0.25 * (1.0 - xi) * (1.0 - eta) * (xi + eta + 1.0),
        0.25 * (1.0 + xi) * (1.0 - eta) * (xi - eta - 1.0),
        0.25 * (1.0 + xi) * (1.0 + eta) * (xi + eta - 1.0),
        -0.25 * (1.0 - xi) * (1.0 + eta) * (xi - eta + 1.0),
        0.5 * (1.0 - xi * xi) * (1.0 - eta),
		0.5 * (1.0 + xi) * (1.0 - eta * eta),
		0.5 * (1.0 - xi * xi) * (1.0 + eta),
		0.5 * (1.0 - xi) * (1.0 - eta * eta),
    };

    for (unsigned i = 0; i < 2; i++)
        for (unsigned j = 0; j < 8; j++)
            N[i][(2 * j + i)] = n[j];
}

//!	Calculate derivative of Q8 element shape function matrix B and Jacobian determination
void CQ8::ElementStrainFunction(double(&B)[3][16], double* det, double xi, double eta)
{
    // Calculate the Grad(N) matrix, natural coordinate (scalar=4)
	double GN[2][8] = {
		{-(eta + 2.0 * xi) * (eta - 1), (eta - 2.0 * xi) * (eta - 1), (eta + 2.0 * xi) * (1 + eta), (eta - 2.0 * xi) * (-eta - 1),4.0 * xi * (eta - 1),2.0 * (1 - eta * eta),-4.0 * xi * (eta + 1),2.0 * (eta * eta - 1)},
		{-(2.0 * eta + xi) * (xi - 1), (2.0 * eta - xi) * (xi + 1), (2.0 * eta + xi) * (1 + xi), (2.0 * eta - xi) * (1 - xi),2.0 * (xi * xi - 1),-4.0 * eta * (xi + 1),2.0 * (1 - xi * xi),4.0 * eta * (xi - 1)},
	};
    
    // x and y coordinates
    double C[8][2];
    for (unsigned int i = 0; i < 8; i++)
    {
        C[i][0] = nodes_[i]->XYZ[0];
        C[i][1] = nodes_[i]->XYZ[1];
    };

    // Calculate Jacobian matrix element (scalar=4) 
    double J[2][2] = { 0 };
    
    // calculate dx & dy
    double x21 = C[1][0] - C[0][0]; // x2 - x1
    double x43 = C[3][0] - C[2][0]; // x4 - x3
    double x68 = C[5][0] - C[7][0]; // x6 - x8
    double x57 = C[4][0] - C[6][0]; // x5 - x7
    double x41 = C[3][0] - C[0][0]; // x4 - x1
    double x32 = C[2][0] - C[1][0]; // x3 - x2
    double x86 = C[7][0] - C[5][0]; // x8 - x6 

    double y21 = C[1][1] - C[0][1]; // y2 - y1
    double y43 = C[3][1] - C[2][1]; // y4 - y3
    double y68 = C[5][1] - C[7][1]; // y6 - y8
    double y57 = C[4][1] - C[6][1]; // y5 - y7
    double y41 = C[3][1] - C[0][1]; // y4 - y1
    double y32 = C[2][1] - C[1][1]; // y3 - y2
    double y86 = C[7][1] - C[5][1]; // y8 - y6 


    J[0][0] = x21 * eta * (eta - 1.0) - x43 * eta * (eta + 1)
        - 2.0 * xi * (eta - 1) * (C[0][0] + C[1][0])
        + 2.0 * xi * (eta + 1) * (C[2][0] + C[3][0])
        + 2.0 * (1 - eta * eta) * x68
        + 4.0 * xi * eta * x57
        - 4.0 * xi * (C[4][0] + C[6][0]);


    J[0][1] = y21 * eta * (eta - 1.0) - y43 * eta * (eta + 1)
        - 2.0 * xi * (eta - 1) * (C[0][1] + C[1][1])
        + 2.0 * xi * (eta + 1) * (C[2][1] + C[3][1])
        + 2.0 * (1 - eta * eta) * y68
        + 4.0 * xi * eta * y57
        - 4.0 * xi * (C[4][1] + C[6][1]);


    J[1][0] = x41 * xi * (xi - 1.0) + x32 * xi * (xi + 1)
        - 2.0 * eta * (xi - 1) * (C[0][0] + C[3][0])
        + 2.0 * eta * (xi + 1) * (C[1][0] + C[2][0])
        + 2.0 * (xi * xi - 1) * x57
        + 4.0 * xi * eta * x86
        - 4.0 * eta * (C[5][0] + C[7][0]);


    J[1][1] = y41 * xi * (xi - 1.0) + y32 * xi * (xi + 1)
        - 2.0 * eta * (xi - 1) * (C[0][1] + C[3][1])
        + 2.0 * eta * (xi + 1) * (C[1][1] + C[2][1])
        + 2.0 * (xi * xi - 1) * y57
        + 4.0 * xi * eta * y86
        - 4.0 * eta * (C[5][1] + C[7][1]);
    // Calculate determinant of Jacobian
	*det = 1.0 / 16.0 * (J[0][0] * J[1][1] - J[0][1] * J[1][0]); // exact

    // Calculate inverse matrix of Jacobian
    double invJ[2][2];
	double invDet = 1.0 / (16.0 * (*det));
    invJ[0][0] = J[1][1] * invDet;
    invJ[0][1] = -J[0][1] * invDet;
    invJ[1][0] = -J[1][0] * invDet;
    invJ[1][1] = J[0][0] * invDet;

    // Calculate the Grad matrix, physical coordinate
    double G[2][8] = { 0 };
    for (unsigned int i = 0; i < 2; i++)
        for (unsigned int j = 0; j < 8; j++)
            for (unsigned int k = 0; k < 2; k++)
                G[i][j] += invJ[i][k] * GN[k][j];

    // Derivative of Q8 element shape function matrix B
    for (int i = 0; i < 8; i++) {
        // 1st row��N_{i,x} 0 ...
        B[0][2 * i] = G[0][i];
        B[0][2 * i + 1] = 0.0;

        // 2nd row��0 N_{i,y} ...
        B[1][2 * i] = 0.0;
        B[1][2 * i + 1] = G[1][i];

        // 3rd row��N_{i,y} N_{i,x} ...
        B[2][2 * i] = G[1][i];
        B[2][2 * i + 1] = G[0][i];
    }
}
