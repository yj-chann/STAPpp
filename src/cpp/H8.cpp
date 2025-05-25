/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "H8.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CH8::CH8()
{
    NEN_ = 8;	// Each element has 8 nodes
    nodes_ = new CNode * [NEN_]; // pointer array to CNode objs 

    ND_ = 24; //    LM Height
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//	Desconstructor
CH8::~CH8()
{
}

//	Read element data from stream Input
bool CH8::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;	// Material property set number
    unsigned int N1, N2, N3, N4, N5, N6, N7, N8;;	// 1~8 node number, Counterclockwise

    Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;
    ElementMaterial_ = dynamic_cast<CH8Material*>(MaterialSets) + MSet - 1;
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
void CH8::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
        << setw(11) << nodes_[1]->NodeNumber
        << setw(11) << nodes_[2]->NodeNumber
        << setw(11) << nodes_[3]->NodeNumber
        << setw(11) << nodes_[4]->NodeNumber
        << setw(11) << nodes_[5]->NodeNumber
        << setw(11) << nodes_[6]->NodeNumber
        << setw(11) << nodes_[7]->NodeNumber
        << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CH8::ElementStiffness(double* Matrix)  // use full integration (2x2x2)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    //	Calculate element stiffness matrix

    CH8Material* material_ = dynamic_cast<CH8Material*>(ElementMaterial_);	// Pointer to material of the element

    double E = material_->E;
    double nu = material_->nu;

	double k = E / (1.0 + nu) / (1.0 - 2.0 * nu);

    double D[6][6] = { 0 };	// Build the constitutive matrix D (6x6) 
	D[0][0] = k * (1.0 - nu);
    D[0][1] = k * nu;
    D[0][2] = k * nu;
    D[1][0] = k * nu;
    D[1][1] = k * (1.0 - nu);
    D[1][2] = k * nu;
    D[2][0] = k * nu;
    D[2][1] = k * nu;
    D[2][2] = k * (1.0 - nu);
	D[3][3] = k * (1.0 - 2.0 * nu) / 2.0;
    D[4][4] = k * (1.0 - 2.0 * nu) / 2.0;
    D[5][5] = k * (1.0 - 2.0 * nu) / 2.0;

    double GaussPoints[2] = {
        -sqrt(3) / 3.0,
        sqrt(3) / 3.0,
    };
    double Ke[24][24] = { 0 };
    for (unsigned int i = 0; i < 2; i++)
    {
        for (unsigned int j = 0; j < 2; j++)
        {
            for (unsigned int k = 0; k < 2; k++)
            {
                double B[6][24] = { 0 };
                double det;
                ElementStrainFunction(B, &det, GaussPoints[i], GaussPoints[j],GaussPoints[k]);
                for (int loop1 = 0; loop1 < 24; loop1++)
                    for (int loop2 = 0; loop2 < 24; loop2++)
                        for (int loop3 = 0; loop3 < 6; loop3++)
                            for (int loop4 = 0; loop4 < 6; loop4++)
                                Ke[loop1][loop2] += B[loop3][loop1] * D[loop3][loop4] * B[loop4][loop2] * det;  // Weight is 1
            }
        }
    }

    int index = 0;
    for (int j = 0; j < 24; j++)	// Column-major order
    {
        for (int i = j; i >= 0; i--) {    // Upper triangle (i ¡Ü j)
            Matrix[index++] = Ke[i][j];   // Pack into one dimensional element stiffness matrix
        }
    }
}

//	Calculate element stress
void CH8::ElementStress(double* stress, double* Displacement)
{
    CH8Material* material_ = dynamic_cast<CH8Material*>(ElementMaterial_);	// Pointer to material of the element

    double E = material_->E;
    double nu = material_->nu;


    double k = E / (1.0 + nu) / (1.0 - 2.0 * nu);

    double D[6][6] = { 0 };	// Build the constitutive matrix D (6x6) 
    D[0][0] = k * (1.0 - nu);
    D[0][1] = k * nu;
    D[0][2] = k * nu;
    D[1][0] = k * nu;
    D[1][1] = k * (1.0 - nu);
    D[1][2] = k * nu;
    D[2][0] = k * nu;
    D[2][1] = k * nu;
    D[2][2] = k * (1.0 - nu);
    D[3][3] = k * (1.0 - 2.0 * nu) / 2.0;
    D[4][4] = k * (1.0 - 2.0 * nu) / 2.0;
    D[5][5] = k * (1.0 - 2.0 * nu) / 2.0;

    //  Get element displacement vector
    double d[24] = { 0 };
    for (unsigned int i = 0; i < 24; i++)
        if (LocationMatrix_[i])
            d[i] = Displacement[LocationMatrix_[i] - 1];

    // x and y coordinates
    double C[24];
    for (unsigned int i = 0; i < 8; i++)
    {
        C[3 * i] = nodes_[i]->XYZ[0];
        C[3 * i + 1] = nodes_[i]->XYZ[1];
        C[3 * i + 2] = nodes_[i]->XYZ[2];
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
            for (unsigned int g = 0; g< 2; g++)
            {
                double B[6][24] = { 0 };
                double det;
                ElementStrainFunction(B, &det, GaussPoints[i], GaussPoints[j],GaussPoints[g]);
                for (unsigned int k = 0; k < 6; k++)
                    for (unsigned int l = 0; l < 6; l++)
                        for (unsigned int m = 0; m < 24; m++)
                            stress[9 * index + 3 + k] += D[k][l] * B[l][m] * d[m];
                double N[3][24] = { 0 };
                ElementShapeFunction(N, GaussPoints[i], GaussPoints[j], GaussPoints[g]);
                for (unsigned int k = 0; k < 3; k++)
                    for (unsigned int l = 0; l < 24; l++)
                        stress[9 * index + k] += N[k][l] * C[l];
                index++;
            }
        }
    }
}

// Calculate H8 element shape function matrix N
void CH8::ElementShapeFunction(double(&N)[3][24], double xi, double eta,double zeta)
{
	double n[8] = {
	   (1 - xi) * (1 - eta) * (1 - zeta) / 8.0,
	   (1 + xi) * (1 - eta) * (1 - zeta) / 8.0,
	   (1 + xi) * (1 + eta) * (1 - zeta) / 8.0,
	   (1 - xi) * (1 + eta) * (1 - zeta) / 8.0,
	   (1 - xi) * (1 - eta) * (1 + zeta) / 8.0,
	   (1 + xi) * (1 - eta) * (1 + zeta) / 8.0,
	   (1 + xi) * (1 + eta) * (1 + zeta) / 8.0,
	   (1 - xi) * (1 + eta) * (1 + zeta) / 8.0,
	};

    for (unsigned i = 0; i < 3; i++)
        for (unsigned j = 0; j < 8; j++)
            N[i][(3 * j + i)] = n[j];
}

//!	Calculate derivative of H8 element shape function matrix B and Jacobian determination
void CH8::ElementStrainFunction(double(&B)[6][24], double* det, double xi, double eta,double zeta)
{
    // Calculate the Grad(N) matrix, natural coordinate (scalar=8)
    double GN[3][8] = {
        { -(1.0 - eta) * (1.0 - zeta), (1.0 - eta) * (1.0 - zeta), (1.0 + eta) * (1.0 - zeta), -(1.0 + eta) * (1.0 - zeta), -(1.0 - eta) * (1.0 + zeta), (1.0 - eta) * (1.0 + zeta), (1.0 + eta) * (1.0 + zeta), -(1.0 + eta) * (1.0 + zeta)},
        { -(1.0 - xi) * (1.0 - zeta), -(1.0 + xi) * (1.0 - zeta), (1.0 + xi) * (1.0 - zeta), (1.0 - xi) * (1.0 - zeta), -(1.0 - xi) * (1.0 + zeta), -(1.0 + xi) * (1.0 + zeta), (1.0 + xi) * (1.0 + zeta), (1.0 - xi) * (1.0 + zeta)},
        { -(1.0 - xi) * (1.0 - eta), -(1.0 + xi) * (1.0 - eta), -(1.0 + xi) * (1.0 + eta), -(1.0 - xi) * (1.0 + eta), (1.0 - xi) * (1.0 - eta), (1.0 + xi) * (1.0 - eta), (1.0 + xi) * (1.0 + eta), (1.0 - xi) * (1.0 + eta)}
    };

    // x and y and z coordinates
    double C[8][3];
    for (unsigned int i = 0; i < 8; i++)
    {
        C[i][0] = nodes_[i]->XYZ[0];
        C[i][1] = nodes_[i]->XYZ[1];
        C[i][2] = nodes_[i]->XYZ[2];
    };

    // Calculate Jacobian matrix element (scalar=8) 
    double J[3][3] = { 0 };

    // calculate dx & dy & dz
    double x21 = C[1][0] - C[0][0];
    double y21 = C[1][1] - C[0][1];
    double z21 = C[1][2] - C[0][2];

    double x43 = C[3][0] - C[2][0];
    double y43 = C[3][1] - C[2][1];
    double z43 = C[3][2] - C[2][2];

    double x65 = C[5][0] - C[4][0];
    double y65 = C[5][1] - C[4][1];
    double z65 = C[5][2] - C[4][2];

    double x87 = C[7][0] - C[6][0];
    double y87 = C[7][1] - C[6][1];
    double z87 = C[7][2] - C[6][2];

    double x41 = C[3][0] - C[0][0];
    double y41 = C[3][1] - C[0][1];
    double z41 = C[3][2] - C[0][2];

    double x32 = C[2][0] - C[1][0];
    double y32 = C[2][1] - C[1][1];
    double z32 = C[2][2] - C[1][2];

    double x85 = C[7][0] - C[4][0];
    double y85 = C[7][1] - C[4][1];
    double z85 = C[7][2] - C[4][2];

    double x76 = C[6][0] - C[5][0];
    double y76 = C[6][1] - C[5][1];
    double z76 = C[6][2] - C[5][2];

    double x51 = C[4][0] - C[0][0];
    double y51 = C[4][1] - C[0][1];
    double z51 = C[4][2] - C[0][2];

    double x62 = C[5][0] - C[1][0];
    double y62 = C[5][1] - C[1][1];
    double z62 = C[5][2] - C[1][2];

    double x73 = C[6][0] - C[2][0];
    double y73 = C[6][1] - C[2][1];
    double z73 = C[6][2] - C[2][2];

    double x84 = C[7][0] - C[3][0];
    double y84 = C[7][1] - C[3][1];
    double z84 = C[7][2] - C[3][2];

  
    J[0][0] = x21 * (1.0 - eta) * (1.0 - zeta) - x43 * (1.0 + eta) * (1.0 - zeta)
        + x65 * (1.0 - eta) * (1.0 + zeta) - x87 * (1.0 + eta) * (1.0 + zeta);

    J[0][1] = y21 * (1.0 - eta) * (1.0 - zeta) - y43 * (1.0 + eta) * (1.0 - zeta)
        + y65 * (1.0 - eta) * (1.0 + zeta) - y87 * (1.0 + eta) * (1.0 + zeta);

    J[0][2] = z21 * (1.0 - eta) * (1.0 - zeta) - z43 * (1.0 + eta) * (1.0 - zeta)
        + z65 * (1.0 - eta) * (1.0 + zeta) - z87 * (1.0 + eta) * (1.0 + zeta);

    J[1][0] = x41 * (1.0 - xi) * (1.0 - zeta) + x32 * (1.0 + xi) * (1.0 - zeta)
        + x85 * (1.0 - xi) * (1.0 + zeta) + x76 * (1.0 + xi) * (1.0 + zeta);

    J[1][1] = y41 * (1.0 - xi) * (1.0 - zeta) + y32 * (1.0 + xi) * (1.0 - zeta)
        + y85 * (1.0 - xi) * (1.0 + zeta) + y76 * (1.0 + xi) * (1.0 + zeta);

    J[1][2] = z41 * (1.0 - xi) * (1.0 - zeta) + z32 * (1.0 + xi) * (1.0 - zeta)
        + z85 * (1.0 - xi) * (1.0 + zeta) + z76 * (1.0 + xi) * (1.0 + zeta);

    J[2][0] = x51 * (1.0 - xi) * (1.0 - eta) + x62 * (1.0 + xi) * (1.0 - eta)
        + x73 * (1.0 + xi) * (1.0 + eta) + x84 * (1.0 - xi) * (1.0 + eta);

    J[2][1] = y51 * (1.0 - xi) * (1.0 - eta) + y62 * (1.0 + xi) * (1.0 - eta)
        + y73 * (1.0 + xi) * (1.0 + eta) + y84 * (1.0 - xi) * (1.0 + eta);

    J[2][2] = z51 * (1.0 - xi) * (1.0 - eta) + z62 * (1.0 + xi) * (1.0 - eta)
        + z73 * (1.0 + xi) * (1.0 + eta) + z84 * (1.0 - xi) * (1.0 + eta);

    // Calculate determinant of Jacobian
    *det = 1.0 / 512.0 * (J[0][0] * J[1][1] * J[2][2] + J[0][1] * J[1][2] * J[2][0] + J[0][2] * J[1][0] * J[2][1]
        - J[0][2] * J[1][1] * J[2][0] - J[0][1] * J[1][0] * J[2][2] - J[0][0] * J[1][2] * J[2][1]); // exact

    // Calculate inverse matrix of Jacobian
    double invJ[3][3];
    double invDet = 1.0 / (512.0 * (*det));
	invJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * invDet;
	invJ[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) * invDet;
	invJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * invDet;
	invJ[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) * invDet;
	invJ[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * invDet;
	invJ[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) * invDet;
	invJ[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * invDet;
	invJ[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) * invDet;
	invJ[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * invDet;

    // Calculate the Grad matrix, physical coordinate
    double G[3][8] = { 0 };
    for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 8; j++)
            for (unsigned int k = 0; k < 3; k++)
                G[i][j] += invJ[i][k] * GN[k][j];

    // Derivative of H8 element shape function matrix B
    for (int i = 0; i < 8; i++) {
        // 1st row£ºN_{i,x} 0  0...
        B[0][3 * i] = G[0][i];
        B[0][3 * i + 1] = 0.0;
        B[0][3 * i + 2] = 0.0;

        // 2nd row£º0 N_{i,y} 0...
        B[1][3 * i] = 0.0;
        B[1][3 * i + 1] = G[1][i];
		B[1][3 * i + 2] = 0.0;

        // 3rd row£º0 0 N_{i,y} ...
        B[2][3 * i] = 0.0;
		B[2][3 * i + 1] = 0.0;
        B[2][3 * i + 2] = G[2][i];

        // 4th row£ºN_{i,y} N_{i,x},0 ...
        B[3][3 * i] = G[1][i];
        B[3][3 * i + 1] = G[0][i];
        B[3][3 * i + 2] = 0.0;

        // 5th row£ºN_{i,z} 0 N_{i,x} ...
        B[3][3 * i] = G[2][i];
        B[3][3 * i + 1] = 0.0;
        B[3][3 * i + 2] = G[0][i];

        // 6th row£º0 N_{i,z} N_{i,y} ...
        B[3][3 * i] = 0.0;
        B[3][3 * i + 1] = G[2][i];
        B[3][3 * i + 2] = G[1][i];

    }
}
