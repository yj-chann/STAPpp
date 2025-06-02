/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Tet4.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CTet4::CTet4()
{
    NEN_ = 4;	// Each element has 4 nodes
    nodes_ = new CNode * [NEN_]; // pointer array to CNode objs 

    ND_ = 12; //    LM Height
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//	Desconstructor
CTet4::~CTet4()
{
}

//	Read element data from stream Input
bool CTet4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;	// Material property set number
    unsigned int N1, N2, N3, N4;	// 1~4 node number

    Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CTet4Material*>(MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];

    return true;
}

//	Write element data to stream
void CTet4::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
        << setw(11) << nodes_[1]->NodeNumber
        << setw(11) << nodes_[2]->NodeNumber
        << setw(11) << nodes_[3]->NodeNumber
        << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CTet4::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    //	Calculate Tet4 element volume
	double X[3];
    double Y[3];
    double Z[3];		//	x, y and z coordinates
	for (unsigned int i = 0; i < 4; i++)
	{
		X[i] = nodes_[i]->XYZ[0];
        Y[i] = nodes_[i]->XYZ[1];
        Z[i] = nodes_[i]->XYZ[2];
    }

    double a1 = X[1]*(Y[2]*Z[3]-Y[3]*Z[2]) - X[2]*(Y[1]*Z[3]-Y[3]*Z[1]) + X[3]*(Y[1]*Z[2]-Y[2]*Z[1]);
    double a2 = - X[0]*(Y[2]*Z[3]-Y[3]*Z[2]) + X[2]*(Y[0]*Z[3]-Y[3]*Z[0]) - X[3]*(Y[0]*Z[2]-Y[2]*Z[0]);
    double a3 = X[0]*(Y[1]*Z[3]-Y[3]*Z[1]) - X[1]*(Y[0]*Z[3]-Y[3]*Z[0]) + X[3]*(Y[0]*Z[1]-Y[1]*Z[0]);
    double a4 = - X[0]*(Y[1]*Z[2]-Y[2]*Z[1]) + X[1]*(Y[0]*Z[2]-Y[2]*Z[0]) - X[2]*(Y[0]*Z[1]-Y[1]*Z[0]);

    double b1 = -(Y[2]*Z[3]-Y[3]*Z[2]) +(Y[1]*Z[3]-Y[3]*Z[1]) -(Y[1]*Z[2]-Y[2]*Z[1]);
    double b2 = (Y[2]*Z[3]-Y[3]*Z[2]) -(Y[0]*Z[3]-Y[3]*Z[0]) +(Y[0]*Z[2]-Y[2]*Z[0]);
    double b3 = -(Y[1]*Z[3]-Y[3]*Z[1]) +(Y[0]*Z[3]-Y[3]*Z[0]) -(Y[0]*Z[1]-Y[1]*Z[0]);
    double b4 = (Y[1]*Z[2]-Y[2]*Z[1]) -(Y[0]*Z[2]-Y[2]*Z[0]) +(Y[0]*Z[1]-Y[1]*Z[0]);

    double c1 = (X[2]*Z[3]-X[3]*Z[2]) -(X[1]*Z[3]-X[3]*Z[1]) +(X[1]*Z[2]-X[2]*Z[1]);
    double c2 = -(X[2]*Z[3]-X[3]*Z[2]) +(X[0]*Z[3]-X[3]*Z[0]) -(X[0]*Z[2]-X[2]*Z[0]);
    double c3 = (X[1]*Z[3]-X[3]*Z[1]) -(X[0]*Z[3]-X[3]*Z[0]) +(X[0]*Z[1]-X[1]*Z[0]);
    double c4 = -(X[1]*Z[2]-X[2]*Z[1]) +(X[0]*Z[2]-X[2]*Z[0]) -(X[0]*Z[1]-X[1]*Z[0]);

    double d1 = -(X[2]*Y[3]-X[3]*Y[2]) +(X[1]*Y[3]-X[3]*Y[1]) -(X[1]*Y[2]-X[2]*Y[1]);
    double d2 = (X[2]*Y[3]-X[3]*Y[2]) -(X[0]*Y[3]-X[3]*Y[0]) +(X[0]*Y[2]-X[2]*Y[0]);
    double d3 = -(X[1]*Y[3]-X[3]*Y[1]) +(X[0]*Y[3]-X[3]*Y[0]) -(X[0]*Y[1]-X[1]*Y[0]);
    double d4 = (X[1]*Y[2]-X[2]*Y[1]) -(X[0]*Y[2]-X[2]*Y[0]) +(X[0]*Y[1]-X[1]*Y[0]);

    double V = fabs(a1+a2+a3+a4) / 6.0;

    //	Calculate element strain matrix
	double B[6][12] = {0};
	B[0][0] = b1/(6.0*V);
    B[0][3] = b2/(6.0*V);
    B[0][6] = b3/(6.0*V);
    B[0][9] = b4/(6.0*V);
    B[1][1] = c1/(6.0*V);
    B[1][4] = c2/(6.0*V);
    B[1][7] = c3/(6.0*V);
    B[1][10] = c4/(6.0*V);
    B[2][2] = d1/(6.0*V);
    B[2][5] = d2/(6.0*V);
    B[2][8] = d3/(6.0*V);
    B[2][11] = d4/(6.0*V);
    B[3][0] = c1/(6.0*V);
    B[3][1] = b1/(6.0*V);
    B[3][3] = c2/(6.0*V);
    B[3][4] = b2/(6.0*V);
    B[3][6] = c3/(6.0*V);
    B[3][7] = b3/(6.0*V);
    B[3][9] = c4/(6.0*V);
    B[3][10] = b4/(6.0*V);
    B[4][0] = d1/(6.0*V);
    B[4][2] = b1/(6.0*V);
    B[4][3] = d2/(6.0*V);
    B[4][5] = b2/(6.0*V);
    B[4][6] = d3/(6.0*V);
    B[4][8] = b3/(6.0*V);
    B[4][9] = d4/(6.0*V);
    B[4][11] = b4/(6.0*V);
    B[5][1] = d1/(6.0*V);
    B[5][2] = c1/(6.0*V);
    B[5][4] = d2/(6.0*V);
    B[5][5] = c2/(6.0*V);
    B[5][7] = d3/(6.0*V);
    B[5][8] = c3/(6.0*V);
    B[5][10] = d4/(6.0*V);
    B[5][11] = c4/(6.0*V);

    //	Calculate element stiffness matrix

    CTet4Material* material_ = dynamic_cast<CTet4Material*>(ElementMaterial_);	// Pointer to material of the element

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

	double Ke[12][12] = {0};
    for (int i = 0; i < 12; i++)
		for (int j = 0; j < 12; j++)
			for (int k = 0; k < 6; k++)
				for (int l = 0; l < 6; l++)
					Ke[i][j] += V * B[k][i] * D[k][l] * B[l][j];

    int index = 0;
    for (int j = 0; j < 12; j++)	// Column-major order
    {
        for (int i = j; i >= 0; i--) {    // Upper triangle (i ≤ j)
            Matrix[index++] = Ke[i][j];   // Pack into one dimensional element stiffness matrix
        }
    }
}

//	Calculate element stress
void CTet4::ElementStress(double* stress, double* Displacement)
{
    //	Calculate Tet4 element volume
	double X[3];
    double Y[3];
    double Z[3];		//	x, y and z coordinates
	for (unsigned int i = 0; i < 4; i++)
	{
		X[i] = nodes_[i]->XYZ[0];
        Y[i] = nodes_[i]->XYZ[1];
        Z[i] = nodes_[i]->XYZ[2];
    }

    double a1 = X[1]*(Y[2]*Z[3]-Y[3]*Z[2]) - X[2]*(Y[1]*Z[3]-Y[3]*Z[1]) + X[3]*(Y[1]*Z[2]-Y[2]*Z[1]);
    double a2 = - X[0]*(Y[2]*Z[3]-Y[3]*Z[2]) + X[2]*(Y[0]*Z[3]-Y[3]*Z[0]) - X[3]*(Y[0]*Z[2]-Y[2]*Z[0]);
    double a3 = X[0]*(Y[1]*Z[3]-Y[3]*Z[1]) - X[1]*(Y[0]*Z[3]-Y[3]*Z[0]) + X[3]*(Y[0]*Z[1]-Y[1]*Z[0]);
    double a4 = - X[0]*(Y[1]*Z[2]-Y[2]*Z[1]) + X[1]*(Y[0]*Z[2]-Y[2]*Z[0]) - X[2]*(Y[0]*Z[1]-Y[1]*Z[0]);

    double b1 = -(Y[2]*Z[3]-Y[3]*Z[2]) +(Y[1]*Z[3]-Y[3]*Z[1]) -(Y[1]*Z[2]-Y[2]*Z[1]);
    double b2 = (Y[2]*Z[3]-Y[3]*Z[2]) -(Y[0]*Z[3]-Y[3]*Z[0]) +(Y[0]*Z[2]-Y[2]*Z[0]);
    double b3 = -(Y[1]*Z[3]-Y[3]*Z[1]) +(Y[0]*Z[3]-Y[3]*Z[0]) -(Y[0]*Z[1]-Y[1]*Z[0]);
    double b4 = (Y[1]*Z[2]-Y[2]*Z[1]) -(Y[0]*Z[2]-Y[2]*Z[0]) +(Y[0]*Z[1]-Y[1]*Z[0]);

    double c1 = (X[2]*Z[3]-X[3]*Z[2]) -(X[1]*Z[3]-X[3]*Z[1]) +(X[1]*Z[2]-X[2]*Z[1]);
    double c2 = -(X[2]*Z[3]-X[3]*Z[2]) +(X[0]*Z[3]-X[3]*Z[0]) -(X[0]*Z[2]-X[2]*Z[0]);
    double c3 = (X[1]*Z[3]-X[3]*Z[1]) -(X[0]*Z[3]-X[3]*Z[0]) +(X[0]*Z[1]-X[1]*Z[0]);
    double c4 = -(X[1]*Z[2]-X[2]*Z[1]) +(X[0]*Z[2]-X[2]*Z[0]) -(X[0]*Z[1]-X[1]*Z[0]);

    double d1 = -(X[2]*Y[3]-X[3]*Y[2]) +(X[1]*Y[3]-X[3]*Y[1]) -(X[1]*Y[2]-X[2]*Y[1]);
    double d2 = (X[2]*Y[3]-X[3]*Y[2]) -(X[0]*Y[3]-X[3]*Y[0]) +(X[0]*Y[2]-X[2]*Y[0]);
    double d3 = -(X[1]*Y[3]-X[3]*Y[1]) +(X[0]*Y[3]-X[3]*Y[0]) -(X[0]*Y[1]-X[1]*Y[0]);
    double d4 = (X[1]*Y[2]-X[2]*Y[1]) -(X[0]*Y[2]-X[2]*Y[0]) +(X[0]*Y[1]-X[1]*Y[0]);

    double V = fabs(a1+a2+a3+a4) / 6.0;

    //	Calculate element strain matrix
	double B[6][12] = {0};
	B[0][0] = b1/(6.0*V);
    B[0][3] = b2/(6.0*V);
    B[0][6] = b3/(6.0*V);
    B[0][9] = b4/(6.0*V);
    B[1][1] = c1/(6.0*V);
    B[1][4] = c2/(6.0*V);
    B[1][7] = c3/(6.0*V);
    B[1][10] = c4/(6.0*V);
    B[2][2] = d1/(6.0*V);
    B[2][5] = d2/(6.0*V);
    B[2][8] = d3/(6.0*V);
    B[2][11] = d4/(6.0*V);
    B[3][0] = c1/(6.0*V);
    B[3][1] = b1/(6.0*V);
    B[3][3] = c2/(6.0*V);
    B[3][4] = b2/(6.0*V);
    B[3][6] = c3/(6.0*V);
    B[3][7] = b3/(6.0*V);
    B[3][9] = c4/(6.0*V);
    B[3][10] = b4/(6.0*V);
    B[4][0] = d1/(6.0*V);
    B[4][2] = b1/(6.0*V);
    B[4][3] = d2/(6.0*V);
    B[4][5] = b2/(6.0*V);
    B[4][6] = d3/(6.0*V);
    B[4][8] = b3/(6.0*V);
    B[4][9] = d4/(6.0*V);
    B[4][11] = b4/(6.0*V);
    B[5][1] = d1/(6.0*V);
    B[5][2] = c1/(6.0*V);
    B[5][4] = d2/(6.0*V);
    B[5][5] = c2/(6.0*V);
    B[5][7] = d3/(6.0*V);
    B[5][8] = c3/(6.0*V);
    B[5][10] = d4/(6.0*V);
    B[5][11] = c4/(6.0*V);

    CTet4Material* material_ = dynamic_cast<CTet4Material*>(ElementMaterial_);	// Pointer to material of the element

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
    double d[12] = { 0 };
    for (unsigned int i = 0; i < 12; i++)
        if (LocationMatrix_[i])
            d[i] = Displacement[LocationMatrix_[i] - 1];

    unsigned int Node_NDF = ND_ / NEN_;
    // Resolving non-homogeneous essential boundary conditions
    if (NonHomo_)
        for (unsigned int i = 0; i < ND_; i++)
            if (!LocationMatrix_[i] && nodes_[i/Node_NDF]->BC[i%Node_NDF] != 0)
                d[i] = nodes_[i/Node_NDF]->BC[i%Node_NDF];

    // Calculate element stress (constant), sigma = D * B * d
    double Bd[6] = {0};
    for (unsigned int i = 0; i < 6; i++)
        for (unsigned int j = 0; j < 12; j++)
            Bd[i] += B[i][j] * d[j];

	*stress = *(stress + 1) = *(stress + 2) = *(stress + 3) = *(stress + 4) = *(stress + 5) = 0.0;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            stress[i] += D[i][j] * Bd[j];
}

//	Calculate element non-homogeneous essential boundary conditions
void CTet4::ElementNonHomo(double* Matrix, double* NonForce)
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
        d[3*i] = nodes_[i]->BC[0];
        d[3*i+1] = nodes_[i]->BC[1];
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
