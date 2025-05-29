/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "S8R5.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CS8R5::CS8R5()
{
	NEN_ = 8;	// Each element has 8 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 40;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CS8R5::~CS8R5()
{
}

//	Read element data from stream Input
bool CS8R5::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8;	// 1~8 node number, Counterclockwise

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;
    ElementMaterial_ = dynamic_cast<CS8R5Material*>(MaterialSets) + MSet - 1;
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
void CS8R5::Write(COutputter& output)
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
//  Use reduced integration
void CS8R5::ElementStiffness(double* Matrix)  // Local stiffness matrix can be derived analytically, then use DCM matrix
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate element stiffness matrix

	CS8R5Material* material_ = dynamic_cast<CS8R5Material*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E;
    double nu = material_->nu;
    double shcof = material_->k;  // Correction factor

    double k = E / (1.0 - nu * nu);
	
	double D[5][5] = {0};	// Build the constitutive matrix D (5x5) in local frame
	D[0][0] = k;
    D[0][1] = k * nu;
    D[1][0] = k * nu;
    D[1][1] = k;
    D[2][2] = k * (1.0 - nu) / 2.0;
    D[3][3] = k * shcof * (1.0 - nu) / 2.0;
    D[4][4] = k * shcof * (1.0 - nu) / 2.0;

    double GaussPoints[2] = {
        -sqrt(3) / 3.0,
        sqrt(3) / 3.0,
    };

    double Ke[40][40] = {0};
    for (unsigned int gp1 = 0; gp1 < 2; gp1++)
    {
        for (unsigned int gp2 = 0; gp2 < 2; gp2++)
        {   
            for (unsigned int gp3 = 0; gp3 < 2; gp3++)
            {
                double B[5][40] = {0};
                double det;
                ElementStrainFunction(B, &det, GaussPoints[gp1], GaussPoints[gp2], GaussPoints[gp3]);
                for (int i = 0; i < 40; i++)
                    for (int j = 0; j < 40; j++)
                        for (int k = 0; k < 5; k++)
                            for (int l = 0; l < 5; l++)
                                Ke[i][j] += B[k][i] * D[k][l] * B[l][j] * det;  // Weight is 1
            }            
        }
    }

	int index = 0;
    for (int j = 0; j < 40; j++)	// Column-major order
	{         
        for (int i = j; i >= 0; i--) {    // Upper triangle (i ≤ j)
            Matrix[index++] = Ke[i][j];   // Pack into one dimensional element stiffness matrix
        }
    }
    
}

//	Calculate element stress, for S8R5 element
void CS8R5::ElementStress(double* stress, double* Displacement)
{	    
	
}

//	Calculate element non-homogeneous essential boundary conditions
void CS8R5::ElementNonHomo(double* Matrix, double* NonForce)
{   
    double Ke[40][40] = {0};
    unsigned int index = 0;
    for (int j = 0; j < 40; j++) 
        for (int i = j; i >= 0; i--) 
            Ke[i][j] = Matrix[index++];

    for (unsigned int i = 0; i < 40; i++) 
        for (unsigned int j = 0; j < i; j++) 
            Ke[i][j] = Ke[j][i];        

    double d[40] = {0};
    for (unsigned int i = 0; i < 8; i++)
    {
        d[5*i] = nodes_[i]->BC[0];
        d[5*i+1] = nodes_[i]->BC[1];
		d[5*i+2] = nodes_[i]->BC[2];
        d[5*i+1] = nodes_[i]->BC[3];
		d[5*i+2] = nodes_[i]->BC[4];
    }

    for (unsigned int i = 0; i < 40; i++)
    {   
        if (LocationMatrix_[i] == 0)
            continue;
        for (unsigned int j = 0; j < 40; j++)
        {
            if (LocationMatrix_[j] != 0)
                continue;
            NonForce[i] += Ke[i][j] * d[j];
        }
    }
}

//	Calculate element shape function matrix N at parent coordinate xi, eta, zeta
void CS8R5::ElementShapeFunction(double (&N)[3][40], double xi, double eta, double zeta)
{
    double N_I[8] = {
        -0.25 * (1.0 - xi) * (1.0 - eta) * (xi + eta + 1.0),
        0.25 * (1.0 + xi) * (1.0 - eta) * (xi - eta - 1.0),
        0.25 * (1.0 + xi) * (1.0 + eta) * (xi + eta - 1.0),
        -0.25 * (1.0 - xi) * (1.0 + eta) * (xi - eta + 1.0),
        0.5 * (1.0 - xi * xi) * (1.0 - eta),
		0.5 * (1.0 + xi) * (1.0 - eta * eta),
		0.5 * (1.0 - xi * xi) * (1.0 + eta),
		0.5 * (1.0 - xi) * (1.0 - eta * eta),
    };

    CS8R5Material* material_ = dynamic_cast<CS8R5Material*>(ElementMaterial_);	// Pointer to material of the element
    double t = material_->t;
    double v3[3] = {0};
    v3[0] = material_->v3[0];
    v3[1] = material_->v3[1];
    v3[2] = material_->v3[2];

    // unsigned int i = 0;
    // for (; i < 3; i++)
    // {
    //     if (v3[i])
    //         break;
    // }
    // double v1[3] = {0}; // unit vector normal to v3, right-hand system
    // double v2[3] = {0}; // unit vector normal to v3, right-hand system
    // v1[(i+1)%3] = 1.0;
    // v2[(i+2)%3] = 1.0;

    double v1[3] = {0};
    double e1[3] = {1.0, 0.0, 0.0};
    for (unsigned int i = 0; i < 3; i++)
        v1[i] = e1[(i+1)%3] * v3[(i+2)%3] - e1[(i+2)%3] * v3[(i+1)%3];

    double v2[3] = {0};
    for (unsigned int i = 0; i < 3; i++)
        v2[i] = v3[(i+1)%3] * v1[(i+2)%3] - v3[(i+2)%3] * v1[(i+1)%3];

    for (unsigned i = 0; i < 8; i++)
    {
        N[0][5*i] = N_I[i];
        N[0][5*i+3] = zeta / 2.0 * t * N_I[i]*v1[0];
        N[0][5*i+4] = -zeta / 2.0 * t * N_I[i]*v2[0];

        N[1][5*i+1] = N_I[i];
        N[1][5*i+3] = zeta / 2.0 * t * N_I[i]*v1[1];
        N[1][5*i+4] = -zeta / 2.0 * t * N_I[i]*v2[1];

        N[2][5*i+2] = N_I[i];
        N[2][5*i+3] = zeta / 2.0 * t * N_I[i]*v1[2];
        N[2][5*i+4] = -zeta / 2.0 * t * N_I[i]*v2[2];
    }
}

//	Calculate derivative of element shape function matrix B(global frame) at parent coordinate xi, eta, zeta
void CS8R5::ElementStrainFunction(double (&B)[5][40], double* det, double xi, double eta, double zeta)
{   
    CS8R5Material* material_ = dynamic_cast<CS8R5Material*>(ElementMaterial_);	// Pointer to material of the element
    double t = material_->t;
    double v3[3] = {0};
    v3[0] = material_->v3[0];
    v3[1] = material_->v3[1];
    v3[2] = material_->v3[2];

    // unsigned int i = 0;
    // for (; i < 3; i++)
    // {
    //     if (v3[i])
    //         break;
    // }
    // double v1[3] = {0}; // unit vector normal to v3, right-hand system
    // double v2[3] = {0}; // unit vector normal to v3, right-hand system
    // v1[(i+1)%3] = 1.0;
    // v2[(i+2)%3] = 1.0;

    double v1[3] = {0};
    double e1[3] = {1.0, 0.0, 0.0};
    for (unsigned int i = 0; i < 3; i++)
        v1[i] = e1[(i+1)%3] * v3[(i+2)%3] - e1[(i+2)%3] * v3[(i+1)%3];

    double v2[3] = {0};
    for (unsigned int i = 0; i < 3; i++)
        v2[i] = v3[(i+1)%3] * v1[(i+2)%3] - v3[(i+2)%3] * v1[(i+1)%3];

    double V3[3] = {0}; // the vetor pointing from the bottom nodes to the top nodes
    V3[0] = t * v3[0];
    V3[1] = t * v3[1];
    V3[2] = t * v3[2];

    double N_I[8] = {
        -0.25 * (1.0 - xi) * (1.0 - eta) * (xi + eta + 1.0),
        0.25 * (1.0 + xi) * (1.0 - eta) * (xi - eta - 1.0),
        0.25 * (1.0 + xi) * (1.0 + eta) * (xi + eta - 1.0),
        -0.25 * (1.0 - xi) * (1.0 + eta) * (xi - eta + 1.0),
        0.5 * (1.0 - xi * xi) * (1.0 - eta),
		0.5 * (1.0 + xi) * (1.0 - eta * eta),
		0.5 * (1.0 - xi * xi) * (1.0 + eta),
		0.5 * (1.0 - xi) * (1.0 - eta * eta),
    };

    // derivative of shape function NI in parent coordinates, relevant to xi and eta
	double DNDxi[8] = {0};
	DNDxi[0] = 0.25 * (2.0 * xi + eta) * (1.0 - eta);
	DNDxi[1] = 0.25 * (2.0 * xi - eta) * (1.0 - eta);
	DNDxi[2] = 0.25 * (2.0 * xi + eta) * (1.0 + eta);
	DNDxi[3] = 0.25 * (2.0 * xi - eta) * (1.0 + eta);
	DNDxi[4] = -xi * (1.0 - eta);
	DNDxi[5] = 0.5 * (1.0 - eta * eta);
	DNDxi[6] = -xi * (1.0 + eta);
	DNDxi[7] = -0.5 * (1.0 - eta * eta);
	
	double DNDeta[8] = {0};
	DNDeta[0] = 0.25 * (xi + 2.0 * eta) * (1.0 - xi);
	DNDeta[1] = 0.25 * (-xi + 2.0 * eta) * (1.0 + xi);
	DNDeta[2] = 0.25 * (xi + 2.0 * eta) * (1.0 + xi);
	DNDeta[3] = 0.25 * (-xi + 2.0 * eta) * (1.0 - xi);
	DNDeta[4] = -0.5 * (1.0 - xi * xi);
	DNDeta[5] = (1.0 + xi) * (-eta);
	DNDeta[6] = 0.5 * (1.0 - xi * xi);
	DNDeta[7] = (1.0 - xi) * (-eta);
	
    // compute Jacobian matrix J, x(xi, eta, zeta), y(xi, eta, zeta), z(xi, eta, zeta)
    double DXDxi[3] = {0};
    double DXDeta[3] = {0};
    double DXDzeta[3] = {0};

    // x, y and z coordinates
    double C[8][3];
	for (unsigned int i = 0; i < 8; i++)
	{
		C[i][0] = nodes_[i]->XYZ[0];
        C[i][1] = nodes_[i]->XYZ[1];
        C[i][2] = nodes_[i]->XYZ[2];
    };   
    
    for (unsigned int i = 0; i < 8; i++)
    {
        DXDxi[0] += DNDxi[i] * (C[i][0] + 0.5 * zeta * V3[0]);
		DXDxi[1] += DNDxi[i] * (C[i][1] + 0.5 * zeta * V3[1]);
		DXDxi[2] += DNDxi[i] * (C[i][2] + 0.5 * zeta * V3[2]);
		DXDeta[0] += DNDeta[i] * (C[i][0] + 0.5 * zeta * V3[0]);
		DXDeta[1] += DNDeta[i] * (C[i][1] + 0.5 * zeta * V3[1]);
		DXDeta[2] += DNDeta[i] * (C[i][2] + 0.5 * zeta * V3[2]);
		DXDzeta[0] += N_I[i] * 0.5 * V3[0];
		DXDzeta[1] += N_I[i] * 0.5 * V3[1];
		DXDzeta[2] += N_I[i] * 0.5 * V3[2];
    }

    double J[3][3] = {
        {DXDxi[0], DXDxi[1], DXDxi[2]},
        {DXDeta[0], DXDeta[1], DXDeta[2]},
        {DXDzeta[0], DXDzeta[1], DXDzeta[2]}
    };

    *det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
         - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) 
         + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    // Calculate adjoint matrix of J
    double adj[3][3];
    adj[0][0] = +(J[1][1] * J[2][2] - J[1][2] * J[2][1]);
    adj[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]);
    adj[0][2] = +(J[0][1] * J[1][2] - J[0][2] * J[1][1]);

    adj[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]);
    adj[1][1] = +(J[0][0] * J[2][2] - J[0][2] * J[2][0]);
    adj[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]);

    adj[2][0] = +(J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    adj[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]);
    adj[2][2] = +(J[0][0] * J[1][1] - J[0][1] * J[1][0]);
    // Calculate inverse matrix of J
    double invJ[3][3];
    for (int i = 0; i < 3; i++) 
        for (int j = 0; j < 3; j++) 
            invJ[i][j] = adj[i][j] / *det;
    
    double J_9[9][9] = {0}; // block diagnoal matrix with 3 invJ
    for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
            for (unsigned int k = 0; k < 3; k++)
                J_9[3*i+j][3*i+k] = invJ[j][k];

    // derivative of shape function matrix N in parent coordinates xi, eta and zeta
    double DuDp[9][40] = {0};   // p: parent coordinates
    for (unsigned int i = 0; i < 8; i++)
    {   
        DuDp[0][5*i] = DNDxi[i];
        DuDp[0][5*i+3] = 0.5 * zeta * t * DNDxi[i] * v1[0];
        DuDp[0][5*i+4] = -0.5 * zeta * t * DNDxi[i] * v2[0];

        DuDp[1][5*i] = DNDeta[i];
        DuDp[1][5*i+3] = 0.5 * zeta * t * DNDeta[i] * v1[0];
        DuDp[1][5*i+4] = -0.5 * zeta * t * DNDeta[i] * v2[0];

        DuDp[2][5*i+3] = 0.5 * t * N_I[i] * v1[0];
        DuDp[2][5*i+4] = -0.5 * t * N_I[i] * v2[0];
        
        DuDp[3][5*i+1] = DNDxi[i];
        DuDp[3][5*i+3] = 0.5 * zeta * t * DNDxi[i] * v1[1];
        DuDp[3][5*i+4] = -0.5 * zeta * t * DNDxi[i] * v2[1];

        DuDp[4][5*i+1] = DNDeta[i];
        DuDp[4][5*i+3] = 0.5 * zeta * t * DNDeta[i] * v1[1];
        DuDp[4][5*i+4] = -0.5 * zeta * t * DNDeta[i] * v2[1];

        DuDp[5][5*i+3] = 0.5 * t * N_I[i] * v1[1];
        DuDp[5][5*i+4] = -0.5 * t * N_I[i] * v2[1];

        DuDp[6][5*i+2] = DNDxi[i];
        DuDp[6][5*i+3] = 0.5 * zeta * t * DNDxi[i] * v1[2];
        DuDp[6][5*i+4] = -0.5 * zeta * t * DNDxi[i] * v2[2];

        DuDp[7][5*i+2] = DNDeta[i];
        DuDp[7][5*i+3] = 0.5 * zeta * t * DNDeta[i] * v1[2];
        DuDp[7][5*i+4] = -0.5 * zeta * t * DNDeta[i] * v2[2];

        DuDp[8][5*i+3] = 0.5 * t * N_I[i] * v1[2];
        DuDp[8][5*i+4] = -0.5 * t * N_I[i] * v2[2];
    }

    // Calculate derivative of shape function matrix N in global coordinates
    double DuDg[9][40] = {0};   // g: global coordinates
    for (unsigned int i = 0; i < 9; i++)
        for (unsigned int j = 0; j < 40; j++)
            for (unsigned int k = 0; k < 9; k++)
                DuDg[i][j] += J_9[i][k] * DuDp[k][j];

    // transformation matrix H between strain-displacement matrix B and displacement derivative matrix
    double H[6][9] = {0};
    H[0][0] = 1;    // normal strain, epsilon_xx

    H[1][4] = 1;    // normal strain, epsilon_yy

    H[2][8] = 1;    // normal strain, epsilon_zz

    H[3][1] = 1;    // shear strain, epsilon_xy
    H[3][3] = 1;
    
    H[4][5] = 1;    // shear strain, epsilon_yz
    H[4][7] = 1;

    H[5][2] = 1;    // shear strain, epsilon_zx
    H[5][6] = 1;

    // Calculate strain-displacement matrix in global coordinates
    double B1[6][40] = {0};
    for (unsigned int i = 0; i < 6; i++)
        for (unsigned int j = 0; j < 40; j++)
            for (unsigned int k = 0; k < 9; k++)
                B1[i][j] += H[i][k] * DuDg[k][j];

    // Calculate strain transformation matrix T between global and local coordinates
    double vector[3] = {0};   // cross product
    for (unsigned int i = 0; i < 3; i++)
        vector[i] = DXDxi[(i+1)%3] * DXDeta[(i+2)%3] - DXDxi[(i+2)%3] * DXDeta[(i+1)%3];
    double vector_norm2 = vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
    double vector_norm = sqrt(vector_norm2);

	double v3_[3] = {0};
    v3_[0] = vector[0] / vector_norm;
    v3_[1] = vector[1] / vector_norm;
    v3_[2] = vector[2] / vector_norm;

    double v1_[3] = {0};
    // double e1[3] = {1,0,0};
    for (unsigned int i = 0; i < 3; i++)
        v1_[i] = e1[(i+1)%3] * v3_[(i+2)%3] - e1[(i+2)%3] * v3_[(i+1)%3];

    double v2_[3] = {0};
    for (unsigned int i = 0; i < 3; i++)
        v2_[i] = v3_[(i+1)%3] * v1_[(i+2)%3] - v3_[(i+2)%3] * v1_[(i+1)%3];

    double T[6][6] = {0};
    T[0][0] = v1_[0] * v1_[0];
    T[0][1] = v1_[1] * v1_[1];
    T[0][2] = v1_[2] * v1_[2];
    T[0][3] = v1_[0] * v1_[1];
    T[0][4] = v1_[1] * v1_[2];
    T[0][5] = v1_[2] * v1_[0];

    T[1][0] = v2_[0] * v2_[0];
    T[1][1] = v2_[1] * v2_[1];
    T[1][2] = v2_[2] * v2_[2];
    T[1][3] = v2_[0] * v2_[1];
    T[1][4] = v2_[1] * v2_[2];
    T[1][5] = v2_[2] * v2_[0];

    T[2][0] = v3_[0] * v3_[0];
    T[2][1] = v3_[1] * v3_[1];
    T[2][2] = v3_[2] * v3_[2];
    T[2][3] = v3_[0] * v3_[1];
    T[2][4] = v3_[1] * v3_[2];
    T[2][5] = v3_[2] * v3_[0];

    T[3][0] = 2.0 * v1_[0] * v2_[0];
    T[3][1] = 2.0 * v1_[1] * v2_[1];
    T[3][2] = 2.0 * v1_[2] * v2_[2];
    T[3][3] = v1_[0] * v2_[1] + v1_[1] * v2_[0];
    T[3][4] = v1_[1] * v2_[2] + v1_[2] * v2_[1];
    T[3][5] = v1_[2] * v2_[0] + v1_[0] * v2_[2];

    T[4][0] = 2.0 * v2_[0] * v3_[0];
    T[4][1] = 2.0 * v2_[1] * v3_[1];
    T[4][2] = 2.0 * v2_[2] * v3_[2];
    T[4][3] = v2_[0] * v3_[1] + v2_[1] * v3_[0];
    T[4][4] = v2_[1] * v3_[2] + v2_[2] * v3_[1];
    T[4][5] = v2_[2] * v3_[0] + v2_[0] * v3_[2];

    T[5][0] = 2.0 * v3_[0] * v1_[0];
    T[5][1] = 2.0 * v3_[1] * v1_[1];
    T[5][2] = 2.0 * v3_[2] * v1_[2];
    T[5][3] = v3_[0] * v1_[1] + v3_[1] * v1_[0];
    T[5][4] = v3_[1] * v1_[2] + v3_[2] * v1_[1];
    T[5][5] = v3_[2] * v1_[0] + v3_[0] * v1_[2];

    // Calculate strain-displacement matrix in local coordinates
    double B2[6][40] = {0};
    for (unsigned int i = 0; i < 6; i++)
        for (unsigned int j = 0; j < 40; j++)
            for (unsigned int k = 0; k < 6; k++)
                B2[i][j] += T[i][k] * B1[k][j];

    // noticing that normal strain in z_ direction is 0
    // B is in the local frame
    for (unsigned int i = 0; i < 5; i++){
        for (unsigned int j = 0; j < 40; j++){
            if (i < 2)
                B[i][j] = B2[i][j];
            else
                B[i][j] = B2[i+1][j];
        }   
    }        
}

// Calcualte detJ at parent coordinate xi, eta, zeta = 1
void CS8R5::ElementDetTop(double* det, double xi, double eta)
{
    // x and y coordinates
    double C[8][2];
	for (unsigned int i = 0; i < 8; i++)
	{
		C[i][0] = nodes_[i]->XYZ[0];
        C[i][1] = nodes_[i]->XYZ[1];
    };

    // derivative of shape function NI in parent coordinates, relevant to xi and eta
	double DNDp[2][8] = {0};   // p: parent coordinates
    // DNDxi
	DNDp[0][0] = 0.25 * (2.0 * xi + eta) * (1.0 - eta);
	DNDp[0][1] = 0.25 * (2.0 * xi - eta) * (1.0 - eta);
	DNDp[0][2] = 0.25 * (2.0 * xi + eta) * (1.0 + eta);
	DNDp[0][3] = 0.25 * (2.0 * xi - eta) * (1.0 + eta);
	DNDp[0][4] = -xi * (1.0 - eta);
	DNDp[0][5] = 0.5 * (1.0 - eta * eta);
	DNDp[0][6] = -xi * (1.0 + eta);
	DNDp[0][7] = -0.5 * (1.0 - eta * eta);
	// DNDeta
	DNDp[1][0] = 0.25 * (xi + 2.0 * eta) * (1.0 - xi);
	DNDp[1][1] = 0.25 * (-xi + 2.0 * eta) * (1.0 + xi);
	DNDp[1][2] = 0.25 * (xi + 2.0 * eta) * (1.0 + xi);
	DNDp[1][3] = 0.25 * (-xi + 2.0 * eta) * (1.0 - xi);
	DNDp[1][4] = -0.5 * (1.0 - xi * xi);
	DNDp[1][5] = (1.0 + xi) * (-eta);
	DNDp[1][6] = 0.5 * (1.0 - xi * xi);
	DNDp[1][7] = (1.0 - xi) * (-eta);

    double J[2][2] = {0};
    for (unsigned int i = 0; i < 2; i++)
        for (unsigned int j = 0; j < 2; j++)
            for (unsigned int k = 0; k < 8; k++)    
                J[i][j] += DNDp[i][k] * C[k][j];

    *det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

