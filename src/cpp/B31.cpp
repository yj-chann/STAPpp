/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "B31.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CB31::CB31()
{
	NEN_ = 2;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 12;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CB31::~CB31()
{
}

//	Read element data from stream Input
bool CB31::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// first node number, second node number

	Input >> N1 >> N2 >> MSet;
    ElementMaterial_ = dynamic_cast<CB31Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CB31::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
//  Use selective reduced integration
void CB31::ElementStiffness(double* Matrix)  // Local stiffness matrix can be derived analytically, then use DCM matrix
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate element stiffness matrix

	CB31Material* material_ = dynamic_cast<CB31Material*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E;
    double G = material_->G;
    double Area = material_->Area;
	double Iy = material_->Iy;
    double Iz = material_->Iz;
    double J = material_->J;
    double k = material_->k;
    unsigned int direction = material_->direction;

//	Calculate beam length
	double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[6];	//  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[2] * DX[2];
	DX2[3] = DX[0] * DX[1];
	DX2[4] = DX[1] * DX[2];
	DX2[5] = DX[0] * DX[2];

	double L2 = DX2[0] + DX2[1] + DX2[2];
	double L = sqrt(L2);

    double Ke[12][12] = {0};

    switch (direction)
    {
        case 1: // Beam in the YOZ plane
        {
            Ke[1][1] = E * Area / L;
            Ke[1][7] = -E * Area / L;

            Ke[2][2] = k * G * Area / L;
            Ke[2][6] = k * G * Area / 2.0;
            Ke[2][8] = -k * G * Area / L;
            Ke[2][0] = k * G * Area / 2.0;

            Ke[3][3] = k * G * Area / L;
            Ke[3][5] = -k * G * Area / 2.0;
            Ke[3][9] = -k * G * Area / L;
            Ke[3][11] = -k * G * Area / 2.0;

            Ke[4][4] = G * J / L;
            Ke[4][10] = -G * J / L;

            Ke[5][3] = -k * G * Area / 2.0;
            Ke[5][5] = k * G * Area * L / 4.0 + E * Iy / L;
            Ke[5][9] = k * G * Area / 2.0;
            Ke[5][11] = k * G * Area * L / 4.0 - E * Iy / L;

            Ke[6][2] = k * G * Area / 2.0;
            Ke[6][6] = k * G * Area * L / 4.0 + E * Iz / L;
            Ke[6][8] = -k * G * Area / 2.0;
            Ke[6][0] = k * G * Area * L / 4.0 - E * Iz / L;

            Ke[7][1] = -E * Area / L;
            Ke[7][7] = E * Area / L;

            Ke[8][2] = -k * G * Area / L;
            Ke[8][6] = -k * G * Area / 2.0;
            Ke[8][8] = k * G * Area / L;
            Ke[8][0] = -k * G * Area / 2.0;

            Ke[9][3] = -k * G * Area / L;
            Ke[9][5] = k * G * Area / 2.0;
            Ke[9][9] = k * G * Area / L;
            Ke[9][11] = k * G * Area / 2.0;

            Ke[10][4] = -G * J / L;
            Ke[10][10] = G * J / L;

            Ke[11][3] = -k * G * Area / 2.0;
            Ke[11][5] = k * G * Area * L / 4.0 - E * Iy / L;
            Ke[11][9] = k * G * Area / 2.0;
            Ke[11][11] = k * G * Area * L / 4.0 + E * Iy / L;

            Ke[0][2] = k * G * Area / 2.0;
            Ke[0][6] = k * G * Area * L / 4.0 - E * Iz / L;
            Ke[0][8] = -k * G * Area / 2.0;
            Ke[0][0] = k * G * Area * L / 4.0 + E * Iz / L;
            break;
        }

        case 2: // Beam in the ZOX plane
        {
            Ke[2][2] = E * Area / L;
            Ke[2][8] = -E * Area / L;

            Ke[3][3] = k * G * Area / L;
            Ke[3][7] = k * G * Area / 2.0;
            Ke[3][9] = -k * G * Area / L;
            Ke[3][1] = k * G * Area / 2.0;

            Ke[4][4] = k * G * Area / L;
            Ke[4][6] = -k * G * Area / 2.0;
            Ke[4][10] = -k * G * Area / L;
            Ke[4][0] = -k * G * Area / 2.0;

            Ke[5][5] = G * J / L;
            Ke[5][11] = -G * J / L;

            Ke[6][4] = -k * G * Area / 2.0;
            Ke[6][6] = k * G * Area * L / 4.0 + E * Iy / L;
            Ke[6][10] = k * G * Area / 2.0;
            Ke[6][0] = k * G * Area * L / 4.0 - E * Iy / L;

            Ke[7][3] = k * G * Area / 2.0;
            Ke[7][7] = k * G * Area * L / 4.0 + E * Iz / L;
            Ke[7][9] = -k * G * Area / 2.0;
            Ke[7][1] = k * G * Area * L / 4.0 - E * Iz / L;

            Ke[8][2] = -E * Area / L;
            Ke[8][8] = E * Area / L;

            Ke[9][3] = -k * G * Area / L;
            Ke[9][7] = -k * G * Area / 2.0;
            Ke[9][9] = k * G * Area / L;
            Ke[9][1] = -k * G * Area / 2.0;

            Ke[10][4] = -k * G * Area / L;
            Ke[10][6] = k * G * Area / 2.0;
            Ke[10][10] = k * G * Area / L;
            Ke[10][0] = k * G * Area / 2.0;

            Ke[11][5] = -G * J / L;
            Ke[11][11] = G * J / L;

            Ke[0][4] = -k * G * Area / 2.0;
            Ke[0][6] = k * G * Area * L / 4.0 - E * Iy / L;
            Ke[0][10] = k * G * Area / 2.0;
            Ke[0][0] = k * G * Area * L / 4.0 + E * Iy / L;

            Ke[1][3] = k * G * Area / 2.0;
            Ke[1][7] = k * G * Area * L / 4.0 - E * Iz / L;
            Ke[1][9] = -k * G * Area / 2.0;
            Ke[1][1] = k * G * Area * L / 4.0 + E * Iz / L;
            break;
        }

        case 3: // Beam in the XOY plane
        {
            Ke[0][0] = E * Area / L;
            Ke[0][6] = -E * Area / L;

            Ke[1][1] = k * G * Area / L;
            Ke[1][5] = k * G * Area / 2.0;
            Ke[1][7] = -k * G * Area / L;
            Ke[1][11] = k * G * Area / 2.0;

            Ke[2][2] = k * G * Area / L;
            Ke[2][4] = -k * G * Area / 2.0;
            Ke[2][8] = -k * G * Area / L;
            Ke[2][10] = -k * G * Area / 2.0;

            Ke[3][3] = G * J / L;
            Ke[3][9] = -G * J / L;

            Ke[4][2] = -k * G * Area / 2.0;
            Ke[4][4] = k * G * Area * L / 4.0 + E * Iy / L;
            Ke[4][8] = k * G * Area / 2.0;
            Ke[4][10] = k * G * Area * L / 4.0 - E * Iy / L;

            Ke[5][1] = k * G * Area / 2.0;
            Ke[5][5] = k * G * Area * L / 4.0 + E * Iz / L;
            Ke[5][7] = -k * G * Area / 2.0;
            Ke[5][11] = k * G * Area * L / 4.0 - E * Iz / L;

            Ke[6][0] = -E * Area / L;
            Ke[6][6] = E * Area / L;

            Ke[7][1] = -k * G * Area / L;
            Ke[7][5] = -k * G * Area / 2.0;
            Ke[7][7] = k * G * Area / L;
            Ke[7][11] = -k * G * Area / 2.0;

            Ke[8][2] = -k * G * Area / L;
            Ke[8][4] = k * G * Area / 2.0;
            Ke[8][8] = k * G * Area / L;
            Ke[8][10] = k * G * Area / 2.0;

            Ke[9][3] = -G * J / L;
            Ke[9][9] = G * J / L;

            Ke[10][2] = -k * G * Area / 2.0;
            Ke[10][4] = k * G * Area * L / 4.0 - E * Iy / L;
            Ke[10][8] = k * G * Area / 2.0;
            Ke[10][10] = k * G * Area * L / 4.0 + E * Iy / L;

            Ke[11][1] = k * G * Area / 2.0;
            Ke[11][5] = k * G * Area * L / 4.0 - E * Iz / L;
            Ke[11][7] = -k * G * Area / 2.0;
            Ke[11][11] = k * G * Area * L / 4.0 + E * Iz / L;
            break;
        }
    }

// Calculate transformation Matrix
    double Te[12][12] = {0};

    switch (direction)
    {
        case 1: // Beam in the YOZ plane
        {
            // Calculate direction cosine matrix
            //! In this case, since the x-axis of the local coordinate system of the B31 element coincides with the x-axis of the global coordinate system,
            //! the analysis is simplified. Otherwise, we would need to determine the local y-axis and z-axis based on the cross-section shape, 
            //! which is a time-consuming task. It should be noted that, according to the inp file, the cross-section in this problem is of a box-type (rectangular hollow) shape.
                double lypy = DX[1] / L;    // Cosine of angle between beam and y-axis
                double lypz = DX[2] / L;    // Cosine of angle between beam and z-axis
            
            Te[0][0] = 1.0;
            Te[1][1] = lypy;
            Te[1][2] = lypz;
            Te[2][1] = -lypz;
            Te[2][2] = lypy;

            Te[3][3] = 1.0;
            Te[4][4] = lypy;
            Te[4][5] = lypz;
            Te[5][4] = -lypz;
            Te[5][5] = lypy;

            Te[6][6] = 1.0;
            Te[7][7] = lypy;
            Te[7][8] = lypz;
            Te[8][7] = -lypz;
            Te[8][8] = lypy;

            Te[9][9] = 1.0;
            Te[10][10] = lypy;
            Te[10][11] = lypz;
            Te[11][10] = -lypz;
            Te[11][11] = lypy;
            break;
        }

        case 2: // Beam in the ZOX plane
        {
            // Calculate direction cosine matrix
            //! In this case, since the y-axis of the local coordinate system of the B31 element coincides with the y-axis of the global coordinate system,
            //! the analysis is simplified. Otherwise, we would need to determine the local z-axis and x-axis based on the cross-section shape, 
            //! which is a time-consuming task. It should be noted that, according to the inp file, the cross-section in this problem is of a box-type (rectangular hollow) shape.
                double lzpz = DX[2] / L;    // Cosine of angle between beam and z-axis
                double lzpx = DX[0] / L;    // Cosine of angle between beam and x-axis
            
            Te[0][0] = lzpz;
            Te[0][2] = -lzpx;
            Te[1][1] = 1.0;
            Te[2][0] = lzpx;
            Te[2][2] = lzpz;

            Te[3][3] = lzpz;
            Te[3][5] = -lzpx;
            Te[4][4] = 1.0;
            Te[5][3] = lzpx;
            Te[5][5] = lzpz;

            Te[6][6] = lzpz;
            Te[6][8] = -lzpx;
            Te[7][7] = 1.0;
            Te[8][6] = lzpx;
            Te[8][8] = lzpz;

            Te[9][9] = lzpz;
            Te[9][11] = -lzpx;
            Te[10][10] = 1.0;
            Te[11][9] = lzpx;
            Te[11][11] = lzpz;
            break;
        }

        case 3: // Beam in the XOY plane
        {
            // Calculate direction cosine matrix
            //! In this case, since the z-axis of the local coordinate system of the B31 element coincides with the z-axis of the global coordinate system,
            //! the analysis is simplified. Otherwise, we would need to determine the local x-axis and y-axis based on the cross-section shape, 
            //! which is a time-consuming task. It should be noted that, according to the inp file, the cross-section in this problem is of a box-type (rectangular hollow) shape.
                double lxpx = DX[0] / L;    // Cosine of angle between beam and x-axis
                double lxpy = DX[1] / L;    // Cosine of angle between beam and y-axis

            Te[0][0] = lxpx;
            Te[0][1] = lxpy;
            Te[1][0] = -lxpy;
            Te[1][1] = lxpx;
            Te[2][2] = 1.0;

            Te[3][3] = lxpx;
            Te[3][4] = lxpy;
            Te[4][3] = -lxpy;
            Te[4][4] = lxpx;
            Te[5][5] = 1.0;

            Te[6][6] = lxpx;
            Te[6][7] = lxpy;
            Te[7][6] = -lxpy;
            Te[7][7] = lxpx;
            Te[8][8] = 1.0;

            Te[9][9] = lxpx;
            Te[9][10] = lxpy;
            Te[10][9] = -lxpy;
            Te[10][10] = lxpx;
            Te[11][11] = 1.0;
            break;
        }
    }

    double K[12][12] = {0};
    for (int i = 0; i < 12; i++)
        for (int j = 0; j < 12; j++)
            for (int k = 0; k < 12; k++)
                for (int l = 0; l < 12; l++)
                    K[i][j] += Te[k][i] * Ke[k][l] * Te[l][j];

	int index = 0;
    for (int j = 0; j < 12; j++)	// Column-major order
	{         
        for (int i = j; i >= 0; i--) {    // Upper triangle (i ≤ j)
            Matrix[index++] = K[i][j];   // Pack into one dimensional element stiffness matrix
        }
    }
}

//	Calculate element stress, for Beam element, they are moment and shear force in Gausspoints
void CB31::ElementStress(double* stress, double* Displacement)
{	    
	
}

//	Calculate element non-homogeneous essential boundary conditions
void CB31::ElementNonHomo(double* Matrix, double* NonForce)
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
    for (unsigned int i = 0; i < 2; i++)
    {
        d[6*i] = nodes_[i]->BC[0];
        d[6*i+1] = nodes_[i]->BC[1];
		d[6*i+2] = nodes_[i]->BC[2];
        d[6*i+3] = nodes_[i]->BC[3];
        d[6*i+4] = nodes_[i]->BC[4];
		d[6*i+5] = nodes_[i]->BC[5];
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

//	Calculate element shape function matrix N at parent coordinate xi
void CB31::ElementShapeFunction(double* N, double xi)
{
    
}

//	Calculate derivative of element shape function matrix B at parent coordinate xi(m = 2, actually second derivative)
void CB31::ElementStrainFunction(double* B, double xi)
{   
    
}

//	Calculate second derivative of element shape function matrix S at parent coordinate xi
void CB31::ElementDerivativeStrainFunction(double *S, double xi)
{
    
}
