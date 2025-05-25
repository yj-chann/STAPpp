/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}

//	Read material data from stream Input
bool CCSTMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu >> t >> plane_strain;	// Young's modulus, Poisson's ratio, element thickness and Plane Strain or not

	return true;
}

//	Write material data to Stream
void CCSTMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << setw(16) << t << setw(8) << plane_strain << endl;
}

//	Read material data from stream Input
bool CQ4Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu >> t >> plane_strain;	// Young's modulus, Poisson's ratio, element thickness and Plane Strain or not

	return true;
}

//	Write material data to Stream
void CQ4Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << setw(16) << t << setw(8) << plane_strain << endl;
}

//	Read material data from stream Input
bool CQ8Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu >> t >> plane_strain;	// Young's modulus, Poisson's ratio, element thickness and Plane Strain or not

	return true;
}

//	Write material data to Stream
void CQ8Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << setw(16) << t << setw(8) << plane_strain << endl;
}

//	Read material data from stream Input
bool CB21EBMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area >> I;	// Young's modulus, Cross-section area, Moment of inertia

	return true;
}

//	Write material data to Stream
void CB21EBMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << setw(16) << I << endl;
}

//	Read material data from stream Input
bool CB31Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu >> Area >> Iy >> Iz >> J >> k;	// Young's modulus, Shear modulus, Cross-section area, Moment of inertia, Torsional Moment, Correction factor

	G = E / 2.0 / (1 + nu);

	return true;
}

//	Write material data to Stream
void CB31Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << G << setw(16) << Area << setw(16) << Iy << setw(16) << Iz << setw(16) << J << setw(16) << k << endl;
}