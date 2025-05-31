/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "LoadCaseData.h"

#include <iomanip>
#include <iostream>

using namespace std;

CLoadCaseData :: ~CLoadCaseData()
{
	delete [] node;
	delete [] dof;
	delete [] coordinate;
	delete [] load;
}

void CLoadCaseData :: Allocate(unsigned int LL, unsigned int num) 
{
	switch (LL)
	{
	case 1:	// All concentrated loads in node points
		nloads = num;
		node = new unsigned int[nloads];
		dof = new unsigned int[nloads];
		coordinate = nullptr;	// In case 1, coordinate is useless
		load = new double[nloads];
		break;	
	case 2:	// All concentrated loads in inner point of CST element(Dirac delta function)
		nloads = num;
		node = nullptr;	// In case 2, node is useless
		dof = nullptr;  // In case 2, dof is useless
		coordinate = new double[2*nloads];	// In case 2, coordinate means x and y coordinates
		load = new double[2*nloads];
		break;
	case 3:	// All body forces of CST element
		nloads = num;
		node = new unsigned int[nloads];	// In case 3, node is element number(CST elements with body forces)
		dof = nullptr;  // In case 3, dof is useless
		coordinate = nullptr;	// In case 3, coordinate is useless
		load = new double[6*nloads];	// In case 3, every element body force must be input in order of element node number, b_x1 b_y1 b_x2 b_y2 b_x3 b_y3
		break;
	case 4:	// All surface forces of CST element
		nloads = num;
		node = new unsigned int[2*nloads];
		dof = new unsigned int[nloads];  // In case 4, dof is element number(CST elements with surface forces)
		coordinate = nullptr;	// In case 4, coordinate is useless
		load = new double[4*nloads];	// In case 4, every element surface force must be input in order of element node number, t_x1 t_y1 t_x2 t_y2
		break;
	case 5:	// All body forces of Q4 element
		nloads = num;
		node = new unsigned int[nloads];	// In case 5, node is element number(CST elements with body forces)
		dof = nullptr;  // In case 5, dof is useless
		coordinate = nullptr;	// In case 5, coordinate is useless
		load = new double[8*nloads];	// In case 5, every element body force must be input in order of element node number, b_x1 b_y1 b_x2 b_y2 b_x3 b_y3 b_x4 b_y4
		break;
	case 6:	// All surface forces of Q4 element
		nloads = num;
		node = new unsigned int[2*nloads];
		dof = new unsigned int[nloads];  // In case 6, dof is element number(CST elements with surface forces)
		coordinate = nullptr;	// In case 6, coordinate is useless
		load = new double[4*nloads];	// In case 6, every element surface force must be input in order of element node number, t_x1 t_y1 t_x2 t_y2
		break;
		// Newly Added Q8
	case 7:	// All body forces of Q8 element
		nloads = num;
		node = new unsigned int[nloads];	// In case 7, node is element number
		dof = nullptr;  // In case 7, dof is useless
		coordinate = nullptr;	// In case 7, coordinate is useless
		load = new double[16*nloads];	// In case 7, every element body force must be input in order of element node number, b_x1 b_y1 ...
		break;
	case 8:	// All surface forces of Q8 element
		nloads = num;
		node = new unsigned int[3*nloads];
		dof = new unsigned int[nloads];  // In case 8, dof is element number
		coordinate = nullptr;	// In case 8, coordinate is useless
		load = new double[6*nloads];	// In case 8, every element surface force must be input in order of element node number, t_x1 t_y1 ...
		break;
	case 9:
		break;
	case 10:
		break;
	case 11:	// All body forces of H8 element
		nloads = num;
		node = new unsigned int[nloads];	// In case 11, node is element number
		dof = nullptr;  // In case 11, dof is useless
		coordinate = nullptr;	// In case 11, coordinate is useless
		load = new double[24*nloads];	// In case 11, every element body force must be input in order of element node number, b_x1 b_y1 b_z1 ...
		break;
	case 12:	// All surface forces of H8 element
		nloads = num;
		node = new unsigned int[4*nloads];
		dof = new unsigned int[nloads];  // In case 12, dof is element number
		coordinate = nullptr;	// In case 12, coordinate is useless
		load = new double[12*nloads];	// In case 12, every element surface force must be input in order of element node number, t_x1 t_y1 t_z1 ...
		break;
	case 13:	// All surface forces of S8R5 element
		nloads = num;
		node = new unsigned int[nloads];	// In case 13, node is element number
		dof = nullptr;  // In case 13, dof is useless
		coordinate = nullptr;	// In case 13, coordinate is useless
		load = new double[24*nloads];	// In case 13, every element body force must be input in order of element node number, t_x1 t_y1 t_z1 ...
		break;
	case 14:	// All body forces of S8R5 element
		nloads = num;
		node = new unsigned int[nloads];	// In case 14, node is element number
		dof = nullptr;  // In case 14, dof is useless
		coordinate = nullptr;	// In case 14, coordinate is useless
		load = new double[24*nloads];	// In case 14, every element body force must be input in order of element node number, b_x1 b_y1 b_z1 ...
		break;
	case 15:	// All body forces Mindlin-Reissner Plate element, except body moment
		nloads = num;
		node = new unsigned int[nloads];	// In case 15, node is element number
		dof = nullptr;  // In case 15, dof is useless
		coordinate = nullptr;	// In case 15, coordinate is useless
		load = new double[4*nloads];	// In case 15, every element body force must be input in order of element node number, b_z1 b_z2 b_z3 ...
		break;
	case 16:	// All body forces basic Plate element
		nloads = num;
		node = new unsigned int[nloads];	// In case 16, node is element number
		dof = nullptr;  // In case 16, dof is useless
		coordinate = nullptr;	// In case 16, coordinate is useless
		load = new double[4*nloads];	// In case 16, every element body force must be input in order of element node number, b_z1 b_z2 b_z3 ...
		break;
	default:
		std::cerr << "LodaCase " << LL << " not available. See CLoadCaseData::Allocate." << std::endl;
		exit(5);
		break;
	}	
}; 

//	Read load case data from stream Input
bool CLoadCaseData :: Read(unsigned int LL, ifstream& Input)
{
//	Load case number (LL) and number of concentrated loads/ elements/ boundaries in this load case(NL)
	
	unsigned int NL;

	Input >> NL;

	LoadCaseType_ = LL;
	Allocate(LL, NL);

	switch (LL)
	{		
	case 1:	// All concentrated loads in node points
		for (unsigned int i = 0; i < NL; i++)
			Input >> node[i] >> dof[i] >> load[i];
		break;	
	case 2:	// All concentrated loads in inner point of CST element(Dirac delta function)
		for (unsigned int i = 0; i < NL; i++)
			Input >> coordinate[2*i] >> coordinate[2*i+1] >> load[2*i] >> load[2*i+1];
		break;
	case 3:	// All body forces of CST element
		for (unsigned int i = 0; i < NL; i++)
			Input >> node[i] >> load[6*i] >> load[6*i+1] >> load[6*i+2] >> load[6*i+3] >> load[6*i+4] >> load[6*i+5];
		break;
	case 4:	// All surface forces of CST element
		for (unsigned int i = 0; i < NL; i++)
			Input >> dof[i] >> node[2*i] >> node[2*i+1] >> load[4*i] >> load[4*i+1] >> load[4*i+2] >> load[4*i+3];
		break;
	case 5:	// All body forces of Q4 element
		for (unsigned int i = 0; i < NL; i++)
			Input >> node[i] >> load[8*i] >> load[8*i+1] >> load[8*i+2] >> load[8*i+3] >> load[8*i+4] >> load[8*i+5] >> load[8*i+6] >> load[8*i+7];
		break;
	case 6:	// All surface forces of Q4 element
		for (unsigned int i = 0; i < NL; i++)
			Input >> dof[i] >> node[2*i] >> node[2*i+1] >> load[4*i] >> load[4*i+1] >> load[4*i+2] >> load[4*i+3];
		break;
	case 7:	// All body forces of Q8 element
		for (unsigned int i = 0; i < NL; i++)
			Input >> node[i] >> load[16*i] >> load[16*i+1] >> load[16*i+2] >> load[16*i+3] >> load[16*i+4] >> load[16*i+5] >> load[16*i+6] >> load[16*i+7] >> load[16*i+8] >> load[16*i+9] >> load[16*i+10] >> load[16*i+11] >> load[16*i+12] >> load[16*i+13] >> load[16*i+14] >> load[16*i+15];
		break;
	case 8:	// All surface forces of Q8 element
		for (unsigned int i = 0; i < NL; i++)
			Input >> dof[i] >> node[3*i] >> node[3*i+1] >> node[3*i+2] >> load[6*i] >> load[6*i+1] >> load[6*i+2] >> load[6*i+3] >> load[6*i+4] >> load[6*i+5];
		break;
	case 9:
		break;
	case 10:
		break;
	case 11:	// All body forces of H8 element
		for (unsigned int i = 0; i < NL; i++)
			Input >> node[i] >> load[24*i] >> load[24*i+1] >> load[24*i+2] >> load[24*i+3] >> load[24*i+4] >> load[24*i+5] >> load[24*i+6] >> load[24*i+7] >> load[24*i+8] >> load[24*i+9] >> load[24*i+10] >> load[24*i+11] >> load[24*i+12] >> load[24*i+13] >> load[24*i+14] >> load[24*i+15] >> load[24*i+16] >> load[24*i+17] >> load[24*i+18] >> load[24*i+19] >> load[24*i+20] >> load[24*i+21] >> load[24*i+22] >> load[24*i+23];
		break;
	case 12:	// All surface forces of H8 element
		for (unsigned int i = 0; i < NL; i++)
			Input >> dof[i] >> node[4*i] >> node[4*i+1] >> node[4*i+2] >> node[4*i+3] >> load[12*i] >> load[12*i+1] >> load[12*i+2] >> load[12*i+3] >> load[12*i+4] >> load[12*i+5] >> load[12*i+6] >> load[12*i+7] >> load[12*i+8] >> load[12*i+9] >> load[12*i+10] >> load[12*i+11];
		break;
	case 13:	// All surface forces of S8R5 element
		for (unsigned int i = 0; i < NL; i++)
			Input >> node[i] >> load[24*i] >> load[24*i+1] >> load[24*i+2] >> load[24*i+3] >> load[24*i+4] >> load[24*i+5] >> load[24*i+6] >> load[24*i+7] >> load[24*i+8] >> load[24*i+9] >> load[24*i+10] >> load[24*i+11] >> load[24*i+12] >> load[24*i+13] >> load[24*i+14] >> load[24*i+15] >> load[24*i+16] >> load[24*i+17] >> load[24*i+18] >> load[24*i+19] >> load[24*i+20] >> load[24*i+21] >> load[24*i+22] >> load[24*i+23];
		break;
	case 14:	// All body forces of S8R5 element
		for (unsigned int i = 0; i < NL; i++)
			Input >> node[i] >> load[24*i] >> load[24*i+1] >> load[24*i+2] >> load[24*i+3] >> load[24*i+4] >> load[24*i+5] >> load[24*i+6] >> load[24*i+7] >> load[24*i+8] >> load[24*i+9] >> load[24*i+10] >> load[24*i+11] >> load[24*i+12] >> load[24*i+13] >> load[24*i+14] >> load[24*i+15] >> load[24*i+16] >> load[24*i+17] >> load[24*i+18] >> load[24*i+19] >> load[24*i+20] >> load[24*i+21] >> load[24*i+22] >> load[24*i+23];
		break;
	case 15:	// All body forces of Mindlin-Reissner Plate element
		for (unsigned int i = 0; i < NL; i++)
			Input >> node[i] >> load[4*i] >> load[4*i+1] >> load[4*i+2] >> load[4*i+3];
		break;
	case 16:	// All body forces of basic Plate element
		for (unsigned int i = 0; i < NL; i++)
			Input >> node[i] >> load[4*i] >> load[4*i+1] >> load[4*i+2] >> load[4*i+3];
		break;
	default:
		std::cerr << "LodaCase " << LL << " not available. See CLoadCaseData::Read." << std::endl;
		exit(5);
		break;
	}

	return true;
}

//	Write load case data to stream
void CLoadCaseData::Write(unsigned int lcase, COutputter& output)
{
	switch (lcase)
	{
	case 1:	// All concentrated loads in node points
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(7) << node[i] << setw(13) << dof[i] << setw(19) << load[i] << endl;
		break;	
	case 2:	// All concentrated loads in inner point of CST element(Dirac delta function)
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(13) << coordinate[2*i] << setw(13) << coordinate[2*i+1] << setw(19) << load[2*i] << setw(19) << load[2*i+1] << endl;	
		break;
	case 3:	// All body forces of CST element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << node[i] << setw(18) << load[6*i] << setw(18) << load[6*i+1] << setw(18) << load[6*i+2] << setw(18) << load[6*i+3] << setw(18) << load[6*i+4] << setw(18) << load[6*i+5] << endl;	
		break;
	case 4:	// All surface forces of CST element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << dof[i] << setw(13) << node[2*i] << setw(13) << node[2*i+1] << setw(18) << load[4*i] << setw(18) << load[4*i+1] << setw(18) << load[4*i+2] << setw(18) << load[4*i+3] << endl;	
		break;
	case 5:	// All body forces of Q4 element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << node[i] << setw(18) << load[8*i] << setw(18) << load[8*i+1] << setw(18) << load[8*i+2] << setw(18) << load[8*i+3] << setw(18) << load[8*i+4] << setw(18) << load[8*i+5] << setw(18) << load[8*i+6] << setw(18) << load[8*i+7] << endl;	
		break;
	case 6:	// All surface forces of Q4 element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << dof[i] << setw(13) << node[2*i] << setw(13) << node[2*i+1] << setw(18) << load[4*i] << setw(18) << load[4*i+1] << setw(18) << load[4*i+2] << setw(18) << load[4*i+3] << endl;	
		break;
	case 7:	// All body forces of Q8 element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << node[i] << setw(18) << load[16*i] << setw(18) << load[16*i+1] << setw(18) << load[16*i+2] << setw(18) << load[16*i+3] << setw(18) << load[16*i+4] << setw(18) << load[16*i+5] << setw(18) << load[16*i+6] << setw(18) << load[16*i+7] << setw(18) << load[16*i+8] << setw(18) << load[16*i+9] << setw(18) << load[16*i+10] << setw(18) << load[16*i+11] << setw(18) << load[16*i+12] << setw(18) << load[16*i+13] << setw(18) << load[16*i+14] << setw(18) << load[16*i+15] << endl;
		break;
	case 8:	// All surface forces of Q8 element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << dof[i] << setw(13) << node[3*i] << setw(13) << node[3*i+1] << setw(13) << node[3*i+2] << setw(18) << load[6*i] << setw(18) << load[6*i+1] << setw(18) << load[6*i+2] << setw(18) << load[6*i+3] << setw(18) << load[6*i+4] << setw(18) << load[6*i+5] << endl;
		break;
	case 9:
		break;
	case 10:
		break;
	case 11:	// All body forces of H8 element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << node[i] << setw(18) << load[24*i] << setw(18) << load[24*i+1] << setw(18) << load[24*i+2] << setw(18) << load[24*i+3] << setw(18) << load[24*i+4] << setw(18) << load[24*i+5] << setw(18) << load[24*i+6] << setw(18) << load[24*i+7] << setw(18) << load[24*i+8] << setw(18) << load[24*i+9] << setw(18) << load[24*i+10] << setw(18) << load[24*i+11] << setw(18) << load[24*i+12] << setw(18) << load[24*i+13] << setw(18) << load[24*i+14] << setw(18) << load[24*i+15] << setw(18) << load[24*i+16] << setw(18) << load[24*i+17] << setw(18) << load[24*i+18] << setw(18) << load[24*i+19] << setw(18) << load[24*i+20] << setw(18) << load[24*i+21] << setw(18) << load[24*i+22] << setw(18) << load[24*i+23] << endl;
		break;
	case 12:	// All surface forces of H8 element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << dof[i] << setw(13) << node[4*i] << setw(13) << node[4*i+1] << setw(13) << node[4*i+2] << setw(13) << node[4*i+3] << setw(18) << load[12*i] << setw(18) << load[12*i+1] << setw(18) << load[12*i+2] << setw(18) << load[12*i+3] << setw(18) << load[12*i+4] << setw(18) << load[12*i+5] << setw(18) << load[12*i+6] << setw(18) << load[12*i+7] << setw(18) << load[12*i+8] << setw(18) << load[12*i+9] << setw(18) << load[12*i+10] << setw(18) << load[12*i+11] << endl;
		break;
	case 13:	// All surface forces of S8R5 element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << node[i] << setw(18) << load[24*i] << setw(18) << load[24*i+1] << setw(18) << load[24*i+2] << setw(18) << load[24*i+3] << setw(18) << load[24*i+4] << setw(18) << load[24*i+5] << setw(18) << load[24*i+6] << setw(18) << load[24*i+7] << setw(18) << load[24*i+8] << setw(18) << load[24*i+9] << setw(18) << load[24*i+10] << setw(18) << load[24*i+11] << setw(18) << load[24*i+12] << setw(18) << load[24*i+13] << setw(18) << load[24*i+14] << setw(18) << load[24*i+15] << setw(18) << load[24*i+16] << setw(18) << load[24*i+17] << setw(18) << load[24*i+18] << setw(18) << load[24*i+19] << setw(18) << load[24*i+20] << setw(18) << load[24*i+21] << setw(18) << load[24*i+22] << setw(18) << load[24*i+23] << endl;
		break;
	case 14:	// All body forces of S8R5 element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << node[i] << setw(18) << load[24*i] << setw(18) << load[24*i+1] << setw(18) << load[24*i+2] << setw(18) << load[24*i+3] << setw(18) << load[24*i+4] << setw(18) << load[24*i+5] << setw(18) << load[24*i+6] << setw(18) << load[24*i+7] << setw(18) << load[24*i+8] << setw(18) << load[24*i+9] << setw(18) << load[24*i+10] << setw(18) << load[24*i+11] << setw(18) << load[24*i+12] << setw(18) << load[24*i+13] << setw(18) << load[24*i+14] << setw(18) << load[24*i+15] << setw(18) << load[24*i+16] << setw(18) << load[24*i+17] << setw(18) << load[24*i+18] << setw(18) << load[24*i+19] << setw(18) << load[24*i+20] << setw(18) << load[24*i+21] << setw(18) << load[24*i+22] << setw(18) << load[24*i+23] << endl;
		break;
	case 15:	// All body forces of Mindlin-Reissner Plate element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << node[i] << setw(18) << load[4*i] << setw(18) << load[4*i+1] << setw(18) << load[4*i+2] << setw(18) << load[4*i+3] << endl;	
		break;
	case 16:	// All body forces of basic Plate element
		for (unsigned int i = 0; i < nloads; i++)
			output << setw(4) << node[i] << setw(18) << load[4*i] << setw(18) << load[4*i+1] << setw(18) << load[4*i+2] << setw(18) << load[4*i+3] << endl;	
		break;
	default:
		std::cerr << "LodaCase " << lcase << " not available. See CLoadCaseData::Write." << std::endl;
		exit(5);
		break;
	}	
}
