/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

using namespace std;

//! Class LoadData is used to store load data
class CLoadCaseData
{
public:

	unsigned int LoadCaseType_;	//!< Load case 
	unsigned int nloads;		//!< Number of concentrated loads in this load case
	unsigned int* node;			//!< Node number to which this load is applied
	unsigned int* dof;			//!< Degree of freedom number for this load component
	double* coordinate;			//!< X Y Z coordinates of concentrated loads in inner point
	double* load;				//!< Magnitude of load

public:

	CLoadCaseData() : nloads(0), node(NULL), dof(NULL), coordinate(NULL), load(NULL) {};
	~CLoadCaseData();

//!	Set nloads, and new array node, dof and load
	void Allocate(unsigned int LL, unsigned int num);

//!	Read load case data from stream Input
	bool Read(unsigned int LL, ifstream& Input);

//!	Write load case data to stream
	void Write(unsigned int lcase, COutputter& output);
};
