/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Material.h"

using namespace std;

//	Clear an array
template <class type> void clear(type* a, unsigned int N)
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;

	NUMEG = 0;
	EleGrpList = nullptr;

	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;

	NEQ = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;
}

//	Desconstructor
CDomain::~CDomain()
{
	delete[] NodeList;

	delete[] EleGrpList;

	delete[] NLOAD;
	delete[] LoadCases;

	delete[] Force;
	delete StiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::GetInstance()
{
	if (!_instance)
		_instance = new CDomain();

	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::GetInstance(OutFile);

	//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

	//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

	//	Read nodal point data
	if (ReadNodalPoints())
		Output->OutputNodeInfo();
	else
		return false;

	//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

	//	Read load data
	if (ReadLoadCases())
		Output->OutputLoadInfo();
	else
		return false;

	//	Read element data
	if (ReadElements())
		Output->OutputElementInfo();
	else
		return false;

	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{

	//	Read nodal point data lines
	NodeList = new CNode[NUMNP];

	//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
	{
		if (!NodeList[np].Read(Input))
			return false;

		if (NodeList[np].NodeNumber != np + 1)
		{
			cerr << "*** Error *** Nodes must be inputted in order !" << endl
				<< "   Expected node number : " << np + 1 << endl
				<< "   Provided node number : " << NodeList[np].NodeNumber << endl;

			return false;
		}
	}

	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof])
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
	//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

	//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
	{
		unsigned int LL;
		Input >> LL;

		// if (LL != lcase + 1)
		// {
		//     cerr << "*** Error *** Load case must be inputted in order !" << endl
		//     << "   Expected load case : " << lcase + 1 << endl
		//     << "   Provided load case : " << LL << endl;

		//     return false;
		// }

		LoadCases[lcase].Read(LL, Input);
	}

	return true;
}

// Read element data
bool CDomain::ReadElements()
{
	EleGrpList = new CElementGroup[NUMEG];

	//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
		if (!EleGrpList[EleGrp].Read(Input))
			return false;

	return true;
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
#ifdef _DEBUG_
	COutputter* Output = COutputter::GetInstance();
	*Output << setw(9) << "Ele = " << setw(22) << "Location Matrix" << endl;
#endif

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
	{
		CElementGroup& ElementGrp = EleGrpList[EleGrp];
		unsigned int NUME = ElementGrp.GetNUME();

		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
		{
			CElement& Element = ElementGrp[Ele];

			// Generate location matrix
			Element.GenerateLocationMatrix();

#ifdef _DEBUG_
			unsigned int* LocationMatrix = Element.GetLocationMatrix();

			*Output << setw(9) << Ele + 1;
			for (int i = 0; i < Element.GetND(); i++)
				*Output << setw(5) << LocationMatrix[i];
			*Output << endl;
#endif

			StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
		}
	}

	StiffnessMatrix->CalculateMaximumHalfBandwidth();

#ifdef _DEBUG_
	* Output << endl;
	Output->PrintColumnHeights();
#endif

}

//    Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//    and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
	//    Allocate for global force/displacement vector
	Force = new double[NEQ];

	//  Create the banded stiffness matrix
	StiffnessMatrix = new CSkylineMatrix<double>(NEQ);

	//    Calculate column heights
	CalculateColumnHeights();

	//    Calculate address of diagonal elements in banded matrix
	StiffnessMatrix->CalculateDiagnoalAddress();

	//    Allocate for banded global stiffness matrix
	StiffnessMatrix->Allocate();

	COutputter* Output = COutputter::GetInstance();
	Output->OutputTotalSystemData();
}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
	//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		CElementGroup& ElementGrp = EleGrpList[EleGrp];
		unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

		//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = ElementGrp[Ele];
			Element.ElementStiffness(Matrix);
			StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
		}

		delete[] Matrix;
		Matrix = nullptr;
	}

#ifdef _DEBUG_
	COutputter* Output = COutputter::GetInstance();
	Output->PrintStiffnessMatrix();
#endif

}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE)	// LoadCase only means the order of LoadCase, not Type!!!
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

	if (LoadCase == 1)
		clear(Force, NEQ);

	switch (LoadData->LoadCaseType_)
	{
	case 1:	//	Loop over for all concentrated loads in load case LoadCase	
	{
		for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
		{
			unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];

			if (dof) // The DOF is activated
				Force[dof - 1] += LoadData->load[lnum];
		}
		break;
	}
	case 2:	// All concentrated loads in inner point of CST element(Dirac delta function)
	{
		break;	// Given that concentrated forces are scarcely found in natural systems, 
		// this component is temporarily omitted and will be addressed subsequently.
	}
	case 3:	// All body forces of CST element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 2)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->node[LoadOrder];
			double f_Omega[6] = { 0 };
			double b_x1 = LoadData->load[6 * LoadOrder];
			double b_y1 = LoadData->load[6 * LoadOrder + 1];
			double b_x2 = LoadData->load[6 * LoadOrder + 2];
			double b_y2 = LoadData->load[6 * LoadOrder + 3];
			double b_x3 = LoadData->load[6 * LoadOrder + 4];
			double b_y3 = LoadData->load[6 * LoadOrder + 5];

			// Calculate area of the CST element
			CElement& Element = ElementGrp[Ele - 1];
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CCSTMaterial* material_ = dynamic_cast<CCSTMaterial*>(ElementMaterial);
			double t = material_->t;

			CNode** nodes = Element.GetNodes();

			double X[3] = {
				nodes[0]->XYZ[0],
				nodes[1]->XYZ[0],
				nodes[2]->XYZ[0]
			};
			double Y[3] = {
				nodes[0]->XYZ[1],
				nodes[2]->XYZ[1],
				nodes[2]->XYZ[1],
			};
			double A = 0.5 * (X[0] * Y[1] - X[1] * Y[0] + X[2] * Y[0] - X[0] * Y[2] + X[1] * Y[2] - X[2] * Y[1]);

			f_Omega[0] = (2 * b_x1 + b_x2 + b_x3) * A * t / 12.0;
			f_Omega[1] = (2 * b_y1 + b_y2 + b_y3) * A * t / 12.0;
			f_Omega[2] = (b_x1 + 2 * b_x2 + b_x3) * A * t / 12.0;
			f_Omega[3] = (b_y1 + 2 * b_y2 + b_y3) * A * t / 12.0;
			f_Omega[4] = (b_x1 + b_x2 + 2 * b_x3) * A * t / 12.0;
			f_Omega[5] = (b_y1 + b_y2 + 2 * b_y3) * A * t / 12.0;

			unsigned int dof[6] = {
				nodes[0]->bcode[0],
				nodes[0]->bcode[1],
				nodes[1]->bcode[0],
				nodes[1]->bcode[1],
				nodes[2]->bcode[0],
				nodes[2]->bcode[1],
			};
			for (int i = 0; i < 6; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Omega[i];
		}
		break;
	}
	case 4:	// All surface forces of CST element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 2)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->dof[LoadOrder];
			double f_Gamma[4] = { 0 };
			unsigned int ElementNode[2] = {
				LoadData->node[2 * LoadOrder],
				LoadData->node[2 * LoadOrder + 1],
			};
			double t_x1 = LoadData->load[4 * LoadOrder];
			double t_y1 = LoadData->load[4 * LoadOrder + 1];
			double t_x2 = LoadData->load[4 * LoadOrder + 2];
			double t_y2 = LoadData->load[4 * LoadOrder + 3];

			// Calculate length of the line with surface force
			CElement& Element = ElementGrp[Ele - 1];
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CCSTMaterial* material_ = dynamic_cast<CCSTMaterial*>(ElementMaterial);
			double t = material_->t;

			CNode** nodes = Element.GetNodes();

			double X[2] = {
				nodes[ElementNode[0] - 1]->XYZ[0],
				nodes[ElementNode[1] - 1]->XYZ[0],
			};
			double Y[2] = {
				nodes[ElementNode[0] - 1]->XYZ[1],
				nodes[ElementNode[1] - 1]->XYZ[1],
			};
			double L2 = (X[1] - X[0]) * (X[1] - X[0]) + (Y[1] - Y[0]) * (Y[1] - Y[0]);
			double L = sqrt(L2);

			f_Gamma[0] = (2 * t_x1 + t_x2) * L * t / 6.0;
			f_Gamma[1] = (2 * t_y1 + t_y2) * L * t / 6.0;
			f_Gamma[2] = (t_x1 + 2 * t_x2) * L * t / 6.0;
			f_Gamma[3] = (t_y1 + 2 * t_y2) * L * t / 6.0;

			unsigned int dof[4] = {
				nodes[ElementNode[0] - 1]->bcode[0],
				nodes[ElementNode[0] - 1]->bcode[1],
				nodes[ElementNode[1] - 1]->bcode[0],
				nodes[ElementNode[1] - 1]->bcode[1],
			};
			for (int i = 0; i < 4; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Gamma[i];
		}
		break;
	}
	case 5:	// All body forces of Q4 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 3)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->node[LoadOrder];
			double b[8] = {
				LoadData->load[8 * LoadOrder],
				LoadData->load[8 * LoadOrder + 1],
				LoadData->load[8 * LoadOrder + 2],
				LoadData->load[8 * LoadOrder + 3],
				LoadData->load[8 * LoadOrder + 4],
				LoadData->load[8 * LoadOrder + 5],
				LoadData->load[8 * LoadOrder + 6],
				LoadData->load[8 * LoadOrder + 7],
			};


			CElement& Element = ElementGrp[Ele - 1];
			CQ4* Element_ = dynamic_cast<CQ4*>(&Element);
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial);
			double t = material_->t;
			CNode** nodes = Element.GetNodes();

			double f_Omega[8] = { 0 };
			double GaussPoints[2] = {
				-sqrt(3) / 3.0,
				sqrt(3) / 3.0,
			};
			for (unsigned int i = 0; i < 2; i++)	// use full integration here, then is reduced integration
			{
				for (unsigned int j = 0; j < 2; j++)
				{
					double B[3][8] = { 0 };
					double det;
					Element_->ElementStrainFunction(B, &det, GaussPoints[i], GaussPoints[j]);
					double N[2][8] = { 0 };
					Element_->ElementShapeFunction(N, GaussPoints[i], GaussPoints[j]);
					for (int k = 0; k < 8; k++)
						for (int l = 0; l < 2; l++)
							for (int m = 0; m < 8; m++)
								f_Omega[k] += t * N[l][k] * N[l][m] * b[m] * det;  // Weight is 1
				}
			}

			unsigned int dof[8] = {
				nodes[0]->bcode[0],
				nodes[0]->bcode[1],
				nodes[1]->bcode[0],
				nodes[1]->bcode[1],
				nodes[2]->bcode[0],
				nodes[2]->bcode[1],
				nodes[3]->bcode[0],
				nodes[3]->bcode[1],
			};

			for (int i = 0; i < 8; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Omega[i];
		}
		break;
	}
	case 6:	// All surface forces of Q4 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 3)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->dof[LoadOrder];
			double f_Gamma[4] = { 0 };
			unsigned int ElementNode[2] = {
				LoadData->node[2 * LoadOrder],
				LoadData->node[2 * LoadOrder + 1],
			};
			double t_x1 = LoadData->load[4 * LoadOrder];
			double t_y1 = LoadData->load[4 * LoadOrder + 1];
			double t_x2 = LoadData->load[4 * LoadOrder + 2];
			double t_y2 = LoadData->load[4 * LoadOrder + 3];

			// Calculate length of the line with surface force
			CElement& Element = ElementGrp[Ele - 1];
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial);
			double t = material_->t;

			CNode** nodes = Element.GetNodes();

			double X[2] = {
				nodes[ElementNode[0] - 1]->XYZ[0],
				nodes[ElementNode[1] - 1]->XYZ[0],
			};
			double Y[2] = {
				nodes[ElementNode[0] - 1]->XYZ[1],
				nodes[ElementNode[1] - 1]->XYZ[1],
			};
			double L2 = (X[1] - X[0]) * (X[1] - X[0]) + (Y[1] - Y[0]) * (Y[1] - Y[0]);
			double L = sqrt(L2);

			f_Gamma[0] = (2 * t_x1 + t_x2) * L * t / 6.0;
			f_Gamma[1] = (2 * t_y1 + t_y2) * L * t / 6.0;
			f_Gamma[2] = (t_x1 + 2 * t_x2) * L * t / 6.0;
			f_Gamma[3] = (t_y1 + 2 * t_y2) * L * t / 6.0;

			unsigned int dof[4] = {
				nodes[ElementNode[0] - 1]->bcode[0],
				nodes[ElementNode[0] - 1]->bcode[1],
				nodes[ElementNode[1] - 1]->bcode[0],
				nodes[ElementNode[1] - 1]->bcode[1],
			};
			for (int i = 0; i < 4; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Gamma[i];
		}
		break;
	}
	case 7:	// All body forces of Q8 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 4)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->node[LoadOrder];
			double b[16] = {
				LoadData->load[16 * LoadOrder],
				LoadData->load[16 * LoadOrder + 1],
				LoadData->load[16 * LoadOrder + 2],
				LoadData->load[16 * LoadOrder + 3],
				LoadData->load[16 * LoadOrder + 4],
				LoadData->load[16 * LoadOrder + 5],
				LoadData->load[16 * LoadOrder + 6],
				LoadData->load[16 * LoadOrder + 7],
				LoadData->load[16 * LoadOrder + 8],
				LoadData->load[16 * LoadOrder + 9],
				LoadData->load[16 * LoadOrder + 10],
				LoadData->load[16 * LoadOrder + 11],
				LoadData->load[16 * LoadOrder + 12],
				LoadData->load[16 * LoadOrder + 13],
				LoadData->load[16 * LoadOrder + 14],
				LoadData->load[16 * LoadOrder + 15],

			};

			CElement& Element = ElementGrp[Ele - 1];
			CQ8* Element_ = dynamic_cast<CQ8*>(&Element);
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CQ8Material* material_ = dynamic_cast<CQ8Material*>(ElementMaterial);
			double t = material_->t;
			CNode** nodes = Element.GetNodes();

			double f_Omega[16] = { 0 };
			double GaussPoints[2] = {
				-sqrt(3) / 3.0,
				sqrt(3) / 3.0,
			};
			for (unsigned int i = 0; i < 2; i++)	//  reduced integration
			{
				for (unsigned int j = 0; j < 2; j++)
				{
					double B[3][16] = { 0 };
					double det;
					Element_->ElementStrainFunction(B, &det, GaussPoints[i], GaussPoints[j]);
					double N[2][16] = { 0 };
					Element_->ElementShapeFunction(N, GaussPoints[i], GaussPoints[j]);
					for (int k = 0; k < 16; k++)
						for (int l = 0; l < 2; l++)
							for (int m = 0; m < 16; m++)
								f_Omega[k] += t * N[l][k] * N[l][m] * b[m] * det;  // Weight is 1
				}
			}

			unsigned int dof[16] = {
				nodes[0]->bcode[0],
				nodes[0]->bcode[1],
				nodes[1]->bcode[0],
				nodes[1]->bcode[1],
				nodes[2]->bcode[0],
				nodes[2]->bcode[1],
				nodes[3]->bcode[0],
				nodes[3]->bcode[1],
				nodes[4]->bcode[0],
				nodes[4]->bcode[1],
				nodes[5]->bcode[0],
				nodes[5]->bcode[1],
				nodes[6]->bcode[0],
				nodes[6]->bcode[1],
				nodes[7]->bcode[0],
				nodes[7]->bcode[1],
			};

			for (int i = 0; i < 16; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Omega[i];
		}
		break;
	}
	case 8:	// All surface forces of Q8 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 4)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->dof[LoadOrder];
			double f_Gamma[6] = { 0 };
			unsigned int ElementNode[3] = {
				LoadData->node[3 * LoadOrder],
				LoadData->node[3 * LoadOrder + 1],
				LoadData->node[3 * LoadOrder + 2],
			};
			double t_x1 = LoadData->load[6 * LoadOrder];
			double t_y1 = LoadData->load[6 * LoadOrder + 1];
			double t_x2 = LoadData->load[6 * LoadOrder + 2];
			double t_y2 = LoadData->load[6 * LoadOrder + 3];
			double t_x3 = LoadData->load[6 * LoadOrder + 4];
			double t_y3 = LoadData->load[6 * LoadOrder + 5];


			CElement& Element = ElementGrp[Ele - 1];
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CQ8Material* material_ = dynamic_cast<CQ8Material*>(ElementMaterial);
			double t = material_->t;

			CNode** nodes = Element.GetNodes();

			double X[3] = {
				nodes[ElementNode[0] - 1]->XYZ[0],
				nodes[ElementNode[1] - 1]->XYZ[0],
				nodes[ElementNode[2] - 1]->XYZ[0],
			};
			double Y[3] = {
				nodes[ElementNode[0] - 1]->XYZ[1],
				nodes[ElementNode[1] - 1]->XYZ[1],
				nodes[ElementNode[2] - 1]->XYZ[1],
			};

			double GaussPoints[2] = {
				-sqrt(3) / 3.0,
				 sqrt(3) / 3.0,
			};
			// Scaling Factor
			double f12 = ((X[1] + X[0] - 2 * X[2]) * GaussPoints[0] + 0.5 * (X[1] - X[0])) * ((X[1] + X[0] - 2 * X[2]) * GaussPoints[0] + 0.5 * (X[1] - X[0])) + ((Y[1] + Y[0] - 2 * Y[2]) * GaussPoints[0] + 0.5 * (Y[1] - Y[0])) * ((Y[1] + Y[0] - 2 * Y[2]) * GaussPoints[0] + 0.5 * (Y[1] - Y[0]));
			double f1 = sqrt(f12);
			// Scaling Factor
			double f22 = ((X[1] + X[0] - 2 * X[2]) * GaussPoints[1] + 0.5 * (X[1] - X[0])) * ((X[1] + X[0] - 2 * X[2]) * GaussPoints[1] + 0.5 * (X[1] - X[0])) + ((Y[1] + Y[0] - 2 * Y[2]) * GaussPoints[1] + 0.5 * (Y[1] - Y[0])) * ((Y[1] + Y[0] - 2 * Y[2]) * GaussPoints[1] + 0.5 * (Y[1] - Y[0]));
			double f2 = sqrt(f22);

			// Nitxi Nityi
			double Nitxi1 = 0.5 * GaussPoints[0] * (GaussPoints[0] - 1) * t_x1 + (1 - GaussPoints[0] * GaussPoints[0]) * t_x3 + 0.5 * GaussPoints[0] * (GaussPoints[0] + 1) * t_x2;
			double Nitxi2 = 0.5 * GaussPoints[1] * (GaussPoints[1] - 1) * t_x1 + (1 - GaussPoints[1] * GaussPoints[1]) * t_x3 + 0.5 * GaussPoints[1] * (GaussPoints[1] + 1) * t_x2;
			double Nityi1 = 0.5 * GaussPoints[0] * (GaussPoints[0] - 1) * t_y1 + (1 - GaussPoints[0] * GaussPoints[0]) * t_y3 + 0.5 * GaussPoints[0] * (GaussPoints[0] + 1) * t_y2;
			double Nityi2 = 0.5 * GaussPoints[1] * (GaussPoints[1] - 1) * t_y1 + (1 - GaussPoints[1] * GaussPoints[1]) * t_y3 + 0.5 * GaussPoints[1] * (GaussPoints[1] + 1) * t_y2;
			f_Gamma[0] = 0.5 * GaussPoints[0] * (GaussPoints[0] - 1) * Nitxi1 * f1 + 0.5 * GaussPoints[1] * (GaussPoints[1] - 1) * Nitxi2 * f2;
			f_Gamma[1] = 0.5 * GaussPoints[0] * (GaussPoints[0] - 1) * Nityi1 * f1 + 0.5 * GaussPoints[1] * (GaussPoints[1] - 1) * Nityi2 * f2;
			f_Gamma[2] = 0.5 * GaussPoints[0] * (GaussPoints[0] + 1) * Nitxi1 * f1 + 0.5 * GaussPoints[1] * (GaussPoints[1] + 1) * Nitxi2 * f2;
			f_Gamma[3] = 0.5 * GaussPoints[0] * (GaussPoints[0] + 1) * Nityi1 * f1 + 0.5 * GaussPoints[1] * (GaussPoints[1] + 1) * Nityi2 * f2;
			f_Gamma[4] = (1 - GaussPoints[0] * GaussPoints[0]) * Nitxi1 * f1 + (1 - GaussPoints[1] * GaussPoints[1]) * Nitxi2 * f2;
			f_Gamma[5] = (1 - GaussPoints[0] * GaussPoints[0]) * Nityi1 * f1 + (1 - GaussPoints[1] * GaussPoints[1]) * Nityi2 * f2;

			unsigned int dof[6] = {
				nodes[ElementNode[0] - 1]->bcode[0],
				nodes[ElementNode[0] - 1]->bcode[1],
				nodes[ElementNode[1] - 1]->bcode[0],
				nodes[ElementNode[1] - 1]->bcode[1],
				nodes[ElementNode[2] - 1]->bcode[0],
				nodes[ElementNode[2] - 1]->bcode[1],
			};
			for (int i = 0; i < 6; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Gamma[i];
		}
		break;
	}


	case 11:	// All body forces of H8 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 7)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->node[LoadOrder];
			double b[24] = {
				LoadData->load[24 * LoadOrder],
				LoadData->load[24 * LoadOrder + 1],
				LoadData->load[24 * LoadOrder + 2],
				LoadData->load[24 * LoadOrder + 3],
				LoadData->load[24 * LoadOrder + 4],
				LoadData->load[24 * LoadOrder + 5],
				LoadData->load[24 * LoadOrder + 6],
				LoadData->load[24 * LoadOrder + 7],
				LoadData->load[24 * LoadOrder + 8],
				LoadData->load[24 * LoadOrder + 9],
				LoadData->load[24 * LoadOrder + 10],
				LoadData->load[24 * LoadOrder + 11],
				LoadData->load[24 * LoadOrder + 12],
				LoadData->load[24 * LoadOrder + 13],
				LoadData->load[24 * LoadOrder + 14],
				LoadData->load[24 * LoadOrder + 15],
				LoadData->load[24 * LoadOrder + 16],
				LoadData->load[24 * LoadOrder + 17],
				LoadData->load[24 * LoadOrder + 18],
				LoadData->load[24 * LoadOrder + 19],
				LoadData->load[24 * LoadOrder + 20],
				LoadData->load[24 * LoadOrder + 21],
				LoadData->load[24 * LoadOrder + 22],
				LoadData->load[24 * LoadOrder + 23],


			};

			CElement& Element = ElementGrp[Ele - 1];
			CH8* Element_ = dynamic_cast<CH8*>(&Element);
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CH8Material* material_ = dynamic_cast<CH8Material*>(ElementMaterial);
			CNode** nodes = Element.GetNodes();

			double f_Omega[24] = { 0 };
			double GaussPoints[2] = {
				-sqrt(3) / 3.0,
				sqrt(3) / 3.0,
			};
			for (unsigned int i = 0; i < 2; i++)	//  full integration
			{
				for (unsigned int j = 0; j < 2; j++)
				{
					for (unsigned int g = 0; g < 2; g++)
					{
						double B[6][24] = { 0 };
						double det;
						Element_->ElementStrainFunction(B, &det, GaussPoints[i], GaussPoints[j], GaussPoints[g]);
						double N[3][24] = { 0 };
						Element_->ElementShapeFunction(N, GaussPoints[i], GaussPoints[j], GaussPoints[g]);
						for (int k = 0; k < 24; k++)
							for (int l = 0; l < 3; l++)
								for (int m = 0; m < 24; m++)
									f_Omega[k] += N[l][k] * N[l][m] * b[m] * det;  // Weight is 1
					}
				}
			}

			unsigned int dof[24] = {
				nodes[0]->bcode[0],
				nodes[0]->bcode[1],
				nodes[0]->bcode[2],
				nodes[1]->bcode[0],
				nodes[1]->bcode[1],
				nodes[1]->bcode[2],
				nodes[2]->bcode[0],
				nodes[2]->bcode[1],
				nodes[2]->bcode[2],
				nodes[3]->bcode[0],
				nodes[3]->bcode[1],
				nodes[3]->bcode[2],
				nodes[4]->bcode[0],
				nodes[4]->bcode[1],
				nodes[4]->bcode[2],
				nodes[5]->bcode[0],
				nodes[5]->bcode[1],
				nodes[5]->bcode[2],
				nodes[6]->bcode[0],
				nodes[6]->bcode[1],
				nodes[6]->bcode[2],
				nodes[7]->bcode[0],
				nodes[7]->bcode[1],
				nodes[7]->bcode[2],
			};

			for (int i = 0; i < 24; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Omega[i];
		}
		break;
	}
	case 12:	// All surface forces of H8 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 7)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->dof[LoadOrder];
			double f_Gamma[12] = { 0 };
			unsigned int ElementNode[4] = {
				LoadData->node[4 * LoadOrder],
				LoadData->node[4 * LoadOrder + 1],
				LoadData->node[4 * LoadOrder + 2],
				LoadData->node[4 * LoadOrder + 3],
			};
			double t_x1 = LoadData->load[12 * LoadOrder];
			double t_y1 = LoadData->load[12 * LoadOrder + 1];
			double t_z1 = LoadData->load[12 * LoadOrder + 2];
			double t_x2 = LoadData->load[12 * LoadOrder + 3];
			double t_y2 = LoadData->load[12 * LoadOrder + 4];
			double t_z2 = LoadData->load[12 * LoadOrder + 5];
			double t_x3 = LoadData->load[12 * LoadOrder + 6];
			double t_y3 = LoadData->load[12 * LoadOrder + 7];
			double t_z3 = LoadData->load[12 * LoadOrder + 8];
			double t_x4 = LoadData->load[12 * LoadOrder + 9];
			double t_y4 = LoadData->load[12 * LoadOrder + 10];
			double t_z4 = LoadData->load[12 * LoadOrder + 11];


			CElement& Element = ElementGrp[Ele - 1];
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CH8Material* material_ = dynamic_cast<CH8Material*>(ElementMaterial);

			CNode** nodes = Element.GetNodes();

			double X[4] = {
				nodes[ElementNode[0] - 1]->XYZ[0],
				nodes[ElementNode[1] - 1]->XYZ[0],
				nodes[ElementNode[2] - 1]->XYZ[0],
				nodes[ElementNode[3] - 1]->XYZ[0],
			};
			double Y[4] = {
				nodes[ElementNode[0] - 1]->XYZ[1],
				nodes[ElementNode[1] - 1]->XYZ[1],
				nodes[ElementNode[2] - 1]->XYZ[1],
				nodes[ElementNode[3] - 1]->XYZ[1],
			};
			double Z[4] = {
				nodes[ElementNode[0] - 1]->XYZ[2],
				nodes[ElementNode[1] - 1]->XYZ[2],
				nodes[ElementNode[2] - 1]->XYZ[2],
				nodes[ElementNode[3] - 1]->XYZ[2],
			};

			double GaussPoints[2] = {
				-sqrt(3) / 3.0,
				 sqrt(3) / 3.0,
			};

			double x_21 = X[1] - X[0], y_21 = Y[1] - Y[0], z_21 = Z[1] - Z[0];
			double x_43 = X[3] - X[2], y_43 = Y[3] - Y[2], z_43 = Z[3] - Z[2];
			double x_41 = X[3] - X[0], y_41 = Y[3] - Y[0], z_41 = Z[3] - Z[0];
			double x_32 = X[2] - X[1], y_32 = Y[2] - Y[1], z_32 = Z[2] - Z[1];

			// Jij_kl
			double J11_11 = 2 * x_21 * (1 - GaussPoints[0]) - 2 * x_43 * (1 + GaussPoints[0]);
			double J12_11 = 2 * y_21 * (1 - GaussPoints[0]) - 2 * y_43 * (1 + GaussPoints[0]);
			double J13_11 = 2 * z_21 * (1 - GaussPoints[0]) - 2 * z_43 * (1 + GaussPoints[0]);
			double J21_11 = 2 * x_41 * (1 - GaussPoints[0]) + 2 * x_32 * (1 + GaussPoints[0]);
			double J22_11 = 2 * y_41 * (1 - GaussPoints[0]) + 2 * y_32 * (1 + GaussPoints[0]);
			double J23_11 = 2 * z_41 * (1 - GaussPoints[0]) + 2 * z_32 * (1 + GaussPoints[0]);

			double J11_12 = 2 * x_21 * (1 - GaussPoints[1]) - 2 * x_43 * (1 + GaussPoints[1]);
			double J12_12 = 2 * y_21 * (1 - GaussPoints[1]) - 2 * y_43 * (1 + GaussPoints[1]);
			double J13_12 = 2 * z_21 * (1 - GaussPoints[1]) - 2 * z_43 * (1 + GaussPoints[1]);
			double J21_12 = 2 * x_41 * (1 - GaussPoints[0]) + 2 * x_32 * (1 + GaussPoints[0]);
			double J22_12 = 2 * y_41 * (1 - GaussPoints[0]) + 2 * y_32 * (1 + GaussPoints[0]);
			double J23_12 = 2 * z_41 * (1 - GaussPoints[0]) + 2 * z_32 * (1 + GaussPoints[0]);

			double J11_21 = 2 * x_21 * (1 - GaussPoints[0]) - 2 * x_43 * (1 + GaussPoints[0]);
			double J12_21 = 2 * y_21 * (1 - GaussPoints[0]) - 2 * y_43 * (1 + GaussPoints[0]);
			double J13_21 = 2 * z_21 * (1 - GaussPoints[0]) - 2 * z_43 * (1 + GaussPoints[0]);
			double J21_21 = 2 * x_41 * (1 - GaussPoints[1]) + 2 * x_32 * (1 + GaussPoints[1]);
			double J22_21 = 2 * y_41 * (1 - GaussPoints[1]) + 2 * y_32 * (1 + GaussPoints[1]);
			double J23_21 = 2 * z_41 * (1 - GaussPoints[1]) + 2 * z_32 * (1 + GaussPoints[1]);

			double J11_22 = 2 * x_21 * (1 - GaussPoints[1]) - 2 * x_43 * (1 + GaussPoints[1]);
			double J12_22 = 2 * y_21 * (1 - GaussPoints[1]) - 2 * y_43 * (1 + GaussPoints[1]);
			double J13_22 = 2 * z_21 * (1 - GaussPoints[1]) - 2 * z_43 * (1 + GaussPoints[1]);
			double J21_22 = 2 * x_41 * (1 - GaussPoints[1]) + 2 * x_32 * (1 + GaussPoints[1]);
			double J22_22 = 2 * y_41 * (1 - GaussPoints[1]) + 2 * y_32 * (1 + GaussPoints[1]);
			double J23_22 = 2 * z_41 * (1 - GaussPoints[1]) + 2 * z_32 * (1 + GaussPoints[1]);


			double E_11 = J11_11 * J11_11 + J12_11 * J12_11 + J13_11 * J13_11;
			double F_11 = J21_11 * J21_11 + J22_11 * J22_11 + J23_11 * J23_11;
			double G_11 = J11_11 * J21_11 + J12_11 * J22_11 + J13_11 * J23_11;

			double E_12 = J11_12 * J11_12 + J12_12 * J12_12 + J13_12 * J13_12;
			double F_12 = J21_12 * J21_12 + J22_12 * J22_12 + J23_12 * J23_12;
			double G_12 = J11_12 * J21_12 + J12_12 * J22_12 + J13_12 * J23_12;

			double E_21 = J11_21 * J11_21 + J12_21 * J12_21 + J13_21 * J13_21;
			double F_21 = J21_21 * J21_21 + J22_21 * J22_21 + J23_21 * J23_21;
			double G_21 = J11_21 * J21_21 + J12_21 * J22_21 + J13_21 * J23_21;

			double E_22 = J11_22 * J11_22 + J12_22 * J12_22 + J13_22 * J13_22;
			double F_22 = J21_22 * J21_22 + J22_22 * J22_22 + J23_22 * J23_22;
			double G_22 = J11_22 * J21_22 + J12_22 * J22_22 + J13_22 * J23_22;


			// Scaling Factor
			double f11_2 = E_11 * F_11 - G_11 * G_11;
			double f11 = sqrt(f11_2);

			double f12_2 = E_12 * F_12 - G_12 * G_12;
			double f12 = sqrt(f12_2);

			double f21_2 = E_21 * F_21 - G_21 * G_21;
			double f21 = sqrt(f21_2);

			double f22_2 = E_22 * F_22 - G_22 * G_22;
			double f22 = sqrt(f22_2);


			// Ni
			double N1_11 = 0.25 * (1 - GaussPoints[0]) * (1 - GaussPoints[0]);
			double N1_12 = 0.25 * (1 - GaussPoints[0]) * (1 - GaussPoints[1]);
			double N1_21 = 0.25 * (1 - GaussPoints[1]) * (1 - GaussPoints[0]);
			double N1_22 = 0.25 * (1 - GaussPoints[1]) * (1 - GaussPoints[1]);

			double N2_11 = 0.25 * (1 + GaussPoints[0]) * (1 - GaussPoints[0]);
			double N2_12 = 0.25 * (1 + GaussPoints[0]) * (1 - GaussPoints[1]);
			double N2_21 = 0.25 * (1 + GaussPoints[1]) * (1 - GaussPoints[0]);
			double N2_22 = 0.25 * (1 + GaussPoints[1]) * (1 - GaussPoints[1]);


			double N3_11 = 0.25 * (1 + GaussPoints[0]) * (1 + GaussPoints[0]);
			double N3_12 = 0.25 * (1 + GaussPoints[0]) * (1 + GaussPoints[1]);
			double N3_21 = 0.25 * (1 + GaussPoints[1]) * (1 + GaussPoints[0]);
			double N3_22 = 0.25 * (1 + GaussPoints[1]) * (1 + GaussPoints[1]);


			double N4_11 = 0.25 * (1 - GaussPoints[0]) * (1 + GaussPoints[0]);
			double N4_12 = 0.25 * (1 - GaussPoints[0]) * (1 + GaussPoints[1]);
			double N4_21 = 0.25 * (1 - GaussPoints[1]) * (1 + GaussPoints[0]);
			double N4_22 = 0.25 * (1 - GaussPoints[1]) * (1 + GaussPoints[1]);

			// Nitxi Nityi
			double Nitxi_11 = N1_11 * t_x1 + N2_11 * t_x2 + N3_11 * t_x3 + N4_11 * t_x4;
			double Nitxi_12 = N1_12 * t_x1 + N2_12 * t_x2 + N3_12 * t_x3 + N4_12 * t_x4;
			double Nitxi_21 = N1_21 * t_x1 + N2_21 * t_x2 + N3_21 * t_x3 + N4_21 * t_x4;
			double Nitxi_22 = N1_22 * t_x1 + N2_22 * t_x2 + N3_22 * t_x3 + N4_22 * t_x4;

			double Nityi_11 = N1_11 * t_y1 + N2_11 * t_y2 + N3_11 * t_y3 + N4_11 * t_y4;
			double Nityi_12 = N1_12 * t_y1 + N2_12 * t_y2 + N3_12 * t_y3 + N4_12 * t_y4;
			double Nityi_21 = N1_21 * t_y1 + N2_21 * t_y2 + N3_21 * t_y3 + N4_21 * t_y4;
			double Nityi_22 = N1_22 * t_y1 + N2_22 * t_y2 + N3_22 * t_y3 + N4_22 * t_y4;

			double Nitzi_11 = N1_11 * t_z1 + N2_11 * t_z2 + N3_11 * t_z3 + N4_11 * t_z4;
			double Nitzi_12 = N1_12 * t_z1 + N2_12 * t_z2 + N3_12 * t_z3 + N4_12 * t_z4;
			double Nitzi_21 = N1_21 * t_z1 + N2_21 * t_z2 + N3_21 * t_z3 + N4_21 * t_z4;
			double Nitzi_22 = N1_22 * t_z1 + N2_22 * t_z2 + N3_22 * t_z3 + N4_22 * t_z4;

			f_Gamma[0] = N1_11 * Nitxi_11 * f11 + N1_12 * Nitxi_12 * f12 + N1_21 * Nitxi_21 * f21 + N1_22 * Nitxi_22 * f22;
			f_Gamma[1] = N1_11 * Nityi_11 * f11 + N1_12 * Nityi_12 * f12 + N1_21 * Nityi_21 * f21 + N1_22 * Nityi_22 * f22;
			f_Gamma[2] = N1_11 * Nitzi_11 * f11 + N1_12 * Nitzi_12 * f12 + N1_21 * Nitzi_21 * f21 + N1_22 * Nitzi_22 * f22;
			f_Gamma[3] = N2_11 * Nitxi_11 * f11 + N2_12 * Nitxi_12 * f12 + N2_21 * Nitxi_21 * f21 + N2_22 * Nitxi_22 * f22;
			f_Gamma[4] = N2_11 * Nityi_11 * f11 + N2_12 * Nityi_12 * f12 + N2_21 * Nityi_21 * f21 + N2_22 * Nityi_22 * f22;
			f_Gamma[5] = N2_11 * Nitzi_11 * f11 + N2_12 * Nitzi_12 * f12 + N2_21 * Nitzi_21 * f21 + N2_22 * Nitzi_22 * f22;
			f_Gamma[6] = N3_11 * Nitxi_11 * f11 + N3_12 * Nitxi_12 * f12 + N3_21 * Nitxi_21 * f21 + N3_22 * Nitxi_22 * f22;
			f_Gamma[7] = N3_11 * Nityi_11 * f11 + N3_12 * Nityi_12 * f12 + N3_21 * Nityi_21 * f21 + N3_22 * Nityi_22 * f22;
			f_Gamma[8] = N3_11 * Nitzi_11 * f11 + N3_12 * Nitzi_12 * f12 + N3_21 * Nitzi_21 * f21 + N3_22 * Nitzi_22 * f22;
			f_Gamma[9] = N4_11 * Nitxi_11 * f11 + N4_12 * Nitxi_12 * f12 + N4_21 * Nitxi_21 * f21 + N4_22 * Nitxi_22 * f22;
			f_Gamma[10] = N4_11 * Nityi_11 * f11 + N4_12 * Nityi_12 * f12 + N4_21 * Nityi_21 * f21 + N4_22 * Nityi_22 * f22;
			f_Gamma[11] = N4_11 * Nitzi_11 * f11 + N4_12 * Nitzi_12 * f12 + N4_21 * Nitzi_21 * f21 + N4_22 * Nitzi_22 * f22;
			unsigned int dof[12] = {
				nodes[ElementNode[0] - 1]->bcode[0],
				nodes[ElementNode[0] - 1]->bcode[1],
				nodes[ElementNode[0] - 1]->bcode[2],
				nodes[ElementNode[1] - 1]->bcode[0],
				nodes[ElementNode[1] - 1]->bcode[1],
				nodes[ElementNode[1] - 1]->bcode[2],
				nodes[ElementNode[2] - 1]->bcode[0],
				nodes[ElementNode[2] - 1]->bcode[1],
				nodes[ElementNode[2] - 1]->bcode[2],
				nodes[ElementNode[3] - 1]->bcode[0],
				nodes[ElementNode[3] - 1]->bcode[1],
				nodes[ElementNode[3] - 1]->bcode[2],
			};
			for (int i = 0; i < 12; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Gamma[i];
		}
		break;
	}
	default:
		std::cerr << "LodaCase " << LoadData->LoadCaseType_ << " not available. See CDomain::AssembleForce." << std::endl;
		exit(5);
		break;
	}

	return true;
}
