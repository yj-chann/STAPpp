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
	
	clear(Force, NEQ);	// For resolving the non-homogeneous essential boundary conditions

//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		CElementGroup& ElementGrp = EleGrpList[EleGrp];
		unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();		
		double* Matrix = new double[size];

		unsigned int ND = ElementGrp[0].GetND();
		double* NonForce = new double[ND];	// Force caused by the non-homogeneous essential boundary conditions
		// double* dE = new double[ND];	// d including non-homogeneous essential boundary conditions
		// unsigned int NEN = ElementGrp[0].GetNEN();
		// unsigned int DOF = ND / NEN;

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
            Element.ElementStiffness(Matrix);
            StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());		

			if (Element.GetNonHomo()) {
				clear(NonForce, Element.GetND());
				unsigned int* LocationMatrix = Element.GetLocationMatrix();
				Element.ElementNonHomo(Matrix, NonForce);
				for (unsigned int i = 0; i < ND; i++)
				{
					if (LocationMatrix[i] == 0)
						continue;
					Force[LocationMatrix[i] - 1] -= NonForce[i];
				}
			}			
						
			// CNode** nodes = Element.GetNodes();
			// for (unsigned int i = 0; i < NEN; i++) {
			// 	if (nodes[i]->NonHomo) {
			// 		for (unsigned j = 0; j < DOF; j++) {
			// 		}
			// 	}					
			// }
        }

		delete[] Matrix;
		Matrix = nullptr;
		delete[] NonForce;
		NonForce = nullptr;
		// delete[] dE;
		// dE = nullptr;
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

	// if (LoadCase == 1)	// Clear the Force vector in CDomain::AssembleStiffnessMatrix
    // 	clear(Force, NEQ);

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
	case 9:
		break;
	case 10:
		break;
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
			double f_Gamma[12] = {0};
			unsigned int ElementNode[4] = {
				LoadData->node[4*LoadOrder],
				LoadData->node[4*LoadOrder+1],
				LoadData->node[4*LoadOrder+2],
				LoadData->node[4*LoadOrder+3],
			};

			double t_x1 = LoadData->load[12*LoadOrder];
			double t_y1 = LoadData->load[12*LoadOrder+1];
			double t_z1 = LoadData->load[12*LoadOrder+2];
			double t_x2 = LoadData->load[12*LoadOrder+3];
			double t_y2 = LoadData->load[12*LoadOrder+4];
			double t_z2 = LoadData->load[12*LoadOrder+5];
			double t_x3 = LoadData->load[12*LoadOrder+6];
			double t_y3 = LoadData->load[12*LoadOrder+7];
			double t_z3 = LoadData->load[12*LoadOrder+8];
			double t_x4 = LoadData->load[12*LoadOrder+9];
			double t_y4 = LoadData->load[12*LoadOrder+10];
			double t_z4 = LoadData->load[12*LoadOrder+11];

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
			double J11_11 = 0.25 * x_21 * (1 - GaussPoints[0]) - 0.25 * x_43 * (1 + GaussPoints[0]);
			double J12_11 = 0.25 * y_21 * (1 - GaussPoints[0]) - 0.25 * y_43 * (1 + GaussPoints[0]);
			double J13_11 = 0.25 * z_21 * (1 - GaussPoints[0]) - 0.25 * z_43 * (1 + GaussPoints[0]);
			double J21_11 = 0.25 * x_41 * (1 - GaussPoints[0]) + 0.25 * x_32 * (1 + GaussPoints[0]);
			double J22_11 = 0.25 * y_41 * (1 - GaussPoints[0]) + 0.25 * y_32 * (1 + GaussPoints[0]);
			double J23_11 = 0.25 * z_41 * (1 - GaussPoints[0]) + 0.25 * z_32 * (1 + GaussPoints[0]);

			double J11_12 = 0.25 * x_21 * (1 - GaussPoints[1]) - 0.25 * x_43 * (1 + GaussPoints[1]);
			double J12_12 = 0.25 * y_21 * (1 - GaussPoints[1]) - 0.25 * y_43 * (1 + GaussPoints[1]);
			double J13_12 = 0.25 * z_21 * (1 - GaussPoints[1]) - 0.25 * z_43 * (1 + GaussPoints[1]);
			double J21_12 = 0.25 * x_41 * (1 - GaussPoints[0]) + 0.25 * x_32 * (1 + GaussPoints[0]);
			double J22_12 = 0.25 * y_41 * (1 - GaussPoints[0]) + 0.25 * y_32 * (1 + GaussPoints[0]);
			double J23_12 = 0.25 * z_41 * (1 - GaussPoints[0]) + 0.25 * z_32 * (1 + GaussPoints[0]);

			double J11_21 = 0.25 * x_21 * (1 - GaussPoints[0]) - 0.25 * x_43 * (1 + GaussPoints[0]);
			double J12_21 = 0.25 * y_21 * (1 - GaussPoints[0]) - 0.25 * y_43 * (1 + GaussPoints[0]);
			double J13_21 = 0.25 * z_21 * (1 - GaussPoints[0]) - 0.25 * z_43 * (1 + GaussPoints[0]);
			double J21_21 = 0.25 * x_41 * (1 - GaussPoints[1]) + 0.25 * x_32 * (1 + GaussPoints[1]);
			double J22_21 = 0.25 * y_41 * (1 - GaussPoints[1]) + 0.25 * y_32 * (1 + GaussPoints[1]);
			double J23_21 = 0.25 * z_41 * (1 - GaussPoints[1]) + 0.25 * z_32 * (1 + GaussPoints[1]);

			double J11_22 = 0.25 * x_21 * (1 - GaussPoints[1]) - 0.25 * x_43 * (1 + GaussPoints[1]);
			double J12_22 = 0.25 * y_21 * (1 - GaussPoints[1]) - 0.25 * y_43 * (1 + GaussPoints[1]);
			double J13_22 = 0.25 * z_21 * (1 - GaussPoints[1]) - 0.25 * z_43 * (1 + GaussPoints[1]);
			double J21_22 = 0.25 * x_41 * (1 - GaussPoints[1]) + 0.25 * x_32 * (1 + GaussPoints[1]);
			double J22_22 = 0.25 * y_41 * (1 - GaussPoints[1]) + 0.25 * y_32 * (1 + GaussPoints[1]);
			double J23_22 = 0.25 * z_41 * (1 - GaussPoints[1]) + 0.25 * z_32 * (1 + GaussPoints[1]);


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
	case 13:	// All surface forces of S8R5 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 8)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->node[LoadOrder];
			double t[40] = {
				LoadData->load[24*LoadOrder],
				LoadData->load[24*LoadOrder+1],
				LoadData->load[24*LoadOrder+2],
				0,
				0,
				LoadData->load[24*LoadOrder+3],
				LoadData->load[24*LoadOrder+4],
				LoadData->load[24*LoadOrder+5],
				0,
				0,
				LoadData->load[24*LoadOrder+6],
				LoadData->load[24*LoadOrder+7],
				LoadData->load[24*LoadOrder+8],
				0,
				0,
				LoadData->load[24*LoadOrder+9],
				LoadData->load[24*LoadOrder+10],
				LoadData->load[24*LoadOrder+11],
				0,
				0,
				LoadData->load[24*LoadOrder+12],
				LoadData->load[24*LoadOrder+13],
				LoadData->load[24*LoadOrder+14],
				0,
				0,
				LoadData->load[24*LoadOrder+15],
				LoadData->load[24*LoadOrder+16],
				LoadData->load[24*LoadOrder+17],
				0,
				0,
				LoadData->load[24*LoadOrder+18],
				LoadData->load[24*LoadOrder+19],
				LoadData->load[24*LoadOrder+20],
				0,
				0,
				LoadData->load[24*LoadOrder+21],
				LoadData->load[24*LoadOrder+22],
				LoadData->load[24*LoadOrder+23],
				0,
				0,
			};

			CElement& Element = ElementGrp[Ele - 1];
			CS8R5* Element_ = dynamic_cast<CS8R5*>(&Element);
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CS8R5Material* material_ = dynamic_cast<CS8R5Material*>(ElementMaterial);
			CNode** nodes = Element.GetNodes();

			double f_Omega[40] = {0};
			double GaussPoints[2] = {
				-sqrt(3) / 3.0,
				sqrt(3) / 3.0,
			};

			for (unsigned int gp1 = 0; gp1 < 2; gp1++)	//  reduced integration
			{
				for (unsigned int gp2 = 0; gp2 < 2; gp2++)
				{
					double det = 0;
					Element_->ElementDetTop(&det, GaussPoints[gp1], GaussPoints[gp2]);	// Top Area
					double N[3][40] = {0};
					Element_->ElementShapeFunction(N, GaussPoints[gp1], GaussPoints[gp2], 1.0);
					for (int i = 0; i < 40; i++)
						for (int j = 0; j < 3; j++)
							for (int k = 0; k < 40; k++)
								f_Omega[i] += N[j][i] * N[j][k] * t[k] * det;  // Weight is 1			
				}
			}

			unsigned int dof[40] = {
				nodes[0]->bcode[0],
				nodes[0]->bcode[1],
				nodes[0]->bcode[2],
				nodes[0]->bcode[3],
				nodes[0]->bcode[4],
				nodes[1]->bcode[0],
				nodes[1]->bcode[1],
				nodes[1]->bcode[2],
				nodes[1]->bcode[3],
				nodes[1]->bcode[4],
				nodes[2]->bcode[0],
				nodes[2]->bcode[1],
				nodes[2]->bcode[2],
				nodes[2]->bcode[3],
				nodes[2]->bcode[4],
				nodes[3]->bcode[0],
				nodes[3]->bcode[1],
				nodes[3]->bcode[2],
				nodes[3]->bcode[3],
				nodes[3]->bcode[4],
				nodes[4]->bcode[0],
				nodes[4]->bcode[1],
				nodes[4]->bcode[2],
				nodes[4]->bcode[3],
				nodes[4]->bcode[4],
				nodes[5]->bcode[0],
				nodes[5]->bcode[1],
				nodes[5]->bcode[2],
				nodes[5]->bcode[3],
				nodes[5]->bcode[4],
				nodes[6]->bcode[0],
				nodes[6]->bcode[1],
				nodes[6]->bcode[2],
				nodes[6]->bcode[3],
				nodes[6]->bcode[4],
				nodes[7]->bcode[0],
				nodes[7]->bcode[1],
				nodes[7]->bcode[2],
				nodes[7]->bcode[3],
				nodes[7]->bcode[4]
			};

			for (int i = 0; i < 40; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Omega[i];
		}
		break;
	}
	case 14:	// All body forces of S8R5 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 8)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->node[LoadOrder];
			double b[40] = {
				LoadData->load[24*LoadOrder],
				LoadData->load[24*LoadOrder+1],
				LoadData->load[24*LoadOrder+2],
				0,
				0,
				LoadData->load[24*LoadOrder+3],
				LoadData->load[24*LoadOrder+4],
				LoadData->load[24*LoadOrder+5],
				0,
				0,
				LoadData->load[24*LoadOrder+6],
				LoadData->load[24*LoadOrder+7],
				LoadData->load[24*LoadOrder+8],
				0,
				0,
				LoadData->load[24*LoadOrder+9],
				LoadData->load[24*LoadOrder+10],
				LoadData->load[24*LoadOrder+11],
				0,
				0,
				LoadData->load[24*LoadOrder+12],
				LoadData->load[24*LoadOrder+13],
				LoadData->load[24*LoadOrder+14],
				0,
				0,
				LoadData->load[24*LoadOrder+15],
				LoadData->load[24*LoadOrder+16],
				LoadData->load[24*LoadOrder+17],
				0,
				0,
				LoadData->load[24*LoadOrder+18],
				LoadData->load[24*LoadOrder+19],
				LoadData->load[24*LoadOrder+20],
				0,
				0,
				LoadData->load[24*LoadOrder+21],
				LoadData->load[24*LoadOrder+22],
				LoadData->load[24*LoadOrder+23],
				0,
				0,
			};

			CElement& Element = ElementGrp[Ele - 1];
			CS8R5* Element_ = dynamic_cast<CS8R5*>(&Element);
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CS8R5Material* material_ = dynamic_cast<CS8R5Material*>(ElementMaterial);
			CNode** nodes = Element.GetNodes();

			double f_Omega[40] = {0};
			double GaussPoints[2] = {
				-sqrt(3) / 3.0,
				sqrt(3) / 3.0,
			};

			for (unsigned int gp1 = 0; gp1 < 2; gp1++)	//  reduced integration
			{
				for (unsigned int gp2 = 0; gp2 < 2; gp2++)
				{
					for (unsigned int gp3 = 0; gp3 < 2; gp3++)
					{
						double B[5][40] = {0};
						double det = 0;
						Element_->ElementStrainFunction(B, &det, GaussPoints[gp1], GaussPoints[gp2], GaussPoints[gp3]);
						double N[3][40] = {0};
						Element_->ElementShapeFunction(N, GaussPoints[gp1], GaussPoints[gp2], GaussPoints[gp3]);
						for (int i = 0; i < 40; i++)
							for (int j = 0; j < 3; j++)
								for (int k = 0; k < 40; k++)
									f_Omega[i] += N[j][i] * N[j][k] * b[k] * det;  // Weight is 1
					}
				}
			}

			unsigned int dof[40] = {
				nodes[0]->bcode[0],
				nodes[0]->bcode[1],
				nodes[0]->bcode[2],
				nodes[0]->bcode[3],
				nodes[0]->bcode[4],
				nodes[1]->bcode[0],
				nodes[1]->bcode[1],
				nodes[1]->bcode[2],
				nodes[1]->bcode[3],
				nodes[1]->bcode[4],
				nodes[2]->bcode[0],
				nodes[2]->bcode[1],
				nodes[2]->bcode[2],
				nodes[2]->bcode[3],
				nodes[2]->bcode[4],
				nodes[3]->bcode[0],
				nodes[3]->bcode[1],
				nodes[3]->bcode[2],
				nodes[3]->bcode[3],
				nodes[3]->bcode[4],
				nodes[4]->bcode[0],
				nodes[4]->bcode[1],
				nodes[4]->bcode[2],
				nodes[4]->bcode[3],
				nodes[4]->bcode[4],
				nodes[5]->bcode[0],
				nodes[5]->bcode[1],
				nodes[5]->bcode[2],
				nodes[5]->bcode[3],
				nodes[5]->bcode[4],
				nodes[6]->bcode[0],
				nodes[6]->bcode[1],
				nodes[6]->bcode[2],
				nodes[6]->bcode[3],
				nodes[6]->bcode[4],
				nodes[7]->bcode[0],
				nodes[7]->bcode[1],
				nodes[7]->bcode[2],
				nodes[7]->bcode[3],
				nodes[7]->bcode[4]
			};

			for (int i = 0; i < 40; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Omega[i];
		}
		break;
	}
	case 15:	// All body forces Mindlin-Reissner Plate element, except body moment
	{	
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 9)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->node[LoadOrder];
			double b[12] = {
				0,
				0,
				LoadData->load[4*LoadOrder],
				0,
				0,
				LoadData->load[4*LoadOrder+1],
				0,
				0,
				LoadData->load[4*LoadOrder+2],
				0,
				0,
				LoadData->load[4*LoadOrder+3],
			};

			CElement& Element = ElementGrp[Ele - 1];
			CMP* Element_ = dynamic_cast<CMP*>(&Element);
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CMPMaterial* material_ = dynamic_cast<CMPMaterial*>(ElementMaterial);
			CNode** nodes = Element.GetNodes();

			double f_Omega[12] = {0};
			double GaussPoints[2] = {
				-sqrt(3) / 3.0,
				sqrt(3) / 3.0,
			};

			for (unsigned int gp1 = 0; gp1 < 2; gp1++)	//  reduced integration
			{
				for (unsigned int gp2 = 0; gp2 < 2; gp2++)
				{					
					double Bb[3][12] = {0};
					double Bs[2][12] = {0};
					double det = 0;
					Element_->ElementStrainFunction(Bb, Bs, &det, GaussPoints[gp1], GaussPoints[gp2]);
					double N[3][12] = {0};
					Element_->ElementShapeFunction(N, GaussPoints[gp1], GaussPoints[gp2]);
					for (int i = 0; i < 12; i++)
						for (int j = 0; j < 3; j++)
							for (int k = 0; k < 12; k++)
								f_Omega[i] += N[j][i] * N[j][k] * b[k] * det;  // Weight is 1					
				}
			}

			unsigned int dof[12] = {
				nodes[0]->bcode[3],
				nodes[0]->bcode[4],
				nodes[0]->bcode[2],
				nodes[1]->bcode[3],
				nodes[1]->bcode[4],
				nodes[1]->bcode[2],
				nodes[2]->bcode[3],
				nodes[2]->bcode[4],
				nodes[2]->bcode[2],
				nodes[3]->bcode[3],
				nodes[3]->bcode[4],
				nodes[3]->bcode[2],
			};

			for (int i = 0; i < 12; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Omega[i];
		}
		break;
	}
	case 16:	// All body forces of basic Plate element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 10)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->node[LoadOrder];
			double b[12] = {
				LoadData->load[4*LoadOrder],
				0,
				0,
				LoadData->load[4*LoadOrder+1],
				0,
				0,
				LoadData->load[4*LoadOrder+2],
				0,
				0,
				LoadData->load[4*LoadOrder+3],
				0,
				0,
			};

			CElement& Element = ElementGrp[Ele - 1];
			CPlate* Element_ = dynamic_cast<CPlate*>(&Element);
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CPlateMaterial* material_ = dynamic_cast<CPlateMaterial*>(ElementMaterial);
			CNode** nodes = Element.GetNodes();

			double f_Omega[12] = {0};
			double GaussPoints[3] = {
				-sqrt(15) / 5.0,
				0.0,
				sqrt(15) / 5.0,
			};

			double Weight[3] = {
				5.0 / 9.0,
				8.0 / 9.0,
				5.0 / 9.0,
			};

			for (unsigned int gp1 = 0; gp1 < 3; gp1++)	//  full integration
			{
				for (unsigned int gp2 = 0; gp2 < 3; gp2++)
				{					
					double B[3][12] = {0};
					double det = 0;
					Element_->ElementStrainFunction(B, &det, GaussPoints[gp1], GaussPoints[gp2]);
					double N[1][12] = {0};
					Element_->ElementShapeFunction(N, GaussPoints[gp1], GaussPoints[gp2]);
					for (int i = 0; i < 12; i++)		
						for (int j = 0; j < 12; j++)
							f_Omega[i] += Weight[gp1] * Weight[gp2] * N[0][i] * N[0][j] * b[j] * det;  // Weight is 1					
				}
			}

			unsigned int dof[12] = {
				nodes[0]->bcode[2],
				nodes[0]->bcode[3],
				nodes[0]->bcode[4],
				nodes[1]->bcode[2],
				nodes[1]->bcode[3],
				nodes[1]->bcode[4],
				nodes[2]->bcode[2],
				nodes[2]->bcode[3],
				nodes[2]->bcode[4],
				nodes[3]->bcode[2],
				nodes[3]->bcode[3],
				nodes[3]->bcode[4],
			};

			for (int i = 0; i < 12; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Omega[i];
		}
		break;
	}
	case 17:	// All body forces of Tet4 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 11)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->node[LoadOrder];
			double f_Omega[12] = {0};
			double b_x1 = LoadData->load[12*LoadOrder];
			double b_y1 = LoadData->load[12*LoadOrder+1];
			double b_z1 = LoadData->load[12*LoadOrder+2];
			double b_x2 = LoadData->load[12*LoadOrder+3];
			double b_y2 = LoadData->load[12*LoadOrder+4];
			double b_z2 = LoadData->load[12*LoadOrder+5];
			double b_x3 = LoadData->load[12*LoadOrder+6];
			double b_y3 = LoadData->load[12*LoadOrder+7];
			double b_z3 = LoadData->load[12*LoadOrder+8];
			double b_x4 = LoadData->load[12*LoadOrder+9];
			double b_y4 = LoadData->load[12*LoadOrder+10];
			double b_z4 = LoadData->load[12*LoadOrder+11];

			// Calculate volume of the Tet4 element
            CElement& Element = ElementGrp[Ele-1];

			CNode** nodes = Element.GetNodes();

			double X[4] = {
				nodes[0]->XYZ[0],
				nodes[1]->XYZ[0],
				nodes[2]->XYZ[0],
				nodes[3]->XYZ[0]
			};
			double Y[4] = {
				nodes[0]->XYZ[1],
				nodes[2]->XYZ[1],
				nodes[2]->XYZ[1],
				nodes[3]->XYZ[1],
			};
			double Z[4] = {
				nodes[0]->XYZ[2],
				nodes[2]->XYZ[2],
				nodes[2]->XYZ[2],
				nodes[3]->XYZ[2],
			};
			double V = fabs(X[1]*(Y[2]*Z[3]-Y[3]*Z[2]) - X[2]*(Y[1]*Z[3]-Y[3]*Z[1]) + X[3]*(Y[1]*Z[2]-Y[2]*Z[1]) - X[0]*(Y[2]*Z[3]-Y[3]*Z[2]) + X[2]*(Y[0]*Z[3]-Y[3]*Z[0]) - X[3]*(Y[0]*Z[2]-Y[2]*Z[0]) + X[0]*(Y[1]*Z[3]-Y[3]*Z[1]) - X[1]*(Y[0]*Z[3]-Y[3]*Z[0]) + X[3]*(Y[0]*Z[1]-Y[1]*Z[0]) - X[0]*(Y[1]*Z[2]-Y[2]*Z[1]) + X[1]*(Y[0]*Z[2]-Y[2]*Z[0]) - X[2]*(Y[0]*Z[1]-Y[1]*Z[0])) / 6.0;

			f_Omega[0] = (2 * b_x1 + b_x2 + b_x3 + b_x4) * V / 20.0;
			f_Omega[1] = (2 * b_y1 + b_y2 + b_y3 + b_y4) * V / 20.0;
			f_Omega[2] = (2 * b_z1 + b_z2 + b_z3 + b_z4) * V / 20.0;
			f_Omega[3] = (b_x1 + 2 * b_x2 + b_x3 + b_x4) * V / 20.0;
			f_Omega[4] = (b_y1 + 2 * b_y2 + b_y3 + b_y4) * V / 20.0;
			f_Omega[5] = (b_z1 + 2 * b_z2 + b_z3 + b_z4) * V / 20.0;
			f_Omega[6] = (b_x1 + b_x2 + 2 * b_x3 + b_x4) * V / 20.0;
			f_Omega[7] = (b_y1 + b_y2 + 2 * b_y3 + b_y4) * V / 20.0;
			f_Omega[8] = (b_z1 + b_z2 + 2 * b_z3 + b_z4) * V / 20.0;
			f_Omega[9] = (b_x1 + b_x2 + b_x3 + 2 * b_x4) * V / 20.0;
			f_Omega[10] = (b_y1 + b_y2 + b_y3 + 2 * b_y4) * V / 20.0;
			f_Omega[11] = (b_z1 + b_z2 + b_z3 + 2 * b_z4) * V / 20.0;

			unsigned int dof[12] = {
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
			};
			for (int i = 0; i < 12; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Omega[i];		
		}
		break;
	}
	case 18:	// All surface forces of Tet4 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup& ElementGrp = EleGrpList[EleGrp];
        	unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 11)
				break;
		}
		CElementGroup& ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
        {	
			unsigned int Ele = LoadData->dof[LoadOrder];
			double f_Gamma[9] = {0};
			unsigned int ElementNode[3] = {
    			LoadData->node[3*LoadOrder], 
    			LoadData->node[3*LoadOrder+1],
				LoadData->node[3*LoadOrder+2],
			};
			double t_x1 = LoadData->load[9*LoadOrder];
			double t_y1 = LoadData->load[9*LoadOrder+1];
			double t_z1 = LoadData->load[9*LoadOrder+2];
			double t_x2 = LoadData->load[9*LoadOrder+3];
			double t_y2 = LoadData->load[9*LoadOrder+4];
			double t_z2 = LoadData->load[9*LoadOrder+5];
			double t_x3 = LoadData->load[9*LoadOrder+6];
			double t_y3 = LoadData->load[9*LoadOrder+7];
			double t_z3 = LoadData->load[9*LoadOrder+8];

			// Calculate area with surface force
            CElement& Element = ElementGrp[Ele-1];

			CNode** nodes = Element.GetNodes();

			double X[3] = {
				nodes[ElementNode[0]-1]->XYZ[0],
				nodes[ElementNode[1]-1]->XYZ[0],
				nodes[ElementNode[2]-1]->XYZ[0],
			};
			double Y[3] = {
				nodes[ElementNode[0]-1]->XYZ[1],
				nodes[ElementNode[1]-1]->XYZ[1],
				nodes[ElementNode[2]-1]->XYZ[1],
			};
			double Z[3] = {
				nodes[ElementNode[0]-1]->XYZ[2],
				nodes[ElementNode[1]-1]->XYZ[2],
				nodes[ElementNode[2]-1]->XYZ[2],
			};
			double A = sqrt(pow((Y[1]-Y[0])*(Z[2]-Z[0])-(Y[2]-Y[0])*(Z[1]-Z[0]),2) + pow((Z[1]-Z[0])*(X[2]-X[0])-(Z[2]-Z[0])*(X[1]-X[0]),2) + pow((X[1]-X[0])*(Y[2]-Y[0])-(X[2]-X[0])*(Y[1]-Y[0]),2)) / 2.0;

			f_Gamma[0] = (2 * t_x1 + t_x2 + t_x3) * A / 12.0;
			f_Gamma[1] = (2 * t_y1 + t_y2 + t_y3) * A / 12.0;
			f_Gamma[2] = (2 * t_z1 + t_z2 + t_z3) * A / 12.0;
			f_Gamma[3] = (t_x1 + 2 * t_x2 + t_x3) * A / 12.0;
			f_Gamma[4] = (t_y1 + 2 * t_y2 + t_y3) * A / 12.0;
			f_Gamma[5] = (t_z1 + 2 * t_z2 + t_z3) * A / 12.0;
			f_Gamma[6] = (t_x1 + t_x2 + 2 * t_x3) * A / 12.0;
			f_Gamma[7] = (t_y1 + t_y2 + 2 * t_y3) * A / 12.0;
			f_Gamma[8] = (t_z1 + t_z2 + 2 * t_z3) * A / 12.0;

			unsigned int dof[9] = {
				nodes[ElementNode[0]-1]->bcode[0],
				nodes[ElementNode[0]-1]->bcode[1],
				nodes[ElementNode[0]-1]->bcode[2],
				nodes[ElementNode[1]-1]->bcode[0],
				nodes[ElementNode[1]-1]->bcode[1],
				nodes[ElementNode[1]-1]->bcode[2],
				nodes[ElementNode[2]-1]->bcode[0],
				nodes[ElementNode[2]-1]->bcode[1],
				nodes[ElementNode[2]-1]->bcode[2],
			};
			for (int i = 0; i < 9; i++)
				if (dof[i])	// The DOF is activated
					Force[dof[i] - 1] += f_Gamma[i];			
		}	
		break;
	}
	case 20: // All surface forces of T6 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup &ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 12)
				break;
		}
		CElementGroup &ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->dof[LoadOrder];
			double f_Gamma[6] = {0}; // 3节点×2自由度

			// T6单元的面力施加在3节点边上
			unsigned int ElementNode[3] = {
				LoadData->node[3 * LoadOrder],
				LoadData->node[3 * LoadOrder + 1],
				LoadData->node[3 * LoadOrder + 2]};

			// 每个节点的面力分量
			double t_x1 = LoadData->load[6 * LoadOrder];
			double t_y1 = LoadData->load[6 * LoadOrder + 1];
			double t_x2 = LoadData->load[6 * LoadOrder + 2];
			double t_y2 = LoadData->load[6 * LoadOrder + 3];
			double t_x3 = LoadData->load[6 * LoadOrder + 4];
			double t_y3 = LoadData->load[6 * LoadOrder + 5];

			// 获取单元信息
			CElement &Element = ElementGrp[Ele - 1];
			CMaterial *ElementMaterial = Element.GetElementMaterial();
			CT6Material *material_ = dynamic_cast<CT6Material *>(ElementMaterial);
			double t = material_->t;

			CNode **nodes = Element.GetNodes();

			// 获取节点坐标
			double X[3] = {
				nodes[ElementNode[0] - 1]->XYZ[0],
				nodes[ElementNode[1] - 1]->XYZ[0],
				nodes[ElementNode[2] - 1]->XYZ[0]};
			double Y[3] = {
				nodes[ElementNode[0] - 1]->XYZ[1],
				nodes[ElementNode[1] - 1]->XYZ[1],
				nodes[ElementNode[2] - 1]->XYZ[1]};

			// 使用2点高斯积分计算曲线边上的面力
			double GaussPoints[2] = {-sqrt(3) / 3.0, sqrt(3) / 3.0};
			double weights[2] = {1.0, 1.0};

			for (int gp = 0; gp < 2; gp++)
			{
				double xi = GaussPoints[gp];
				double weight = weights[gp];

				// 计算形函数及其导数
				double N1 = 0.5 * xi * (xi - 1.0);
				double N2 = 0.5 * xi * (xi + 1.0);
				double N3 = 1.0 - xi * xi;

				double N1_xi = xi - 0.5;
				double N2_xi = xi + 0.5;
				double N3_xi = -2.0 * xi;

				// 计算雅可比(曲线边的长度变化率)
				double dx_dxi = X[0] * N1_xi + X[1] * N2_xi + X[2] * N3_xi;
				double dy_dxi = Y[0] * N1_xi + Y[1] * N2_xi + Y[2] * N3_xi;
				double jac = sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);

				// 插值得到当前点的面力
				double tx = N1 * t_x1 + N2 * t_x2 + N3 * t_x3;
				double ty = N1 * t_y1 + N2 * t_y2 + N3 * t_y3;

				// 累加到节点力
				f_Gamma[0] += N1 * tx * jac * t * weight;
				f_Gamma[1] += N1 * ty * jac * t * weight;
				f_Gamma[2] += N2 * tx * jac * t * weight;
				f_Gamma[3] += N2 * ty * jac * t * weight;
				f_Gamma[4] += N3 * tx * jac * t * weight;
				f_Gamma[5] += N3 * ty * jac * t * weight;
			}

			unsigned int dof[6] = {
				nodes[ElementNode[0] - 1]->bcode[0],
				nodes[ElementNode[0] - 1]->bcode[1],
				nodes[ElementNode[1] - 1]->bcode[0],
				nodes[ElementNode[1] - 1]->bcode[1],
				nodes[ElementNode[2] - 1]->bcode[0],
				nodes[ElementNode[2] - 1]->bcode[1]};

			for (int i = 0; i < 6; i++)
				if (dof[i]) // 自由度被激活
					Force[dof[i] - 1] += f_Gamma[i];
		}
		break;
	}
	break;
	case 19: // All body forces of T6 element
	{
		unsigned int EleGrp = 0;
		for (; EleGrp < NUMEG; EleGrp++)
		{
			CElementGroup &ElementGrp = EleGrpList[EleGrp];
			unsigned int ElementType = ElementGrp.GetElementType();
			if (ElementType == 12)
				break;
		}
		CElementGroup &ElementGrp = EleGrpList[EleGrp];

		for (unsigned int LoadOrder = 0; LoadOrder < LoadData->nloads; LoadOrder++)
		{
			unsigned int Ele = LoadData->node[LoadOrder];
			double b[12] = {// 6节点×2自由度
							LoadData->load[12 * LoadOrder],
							LoadData->load[12 * LoadOrder + 1],
							LoadData->load[12 * LoadOrder + 2],
							LoadData->load[12 * LoadOrder + 3],
							LoadData->load[12 * LoadOrder + 4],
							LoadData->load[12 * LoadOrder + 5],
							LoadData->load[12 * LoadOrder + 6],
							LoadData->load[12 * LoadOrder + 7],
							LoadData->load[12 * LoadOrder + 8],
							LoadData->load[12 * LoadOrder + 9],
							LoadData->load[12 * LoadOrder + 10],
							LoadData->load[12 * LoadOrder + 11]};

			CElement &Element = ElementGrp[Ele - 1];
			CT6 *Element_ = dynamic_cast<CT6 *>(&Element);
			CMaterial *ElementMaterial = Element.GetElementMaterial();
			CT6Material *material_ = dynamic_cast<CT6Material *>(ElementMaterial);
			double t = material_->t;
			CNode **nodes = Element.GetNodes();

			double f_Omega[12] = {0};

			// 使用3点减缩积分方案
			double GaussPoints[3][3] = {
				{1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0}, // 权重为1/3
				{2.0 / 3.0, 1.0 / 6.0, 1.0 / 3.0},
				{1.0 / 6.0, 2.0 / 3.0, 1.0 / 3.0}};

			for (unsigned int gp = 0; gp < 3; gp++)
			{
				double xi = GaussPoints[gp][0];
				double eta = GaussPoints[gp][1];
				double weight = GaussPoints[gp][2];

				double B[3][12] = {0};
				double det;
				Element_->ElementStrainFunction(B, &det, xi, eta);

				double N[2][12] = {0};
				Element_->ElementShapeFunction(N, xi, eta);

				for (int k = 0; k < 12; k++)
					for (int l = 0; l < 2; l++)
						for (int m = 0; m < 12; m++)
							f_Omega[k] += t * N[l][k] * N[l][m] * b[m] * det * weight;
			}

			unsigned int dof[12] = {
				nodes[0]->bcode[0], nodes[0]->bcode[1],
				nodes[1]->bcode[0], nodes[1]->bcode[1],
				nodes[2]->bcode[0], nodes[2]->bcode[1],
				nodes[3]->bcode[0], nodes[3]->bcode[1],
				nodes[4]->bcode[0], nodes[4]->bcode[1],
				nodes[5]->bcode[0], nodes[5]->bcode[1]};

			for (int i = 0; i < 12; i++)
				if (dof[i]) // 自由度被激活
					Force[dof[i] - 1] += f_Omega[i];
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
