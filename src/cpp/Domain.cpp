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
template <class type> void clear( type* a, unsigned int N )
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
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
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
            
            *Output << setw(9) << Ele+1;
            for (int i=0; i<Element.GetND(); i++)
                *Output << setw(5) << LocationMatrix[i];
            *Output << endl;
#endif

            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }
    
    StiffnessMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
    *Output << endl;
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
			
			if(dof) // The DOF is activated
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
			double f_Omega[6] = {0};
			double b_x1 = LoadData->load[6*LoadOrder];
			double b_y1 = LoadData->load[6*LoadOrder+1];
			double b_x2 = LoadData->load[6*LoadOrder+2];
			double b_y2 = LoadData->load[6*LoadOrder+3];
			double b_x3 = LoadData->load[6*LoadOrder+4];
			double b_y3 = LoadData->load[6*LoadOrder+5];

			// Calculate area of the CST element
            CElement& Element = ElementGrp[Ele-1];
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
			double A = 0.5*(X[0]*Y[1] - X[1]*Y[0] + X[2]*Y[0] - X[0]*Y[2] + X[1]*Y[2] - X[2]*Y[1]);

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
			double f_Gamma[4] = {0};
			unsigned int ElementNode[2] = {
    			LoadData->node[2*LoadOrder], 
    			LoadData->node[2*LoadOrder+1],
			};
			double t_x1 = LoadData->load[4*LoadOrder];
			double t_y1 = LoadData->load[4*LoadOrder+1];
			double t_x2 = LoadData->load[4*LoadOrder+2];
			double t_y2 = LoadData->load[4*LoadOrder+3];

			// Calculate length of the line with surface force
            CElement& Element = ElementGrp[Ele-1];
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CCSTMaterial* material_ = dynamic_cast<CCSTMaterial*>(ElementMaterial);
			double t = material_->t;

			CNode** nodes = Element.GetNodes();

			double X[2] = {
				nodes[ElementNode[0]-1]->XYZ[0],
				nodes[ElementNode[1]-1]->XYZ[0],
			};
			double Y[2] = {
				nodes[ElementNode[0]-1]->XYZ[1],
				nodes[ElementNode[1]-1]->XYZ[1],
			};
			double L2 = (X[1] - X[0]) * (X[1] - X[0]) + (Y[1] - Y[0]) * (Y[1] - Y[0]);
			double L = sqrt(L2);

			f_Gamma[0] = (2 * t_x1 + t_x2) * L * t / 6.0;
			f_Gamma[1] = (2 * t_y1 + t_y2) * L * t / 6.0;
			f_Gamma[2] = (t_x1 + 2 * t_x2) * L * t / 6.0;
			f_Gamma[3] = (t_y1 + 2 * t_y2) * L * t / 6.0;

			unsigned int dof[4] = {
				nodes[ElementNode[0]-1]->bcode[0],
				nodes[ElementNode[0]-1]->bcode[1],
				nodes[ElementNode[1]-1]->bcode[0],
				nodes[ElementNode[1]-1]->bcode[1],
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
				LoadData->load[8*LoadOrder],
				LoadData->load[8*LoadOrder+1],
				LoadData->load[8*LoadOrder+2],
				LoadData->load[8*LoadOrder+3],
				LoadData->load[8*LoadOrder+4],
				LoadData->load[8*LoadOrder+5],
				LoadData->load[8*LoadOrder+6],
				LoadData->load[8*LoadOrder+7],
			};

			// Calculate area of the CST element
            CElement& Element = ElementGrp[Ele-1];
			CQ4* Element_ = dynamic_cast<CQ4*>(&Element);
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial);
			double t = material_->t;
			CNode** nodes = Element.GetNodes();

			double f_Omega[8] = {0};
			double GaussPoints[2] = {
				-sqrt(3) / 3.0,
				sqrt(3) / 3.0,
			};
			for (unsigned int i = 0; i < 2; i++)	// use full integration here, then is reduced integration
			{
				for (unsigned int j = 0; j < 2; j++)
				{	
					double B[3][8] = {0};
					double det;
					Element_->ElementStrainFunction(B, &det, GaussPoints[i], GaussPoints[j]);          
					double N[2][8] = {0};
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
			double f_Gamma[4] = {0};
			unsigned int ElementNode[2] = {
    			LoadData->node[2*LoadOrder], 
    			LoadData->node[2*LoadOrder+1],
			};
			double t_x1 = LoadData->load[4*LoadOrder];
			double t_y1 = LoadData->load[4*LoadOrder+1];
			double t_x2 = LoadData->load[4*LoadOrder+2];
			double t_y2 = LoadData->load[4*LoadOrder+3];

			// Calculate length of the line with surface force
            CElement& Element = ElementGrp[Ele-1];
			CMaterial* ElementMaterial = Element.GetElementMaterial();
			CQ4Material* material_ = dynamic_cast<CQ4Material*>(ElementMaterial);
			double t = material_->t;

			CNode** nodes = Element.GetNodes();

			double X[2] = {
				nodes[ElementNode[0]-1]->XYZ[0],
				nodes[ElementNode[1]-1]->XYZ[0],
			};
			double Y[2] = {
				nodes[ElementNode[0]-1]->XYZ[1],
				nodes[ElementNode[1]-1]->XYZ[1],
			};
			double L2 = (X[1] - X[0]) * (X[1] - X[0]) + (Y[1] - Y[0]) * (Y[1] - Y[0]);
			double L = sqrt(L2);

			f_Gamma[0] = (2 * t_x1 + t_x2) * L * t / 6.0;
			f_Gamma[1] = (2 * t_y1 + t_y2) * L * t / 6.0;
			f_Gamma[2] = (t_x1 + 2 * t_x2) * L * t / 6.0;
			f_Gamma[3] = (t_y1 + 2 * t_y2) * L * t / 6.0;

			unsigned int dof[4] = {
				nodes[ElementNode[0]-1]->bcode[0],
				nodes[ElementNode[0]-1]->bcode[1],
				nodes[ElementNode[1]-1]->bcode[0],
				nodes[ElementNode[1]-1]->bcode[1],
			};
			for (int i = 0; i < 4; i++)
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

