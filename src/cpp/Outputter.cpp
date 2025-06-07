/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <ctime>

#include "Domain.h"
#include "Outputter.h"
#include "SkylineMatrix.h"

using namespace std;

//	Output current time and date
void COutputter::PrintTime(const struct tm* ptm, COutputter &output)
{
	const char* weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
							 "Thursday", "Friday", "Saturday"};
	const char* month[] = {"January", "February", "March", "April", "May", "June",
						   "July", "August", "September", "October", "November", "December"};

	output << "        (";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " on ";
	output << month[ptm->tm_mon] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << ", "
		   << weekday[ptm->tm_wday] << ")" << endl
		   << endl;
}

COutputter* COutputter::_instance = nullptr;

//	Constructor
COutputter::COutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutputter* COutputter::GetInstance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
    
	return _instance;
}

//	Print program logo
void COutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "TITLE : " << FEMData->GetTitle() << endl;

	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	PrintTime(timeinfo, *this);
}

//	Print nodal data
void COutputter::OutputNodeInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	CNode* NodeList = FEMData->GetNodeList();

	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	*this << "      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =" << setw(6) << NUMNP << endl;
	*this << "      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =" << setw(6) << NUMEG << endl;
	*this << "      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =" << setw(6) << NLCASE << endl;
	*this << "      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =" << setw(6) << MODEX << endl;
	*this << "         EQ.0, DATA CHECK" << endl
		  << "         EQ.1, EXECUTION" << endl
		  << endl;

	*this << " N O D A L   P O I N T   D A T A" << endl << endl;
	*this << "    NODE       BOUNDARY                         NODAL POINT" << endl
		  << "   NUMBER  CONDITION  CODES                     COORDINATES" << endl;

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this);

	*this << endl;
}

//	Output equation numbers
void COutputter::OutputEquationNumber()
{
	CDomain* FEMData = CDomain::GetInstance();
	unsigned int NUMNP = FEMData->GetNUMNP();

	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES 	OF 	FREEDOM" << endl;
	*this << "        N           X    Y    Z    x    y    z" << endl;

	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this);

	*this << endl;
}

//	Output element data
void COutputter::OutputElementInfo()
{
	//	Print element group control line

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;

		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

		*this << " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =" << setw(5)
			  << ElementType << endl;
		*this << "     EQ.1, TRUSS ELEMENTS" << endl
			  << "     EQ.2, CST" << endl
			  << "     EQ.3, Q4" << endl
			  << "     EQ.4, Q8" << endl
		      << "     EQ.7, H8" << endl
			  << "     EQ.11, Tet4" << endl
			  << "     EQ.12, T6" << endl
			  << endl;

		*this << " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =" << setw(5) << NUME
			  << endl
			  << endl;

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				OutputBarElements(EleGrp);
				break;
			case ElementTypes::CST: // CST element
				OutputCSTElements(EleGrp);
				break;
			case ElementTypes::Q4: // Q4 element
				OutputQ4Elements(EleGrp);
				break;
			case ElementTypes::Q8: // Q8 element
				OutputQ8Elements(EleGrp);
				break;
			case ElementTypes::B21EB: // B21EB element, Euler-Bernoulli
				OutputB21EBElements(EleGrp);
				break;
			case ElementTypes::B31: // B31 element
				OutputB31Elements(EleGrp);
				break;
			case ElementTypes::Tet4: // Tet4 element
				OutputTet4Elements(EleGrp);
				break;
			case ElementTypes::H8: // H8 element
				OutputH8Elements(EleGrp);
				break;
			case ElementTypes::T6: // T6 element
				OutputT6Elements(EleGrp);
				break;
		    default:
		        *this << ElementType << " has not been implemented yet." << endl;
		        break;
		}
	}
}

//	Output bar element data
void COutputter::OutputBarElements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL" << endl
		  << " NUMBER     MODULUS          AREA" << endl
		  << "               E              A" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

//	Output CST element data
void COutputter::OutputCSTElements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S        POISSON'S        THICK        PLANE-" << endl
		  << " NUMBER     MODULUS          RATIO          -NESS        STRAIN" << endl
		  << "               E              nu              t           FLAG" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE       NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I          J        K       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

//	Output Q4 element data
void COutputter::OutputQ4Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		<< endl
		<< endl;

	*this << "  SET       YOUNG'S        POISSON'S        THICK        PLANE-" << endl
		<< " NUMBER     MODULUS          RATIO          -NESS        STRAIN" << endl
		<< "               E              nu              t           FLAG" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE       NODE     	NODE     NODE       MATERIAL" << endl
		<< " NUMBER-N      I          J        K        L       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;
}

//	Output Q8 element data
void COutputter::OutputQ8Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << setw(5) << NUMMAT<< endl
		<< endl;

	*this << "  SET       YOUNG'S        POISSON'S        THICK        PLANE-" << endl
		<< " NUMBER     MODULUS          RATIO          -NESS        STRAIN" << endl
		<< "               E              nu              t           FLAG" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE       NODE     	NODE      NODE      NODE       NODE        NODE        NODE        MATERIAL" << endl
		<< " NUMBER-N      I          J         K         L         M          N           O           P          SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;
	
}

//	Output Tet4 element data
void COutputter::OutputTet4Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << setw(5) << NUMMAT << endl
		<< endl;

	*this << "  SET       YOUNG'S        POISSON'S       " << endl
		<< " NUMBER     MODULUS          RATIO          " << endl
		<< "               E              nu            " << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE       NODE     	NODE      NODE        MATERIAL" << endl
		  << " NUMBER-N      I          J         K         L          SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;

}

//	Output H8 element data
void COutputter::OutputH8Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << setw(5) << NUMMAT << endl
		<< endl;

	*this << "  SET       YOUNG'S        POISSON'S       " << endl
		<< " NUMBER     MODULUS          RATIO          " << endl
		<< "               E              nu            " << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE       NODE     	NODE      NODE      NODE       NODE        NODE        NODE        MATERIAL" << endl
		  << " NUMBER-N      I          J         K         L         M          N           O           P          SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;

}

//	Output B21EB element data
void COutputter::OutputB21EBElements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		<< endl
		<< endl;

	*this << "  SET       YOUNG'S        CROSS-SECTION      MOMENT OF " << endl
		<< " NUMBER     MODULUS          	AREA          	INERTIA  " << endl
		<< "               E              	  A              I      " << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE       NODE        MATERIAL" << endl
		<< " NUMBER-N      I          J      SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;
	
}

//	Output B31 element data
void COutputter::OutputB31Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		<< endl
		<< endl;

	*this << "  SET       YOUNG'S       SHEAR       CROSS-SECTION      MOMENT OF      MOMENT OF      MOMENT OF      CORRECTION" << endl
		<< " NUMBER     MODULUS     MODULUS          AREA          	INERTIA-y      INERTIA-z      TORSION      	FACTOR  " << endl
		<< "               E            G             A              Iy		Iz		J		k      " << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE       NODE        MATERIAL" << endl
		<< " NUMBER-N      I          J      SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;
	
}

//	Output T6 element data
void COutputter::OutputT6Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		<< endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << setw(5) << NUMMAT<< endl
		<< endl;

	*this << "  SET       YOUNG'S        POISSON'S        THICK        PLANE-" << endl
		<< " NUMBER     MODULUS          RATIO          -NESS        STRAIN" << endl
		<< "               E              nu              t           FLAG" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset + 1;
		ElementGroup.GetMaterial(mset).Write(*this);
	}

	*this << endl << endl
		<< " E L E M E N T   I N F O R M A T I O N" << endl;

	*this << " ELEMENT     NODE       NODE     	NODE      NODE      NODE       NODE        MATERIAL" << endl
		<< " NUMBER-N      I          J         K         L         M          N          SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele + 1;
		ElementGroup[Ele].Write(*this);
	}

	*this << endl;
	
}


//	Print load data
void COutputter::OutputLoadInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];
		/*cout << lcase << endl;*/
		switch (LoadData->LoadCaseType_)	
		{
		case 1:	// All concentrated loads in node points
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
			*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << "    NODE       DIRECTION      LOAD" << endl
				<< "   NUMBER                   MAGNITUDE" << endl;
			break;	
		case 2:	// All concentrated loads in inner point of CST element(Dirac delta function)
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << LoadData->LoadCaseType_ << endl;
			*this << "     NUMBER OF CONCENTRATED LOADS IN INNER POINT OF CST ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << "     X-       		Y-      	 X-LOAD      	 Y-LOAD" << endl
				<< "   COORDINATE		COORDINATE		MAGNITUDE		MAGNITUDE" << endl;
			break;	
		case 3:	// All body forces of CST element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << LoadData->LoadCaseType_ << endl;
			*this << "     NUMBER OF BODY FORCES OF CST ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT      BX1-LOAD      	BY1-LOAD      	BX2-LOAD      	BY2-LOAD      	BX3-LOAD      	BY3-LOAD" << endl
				<< " NUMBER	MAGNITUDE	  MAGNITUDE	  MAGNITUDE	  MAGNITUDE	  MAGNITUDE	  MAGNITUDE" << endl;
			break;
		case 4:	// All surface forces of CST element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << LoadData->LoadCaseType_ << endl;
			*this << "     NUMBER OF BODY FORCES OF CST ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT       1-NODE       2-NODE     	TX1-LOAD      	TY1-LOAD      	TX2-LOAD      	TY2-LOAD" << endl
				<< " NUMBER	       NUMBER       NUMBER      MAGNITUDE	MAGNITUDE	MAGNITUDE	 MAGNITUDE" << endl;
			break;
		case 5:	// All body forces of Q4 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << LoadData->LoadCaseType_ << endl;
			*this << "     NUMBER OF BODY FORCES OF Q4 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT      BX1-LOAD      	BY1-LOAD      	BX2-LOAD      	BY2-LOAD      	BX3-LOAD      	BY3-LOAD      	BX4-LOAD      	BY4-LOAD" << endl
				<< " NUMBER	MAGNITUDE	  MAGNITUDE	  MAGNITUDE	  MAGNITUDE	  MAGNITUDE	  MAGNITUDE	  MAGNITUDE	  MAGNITUDE" << endl;
			break;
		case 6:	// All surface forces of Q4 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << LoadData->LoadCaseType_ << endl;
			*this << "     NUMBER OF BODY FORCES OF Q4 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT       1-NODE       2-NODE     	TX1-LOAD      	TY1-LOAD      	TX2-LOAD      	TY2-LOAD" << endl
				<< " NUMBER	       NUMBER       NUMBER      MAGNITUDE	MAGNITUDE	MAGNITUDE	 MAGNITUDE" << endl;
			break;
		case 7:	// All body forces of Q8 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << LoadData->LoadCaseType_ << endl;
			*this << "     BODY FORCES OF Q8 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT      BX1-LOAD      BY1-LOAD      BX2-LOAD      BY2-LOAD      BX3-LOAD      BY3-LOAD      BX4-LOAD      BY4-LOAD      BX5-LOAD      BY5-LOAD      BX6-LOAD      BY6-LOAD      BX7-LOAD      BY7-LOAD      BX8-LOAD      BY8-LOAD" << endl
				<< " NUMBER       MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE" << endl;
			break;
		case 8:	// All surface forces of Q8 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << LoadData->LoadCaseType_ << endl;
			*this << "     SURFACE  FORCES OF Q8 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT       1-NODE        2-NODE        3-NODE        TX1-LOAD       TY1-LOAD       TX2-LOAD       TY2-LOAD       TX3-LOAD       TY3-LOAD" << endl
				<< " NUMBER        NUMBER        NUMBER        NUMBER        MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE" << endl;
			break;
		case 9:
			break;
		case 10:
			break;
		case 11:	// All body forces of H8 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << LoadData->LoadCaseType_ << endl;
			*this << "     BODY FORCES OF H8 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT      BX1-LOAD      BY1-LOAD      BZ1-LOAD      BX2-LOAD      BY2-LOAD      BZ2-LOAD      BX3-LOAD      BY3-LOAD      BZ3-LOAD      BX4-LOAD      BY4-LOAD      BZ4-LOAD      BX5-LOAD      BY5-LOAD      BZ5-LOAD      BX6-LOAD      BY6-LOAD      BZ6-LOAD      BX7-LOAD      BY7-LOAD      BZ7-LOAD      BX8-LOAD      BY8-LOAD      BZ8-LOAD " << endl
				<< " NUMBER       MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE" << endl;
			break;
		case 12:	// All surface forces of H8 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << LoadData->LoadCaseType_ << endl;
			*this << "     SURFACE  FORCES OF H8 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT       1-NODE        2-NODE        3-NODE        4-NODE        TX1-LOAD       TY1-LOAD       TZ1-LOAD       TX2-LOAD       TY2-LOAD       TZ2-LOAD       TX3-LOAD       TY3-LOAD       TZ3-LOAD       TX4-LOAD       TY4-LOAD       TZ4-LOAD" << endl
				<< " NUMBER        NUMBER        NUMBER        NUMBER        NUMBER        MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE" << endl;
			break;
		case 13:	// All surface forces of S8R5 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
			*this << "     BODY FORCES OF S8R5 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT      BX1-LOAD      BY1-LOAD      BZ1-LOAD      BX2-LOAD      BY2-LOAD      BZ2-LOAD      BX3-LOAD      BY3-LOAD      BZ3-LOAD      BX4-LOAD      BY4-LOAD      BZ4-LOAD      BX5-LOAD      BY5-LOAD      BZ5-LOAD      BX6-LOAD      BY6-LOAD      BZ6-LOAD      BX7-LOAD      BY7-LOAD      BZ7-LOAD      BX8-LOAD      BY8-LOAD      BZ8-LOAD " << endl
				<< " NUMBER       MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE" << endl;
			break;
		case 14:	// All body forces of S8R5 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
			*this << "     BODY FORCES OF S8R5 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT      BX1-LOAD      BY1-LOAD      BZ1-LOAD      BX2-LOAD      BY2-LOAD      BZ2-LOAD      BX3-LOAD      BY3-LOAD      BZ3-LOAD      BX4-LOAD      BY4-LOAD      BZ4-LOAD      BX5-LOAD      BY5-LOAD      BZ5-LOAD      BX6-LOAD      BY6-LOAD      BZ6-LOAD      BX7-LOAD      BY7-LOAD      BZ7-LOAD      BX8-LOAD      BY8-LOAD      BZ8-LOAD " << endl
				<< " NUMBER       MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE" << endl;
			break;
		case 15:	// All surface forces of Mindlin-Reissner Plate element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
			*this << "     BODY FORCES OF Mindlin-Reissner Plate ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT      BZ1-LOAD      BZ2-LOAD      BZ3-LOAD      BZ4-LOAD " << endl
				<< " NUMBER       MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE" << endl;
			break;
		case 16:	// All surface forces of Plate element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
			*this << "     BODY FORCES OF Plate ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT      BZ1-LOAD      BZ2-LOAD      BZ3-LOAD      BZ4-LOAD " << endl
				<< " NUMBER       MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE" << endl;
			break;
		case 17:	// All body forces of Tet4 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
			*this << "     BODY FORCES OF Tet4 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT      BX1-LOAD      BY1-LOAD      BZ1-LOAD      BX2-LOAD      BY2-LOAD      BZ2-LOAD      BX3-LOAD      BY3-LOAD      BZ3-LOAD      BX4-LOAD      BY4-LOAD      BZ4-LOAD " << endl
				<< " NUMBER       MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE" << endl;
			break;
		case 18:	// All surface forces of Tet4 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
			*this << "     SURFACE  FORCES OF Tet4 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT       1-NODE        2-NODE        3-NODE        TX1-LOAD       TY1-LOAD       TZ1-LOAD       TX2-LOAD       TY2-LOAD       TZ2-LOAD       TX3-LOAD       TY3-LOAD       TZ3-LOAD" << endl
				<< " NUMBER        NUMBER        NUMBER        NUMBER        NUMBER        MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE" << endl;
			break;
		case 19:	// All body forces of T6 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
			*this << "     BODY FORCES OF T6 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT      BX1-LOAD      BY1-LOAD      BX2-LOAD      BY2-LOAD      BX3-LOAD      BY3-LOAD      BX4-LOAD      BY4-LOAD      BX5-LOAD      BY5-LOAD      BX6-LOAD      BY6-LOAD" << endl
				<< " NUMBER       MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE     MAGNITUDE" << endl;
			break;
		case 20:	// All surface forces of T6 element
			*this << setiosflags(ios::scientific);
			*this << " L O A D   C A S E   D A T A" << endl
				<< endl;

			*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
			*this << "     SURFACE  FORCES OF T6 ELEMENTS . =" << setw(6) << LoadData->nloads << endl
				<< endl;
			*this << " ELEMENT       1-NODE        2-NODE        3-NODE        TX1-LOAD       TY1-LOAD       TX2-LOAD       TY2-LOAD       TX3-LOAD       TY3-LOAD" << endl
				<< " NUMBER        NUMBER        NUMBER        NUMBER        MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE      MAGNITUDE" << endl;
			break;
		default:
			std::cerr << "LodaCase " << LoadData->LoadCaseType_ << " not available. See COutputter::OutputLoadInfo." << std::endl;
			exit(5);
			break;
		}		

		LoadData->Write(LoadData->LoadCaseType_, *this);

		*this << endl;
	}
}

//	Print nodal displacement
void COutputter::OutputNodalDisplacement()
{
	CDomain* FEMData = CDomain::GetInstance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	*this << setiosflags(ios::scientific);

	*this << " D I S P L A C E M E N T S" << endl
		  << endl;
	*this << "  NODE            X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT    X-ROTATION    Y-ROTATION    Z-ROTATION" << endl;

	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		NodeList[np].WriteNodalDisplacement(*this, Displacement);

	*this << endl;
}

//	Calculate stresses
void COutputter::OutputElementStress()
{
	CDomain* FEMData = CDomain::GetInstance();

	double* Displacement = FEMData->GetDisplacement();

	unsigned int NUMEG = FEMData->GetNUMEG();

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)
			  << EleGrpIndex + 1 << endl
			  << endl;

		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
			{
				*this << "  ELEMENT             FORCE            STRESS" << endl
					<< "  NUMBER" << endl;

				double stress;

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(&stress, Displacement);

					CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22) << stress * material.Area << setw(18)
						<< stress << endl;
				}

				*this << endl;

				break;
			}
				
			case ElementTypes::CST: // CST element
			{
				*this << "  ELEMENT        Sigma_xx          sigma_yy              Sigma_xy" << endl
					<< "  NUMBER" << endl;

				double stress_CST[3];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(stress_CST, Displacement);

					*this << setw(5) << Ele + 1 << setw(22) << stress_CST[0] << setw(18)
						<< stress_CST[1] << setw(22) << stress_CST[2] << endl;
				}

				*this << endl;

				break;
			}	

			case ElementTypes::Q4: // Q4 element
			{
			    *this << "ELEMENT NUMBER	  	x-coord	    	    y-coord           Sigma_xx             sigma_yy               Sigma_xy" << endl
					<< "GAUSSPOINT NUMBER" << endl;				

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{	
					CElement& Element = EleGrp[Ele];
					double coord_stress_Q4[20] = {0};  // 2x2 Gauss Points
					Element.ElementStress(coord_stress_Q4, Displacement);
					for (unsigned int i = 0; i < 4; i++)
					{
						*this << setw(8) << Ele + 1 << "-" << setw(1) << i+1 << setw(22) << coord_stress_Q4[5*i] << setw(22) << coord_stress_Q4[5*i+1] 
					        << setw(18) << coord_stress_Q4[5*i+2] << setw(22) << coord_stress_Q4[5*i+3] << setw(22) << coord_stress_Q4[5*i+4] << endl;
					}					
				}

				*this << endl;

				break;
			}	

			case ElementTypes::Q8: // Q8 element
			{
				*this << "ELEMENT NUMBER	  	x-coord	    	    y-coord           Sigma_xx             sigma_yy               Sigma_xy" << endl
					<< "GAUSSPOINT NUMBER" << endl;

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					double coord_stress_Q8[20] = {0}; // 2x2 Gauss Points
					Element.ElementStress(coord_stress_Q8, Displacement);
					for (unsigned int i = 0; i < 4; i++)
					{
						*this << setw(8) << Ele + 1 << "-" << setw(1) << i + 1 << setw(22) << coord_stress_Q8[5*i] << setw(22) << coord_stress_Q8[5*i+1]
							<< setw(18) << coord_stress_Q8[5*i+2] << setw(22) << coord_stress_Q8[5*i+3] << setw(22) << coord_stress_Q8[5*i+4] << endl;
					}
				}

				*this << endl;

				break;
			}

			case ElementTypes::B21EB: // B21EB element, Euler-Bernoulli
			{
				*this << "ELEMENT NUMBER	  	x-coord	    	    y-coord           Moment             Shear Force" << endl
					<< "GAUSSPOINT NUMBER" << endl;

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					double coord_stress_B21EB[8] = {0};// 2 Gauss Points
					Element.ElementStress(coord_stress_B21EB, Displacement);
					for (unsigned int i = 0; i < 2; i++)
					{
						*this << setw(8) << Ele + 1 << "-" << setw(1) << i + 1 << setw(22) << coord_stress_B21EB[4*i] << setw(22) << coord_stress_B21EB[4*i+1]
							<< setw(18) << coord_stress_B21EB[4*i+2] << setw(22) << coord_stress_B21EB[4*i+3]  << endl;
					}
				}

				*this << endl;

				break;
			}
			
			case ElementTypes::Tet4: // Tet4 element
			{
				*this << "ELEMENT NUMBER   Sigma_xx           Sigma_yy           Sigma_zz           Sigma_xy           Sigma_xz           Sigma_yz" << endl
					<< "  NUMBER" << endl;
				
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					double stress_Tet4[6] = {0};
					Element.ElementStress(stress_Tet4, Displacement);

					*this << setw(5) << Ele + 1 << setw(22) << stress_Tet4[0] << setw(18)
						<< stress_Tet4[1] << setw(20) << stress_Tet4[2] << setw(18) << stress_Tet4[3] << setw(18) << stress_Tet4[4] << setw(18) << stress_Tet4[5] << endl;
				}

				*this << endl;

				break;
			}

			case ElementTypes::H8: // H8 element
			{
				*this << "ELEMENT NUMBER   x-coord      y-coord      z-coord           Sigma_xx           Sigma_yy           Sigma_zz           Sigma_xy           Sigma_xz           Sigma_yz" << endl
					<< "GAUSSPOINT NUMBER" << endl;
			
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					double coord_stress_H8[72] = {0};// 2x2x2 Gauss Points
					Element.ElementStress(coord_stress_H8, Displacement);
					for (unsigned int i = 0; i < 8; i++)
					{
						*this << setw(3) << Ele + 1 << "-" << setw(1) << i + 1 << setw(18) << coord_stress_H8[9*i] 
							<< setw(15) << coord_stress_H8[9*i+1] << setw(15) << coord_stress_H8[9*i+2] 
							<< setw(18) << coord_stress_H8[9*i+3] << setw(18) << coord_stress_H8[9*i+4] 
							<< setw(18) << coord_stress_H8[9*i+5] << setw(18) << coord_stress_H8[9*i+6] 
							<< setw(18) << coord_stress_H8[9*i+7] << setw(18) << coord_stress_H8[9*i+8] << endl; //<< setw(8)<< 2 << endl;
					}
				}

				*this << endl;

				break;
			}

			case ElementTypes::T6: // T6 element
			{
				*this << "ELEMENT NUMBER          x-coord              y-coord           Sigma_xx             Sigma_yy             Sigma_xy" << endl
				<< "GAUSSPOINT NUMBER" << endl;
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					double coord_stress_T6[15] = {0};  // 3 Gauss Points × (x, y, σ_xx, σ_yy, σ_xy) = 3×5 = 15
					Element.ElementStress(coord_stress_T6, Displacement);
					for (unsigned int i = 0; i < 3; i++)  // 3 Gauss Points
					{
						*this << setw(8) << Ele + 1 << "-" << setw(1) << i + 1
						<< setw(22) << coord_stress_T6[5 * i]     // x
						<< setw(22) << coord_stress_T6[5 * i + 1] // y
						<< setw(18) << coord_stress_T6[5 * i + 2] // σ_xx
						<< setw(22) << coord_stress_T6[5 * i + 3] // σ_yy
						<< setw(22) << coord_stress_T6[5 * i + 4] // σ_xy
						<< endl;
					}
				}

				*this << endl;

				break;
			}

			default: // Invalid element type
				cerr << "*** Error *** Elment type " << ElementType
					<< " has not been implemented.\n\n";
		}
	}
}

//	Print total system data
void COutputter::OutputTotalSystemData()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "	TOTAL SYSTEM DATA" << endl
		  << endl;

	*this << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = " << FEMData->GetNEQ()
		  << endl
		  << "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = " << FEMData->GetStiffnessMatrix()->size()
		  << endl
		  << "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = " << FEMData->GetStiffnessMatrix()->GetMaximumHalfBandwidth()
		  << endl
		  << "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = " << FEMData->GetStiffnessMatrix()->size() / FEMData->GetNEQ() << endl
		  << endl
		  << endl;
}



#ifdef _DEBUG_

//	Print column heights for debuging
void COutputter::PrintColumnHeights()
{
	*this << "*** _Debug_ *** Column Heights" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int col = 0; col < NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << ColumnHeights[col];
	}

	*this << endl
		  << endl;
}

//	Print address of diagonal elements for debuging
void COutputter::PrintDiagonalAddress()
{
	*this << "*** _Debug_ *** Address of Diagonal Element" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	for (unsigned int col = 0; col <= NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << DiagonalAddress[col];
	}

	*this << endl
		  << endl;
}

//	Print banded and full stiffness matrix for debuging
void COutputter::PrintStiffnessMatrix()
{
	*this << "*** _Debug_ *** Banded stiffness matrix" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++)
	{
		*this << setw(14) << (*StiffnessMatrix)(i);

		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}
	}

	*this << endl
		  << endl;

	*this << "*** _Debug_ *** Full stiffness matrix" << endl;

	for (int I = 1; I <= NEQ; I++)
	{
		for (int J = 1; J <= NEQ; J++)
		{
			int J_new = (J > I) ? J : I;
			int I_new = (J > I) ? I : J;
			int H = DiagonalAddress[J_new] - DiagonalAddress[J_new - 1];
			if (J_new - I_new - H >= 0)
			{
				*this << setw(14) << 0.0;
			}
			else
			{
				*this << setw(14) << (*StiffnessMatrix)(I_new, J_new);
			}
		}

		*this << endl;
	}

	*this << endl;
}

//	Print displacement vector for debuging
void COutputter::PrintDisplacement()
{
	*this << "*** _Debug_ *** Displacement vector" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < NEQ; i++)
	{
		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}

		*this << setw(14) << Force[i];
	}

	*this << endl
		  << endl;
}

#endif
