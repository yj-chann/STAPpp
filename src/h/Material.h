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

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

public:

//! Virtual deconstructor
    virtual ~CMaterial() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output) = 0;

};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

//!	Material class for CST element
class CCSTMaterial : public CMaterial
{
public:

	double nu;        //!< Poisson's ratio of a CST element
	double t;	      //!< Element thickness of a CST element
	int plane_strain; //!< Plane Strain or not

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

//!	Material class for Q4 element
class CQ4Material : public CMaterial
{
public:

	double nu;        //!< Poisson's ratio of a Q4 element
	double t;	      //!< Element thickness of a Q4 element
	int plane_strain; //!< Plane Strain or not

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

//!	Material class for Q8 element
class CQ8Material : public CMaterial
{
public:

	double nu;        //!< Poisson's ratio of a Q8 element
	double t;	      //!< Element thickness of a Q8 element
	int plane_strain; //!< Plane Strain or not

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

//!	Material class for H8 element
class CH8Material : public CMaterial
{
public:

	double nu;        //!< Poisson's ratio of a H8 element


public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//!	Write material data to Stream
	virtual void Write(COutputter& output);
};


//!	Material class for B21EB element
class CB21EBMaterial : public CMaterial
{
public:

	double Area;	  //!< Sectional area of a bar element
	double I;	  	  //!< Moment of inertia of cross-section

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

//!	Material class for B31 element
class CB31Material : public CMaterial
{
public:

	double nu;        //!< Poisson's ratio
	double G;  		  //!< Shear modulus
	double Area;	  //!< Sectional area of a bar element
	double Iy;	  	  //!< Moment of inertia of the x'z' cross-section
	double Iz;	  	  //!< Moment of inertia of the x'y' cross-section
	double J;	  	  //!< Torsional Moment of Inertia
	double k;		  //!< Correction factor (shear area ratio)
	unsigned int direction; //!< 1: YOZ plane; 2: ZOX plane; 3: XOY plane

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

//!	Material class for S8R5 element
class CS8R5Material : public CMaterial
{
public:

	double nu;        //!< Poisson's ratio
	double G;  		  //!< Shear modulus
	double t;	      //!< Element thickness of a Q8 element
	double k;		  //!< Correction factor (shear area ratio)
	double v3[3];	  //!< Unit vector normal to the midsurface

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

//!	Material class for Mindlin-Reissner Plate element
class CMPMaterial : public CMaterial
{
public:

	double nu;        //!< Poisson's ratio
	double G;  		  //!< Shear modulus
	double D_0;	      //!< Bending stiffness
	double t;	      //!< Element thickness of a Q8 element
	double k;		  //!< Correction factor (shear area ratio)

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

	//!	Write material data to Stream
	virtual void Write(COutputter& output);
};




