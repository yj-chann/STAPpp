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

#include "Element.h"

using namespace std;

//! T6 element class (6-node quadratic triangular element)
class CT6 : public CElement
{
public:
    //! Constructor
    CT6();

    //! Destructor
    ~CT6();

    //! Read element data from stream Input
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList) override;

    //! Write element data to stream
    virtual void Write(COutputter& output) override;

    //! Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix) override;

    //! Calculate element stress
    virtual void ElementStress(double* stress, double* Displacement) override;

    //! Calculate element non-homogeneous essential boundary conditions
    virtual void ElementNonHomo(double* Matrix, double* NonForce) override;

    //! Calculate T6 element shape functions
    void ShapeFunction(double xi, double eta, double N[6]);

    //! Calculate derivatives of T6 shape functions with respect to natural coordinates
    void ElementShapeFunction(double (&N)[2][12], double xi, double eta);

    //! Calculate derivatives of T6 shape functions with respect to natural coordinates
    void ShapeFunctionDerivative(double xi, double eta, double dNdxi[6][2]);

    //! Calculate strain-displacement matrix B and Jacobian determinant
    void ElementStrainFunction(double (&B)[3][12], double* detJ, double xi, double eta);
};