/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "T6.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// 构造函数
CT6::CT6()
{
    NEN_ = 6;    // 每个单元有6个节点
    nodes_ = new CNode*[NEN_];
    
    ND_ = 12;    // 6节点 × 2自由度/节点
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

// 析构函数
CT6::~CT6()
{
}

// 从输入流读取单元数据
bool CT6::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;    // 材料属性集编号
    unsigned int N[6];    // 6个节点编号(前3个是角节点，后3个是边中点)

    Input >> N[0] >> N[1] >> N[2] >> N[3] >> N[4] >> N[5] >> MSet;
    ElementMaterial_ = dynamic_cast<CT6Material*>(MaterialSets) + MSet - 1;
    
    for(unsigned int i=0; i<6; i++)
        nodes_[i] = &NodeList[N[i] - 1];

    return true;
}

// 将单元数据写入输出流
void CT6::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(11) << nodes_[1]->NodeNumber
           << setw(11) << nodes_[2]->NodeNumber
           << setw(11) << nodes_[3]->NodeNumber
           << setw(11) << nodes_[4]->NodeNumber
           << setw(11) << nodes_[5]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

// 计算形函数
void CT6::ShapeFunction(double xi, double eta, double N[6])
{
    double zeta = 1.0 - xi - eta;
    
    // 角节点形函数
    N[0] = xi*(2.0*xi - 1.0);        // N1 = ξ1(2ξ1-1)
    N[1] = eta*(2.0*eta - 1.0);      // N2 = ξ2(2ξ2-1)
    N[2] = zeta*(2.0*zeta - 1.0);    // N3 = ξ3(2ξ3-1)
    
    // 边中点形函数
    N[3] = 4.0*xi*eta;               // N4 = 4ξ1ξ2
    N[4] = 4.0*eta*zeta;             // N5 = 4ξ2ξ3
    N[5] = 4.0*zeta*xi;              // N6 = 4ξ3ξ1
}

// 计算形函数对自然坐标的导数
void CT6::ShapeFunctionDerivative(double xi, double eta, double dNdxi[6][2])
{
    double zeta = 1.0 - xi - eta;
    
    // ∂N/∂ξ1
    dNdxi[0][0] = 4.0*xi - 1.0;      // ∂N1/∂ξ1
    dNdxi[1][0] = 0.0;               // ∂N2/∂ξ1
    dNdxi[2][0] = 4.0*xi + 4.0*eta - 3.0; // ∂N3/∂ξ1
    dNdxi[3][0] = 4.0*eta;           // ∂N4/∂ξ1
    dNdxi[4][0] = -4.0*eta;          // ∂N5/∂ξ1
    dNdxi[5][0] = 4.0 - 4.0*eta - 8.0*xi; // ∂N6/∂ξ1
    
    // ∂N/∂ξ2
    dNdxi[0][1] = 0.0;               // ∂N1/∂ξ2
    dNdxi[1][1] = 4.0*eta - 1.0;     // ∂N2/∂ξ2
    dNdxi[2][1] = 4.0*xi + 4.0*eta - 3.0; // ∂N3/∂ξ2
    dNdxi[3][1] = 4.0*xi;            // ∂N4/∂ξ2
    dNdxi[4][1] = 4.0 - 8.0*eta - 4.0*xi; // ∂N5/∂ξ2
    dNdxi[5][1] = -4.0*xi;           // ∂N6/∂ξ2
}

void CT6::ElementShapeFunction(double (&N)[2][12], double xi, double eta)
{
    // T6单元的6个节点标准形函数
    double L1 = 1.0 - xi - eta;
    double L2 = xi;
    double L3 = eta;

    double N1 = L1 * (2 * L1 - 1);
    double N2 = L2 * (2 * L2 - 1);
    double N3 = L3 * (2 * L3 - 1);
    double N4 = 4 * L1 * L2;
    double N5 = 4 * L2 * L3;
    double N6 = 4 * L3 * L1;

    // 每个节点对应2个自由度，分别填充到 N[0][*] (x方向) 和 N[1][*] (y方向)
    for (int i = 0; i < 12; ++i)
        N[0][i] = N[1][i] = 0.0;

    N[0][0] = N1;  N[1][1] = N1;
    N[0][2] = N2;  N[1][3] = N2;
    N[0][4] = N3;  N[1][5] = N3;
    N[0][6] = N4;  N[1][7] = N4;
    N[0][8] = N5;  N[1][9] = N5;
    N[0][10] = N6; N[1][11] = N6;
}

// 计算单元刚度矩阵
void CT6::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    CT6Material* material_ = dynamic_cast<CT6Material*>(ElementMaterial_);
    double E = material_->E;
    double nu = material_->nu;
    double t = material_->t;
    int plane_strain = material_->plane_strain;

    if(plane_strain) {
        E = E / (1.0 - nu*nu);
        nu = nu / (1.0 - nu);
    }

    double k = E / (1.0 - nu*nu);
    double D[3][3] = {{k, k*nu, 0}, 
                      {k*nu, k, 0}, 
                      {0, 0, k*(1.0-nu)/2.0}};

    // 使用3点减缩积分方案
    double GaussPoints[3][3] = {
        {1.0/6.0, 1.0/6.0, 1.0/6.0},  // 权重为1/3
        {2.0/3.0, 1.0/6.0, 1.0/6.0},
        {1.0/6.0, 2.0/3.0, 1.0/6.0}
    };

    double Ke[12][12] = {0};
    
    for(unsigned int gp=0; gp<3; gp++) {
        double xi = GaussPoints[gp][0];
        double eta = GaussPoints[gp][1];
        double weight = GaussPoints[gp][2];
        
        double B[3][12] = {0};
        double detJ;
        ElementStrainFunction(B, &detJ, xi, eta);
        
        // 计算刚度矩阵贡献
        for(int i=0; i<12; i++) {
            for(int j=0; j<12; j++) {
                for(int k=0; k<3; k++) {
                    for(int l=0; l<3; l++) {
                        Ke[i][j] += t * B[k][i] * D[k][l] * B[l][j] * detJ * weight;
                    }
                }
            }
        }
    }

    // 将刚度矩阵打包成一维数组(上三角部分)
    int index = 0;
    for(int j=0; j<12; j++) {
        for(int i=j; i>=0; i--) {
            Matrix[index++] = Ke[i][j];
        }
    }
}

// 计算应变矩阵B和雅可比行列式
void CT6::ElementStrainFunction(double (&B)[3][12], double* detJ, double xi, double eta)
{
    double dNdxi[6][2] = {0};
    ShapeFunctionDerivative(xi, eta, dNdxi);
    
    // 节点坐标
    double X[6], Y[6];
    for(int i=0; i<6; i++) {
        X[i] = nodes_[i]->XYZ[0];
        Y[i] = nodes_[i]->XYZ[1];
    }
    
    // 计算雅可比矩阵
    double J[2][2] = {0};
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
            for(int k=0; k<6; k++) {
                J[i][j] += dNdxi[k][i] * (j==0 ? X[k] : Y[k]);
            }
        }
    }
    
    // 计算雅可比行列式
    *detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    
    // 计算雅可比逆矩阵
    double invJ[2][2];
    double invDet = 1.0 / (*detJ);
    invJ[0][0] =  J[1][1] * invDet;
    invJ[0][1] = -J[0][1] * invDet;
    invJ[1][0] = -J[1][0] * invDet;
    invJ[1][1] =  J[0][0] * invDet;
    
    // 计算形函数对物理坐标的导数
    double dNdx[6], dNdy[6];
    for(int i=0; i<6; i++) {
        dNdx[i] = invJ[0][0]*dNdxi[i][0] + invJ[0][1]*dNdxi[i][1];
        dNdy[i] = invJ[1][0]*dNdxi[i][0] + invJ[1][1]*dNdxi[i][1];
    }
    
    // 组装应变矩阵B
    for(int i=0; i<6; i++) {
        B[0][2*i]   = dNdx[i];  // εxx
        B[1][2*i+1] = dNdy[i];  // εyy
        B[2][2*i]   = dNdy[i];  // γxy
        B[2][2*i+1] = dNdx[i];
    }
}

// 计算单元应力
void CT6::ElementStress(double* stress, double* Displacement)
{
    CT6Material* material_ = dynamic_cast<CT6Material*>(ElementMaterial_);
    double E = material_->E;
    double nu = material_->nu;
    double t = material_->t;
    int plane_strain = material_->plane_strain;

    if(plane_strain) {
        E = E / (1.0 - nu*nu);
        nu = nu / (1.0 - nu);
    }

    double k = E / (1.0 - nu*nu);
    double D[3][3] = {{k, k*nu, 0}, 
                      {k*nu, k, 0}, 
                      {0, 0, k*(1.0-nu)/2.0}};

    // 获取单元位移向量
    double d[12] = {0};
    for(unsigned int i=0; i<12; i++) {
        if(LocationMatrix_[i]) {
            d[i] = Displacement[LocationMatrix_[i]-1];
        }
    }

    // 处理非齐次边界条件
    if(NonHomo_) {
        for(unsigned int i=0; i<12; i++) {
            if(!LocationMatrix_[i] && nodes_[i/2]->BC[i%2] != 0) {
                d[i] = nodes_[i/2]->BC[i%2];
            }
        }
    }

    // 使用3点减缩积分方案计算应力
    double GaussPoints[3][3] = {
        {1.0/6.0, 1.0/6.0, 1.0/6.0},
        {2.0/3.0, 1.0/6.0, 1.0/6.0},
        {1.0/6.0, 2.0/3.0, 1.0/6.0}
    };

    int index = 0;
    for(unsigned int gp=0; gp<3; gp++) {
        double xi = GaussPoints[gp][0];
        double eta = GaussPoints[gp][1];
        
        double B[3][12] = {0};
        double detJ;
        ElementStrainFunction(B, &detJ, xi, eta);
        
        // 计算应力 σ = D*B*d
        for(int i=0; i<3; i++) {
            stress[5*index+2+i] = 0.0;  // 初始化应力分量
            for(int j=0; j<3; j++) {
                for(int k=0; k<12; k++) {
                    stress[5*index+2+i] += D[i][j] * B[j][k] * d[k];
                }
            }
        }
        
        // 计算积分点坐标
        double N[6];
        ShapeFunction(xi, eta, N);
        stress[5*index] = 0.0;   // x坐标
        stress[5*index+1] = 0.0; // y坐标
        for(int i=0; i<6; i++) {
            stress[5*index] += N[i] * nodes_[i]->XYZ[0];
            stress[5*index+1] += N[i] * nodes_[i]->XYZ[1];
        }
        
        index++;
    }
}

//	Calculate element non-homogeneous essential boundary conditions
void CT6::ElementNonHomo(double* Matrix, double* NonForce)
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
    for (unsigned int i = 0; i < 6; i++)
    {
        d[2*i] = nodes_[i]->BC[0];
        d[2*i+1] = nodes_[i]->BC[1];
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
