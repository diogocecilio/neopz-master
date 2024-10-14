// SPDX-FileCopyrightText: 2024 <copyright holder> <email>
// SPDX-License-Identifier: Apache-2.0

#include "tpzklinterpolationspace.h"


TPZKLInterpolationSpace::TPZKLInterpolationSpace(): TPZInterpolationSpace()
{
    //fPreferredOrder = -1;
}

TPZKLInterpolationSpace::TPZKLInterpolationSpace(TPZCompMesh &mesh, const TPZKLInterpolationSpace &copy): TPZInterpolationSpace(mesh, copy)
{
    //fPreferredOrder = copy.fPreferredOrder;
}


TPZKLInterpolationSpace::TPZKLInterpolationSpace(TPZCompMesh &mesh, const TPZKLInterpolationSpace &copy, std::map<int64_t,int64_t> &gl2lcElMap): TPZInterpolationSpace(mesh, copy, gl2lcElMap)
{
    //fPreferredOrder = copy.fPreferredOrder;
}

TPZKLInterpolationSpace::TPZKLInterpolationSpace(TPZCompMesh &mesh, TPZGeoEl *gel): TPZInterpolationSpace(mesh,gel)
{
    //fPreferredOrder = mesh.GetDefaultOrder();
}

TPZKLInterpolationSpace::~TPZKLInterpolationSpace()
{

}

void TPZKLInterpolationSpace::CalcStiffC(TPZCompEl *jel, TPZElementMatrixT<REAL> &ce)
{
    TPZKLInterpolationSpace *intspace2  = dynamic_cast <TPZKLInterpolationSpace *>(jel);
    TPZFMatrix<REAL> elmatt,elmat;
    auto* material =dynamic_cast<TPZKarhunenLoeveMat *>(this->Material());
    TPZElementMatrixT<REAL> fe,cet;
    this->InitializeElementMatrix(ce,fe);
    this->InitializeElementMatrix(cet,fe);
    if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic

    TPZMaterialDataT<REAL> data1,data2;
    this->InitMaterialData(data1);
    intspace2->InitMaterialData(data2);


    data1.p = this->MaxOrder();
    data2.p = intspace2->MaxOrder();

    int dim = Dimension();
    TPZManVector<REAL,3> intpoint1(dim,0.);
    TPZManVector<REAL,3> intpoint2(dim,0.);
    REAL weight1 = 0.,weight2=0.;

    TPZAutoPointer<TPZIntPoints> intrule = GetIntegrationRule().Clone();
    TPZAutoPointer<TPZIntPoints> intrule2 = GetIntegrationRule().Clone();

// 	int order = material->IntegrationRuleOrder(data1.p);
// 	int order2 = material->IntegrationRuleOrder(data2.p);
// 	TPZManVector<int,3> intorder(dim,order);
// 	intrule->SetOrder(intorder);
// 	TPZManVector<int,3> intorder2(dim,order2);
// 	intrule2->SetOrder(intorder2);
// 	//std::cout << " intrulepoints " << intrule->NPoints() <<std::endl;

    int intrulepoints = intrule->NPoints();
    for(int int_ind = 0; int_ind < intrulepoints; int_ind++){

        intrule->Point(int_ind,intpoint1,weight1);
        data1.intLocPtIndex = int_ind;
        this->ComputeRequiredData(data1, intpoint1);
        weight1 *= fabs(data1.detjac);
        for(int int_jnd=0;int_jnd<intrulepoints;int_jnd++)
        {
          intrule2->Point(int_jnd,intpoint2,weight2);
          data2.intLocPtIndex = int_jnd;
          //this->ComputeRequiredData(data2, intpoint2);
          intspace2->ComputeRequiredData(data2, intpoint2);
          weight2 *= fabs(data2.detjac);
          material->ContributeC(data1,data2, weight1,weight2, cet.fMat);
          //cet.fMat.Print(std::cout);
          ce.fMat+=cet.fMat;
        }

}//loop over integratin points
}


void TPZKLInterpolationSpace::CalcStiffB(TPZElementMatrixT<REAL> &be)
{
    auto* material =
    dynamic_cast<TPZKarhunenLoeveMat *>(this->Material());
    TPZElementMatrixT<REAL> ef,bet;
    this->InitializeElementMatrix(be,ef);
    this->InitializeElementMatrix(bet,ef);
    if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic


    TPZMaterialDataT<REAL> data;
    //data.p = this->MaxOrder();
    this->InitMaterialData(data);
    data.p = this->MaxOrder();

    int dim = Dimension();
    TPZManVector<REAL,3> intpoint(dim,0.);
    REAL weight = 0.;

    TPZAutoPointer<TPZIntPoints> intrule = GetIntegrationRule().Clone();
//     int order = material->IntegrationRuleOrder(data.p);
// 	TPZManVector<int,3> intorder(dim,order);
// 	intrule->SetOrder(intorder);
    int intrulepoints = intrule->NPoints();
    for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){
        intrule->Point(int_ind,intpoint,weight);
        data.intLocPtIndex = int_ind;
        this->ComputeRequiredData(data, intpoint);
        weight *= fabs(data.detjac);
        material->ContributeB(data, weight,bet.fMat);
        be.fMat+=bet.fMat;
    }//loop over integratin points
}
