//
//  TPZElasticResponse.h
//  pz
//
//  Created by Erick Slis Raggio Santos on 04/07/2009.
//


#ifndef TPZElasticResponse_H
#define TPZElasticResponse_H

#include "TPZTensor.h"
#include "pzreal.h"


class TPZElasticResponse : public TPZSavable {
    
protected:
    /// First Lamé parameter
    REAL m_lambda;

    /// Second Lamé parameter
    REAL m_mu;

    REAL m_E;

    REAL m_nu;

    /// Reference strain at zero stress state.
    TPZTensor<REAL> m_epsilon_star;

    /// Reference stress at zero strain state.
    TPZTensor<REAL> m_sigma_star;
    
public:
    
    /**
     A unique class identifier
     
     @return The class identifier as integer
     */
    int ClassId() const override;
    
    /**
     Default constructor
     */
    TPZElasticResponse();
    
    /**
     Copy constructor
     */
    TPZElasticResponse(const TPZElasticResponse & other);
    
    /**
     Assignment constructor
     */
    TPZElasticResponse & operator=(const TPZElasticResponse & other);
    
    /**
     Write (persistency)
     
     @param buf The TPZStream object
     @param withclassid The class identifier
     */
    void Write(TPZStream &buf, int withclassid) const override;
    
    /**
     Read (persistency)
     
     @param buf The TPZStream object
     @param context pointer to the associated object
     */
    void Read(TPZStream &buf, void *context) override;
    
    
    /**
     Class name
     
     @return constant char with the class name
     */
    const char * Name() const;
    
    /**
     Print
     
     @param out ostream object to write the output
     */
    void Print(std::ostream & out) const;
    


    void ComputeStress(const TPZTensor<STATE> & epsilon, TPZTensor<STATE> & sigma) const;


    void ComputeStrain(const TPZTensor<STATE> & sigma, TPZTensor<STATE> & epsilon) const;

    /**
     Incremental constitutive relation in Voigt notation
     
     @param De Return the De operator
     */
    void De(TPZFMatrix<STATE> & DeMat)const;

    void InverseDe(TPZFMatrix<STATE> & DeMat) const;

    /**
     Set elastic parameters using engineering data, i.e. Young modulus and Poisson ratio
     
     @param Eyoung Young modulus
     @param Poisson Poisson ratio
     */
    void SetEngineeringData(REAL Eyoung, REAL Poisson);
    
    
    /**
     Set elastic parameters using Lamé data
     
     @param lambda First Lamé parameter
     @param mu Second Lamé parameter (Shear modulus)
     */
    void SetLameData();
    
    /**
     Access to the first Lamé parameter
     
     @return The first Lamé parameter
     */
    REAL Lambda() const;
    
    /**
     Access to the bulk modulus
     
     @return The bulk modulus
     */
    REAL K() const;
    
    /**
     Access to the second Lamé parameter
     
     @return The Second Lamé parameter
     */
    REAL Mu() const;
    
    /**
     Access to the shear modulus
     
     @return The shear modulus
     */
    REAL G() const;
    
    /**
     Access to the Young modulus
     
     @return The Young modulus
     */
    REAL E() const;
    
    /**
     Access to the Poisson ratio
     
     @return The Poisson ratio
     */
    REAL Poisson() const;
    
    /// Set the reference strain
    void SetReferenceStrainData(TPZTensor<REAL> & eps_star);
    
    /// Get the reference strain
    TPZTensor<REAL> & ReferenceStrainData();
    
    /// Set the reference stress
    void SetReferenceStressData(TPZTensor<REAL> & sigma_star);
    
    /// Get the reference strain
    TPZTensor<REAL> & ReferenceStressData();
};

#endif /* TPZElasticResponse_h */
