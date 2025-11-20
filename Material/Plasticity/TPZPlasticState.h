// $Id: TPZPlasticState.h,v 1.8 2009-07-19 05:47:18 erick Exp $

#ifndef TPZPLASTICSTATE_H
#define TPZPLASTICSTATE_H

#include "fadType.h"
#include "TPZTensor.h"
#include <iostream>

/**
 * This class holds the complete set of state variables to define a valid elastoplastic strain state
 */
template <class T>
class TPZPlasticState : public TPZSavable {

public:

    /// Tensors representing the total and plastic strain states
    TPZTensor<T> m_eps_t, m_eps_p;

    /// Plastic volumetric hardeing variable
    T m_hardening;

    /// Identifier for the regime of the material behaviour
    int m_m_type;

	TPZVec<T> fmatprop;

    TPZVec<T> fmatpropinit;

    TPZVec<T> fflux;

    T fpressure;

    TPZVec<T> fdPorePressure;

    TPZVec<T> fSolU;

    TPZTensor<T> fGradSolU;
public:

    /// Default constructor - all values set to zero
    TPZPlasticState(): m_eps_t(), m_eps_p(), m_hardening(T(0.)), m_m_type(0),fpressure(T(0.) ),fmatprop(),fflux(),fdPorePressure(),fSolU(),fGradSolU() { }

    /// Constructor enabling predefinition of hardening
    TPZPlasticState(const T & hardening):m_eps_t(T(0.)), m_eps_p(T(0.)), m_hardening(hardening), m_m_type(0),fpressure(T(0.) ),fmatprop(),fflux(),fdPorePressure(),fSolU(),fGradSolU() { }

    /// Copy constructor
    TPZPlasticState(const TPZPlasticState<T> & source):
    m_eps_t(source.m_eps_t), m_eps_p(source.m_eps_p), m_hardening(source.m_hardening), m_m_type(source.m_m_type),fpressure(source.fpressure ),fmatprop(source.fmatprop),fmatpropinit(source.fmatprop),fflux(source.fflux),fdPorePressure(source.fdPorePressure),fSolU(source.fSolU),fGradSolU(source.fGradSolU){ }

    /// Destructor
    ~TPZPlasticState(){ }

    /// Operator =
    const TPZPlasticState<T> & operator=(const TPZPlasticState<T> & source);

    /// Operator -=
    const TPZPlasticState<T> & operator-=(const TPZPlasticState<T> & source);

    /// Operator +=
    const TPZPlasticState<T> & operator+=(const TPZPlasticState<T> & source);

    /// Operator *=
    const TPZPlasticState<T> & operator*=(const TPZPlasticState<T> & source);

    /// Operator<<
    friend std::ostream& operator<<( std::ostream& Out, const TPZPlasticState<T> & s )
    {
        s.Print(Out);
        return Out;
    }

    /// More complete then Operator << because it allows derivatives supression.
    void Print(std::ostream& Out, int fadDerivatives = 1)const;

    // class Access Members (needed for const PlasticState access)
    /// Tensors representing the total and plastic strain states

    const TPZTensor<T> & EpsT() const
    { return m_eps_t; }

    const TPZTensor<T> & EpsP() const
    { return m_eps_p; }

    const T & VolHardening() const
    { return m_hardening; }

    const int & MType() const
    { return m_m_type; }

    const TPZVec<T> & MatProp() const
    { return fmatprop; }

    const TPZVec<T> & MatPropInit() const
    { return fmatpropinit; }

    const TPZVec<T> & Flux() const
    { return fflux; }

    const T & Pressure() const
    { return fpressure; }

    const TPZVec<T> & GradP() const
    { return fdPorePressure; }

    const TPZVec<T> & DisplacementU() const
    { return fSolU; }

    const  TPZTensor<T> & GradU() const
    { return fGradSolU; }


    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override {
        m_eps_t.Read(buf,context);
        m_eps_p.Read(buf,context);
        buf.Read(&m_hardening);
        buf.Read(&m_m_type);
        buf.Read(fmatprop);
        buf.Read(fmatpropinit);
        buf.Read(fflux);
        buf.Read(&fpressure);
        buf.Read(fdPorePressure);
        buf.Read(fSolU);
        fGradSolU.Read(buf,context);


    }

    void Write(TPZStream &buf, int withclassid) const override{
        m_eps_t.Write(buf,withclassid);
        m_eps_p.Write(buf,withclassid);

        buf.Write(&m_hardening);
        buf.Write(&m_m_type);

        buf.Write(fmatprop);
        buf.Write(fmatpropinit);
        buf.Write(fflux);
        buf.Write(&fpressure);
        buf.Write(fdPorePressure);
        buf.Write(fSolU);
        fGradSolU.Write(buf,withclassid);
    }

    void CleanUp() {
        m_eps_t.Zero();
        m_eps_p.Zero();
        m_hardening = T(0.);
        m_m_type = 0;

        for(int i=0;i<fmatprop.size();i++)fmatprop[i]=0.;
        for(int i=0;i<fmatpropinit.size();i++)fmatpropinit[i]=0.;
        for(int i=0;i<fflux.size();i++)fflux[i]=0.;

        fpressure=0.;

        for(int i=0;i<fdPorePressure.size();i++)fdPorePressure[i]=0.;
        for(int i=0;i<fSolU.size();i++)fSolU[i]=0.;

        fGradSolU.Zero();
    }

    /**
     * Similar to Operator=, but allows copies among different template specializations.
     * When using FAD class Types as template classes, NO DERIVATIVES are copied using this functions,
     * enabling copies from REAL to FAD types and also vice-versa.
     */
    template <class T1>
    void CopyTo(TPZPlasticState<T1> & target) const;

};

template <class T>
int TPZPlasticState<T>::ClassId() const{
    return Hash("TPZPlasticState") ^ ClassIdOrHash<T>() << 1;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator=(const TPZPlasticState<T> & source)
{
    m_eps_t = source.EpsT();
    m_eps_p = source.EpsP();
    m_hardening = source.VolHardening();
    m_m_type = source.MType();
	fmatprop=source.MatProp();
    fmatpropinit=source.MatPropInit();
    fflux = source.Flux();
    fpressure =source.Pressure();
    fdPorePressure=source.GradP();
    fSolU=source.DisplacementU();
    fGradSolU=source.GradU();
    return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator+=(const TPZPlasticState<T> & source)
{
    m_eps_t += source.EpsT();
    m_eps_p += source.EpsP();
    m_hardening+= source.VolHardening();
    return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator-=(const TPZPlasticState<T> & source)
{
    m_eps_t -= source.EpsT();
    m_eps_p -= source.EpsP();
    m_hardening-= source.VolHardening();
    return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator*=(const TPZPlasticState<T> & source)
{
    m_eps_t *= source.EpsT();
    m_eps_p *= source.EpsP();
    m_hardening*= source.VolHardening();
    return *this;
}

template <class T>
inline void TPZPlasticState<T>::Print(std::ostream& Out, int fadDerivatives) const
{
    if (fadDerivatives) {
        Out << "\tm_eps_t = ";
        for (int i = 0; i < 6; i++) Out << m_eps_t[i] << " ";
        Out << std::endl;
        Out << "\tm_eps_p = ";
        for (int i = 0; i < 6; i++) Out << m_eps_p[i] << " ";
        Out << std::endl;
        Out << "\tm_hardening = " << m_hardening << std::endl;
        Out << "\tfpressure = "  << fpressure  << std::endl;

        Out << "\tfmatprop = ";
        for (int i=0;i<(int)fmatprop.size();++i) Out << fmatprop[i] << " ";
        Out << std::endl;

        Out << "\tfmatpropinit = ";
        for (int i=0;i<(int)fmatpropinit.size();++i) Out << fmatpropinit[i] << " ";
        Out << std::endl;

        Out << "\tfflux = ";
        for (int i=0;i<(int)fflux.size();++i) Out << fflux[i] << " ";
        Out << std::endl;

        Out << "\tGradP = ";
        for (int i=0;i<(int)fdPorePressure.size();++i) Out << fdPorePressure[i] << " ";
        Out << std::endl;

        Out << "\tSolU = ";
        for (int i=0;i<(int)fSolU.size();++i) Out << fSolU[i] << " ";
        Out << std::endl;

        Out << "\fGradSolU = ";
        for (int i = 0; i < 6; i++) Out << fGradSolU[i] << " ";;
        Out << std::endl;
    } else {
        Out << "\tm_eps_t = ";
        for (int i = 0; i < 6; i++) Out << TPZExtractVal::val(m_eps_t[i]) << " ";
        Out << std::endl;
        Out << "\tm_eps_p = ";
        for (int i = 0; i < 6; i++) Out << TPZExtractVal::val(m_eps_p[i]) << " ";
        Out << std::endl;
        Out << "\tm_hardening = " << TPZExtractVal::val(m_hardening) << std::endl;
        Out << "\tfpressure = "  << shapeFAD::val(fpressure) << std::endl;

        Out << "\tfmatprop = ";
        for (int i=0;i<(int)fmatprop.size();++i) Out << shapeFAD::val(fmatprop[i]) << " ";
        Out << std::endl;

        Out << "\tfmatpropinit = ";
        for (int i=0;i<(int)fmatpropinit.size();++i) Out << shapeFAD::val(fmatpropinit[i]) << " ";
        Out << std::endl;

        Out << "\tfflux = ";
        for (int i=0;i<(int)fflux.size();++i) Out << shapeFAD::val(fflux[i]) << " ";
        Out << std::endl;

        Out << "\tGradP = ";
        for (int i=0;i<(int)fdPorePressure.size();++i) Out << shapeFAD::val(fdPorePressure[i]) << " ";
        Out << std::endl;

        Out << "\tSolU = ";
        for (int i=0;i<(int)fSolU.size();++i) Out << shapeFAD::val(fSolU[i]) << " ";
        Out << std::endl;


        Out << "\tm_eps_t = ";
        for (int i = 0; i < 6; i++) Out << TPZExtractVal::val(fGradSolU[i]) << " ";
        Out << std::endl;
    }
    Out << "\tm_m_type = " << m_m_type << std::endl;
}

template <class T>
template <class T1>
void TPZPlasticState<T>::CopyTo(TPZPlasticState<T1> & target) const
{
    EpsT().CopyTo(target.m_eps_t);
    EpsP().CopyTo(target.m_eps_p);
    target.m_hardening = TPZExtractVal::val( VolHardening() );
    target.m_m_type    = MType();

    target.fpressure = shapeFAD::val( Pressure() );

    // propriedades / fluxos
    target.fmatprop.Resize( (int)fmatprop.size() );
    target.fmatpropinit.Resize( (int)fmatpropinit.size() );
    target.fflux.Resize( (int)fflux.size() );
    for (int i=0;i<(int)fmatprop.size();    ++i) target.fmatprop[i]     = shapeFAD::val(fmatprop[i]);
    for (int i=0;i<(int)fmatpropinit.size();++i) target.fmatpropinit[i] = shapeFAD::val(fmatpropinit[i]);
    for (int i=0;i<(int)fflux.size();       ++i) target.fflux[i]        = shapeFAD::val(fflux[i]);

    // gradiente de pressão, deslocamentos e gradientes de desloc.
    target.fdPorePressure.Resize( (int)fdPorePressure.size() );
    for (int i=0;i<(int)fdPorePressure.size(); ++i) target.fdPorePressure[i] = shapeFAD::val(fdPorePressure[i]);

    target.fSolU.Resize( (int)fSolU.size() );
    for (int i=0;i<(int)fSolU.size(); ++i) target.fSolU[i] = shapeFAD::val(fSolU[i]);
    GradU().CopyTo(target.fGradSolU);

    // target.fGradSolU.Resize( (int)fGradSolU.size() );
    // for (int i=0;i<(int)fGradSolU.size(); ++i) target.fGradSolU[i] = shapeFAD::val(fGradSolU[i]);
}

#endif

