/**
 * @file
 */


#include "TPZPlasticStepPV.h"
#include "TPZElasticResponse.h"

#include "pzsandlerextPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElasticResponse.h"
#include "TPZYCVonMises.h"
#include "pzlog.h"

//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.material.TPZPlasticStepPV"));
//#endif
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.poroelastoplastic2"));
#endif

#ifdef LOG4CXX
static LoggerPtr logger2(Logger::getLogger("plasticity.poroelastoplastic"));
#endif

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma)
{
	TPZTensor<REAL>::TPZDecomposed DecompSig; // It may be SigTr or SigPr Decomposition, dependes on the part of this method
    TPZTensor<REAL> sigtr;

    //
    TPZTensor<REAL> epsTr, epsPN, epsElaNp1;
    epsPN = fN.fEpsP;
    epsTr = epsTotal;
    epsTr -= epsPN; // Porque soh tem implementado o operator -=

    // Compute and Decomposition of SigTrial
    fER.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
    sigtr.EigenSystem(DecompSig);
    TPZManVector<REAL, 3> sigtrvec(DecompSig.fEigenvalues), sigprvec(3, 0.);

    // ReturMap in the principal values
    STATE nextalpha = -6378.;
    fYC.ProjectSigma(sigtrvec, fN.fAlpha, sigprvec, nextalpha);
    fN.fAlpha = nextalpha;
#ifdef LOG4CXX_KEEP
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    // Reconstruction of sigmaprTensor
    DecompSig.fEigenvalues = sigprvec; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
    sigma = TPZTensor<REAL>(DecompSig);

    fER.ComputeDeformation(sigma, epsElaNp1);
    fN.fEpsT = epsTotal;
    epsPN = epsTotal;
    epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
    fN.fEpsP = epsPN;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) {
    TPZTensor<REAL>::TPZDecomposed DecompSig, DecompEps; // It may be SigTr or SigPr Decomposition, dependes on the part of this method
    TPZTensor<REAL> sigtr;

    TPZTensor<REAL> epsTr, epsPN, epsElaNp1;
    epsPN = fN.fEpsP;
    epsTr = epsTotal;
    epsTr -= epsPN; // Porque soh tem implementado o operator -=

    // Compute and Decomposition of SigTrial
    fER.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
    epsTr.EigenSystemJacobi(DecompEps);
    sigtr.EigenSystemJacobi(DecompSig);

    TPZManVector<REAL, 3> sigtrvec(DecompSig.fEigenvalues), sigprvec(3, 0.);


    // ReturnMap in the principal values
    STATE nextalpha = -6378.;
    STATE printPlastic = fN.Alpha();
    TPZFNMatrix<9> GradSigma(3, 3, 0.);
	fYC.ProjectSigmaDep(sigtrvec, fN.fAlpha, sigprvec, nextalpha, GradSigma);

    fN.fAlpha = nextalpha;
	DecompSig.fEigenvalues = sigprvec; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
	sigma = TPZTensor<REAL>(DecompSig);
	fER.ComputeDeformation(sigma, epsElaNp1);


	//TANGENT MATRIX
	TPZFNMatrix<36> dSigDe(6, 6, 0.);
	fYC.ComputeDep(sigma, epsTr, epsElaNp1, dSigDe);

    fN.fEpsT = epsTotal;
    epsPN = epsTotal;
    epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
    fN.fEpsP = epsPN;
    Dep = dSigDe;
}





template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev, TPZVec<REAL> &conv) {
    TPZTensor<REAL> eps1, eps2, SigmaTemp, Sigma1, Sigma2;
    TPZFNMatrix <36> dSigDe(6, 6, 0.);
    TPZStack<REAL> coef;

    fN.fEpsP.Scale(0.);
    fN.fEpsT.Scale(0.);
    fN.fAlpha = kprev;
    this->ApplyStrainComputeDep(EpsIni, SigmaTemp, dSigDe);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "EpsIni " << EpsIni << "\nSigmaTemp " << SigmaTemp << "\ndSidDe " << dSigDe << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fN.fEpsP.Scale(0.);
    fN.fEpsT.Scale(0.);
    fN.fAlpha = kprev;

    REAL scale = 1.;
    REAL alphatable[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    for (int i = 0; i < 6; i++) {
        alphatable[i] *= scale;
    }
    for (int ia = 0; ia < 5; ia++) {
        REAL alpha1 = alphatable[0];
        REAL alpha2 = alphatable[ia + 1];
        eps1.Scale(0.);
        eps2.Scale(0.);
        eps1 = EpsIni;
        eps2 = EpsIni;
        eps1.Add(deps, alpha1);
        eps2.Add(deps, alpha2);

        fN.fEpsT = EpsIni;
        this->ApplyStrainComputeSigma(eps1, Sigma1);
        fN.fEpsP.Scale(0.);
        fN.fEpsT.Scale(0.);
        fN.fAlpha = kprev;

        fN.fEpsT = EpsIni;
        this->ApplyStrainComputeSigma(eps2, Sigma2);
        fN.fEpsP.Scale(0.);
        fN.fEpsT.Scale(0.);
        fN.fAlpha = kprev;

        TPZFNMatrix <6> deps1(6, 1, 0.), deps2(6, 1, 0.);
        TPZFNMatrix <9> depsMat(3, 3, 0.);
        depsMat = deps;
        deps1 = FromMatToVoight(depsMat);
        deps2 = FromMatToVoight(depsMat);

        TPZFNMatrix <6> tanmult1(6, 1, 0.), tanmult2(6, 1, 0.);
        dSigDe.Multiply(deps1, tanmult1);
        dSigDe.Multiply(deps2, tanmult2);

        for (int i = 0; i < 6; i++) {
            tanmult1(i, 0) *= alpha1;
            tanmult2(i, 0) *= alpha2;
        }

        TPZFNMatrix <9> SigMatTemp33(3, 3, 0.);
        TPZFNMatrix <6> sigprMat(6, 1, 0.), sigpr1Mat(6, 1, 0.), sigpr2Mat(6, 1, 0.);
        SigMatTemp33 = SigmaTemp;
        sigprMat = FromMatToVoight(SigMatTemp33);
        SigMatTemp33 = Sigma1;
        sigpr1Mat = FromMatToVoight(SigMatTemp33);
        SigMatTemp33 = Sigma2;
        sigpr2Mat = FromMatToVoight(SigMatTemp33);

        TPZFNMatrix<6> error1(6, 1, 0.), error2(6, 1, 0.);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sigprMat.Print("sigprMat", sout);
            sigpr1Mat.Print("sigpr1Mat", sout);
            tanmult1.Print("tanmult1", sout);
            sigpr2Mat.Print("sigpr2Mat", sout);
            tanmult2.Print("tanmult2", sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        for (int i = 0; i < 6; i++) {
            error1(i, 0) = sigpr1Mat(i, 0) - sigprMat(i, 0) - tanmult1(i, 0);
            error2(i, 0) = sigpr2Mat(i, 0) - sigprMat(i, 0) - tanmult2(i, 0);
        }
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            error1.Print("error1:", sout);
            error2.Print("error2:", sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        REAL n;
        REAL norm1, norm2;
        norm1 = NormVecOfMat(error1);
        norm2 = NormVecOfMat(error2);
        n = (log(norm1) - log(norm2)) / (log(alpha1) - log(alpha2));
        coef.push_back(n);
    }
    conv = coef;
    std::cout << "coef = " << coef << std::endl;
}

template <class YC_t, class ER_t >
REAL TPZPlasticStepPV<YC_t, ER_t>::ComputeNFromTaylorCheck(REAL alpha1, REAL alpha2, TPZFMatrix<REAL> &error1Mat, TPZFMatrix<REAL> &error2Mat)
{
    REAL norm1, norm2, n;
    norm1 = NormVecOfMat(error1Mat);
    norm2 = NormVecOfMat(error2Mat);
    n = log(norm1 / norm2) / log(alpha1 / alpha2);
    return n;
}

REAL NormVecOfMat(TPZFNMatrix <9> mat)
{
    REAL norm = 0.;
    for (int i = 0; i < mat.Rows(); i++) {
        norm += mat(i, 0) * mat(i, 0);
    }
    norm = sqrt(norm);
    return norm;
}

REAL InnerVecOfMat(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2)
{
    REAL dot = 0.;
    for (int i = 0; i < m1.Rows(); i++) {
        dot += m1(i, 0) * m2(i, 0);
    }
    return dot;
}

TPZFMatrix<REAL> ProdT(TPZManVector<REAL,3> &v1, TPZManVector<REAL,3> &v2) {
    TPZFMatrix<REAL> mat(3, 3, 0.);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mat(i, j) = v1[i] * v2[j];
        }
    }
    return mat;
}


TPZFMatrix<REAL> ProdT2(TPZFMatrix<STATE> v1, TPZFMatrix<STATE> v2) {
	int sz = v1.Rows();
	if (v1.Cols() > 1) {
		std::cout << "tensor product must be two vectors " << std::endl;
		DebugStop();
	}
	TPZFMatrix<REAL> mat(sz, sz, 0.);
	for (int i = 0; i < sz; i++) {
		for (int j = 0; j < sz; j++) {
			mat(i, j) = v1(i, 0) * v2(i, 0);
		}
	}
	return mat;
}

TPZFMatrix<STATE> FromTensorToStandardOrder(TPZTensor<STATE> mat)
{
	TPZFMatrix<STATE>  b(6,1,0.);
	b(0, 0) = mat.XX();
	b(1, 0) = mat.YY();
	b(2, 0) = mat.ZZ();
	b(3, 0) = mat.XZ();
	b(4, 0) = mat.YZ();
	b(5, 0) = mat.XY();

	return b;
}

TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat)
{
    TPZFNMatrix <6> voi(6, 1, 0.);
    int k = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
            voi(k++, 0) = mat(i, j);
        }
    }
    return voi;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrain(const TPZTensor<REAL> &epsTotal)
{

    std::cout << " \n this method is not implemented in PlasticStepPV. ";
    DebugStop();

}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyLoad(const TPZTensor<REAL> & GivenStress, TPZTensor<REAL> &epsTotal)
{

    TPZPlasticState<STATE> prevstate = GetState();
    epsTotal = prevstate.fEpsP;
    TPZTensor<STATE> GuessStress, Diff, Diff2, deps;
    TPZFNMatrix<36, STATE> Dep(6, 6);
    TPZFNMatrix<6, STATE> DiffFN(6, 1);

    ApplyStrainComputeDep(epsTotal, GuessStress, Dep);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        Dep.Print("Dep = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    Diff = GivenStress;
    Diff -= GuessStress;

    STATE norm = Norm(Diff), normprev;
    STATE tol = 1.e-7;
    int counter = 0;

    while (norm>tol && counter<30)
    {
        CopyFromTensorToFMatrix(Diff, DiffFN);
        std::list<long> singular;
        Dep.Solve_LU(&DiffFN, singular);
        CopyFromFMatrixToTensor(DiffFN, Diff);
        TPZTensor<STATE> epsprev(epsTotal);
        normprev = norm;
        STATE scale = 1.;
        int counter2 = 0;
        do {
            for (unsigned int i = 0; i < 6; ++i)epsTotal.fData[i] = epsprev.fData[i] + scale * Diff.fData[i];

            ApplyStrainComputeDep(epsTotal, GuessStress, Dep);
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                Dep.Print("Dep = ", sout, EMathematicaInput);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif

            fN = prevstate;
            Diff2 = GivenStress;
            Diff2 -= GuessStress;
            CopyFromTensorToFMatrix(Diff2, DiffFN);
            Dep.Solve_LU(&DiffFN, singular);
            norm = Norm(Diff2);
            //scale*=0.5;
            counter2++;
        } while (norm >= normprev && counter2 < 30);
        Diff = Diff2;
        counter++;
    }
    ApplyStrainComputeDep(epsTotal, GuessStress, Dep);

}

template <class YC_t, class ER_t >
TPZPlasticState<STATE>  TPZPlasticStepPV<YC_t, ER_t>::GetState() const
{
    return fN;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Phi(const TPZTensor<STATE> &eps, TPZVec<REAL> &phi) const
{
    TPZTensor<STATE> sigma;
    fER.Compute(eps, sigma);
    TPZTensor<STATE>::TPZDecomposed DecSig;
    sigma.EigenSystem(DecSig);
    TPZVec<STATE> sigvec(3);
    sigvec[0] = DecSig.fEigenvalues[0];
    sigvec[1] = DecSig.fEigenvalues[1];
    sigvec[2] = DecSig.fEigenvalues[2];
    fYC.Phi(sigvec, fN.Alpha(), phi);
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::SetState(const TPZPlasticState<REAL> &state)
{
    fN = state;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Read(TPZStream &buf)
{
    fYC.Read(buf);
    fER.Read(buf);
    buf.Read(&fResTol);
    buf.Read(&fMaxNewton);
    fN.Read(buf);
}

/** @brief Object which represents the yield criterium */
//YC_t fYC;

/** @brief Object representing the elastic response */
//ER_t fER;

/** @brief Residual tolerance accepted in the plastic loop processes */
//REAL fResTol;

/** @brief Maximum number of Newton interations allowed in the nonlinear solvers */
//int fMaxNewton;	// COLOCAR = 30 (sugestao do erick!)




/** @brief Plastic State Variables (EpsT, EpsP, Alpha) at the current time step */
//TPZPlasticState<STATE> fN;

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Write(TPZStream &buf) const
{
    fYC.Write(buf);
    fER.Write(buf);
    buf.Write(&fResTol);
    buf.Write(&fMaxNewton);
    fN.Write(buf);

}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::CopyFromFMatrixToTensor(TPZFMatrix<STATE> FNM,TPZTensor<STATE> &copy)
{
    FNM.Resize(6, 1);
    copy.XX() = FNM(0, 0);
    copy.XY() = FNM(1, 0);
    copy.XZ() = FNM(2, 0);
    copy.YY() = FNM(3, 0);
    copy.YZ() = FNM(4, 0);
    copy.ZZ() = FNM(5, 0);
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::CopyFromTensorToFMatrix(TPZTensor<STATE> tensor,TPZFMatrix<STATE> &copy)
{

    copy(0, 0) = tensor.XX();
    copy(1, 0) = tensor.XY();
    copy(2, 0) = tensor.XZ();
    copy(3, 0) = tensor.YY();
    copy(4, 0) = tensor.YZ();
    copy(5, 0) = tensor.ZZ();
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::SetElasticResponse(TPZElasticResponse &ER)
{
    fER = ER;
    fYC.SetElasticResponse(ER);
}
#include "TPZYCVonMises.h"
template class TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>;
template class TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>;
template class TPZPlasticStepPV<TPZYCVonMises, TPZElasticResponse>;
/*
 // Correcao do giro rigido
 for (int i = 0; i < 2; i++) {
 for (int j = i+1; j<3 ; j++) {
 deigeneps = DecompEps.fEigenvalues[i]  - DecompEps.fEigenvalues[j];
 deigensig = sigprvec[i] - sigprvec[j];
 TPZFNMatrix<9,REAL> tempMat(3,1,0.);
 depsMat.Multiply(epsegveFromProj[i], tempMat);
 REAL deij = InnerVecOfMat(tempMat,epsegveFromProj[j]);
 REAL factor = 0.;
 if (!IsZero(deigeneps)) {
 factor = deigensig * deij / deigeneps;
 }
 else {
 factor = fER.G() * ( GradSigma(i,i) - GradSigma(i,j) - GradSigma(j,i) + GradSigma(j,j) ) * deij;
 }
 std::cout << "factor = " << factor << std::endl;
 std::cout << "G = " << fER.G() << std::endl;
 GradSigma.Print("GradSigma");
 tempMat.Redim(3, 3);
 tempMat = ProdT(epsegveFromProj[i],epsegveFromProj[j]) + ProdT(epsegveFromProj[j],epsegveFromProj[i]);
 factorMat += tempMat * factor;
 
 }
 }
 */