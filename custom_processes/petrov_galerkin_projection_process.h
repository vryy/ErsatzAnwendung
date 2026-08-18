//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 14, 2026 $
//
//


#if !defined(KRATOS_ERSATZ_ANWENDUNG_PETROV_GALERKIN_PROJECTION_PROCESS_H_INCLUDED )
#define  KRATOS_ERSATZ_ANWENDUNG_PETROV_GALERKIN_PROJECTION_PROCESS_H_INCLUDED


/* System includes */

/* External includes */

/* Project includes */
#include "custom_processes/pod_process.h"
#include "custom_utilities/pod_utils.h"


namespace Kratos
{

/**
 * This process applies POD mode via projection W^T K V, where W and V are the left and right POD basis matrices respectively,
 * and K is the system matrix.
 */
template<class TSparseSpace, class TDenseSpace, class TModelPart>
class PetrovGalerkinProjectionProcess : public PodProcess<TSparseSpace, TDenseSpace, TModelPart>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( PetrovGalerkinProjectionProcess );

    typedef PodProcess<TSparseSpace, TDenseSpace, TModelPart> BaseType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    /**
     * Default constructor
     */
    PetrovGalerkinProjectionProcess(const Matrix& rPhiL, const Matrix& rPhiR)
    : BaseType(), mPhiL(rPhiL), mPhiR(rPhiR)
    {}

    void ApplyProjection(TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb) override
    {
        ApplyProjection(mPhiL, mPhiR, rA, rDx, rb);
    }

private:

    /**
      * Utility function to apply the projection on the linear system
      */
    template<typename TLocalSystemMatrixType>
    static void ApplyProjection(TLocalSystemMatrixType& rPhiL, TLocalSystemMatrixType& rPhiR,
        TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb)
    {
        const std::size_t rsize = rPhiR.size2(); // size of reduced system

        // construct the reduced linear system
        TLocalSystemMatrixType Ared;
        POD_Utils::WtKV(Ared, rPhiL, rA, rPhiR);

        TSystemVectorType bred(rsize);
        noalias(bred) = prod(trans(rPhiL), rb);

        // solve the reduced linear system
        TSystemVectorType xred(rsize);
        int singular = POD_Utils::Solve(Ared, xred, bred);
        if (singular)
            KRATOS_ERROR << "The reduced system matrix is singular";

        // project back the solution to full space
        noalias(rDx) = prod(rPhiR, xred);
    }

    Matrix mPhiL, mPhiR; // projection matrices

}; /* Class PetrovGalerkinProjectionProcess */

}  /* namespace Kratos.*/

#endif /* KRATOS_ERSATZ_ANWENDUNG_PETROV_GALERKIN_PROJECTION_PROCESS_H_INCLUDED  defined */
