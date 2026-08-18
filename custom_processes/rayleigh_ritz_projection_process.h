//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 07, 2026 $
//
//


#if !defined(KRATOS_ERSATZ_ANWENDUNG_RAYLEIGH_RITZ_PROJECTION_PROCESS_H_INCLUDED )
#define  KRATOS_ERSATZ_ANWENDUNG_RAYLEIGH_RITZ_PROJECTION_PROCESS_H_INCLUDED


/* System includes */

/* External includes */

/* Project includes */
#include "custom_processes/pod_process.h"
#include "custom_utilities/pod_utils.h"


namespace Kratos
{

/**
 * This process applies POD mode via projection V^T K V, where V is the POD basis matrix and K is the system matrix.
 */
template<class TSparseSpace, class TDenseSpace, class TModelPart>
class RayleighRitzProjectionProcess : public PodProcess<TSparseSpace, TDenseSpace, TModelPart>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( RayleighRitzProjectionProcess );

    typedef PodProcess<TSparseSpace, TDenseSpace, TModelPart> BaseType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    /**
     * Default constructor
     */
    RayleighRitzProjectionProcess(const Matrix& rPhi)
    : BaseType(), mPhi(rPhi)
    {}

    void ApplyProjection(TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb) override
    {
        ApplyProjection(mPhi, rA, rDx, rb);
    }

private:

    /**
      * Utility function to apply the projection on the linear system
      */
    template<typename TLocalSystemMatrixType>
    static void ApplyProjection(TLocalSystemMatrixType& rPhi,
        TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb)
    {
        const std::size_t rsize = rPhi.size2(); // size of reduced system

        // construct the reduced linear system
        TLocalSystemMatrixType Ared;
        POD_Utils::WtKV(Ared, rPhi, rA, rPhi);

        TSystemVectorType bred(rsize);
        noalias(bred) = prod(trans(rPhi), rb);

        // solve the reduced linear system
        TSystemVectorType xred(rsize);
        int singular = POD_Utils::Solve(Ared, xred, bred);
        if (singular)
            KRATOS_ERROR << "The reduced system matrix is singular";

        // project back the solution to full space
        noalias(rDx) = prod(rPhi, xred);
    }

    Matrix mPhi; // projection matrix

}; /* Class RayleighRitzProjectionProcess */

}  /* namespace Kratos.*/

#endif /* KRATOS_ERSATZ_ANWENDUNG_RAYLEIGH_RITZ_PROJECTION_PROCESS_H_INCLUDED  defined */
