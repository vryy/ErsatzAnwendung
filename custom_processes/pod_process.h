//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 07, 2026 $
//
//


#if !defined(KRATOS_POD_PROCESS_H_INCLUDED )
#define  KRATOS_POD_PROCESS_H_INCLUDED


/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/pod_utils.h"


namespace Kratos
{

/**
 * This base class provides abstract functions for POD operations.
 * A POD process is also a linear solver, since it can apply projection on the linear system.
 */
template<class TSparseSpace, class TDenseSpace, class TModelPart>
class PodProcess : public Process, public LinearSolver<TSparseSpace, TDenseSpace, TModelPart>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( PodProcess );

    typedef Process BaseType;
    typedef typename TModelPart::DofsArrayType DofsArrayType;

    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef LinearSolver<TSparseSpace, TDenseSpace, TModelPart> TLinearSolver;

    /**
     * Default constructor
     */
    PodProcess() : BaseType()
    {
    }

    /**
     * Assign the dof set
     */
    void SetDofSet(const DofsArrayType& rDofSet)
    {
        mpDofSet = &rDofSet;
    }

    /**
     * Get the dof set
     */
    const DofsArrayType& GetDofSet() const
    {
        if (mpDofSet == nullptr)
            KRATOS_ERROR << "Dof set is not assigned";
        return *mpDofSet;
    }

    /**
     * Apply the projection matrix to the linear system
     */
    virtual void ApplyProjection(TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb)
    {
        KRATOS_ERROR << "Calling base class function";
    }

    bool Solve(TSystemMatrixType& rA, TSystemVectorType& rX, TSystemVectorType& rB) override
    {
        this->ApplyProjection(rA, rX, rB);
        return true;
    }

protected:

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
        POD_Utils::VtKV(Ared, rPhi, rA);

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

    /**
     * Compute the principal components and vectors of the list of snapshots by means of SVD.
     * On the output, the principal components are stored in S and the respective principal vectors
     * on the column of matrix U.
     */
    template<typename TLocalSystemMatrixType, typename TLocalSystemVectorType>
    static void ComputePrincipalComponents(const std::vector<TLocalSystemVectorType>& rSnapshot,
        TLocalSystemMatrixType& U, TLocalSystemVectorType& S)
    {
        const std::size_t n = rSnapshot.size();
        if (n > 0)
        {
            const std::size_t m = rSnapshot[0].size();
            TLocalSystemMatrixType Q(m, n);

            for (std::size_t i = 0; i < n; ++i)
            {
                noalias(column(Q, i)) = rSnapshot[i];
            }

            TLocalSystemMatrixType VT;
            POD_Utils::SVD(Q, U, S, VT);
        }
        else
        {
            U.resize(0, 0, false);
            S.resize(0, false);
        }
    }

    /**
     * Record the current values of the dof set into a vector
     */
    Vector TakeSnapshot() const
    {
        const auto& DofSet = this->GetDofSet();
        const std::size_t SystemSize = DofSet.size();
        Vector dof_vector(SystemSize);

        std::size_t i = 0;
        for (auto it = DofSet.begin(); it != DofSet.end(); ++it, ++i)
        {
            dof_vector[i] = it->GetSolutionStepValue();
        }

        return dof_vector;
    }

private:

    const DofsArrayType* mpDofSet = nullptr;

}; /* Class PodProcess */

}  /* namespace Kratos.*/

#endif /* KRATOS_POD_PROCESS_H_INCLUDED  defined */
