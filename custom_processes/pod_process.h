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
    typedef TModelPart ModelPartType;
    typedef typename TModelPart::DofsArrayType DofsArrayType;

    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    typedef LinearSolver<TSparseSpace, TDenseSpace, TModelPart> TLinearSolver;

    /**
     * Default constructor
     */
    PodProcess() : BaseType()
    {
    }

    /**
     * Assign the model part
     */
    void SetModelPart(const ModelPartType& rModelPart)
    {
        mpModelPart = &rModelPart;
    }

    /**
     * Get the model part
     */
    const ModelPartType& GetModelPart() const
    {
        if (mpModelPart == nullptr)
            KRATOS_ERROR << "Model part is not assigned";
        return *mpModelPart;
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
     * Read in the principal vectors stored in the data file
     */
    void ReadPrincipalComponents(Matrix& Phi, const std::string& filename, const std::string& dataset_name) const
    {
#ifdef ERSATZ_APP_USE_MATIO

        Phi = POD_Utils::ReadMat<Matrix>(filename, dataset_name);

#else

        std::ifstream file(filename, std::ios::binary);

        std::size_t number_of_modes, m;
        file.read(reinterpret_cast<char*>(&number_of_modes), sizeof(number_of_modes));
        file.read(reinterpret_cast<char*>(&m), sizeof(m));

        Phi.resize(m, number_of_modes, false);

        for (std::size_t i = 0; i < number_of_modes; ++i)
        {
            Vector T(m);
            file.read(reinterpret_cast<char*>(&T[0]), T.size() * sizeof(double));

            noalias(column(Phi, i)) = T;
        }

        file.close();

#endif
    }

    /**
     * Record the current values of the free dof into a vector
     */
    Vector TakeSnapshot() const
    {
        const auto& DofSet = this->GetDofSet();

        std::vector<double> values;

        for (auto it = DofSet.begin(); it != DofSet.end(); ++it)
        {
            if (it->IsFree())
                values.push_back(it->GetSolutionStepValue());
        }

        Vector dof_vector(values.size());
        std::copy(values.begin(), values.end(), dof_vector.begin());
        return dof_vector;
    }

private:

    const DofsArrayType* mpDofSet = nullptr;
    const ModelPartType* mpModelPart = nullptr;

}; /* Class PodProcess */

}  /* namespace Kratos.*/

#endif /* KRATOS_POD_PROCESS_H_INCLUDED  defined */
