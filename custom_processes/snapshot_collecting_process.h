//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 07, 2026 $
//
//


#if !defined(KRATOS_ERSATZ_ANWENDUNG_SNAPSHOT_COLLECTING_PROCESS_H_INCLUDED )
#define  KRATOS_ERSATZ_ANWENDUNG_SNAPSHOT_COLLECTING_PROCESS_H_INCLUDED


/* System includes */

/* External includes */

/* Project includes */
#include "custom_processes/pod_process.h"
#include "custom_utilities/pod_utils.h"


namespace Kratos
{

/**
 * This process collects the snapshot of the dof set and computing the POD basis after the simulation ends
 */
template<class TSparseSpace, class TDenseSpace, class TModelPart>
class SnapshotCollectingProcess : public PodProcess<TSparseSpace, TDenseSpace, TModelPart>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( SnapshotCollectingProcess );

    typedef PodProcess<TSparseSpace, TDenseSpace, TModelPart> BaseType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::TLinearSolver TLinearSolver;

    SnapshotCollectingProcess(typename TLinearSolver::Pointer pLinearSystemSolver)
    : BaseType(), mpLinearSystemSolver(pLinearSystemSolver)
    {
    }

    void ExecuteInitialize() override
    {
        // record the initial condition of the system
        mSnapshot.push_back(BaseType::TakeSnapshot());
    }

    void ExecuteFinalizeSolutionStep() override
    {
        mSnapshot.push_back(BaseType::TakeSnapshot());
    }

    void ApplyProjection(TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb) override
    {
        mpLinearSystemSolver->Solve(rA, rDx, rb);
    }

    /**
      * Compute the principal values of the snapshot matrix
      */
    Vector GetPrincipalValues() const
    {
        Matrix U;
        Vector S;
        BaseType::ComputePrincipalComponents(mSnapshot, U, S);

        return S;
    }

    /**
      * Compute the principal components and vectors via SVD and save it to the file
      */
    Matrix GetPrincipalComponents(const std::size_t number_of_modes) const
    {
        Matrix U;
        Vector S;
        BaseType::ComputePrincipalComponents(mSnapshot, U, S);

        const std::size_t m = U.size1();
        const std::size_t k = std::min(mSnapshot.size(), number_of_modes);

        Matrix Phi(m, k);
        for (size_t col = 0; col < k; ++col) {
            for (size_t row = 0; row < m; ++row) {
                Phi(row, col) = U(row, col);
            }
        }
        return Phi;
    }

    /**
     * Compute the principal components and vectors via SVD and save it to the file
     */
    void SavePrincipalComponents(const std::string& filename, const std::size_t number_of_modes) const
    {
        Matrix U;
        Vector S;
        BaseType::ComputePrincipalComponents(mSnapshot, U, S);

        const std::size_t m = U.size1();
        const std::size_t k = std::min(mSnapshot.size(), number_of_modes);

        if (k > 0)
        {
            #ifdef ERSATZ_APP_USE_MATIO

            std::vector<double> col_major_buffer;
            col_major_buffer.reserve(m * k);

            for (size_t col = 0; col < k; ++col) {
                for (size_t row = 0; row < m; ++row) {
                    col_major_buffer.push_back(U(row, col));
                }
            }

            const std::string variable_name = "Phi";
            POD_Utils::WriteMat(filename, variable_name, col_major_buffer, m, k);

            #else

            std::ofstream file(filename, std::ios::binary);
            file.write(reinterpret_cast<const char*>(&k), sizeof(k));
            file.write(reinterpret_cast<const char*>(&m), sizeof(m));

            Vector T(m);
            for (std::size_t i = 0; i < k; ++i)
            {
                noalias(T) = column(U, i);
                file.write(reinterpret_cast<const char*>(&T[0]), T.size() * sizeof(double));
            }

            file.close();

            #endif
        }
    }

    /// Accessor for the snapshot collected during the simulation
    const std::vector<Vector>& GetSnapshot() const
    {
        return mSnapshot;
    }

private:

    std::vector<Vector> mSnapshot = {};
    typename TLinearSolver::Pointer mpLinearSystemSolver;

}; /* Class SnapshotCollectingProcess */

}  /* namespace Kratos.*/

#endif /* KRATOS_ERSATZ_ANWENDUNG_SNAPSHOT_COLLECTING_PROCESS_H_INCLUDED  defined */
