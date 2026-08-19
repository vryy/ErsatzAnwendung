//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 15, 2026 $
//
//


#if !defined(KRATOS_ERSATZ_ANWENDUNG_ECSW_SNAPSHOT_COLLECTING_SCHEME_H_INCLUDED )
#define  KRATOS_ERSATZ_ANWENDUNG_ECSW_SNAPSHOT_COLLECTING_SCHEME_H_INCLUDED


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
class EcswSnapshotCollectingProcess : public SnapshotCollectingProcess<TSparseSpace, TDenseSpace, TModelPart>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( EcswSnapshotCollectingProcess );

    typedef SnapshotCollectingProcess<TSparseSpace, TDenseSpace, TModelPart> BaseType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::TLinearSolver TLinearSolver;

    typedef typename BaseType::ModelPartType ModelPartType;
    typedef typename ModelPartType::IndexType IndexType;

    EcswSnapshotCollectingProcess(typename TLinearSolver::Pointer pLinearSystemSolver)
    : BaseType(pLinearSystemSolver)
    {
    }

    void ExecuteInitialize() override
    {
        BaseType::ExecuteInitialize();

        // collect the "unassembled" elemental residuals
        mElementSnapshot.push_back(this->TakeElementalSnapshot(this->GetModelPart()));
    }

    void ExecuteFinalizeSolutionStep() override
    {
        BaseType::ExecuteFinalizeSolutionStep();

        // collect the "unassembled" elemental residuals
        mElementSnapshot.push_back(this->TakeElementalSnapshot(this->GetModelPart()));
    }

    /**
     * Compute the elemental weighting for ECSW method and save it to the file
     */
    void SaveElementalWeighting(const std::string& filename, std::size_t number_of_modes) const
    {
        Matrix G;
        Vector b;
        this->ConstructEcswSystem(G, b, number_of_modes);

        /* Solve the NNLS problem */
        // TODO: Implement the NNLS solver here.
    }

    /// Construct the system data stemming from collecting snapshots for ECSW.
    /// This system is used to compute the elemental weighting for ECSW method
    /// using Non-Negative Least Squares (NNLS) method.
    void ConstructSystem(Matrix& rG, Vector& rb,
            std::map<IndexType, IndexType>& rElementWeightIndex, std::size_t number_of_modes) const
    {
        /* Compute the projection matrix */
        Matrix U;
        Vector S;
        BaseType::ComputePrincipalComponents(BaseType::GetSnapshot(), U, S);

        const std::size_t m = U.size1();
        const std::size_t nt = BaseType::GetSnapshot().size();
        const std::size_t k = std::min(nt, number_of_modes);

        Matrix Phi(m, k);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < k; ++j)
                Phi(i, j) = U(i, j);

        /* Construct the system for Non - Negative Least Squares (NNLS) solution */

        // count the active elements
        std::size_t n_active_elements = 0;
        const auto& r_model_part = this->GetModelPart();
        for (auto it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it)
        {
            IndexType element_id = it->Id();
            if (it->Is(ACTIVE))
                ++n_active_elements;
        }

        if (rG.size1() != nt * k || rG.size2() != n_active_elements)
            rG.resize(nt * k, n_active_elements, false);

        if (rb.size() != nt * k)
            rb.resize(nt * k, false);

        typedef typename ModelPartType::ElementType::EquationIdVectorType EquationIdVectorType;
        EquationIdVectorType EquationId;
        const auto& CurrentProcessInfo = r_model_part.GetProcessInfo();
        rElementWeightIndex.clear();
        for (std::size_t i = 0; i < nt; ++i)
        {
            const auto& single_snapshot = mElementSnapshot[i];
            Vector bifull(m);
            noalias(bifull) = ZeroVector(m);
            std::size_t ie = 0;
            for (auto it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it, ++ie)
            {
                if (!it->Is(ACTIVE))
                    continue;

                // record the index of the element in the weight vector

                rElementWeightIndex[it->Id()] = ie;

                // assemble to the full residual vector for the i-th snapshot

                const Vector& elemental_residual = single_snapshot.at(it->Id());

                it->EquationIdVector(EquationId, CurrentProcessInfo);

                this->AssembleRHS(bifull, elemental_residual, EquationId);

                // compute and assemble to the reduced elemental contributions for the i-th snapshot

                Matrix localV(EquationId.size(), k);
                for (std::size_t j = 0; j < EquationId.size(); ++j)
                {
                    for (std::size_t l = 0; l < k; ++l)
                    {
                        localV(j, l) = Phi(EquationId[j], l);
                    }
                }

                Vector reduced_elemental_residual = prod(trans(localV), elemental_residual);
                // subrange(rG, i * k, (i + 1) * k, ie, ie + 1) = reduced_elemental_residual;
                for (std::size_t j = 0; j < k; ++j)
                    rG(i * k + j, ie) = reduced_elemental_residual(j);
            }

            // compute and assemble the reduced residual vector for the i-th snapshot

            Vector bi = prod(trans(Phi), bifull);
            noalias(subrange(rb, i * k, (i + 1) * k)) = bi;
        }
    }

private:

    /// Collecting the snapshot of the "unassembled" elemental contributions (i.e., the elemental residuals) for ECSW method
    std::vector<std::map<IndexType, Vector> > mElementSnapshot = {};

    /// Collect the "unassembled" elemental contributions (i.e., the elemental residuals) for ECSW method
    std::map<IndexType, Vector> TakeElementalSnapshot(const ModelPartType& r_model_part) const
    {
        std::map<IndexType, Vector> single_snapshot;
        for (auto it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it)
        {
            IndexType element_id = it->Id();
            Vector elemental_residual;
            it->CalculateRightHandSide(elemental_residual, r_model_part.GetProcessInfo());
            single_snapshot[element_id] = elemental_residual;
        }
        return single_snapshot;
    }

    /// Assemble the elemental contributions to the global residual vector
    void AssembleRHS(
        TSystemVectorType& b,
        const LocalSystemVectorType& RHS_Contribution,
        const typename ModelPartType::ElementType::EquationIdVectorType& EquationId
    ) const
    {
        unsigned int local_size = RHS_Contribution.size();
        unsigned int global_size = b.size();

        for (unsigned int i_local = 0; i_local < local_size; ++i_local)
        {
            unsigned int i_global = EquationId[i_local];
            if (i_global < global_size) //on "free" DOFs
            {
                b[i_global] += RHS_Contribution[i_local];
            }
        }
    }

}; /* Class EcswSnapshotCollectingProcess */

}  /* namespace Kratos.*/

#endif /* KRATOS_ERSATZ_ANWENDUNG_ECSW_SNAPSHOT_COLLECTING_SCHEME_H_INCLUDED  defined */
