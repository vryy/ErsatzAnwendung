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
#include "custom_processes/snapshot_collecting_process.h"


namespace Kratos
{

/**
 * This process constructs the least square system for Energy-Conserving Sampling and Weighting method.
 * To avoid the small residual problem in the case the system is loaded by prescribed displacement,
 * the unconstrained residual forces (including the prescribed dofs) are collected in each
 * Newton-Raphson iteration, where the residual norm is larger than a tolerance.
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
    typedef typename ModelPartType::DofsArrayType DofsArrayType;

    EcswSnapshotCollectingProcess(typename TLinearSolver::Pointer pLinearSystemSolver)
    : BaseType(pLinearSystemSolver)
    {
    }

    void SetForceTolerance(double Tol)
    {
        mForceTolerance = Tol;
    }

    void SetNormalize(bool v)
    {
        mNormalize = v;
    }

    void ExecuteInitialize() override
    {
        BaseType::ExecuteInitialize();

        // collect the "unassembled" elemental residuals
        const auto single_snapshot = this->TakeElementalSnapshot(this->GetModelPart());

        // reconstruct the right-hand-side
        Vector bifull = this->AssembleRHS(single_snapshot);

        // check against the tolerance and collect
        if (norm_2(bifull) > mForceTolerance)
        {
            mElementSnapshot.push_back(single_snapshot);
        }
    }

    void ExecuteBuild() override
    {
        BaseType::ExecuteBuild();

        // collect the "unassembled" elemental residuals
        const auto single_snapshot = this->TakeElementalSnapshot(this->GetModelPart());

        // reconstruct the right-hand-side
        Vector bifull = this->AssembleRHS(single_snapshot);

        // check against the tolerance and collect
        if (norm_2(bifull) > mForceTolerance)
        {
            mElementSnapshot.push_back(single_snapshot);
        }
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
        BaseType::ComputePrincipalComponents(this->GetSnapshot(), U, S);

        const std::size_t m = U.size1();    // system size
        const std::size_t nt0 = this->GetSnapshot().size();
        const std::size_t k = std::min(nt0, number_of_modes); // reduced system size

        Matrix Phi(m, k);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < k; ++j)
                Phi(i, j) = U(i, j);

        /* Construct the system for Non - Negative Least Squares (NNLS) solution */

        const std::size_t nt = mElementSnapshot.size();

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
            double norm_factor = 1.0;
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
            }

            // compute and assemble to the reduced elemental contributions for the i-th snapshot

            if (mNormalize)
            {
                norm_factor = 1.0 / norm_2(bifull);
            }

            for (auto it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it)
            {
                if (!it->Is(ACTIVE))
                    continue;

                const Vector& elemental_residual = single_snapshot.at(it->Id());

                it->EquationIdVector(EquationId, CurrentProcessInfo);

                ie = rElementWeightIndex[it->Id()];

                Matrix localV(EquationId.size(), k);
                for (std::size_t j = 0; j < EquationId.size(); ++j)
                {
                    if (EquationId[j] < m)
                    {
                        noalias(row(localV, j)) = row(Phi, EquationId[j]);
                    }
                    else
                    {
                        noalias(row(localV, j)) = ZeroVector(k);
                    }
                }

                Vector reduced_elemental_residual = prod(trans(localV), elemental_residual);
                // subrange(rG, i * k, (i + 1) * k, ie, ie + 1) = reduced_elemental_residual * norm_factor;
                for (std::size_t j = 0; j < k; ++j)
                    rG(i * k + j, ie) = reduced_elemental_residual(j) * norm_factor;
            }

            // compute and assemble the reduced residual vector for the i-th snapshot

            Vector bi = prod(trans(Phi), bifull);
            noalias(subrange(rb, i * k, (i + 1) * k)) = bi * norm_factor;
        }
    }

private:

    /// Container of the snapshot of the "unassembled" elemental contributions (i.e., the elemental residuals) for ECSW method
    std::vector<std::map<IndexType, Vector> > mElementSnapshot = {};

    /// Tolerance for the FOM force, by which it will be recorded in the least square system
    double mForceTolerance = 1e-6;

    /// Flag to normalize the force before recording in the least square system
    bool mNormalize = false;

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

    /// Reconstruct the right hand side from the elemental snapshot
    Vector AssembleRHS(const std::map<IndexType, Vector>& single_snapshot) const
    {
        const auto& r_model_part = this->GetModelPart();

        const auto& dof_set = this->GetDofSet();

        typedef typename ModelPartType::ElementType::EquationIdVectorType EquationIdVectorType;
        EquationIdVectorType EquationId;
        const auto& CurrentProcessInfo = r_model_part.GetProcessInfo();

        // count the free dofs
        std::size_t m = 0;
        for (auto it = dof_set.begin(); it != dof_set.end(); ++it)
            if (it->IsFree())
                ++m;

        // assemble the global residual forces
        Vector bifull(m);
        noalias(bifull) = ZeroVector(m);
        std::size_t ie = 0;
        for (auto it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it, ++ie)
        {
            if (!it->Is(ACTIVE))
                continue;

            // assemble to the full residual vector for the i-th snapshot

            const Vector& elemental_residual = single_snapshot.at(it->Id());

            it->EquationIdVector(EquationId, CurrentProcessInfo);

            this->AssembleRHS(bifull, elemental_residual, EquationId);
        }

        return bifull;
    }

}; /* Class EcswSnapshotCollectingProcess */

}  /* namespace Kratos.*/

#endif /* KRATOS_ERSATZ_ANWENDUNG_ECSW_SNAPSHOT_COLLECTING_SCHEME_H_INCLUDED  defined */
