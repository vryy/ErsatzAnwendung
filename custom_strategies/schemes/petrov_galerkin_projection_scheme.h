//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Oct 5, 2026 $
//
//

#if !defined(KRATOS_ERSATZ_ANWENDUNG_PETROV_GALERKIN_PROJECTION_SCHEME_H_INCLUDED )
#define  KRATOS_ERSATZ_ANWENDUNG_PETROV_GALERKIN_PROJECTION_SCHEME_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/scheme.h"
#include "ersatz_anwendung_variables.h"


namespace Kratos
{

/**
 * TODO
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>
         class TModelPartType = ModelPart
         >
class PetrovGalerkinProjectionScheme : public RayleighRitzProjectionScheme<TSparseSpace, TDenseSpace, TModelPartType>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( PetrovGalerkinProjectionScheme );

    typedef Scheme<TSparseSpace, TDenseSpace, TModelPartType> SchemeType;
    typedef RayleighRitzProjectionScheme<TSparseSpace, TDenseSpace, TModelPartType> BaseType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;
    typedef typename BaseType::TDenseSpaceType TDenseSpaceType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::ModelPartType ModelPartType;
    typedef typename BaseType::ElementType ElementType;
    typedef typename BaseType::ConditionType ConditionType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    PetrovGalerkinProjectionScheme(typename SchemeType::Pointer pScheme)
    : BaseType(pScheme)
    {
    }

    ///@name Operations
    ///@{

    static typename SchemeType::Pointer Create(typename SchemeType::Pointer pScheme, const LocalSystemMatrixType& Phi, const LocalSystemMatrixType& Psi)
    {
        auto pNewScheme = PetrovGalerkinProjectionScheme::Pointer(new PetrovGalerkinProjectionScheme(pScheme));
        pNewScheme->SetLeftProjectionOperator(Phi);
        pNewScheme->SetRightProjectionOperator(Psi);
        return pNewScheme;
    }

    typename SchemeType::Pointer Clone() const override
    {
        if (BaseType::mpScheme != nullptr)
            return typename SchemeType::Pointer(new PetrovGalerkinProjectionScheme(BaseType::mpScheme->Clone()));
        else
            return nullptr;
    }

#ifdef KRATOS_NONSQUARE_SUPPORT
    void CalculateSystemContributions(
        ElementType& rElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        typename ElementType::EquationIdVectorType& rRowEquationIdVector,
        typename ElementType::EquationIdVectorType& rColEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateSystemContributionsImpl(rElement, LHS_Contribution, RHS_Contribution,
            rRowEquationIdVector, rColEquationIdVector, rCurrentProcessInfo);
    }

    void CalculateSystemContributions(
        ConditionType& rCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        typename ConditionType::EquationIdVectorType& rRowEquationIdVector,
        typename ConditionType::EquationIdVectorType& rColEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateSystemContributionsImpl(rCondition, LHS_Contribution, RHS_Contribution,
            rRowEquationIdVector, rColEquationIdVector, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        ElementType& rElement,
        LocalSystemMatrixType& LHS_Contribution,
        typename ElementType::EquationIdVectorType& rRowEquationIdVector,
        typename ElementType::EquationIdVectorType& rColEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        // compute the elemental contribution of FOM
        BaseType::mpScheme->CalculateLHSContribution(rElement, LHS_Contribution,
            rRowEquationIdVector, rColEquationIdVector, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        ConditionType& rCondition,
        LocalSystemMatrixType& LHS_Contribution,
        typename ConditionType::EquationIdVectorType& rRowEquationIdVector,
        typename ConditionType::EquationIdVectorType& rColEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        // compute the elemental contribution of FOM
        BaseType::mpScheme->CalculateLHSContribution(rCondition, LHS_Contribution,
            rRowEquationIdVector, rColEquationIdVector, rCurrentProcessInfo);
    }
#endif

    ///@}
    ///@name Access
    ///@{

    /// Set the global projection matrix
    void SetLeftProjectionOperator(const LocalSystemMatrixType& Phi)
    {
        BaseType::SetProjectionOperator(Phi);
    }

    /// Set the global projection matrix
    void SetRightProjectionOperator(const LocalSystemMatrixType& Psi)
    {
        mpPsi = &Psi;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        if (BaseType::mpScheme == nullptr)
            return "PetrovGalerkinProjectionScheme<nullptr>";
        else
            return "PetrovGalerkinProjectionScheme<" + BaseType::mpScheme->Info() + ">";
    }

    ///@}

private:

    /// pointer to the global projection matrix. This is used to project back the local constribution
    /// to the global system.
    const LocalSystemMatrixType* mpPsi = nullptr;

    template<typename TEntityType>
    void CalculateSystemContributionsImpl(
        TEntityType& rElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        typename TEntityType::EquationIdVectorType& rRowEquationIdVector,
        typename TEntityType::EquationIdVectorType& rColEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        // compute the elemental contribution of FOM
        BaseType::mpScheme->CalculateSystemContributions(rElement, LHS_Contribution, RHS_Contribution,
                rRowEquationIdVector, rColEquationIdVector, rCurrentProcessInfo);

        if (BaseType::mpPhi == nullptr || mpPsi == nullptr)
            KRATOS_ERROR << "The projection operator is not yet set";

        // construct the ROM contribution
        if (rRowEquationIdVector.size() > 0 && rColEquationIdVector.size())
        {
            // assemble the force of ROM
            const auto& Phi = *BaseType::mpPhi;
            const std::size_t full_row_system_size = Phi.size1();
            const std::size_t reduced_row_system_size = Phi.size2();

            const auto& Psi = *mpPsi;
            const std::size_t full_col_system_size = Psi.size1();
            const std::size_t reduced_col_system_size = Psi.size2();

            LocalSystemMatrixType localV(rRowEquationIdVector.size(), reduced_row_system_size);
            for (std::size_t j = 0; j < rRowEquationIdVector.size(); ++j)
            {
                if (rRowEquationIdVector[j] < full_row_system_size)
                    noalias(row(localV, j)) = row(Phi, rRowEquationIdVector[j]);
                else
                    noalias(row(localV, j)) = ZeroVector(reduced_row_system_size);
            }

            LocalSystemMatrixType localW(rColEquationIdVector.size(), reduced_col_system_size);
            for (std::size_t j = 0; j < rColEquationIdVector.size(); ++j)
            {
                if (rColEquationIdVector[j] < full_col_system_size)
                    noalias(row(localW, j)) = row(Psi, rColEquationIdVector[j]);
                else
                    noalias(row(localW, j)) = ZeroVector(reduced_col_system_size);
            }

            LocalSystemVectorType reduced_elemental_residual = prod(trans(localV), RHS_Contribution);
            LocalSystemMatrixType reduced_elemental_stiffness = prod(trans(localV), Matrix(prod(LHS_Contribution, localW)));

            // modify the output
            RHS_Contribution = reduced_elemental_residual;
            LHS_Contribution = reduced_elemental_stiffness;
            rRowEquationIdVector.resize(reduced_row_system_size);
            for (std::size_t i = 0; i < reduced_row_system_size; ++i)
                rRowEquationIdVector[i] = i;
            rColEquationIdVector.resize(reduced_col_system_size);
            for (std::size_t i = 0; i < reduced_col_system_size; ++i)
                rColEquationIdVector[i] = i;
        }
    }

}; /* Class PetrovGalerkinProjectionScheme */

}  /* namespace Kratos.*/

#endif /* KRATOS_ERSATZ_ANWENDUNG_PETROV_GALERKIN_PROJECTION_SCHEME_H_INCLUDED  defined */
