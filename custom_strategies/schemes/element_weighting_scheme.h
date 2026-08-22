//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 17, 2026 $
//
//


#if !defined(KRATOS_ERSATZ_ANWENDUNG_ELEMENT_WEIGHTING_SCHEME_H_INCLUDED )
#define  KRATOS_ERSATZ_ANWENDUNG_ELEMENT_WEIGHTING_SCHEME_H_INCLUDED


/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/scheme.h"


namespace Kratos
{

/**
 * This scheme is an adapter of other scheme which applies weighting factors to the element contributions.
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>
         class TModelPartType = ModelPart
         >
class ElementWeightingScheme : public Scheme<TSparseSpace, TDenseSpace, TModelPartType>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( ElementWeightingScheme );

    typedef Scheme<TSparseSpace, TDenseSpace, TModelPartType> BaseType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;
    typedef typename BaseType::TDenseSpaceType TDenseSpaceType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename TSparseSpaceType::IndexType IndexType;
    typedef typename BaseType::ModelPartType ModelPartType;
    typedef typename BaseType::ElementType ElementType;
    typedef typename BaseType::ConditionType ConditionType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef KRATOS_DOUBLE_TYPE WeightType;

    ElementWeightingScheme(typename BaseType::Pointer pScheme)
    : BaseType(), mpScheme(pScheme)
    {
    }

    ///@name Operations
    ///@{

    typename BaseType::Pointer Clone() const override
    {
        if (mpScheme != nullptr)
            return typename BaseType::Pointer(new ElementWeightingScheme(mpScheme->Clone()));
        else
            return nullptr;
    }

    void Initialize(ModelPartType& rModelPart) override
    {
        mpScheme->Initialize(rModelPart);
    }

    bool SchemeIsInitialized() const override
    {
        return mpScheme->SchemeIsInitialized();
    }

    void SetSchemeIsInitialized(bool SchemeIsInitializedFlag = true) override
    {
        mpScheme->SetSchemeIsInitialized(SchemeIsInitializedFlag);
    }

    bool ElementsAreInitialized() const override
    {
        return mpScheme->ElementsAreInitialized();
    }

    void SetElementsAreInitialized(bool ElementsAreInitializedFlag = true) override
    {
        mpScheme->SetElementsAreInitialized(ElementsAreInitializedFlag);
    }

    bool ConditionsAreInitialized() const override
    {
        return mpScheme->ConditionsAreInitialized();
    }

    void SetConditionsAreInitialized(bool ConditionsAreInitializedFlag = true) override
    {
        mpScheme->SetConditionsAreInitialized(ConditionsAreInitializedFlag);
    }

    void InitializeElements(ModelPartType& rModelPart) override
    {
        mpScheme->InitializeElements(rModelPart);
    }

    void InitializeConditions(ModelPartType& rModelPart) override
    {
        mpScheme->InitializeConditions(rModelPart);
    }

    void InitializeSolutionStep(
        ModelPartType& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        mpScheme->InitializeSolutionStep(rModelPart, A, Dx, b);
    }

    void FinalizeSolutionStep(
        ModelPartType& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        mpScheme->FinalizeSolutionStep(rModelPart, A, Dx, b);
    }

    void InitializeNonLinIteration(
        ModelPartType& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        mpScheme->InitializeNonLinIteration(rModelPart, A, Dx, b);
    }

    void FinalizeNonLinIteration(
        ModelPartType& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        mpScheme->FinalizeNonLinIteration(rModelPart, A, Dx, b);
    }

    void Predict(
        ModelPartType& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        mpScheme->Predict(rModelPart, rDofSet, A, Dx, b);
    }

    void Update(
        ModelPartType& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        mpScheme->Update(rModelPart, rDofSet, A, Dx, b);
    }

    void CalculateOutputData(
        ModelPartType& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        mpScheme->CalculateOutputData(rModelPart, rDofSet, A, Dx, b);
    }

    void CleanOutputData() override
    {
        mpScheme->CleanOutputData();
    }

    void Clean() override
    {
        mpScheme->Clean();
    }

    void CleanMemory(ElementType& rCurrentElement) override
    {
        mpScheme->CleanMemory(rCurrentElement);
    }

    void CleanMemory(ConditionType& rCurrentCondition) override
    {
        mpScheme->CleanMemory(rCurrentCondition);
    }

    void Clear() override
    {
        mpScheme->Clear();
    }

    int Check(const ModelPartType& rModelPart) const override
    {
        return mpScheme->Check(rModelPart);
    }

    void CalculateSystemContributions(
        ElementType& rElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        typename ElementType::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        // apply weighting scheme
        auto it = mElementWeights.find(rElement.Id());
        if (it != mElementWeights.end())
        {
            WeightType weight = it->second;

            if (weight > 0.0)
            {
                mpScheme->CalculateSystemContributions(rElement,
                    LHS_Contribution, RHS_Contribution, rEquationIdVector, rCurrentProcessInfo);

                LHS_Contribution *= weight;
                RHS_Contribution *= weight;

                return;
            }
            //else
            //    std::cout << "Info: Element " << rElement.Id() << " is ignored." << std::endl;
        }
        else
        {
            KRATOS_ERROR << "ElementWeightingScheme::CalculateSystemContributions: "
                         << "No weight found for element ID " << rElement.Id()
                         << std::endl;
        }

        LHS_Contribution.resize(0, 0);
        RHS_Contribution.resize(0);
        rEquationIdVector.resize(0);
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
        mpScheme->CalculateSystemContributions(rElement, LHS_Contribution, RHS_Contribution, rRowEquationIdVector, rColEquationIdVector, rCurrentProcessInfo);

        /*
        // apply weighting scheme
        auto it = mElementWeights.find(rElement.Id());
        if (it != mElementWeights.end())
        {
            WeightType weight = it->second;
            LHS_Contribution *= weight;
            RHS_Contribution *= weight;
        }
        else
        {
            std::cout << "Warning: ElementWeightingScheme::CalculateSystemContributions: No weight found for element ID " << rElement.Id()
                << ". Using default weight of 1.0." << std::endl;
        }
        */
    }
#endif

    void CalculateSystemContributions(
        ConditionType& rCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        typename ConditionType::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpScheme->CalculateSystemContributions(rCondition, LHS_Contribution, RHS_Contribution, rEquationIdVector, rCurrentProcessInfo);
    }

#ifdef KRATOS_NONSQUARE_SUPPORT
    void CalculateSystemContributions(
        ConditionType& rCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        typename ConditionType::EquationIdVectorType& rRowEquationIdVector,
        typename ConditionType::EquationIdVectorType& rColEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpScheme->CalculateSystemContributions(rCondition, LHS_Contribution, RHS_Contribution, rRowEquationIdVector, rColEquationIdVector, rCurrentProcessInfo);
    }
#endif

    void CalculateRHSContribution(
        ElementType& rElement,
        LocalSystemVectorType& RHS_Contribution,
        typename ElementType::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpScheme->CalculateRHSContribution(rElement, RHS_Contribution, rEquationIdVector, rCurrentProcessInfo);

        /*
        // apply weighting scheme
        auto it = mElementWeights.find(rElement.Id());
        if (it != mElementWeights.end())
        {
            WeightType weight = it->second;
            RHS_Contribution *= weight;
        }
        else
        {
            std::cout << "Warning: ElementWeightingScheme::CalculateSystemContributions: No weight found for element ID " << rElement.Id()
                << ". Using default weight of 1.0." << std::endl;
        }
        */
    }

    void CalculateRHSContribution(
        ConditionType& rCondition,
        LocalSystemVectorType& RHS_Contribution,
        typename ConditionType::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpScheme->CalculateRHSContribution(rCondition, RHS_Contribution, rEquationIdVector, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        ElementType& rElement,
        LocalSystemMatrixType& LHS_Contribution,
        typename ElementType::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpScheme->CalculateLHSContribution(rElement, LHS_Contribution, rEquationIdVector, rCurrentProcessInfo);

        /*
        // apply weighting scheme
        auto it = mElementWeights.find(rElement.Id());
        if (it != mElementWeights.end())
        {
            WeightType weight = it->second;
            LHS_Contribution *= weight;
        }
        else
        {
            std::cout << "Warning: ElementWeightingScheme::CalculateSystemContributions: No weight found for element ID " << rElement.Id()
                << ". Using default weight of 1.0." << std::endl;
        }
        */
    }

#ifdef KRATOS_NONSQUARE_SUPPORT
    void CalculateLHSContribution(
        ElementType& rElement,
        LocalSystemMatrixType& LHS_Contribution,
        typename ElementType::EquationIdVectorType& rRowEquationIdVector,
        typename ElementType::EquationIdVectorType& rColEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpScheme->CalculateLHSContribution(rElement, LHS_Contribution, rRowEquationIdVector, rColEquationIdVector, rCurrentProcessInfo);

        /*
        // apply weighting scheme
        auto it = mElementWeights.find(rElement.Id());
        if (it != mElementWeights.end())
        {
            WeightType weight = it->second;
            LHS_Contribution *= weight;
        }
        else
        {
            std::cout << "Warning: ElementWeightingScheme::CalculateSystemContributions: No weight found for element ID " << rElement.Id()
                << ". Using default weight of 1.0." << std::endl;
        }
        */
    }
#endif

    void CalculateLHSContribution(
        ConditionType& rCondition,
        LocalSystemMatrixType& LHS_Contribution,
        typename ConditionType::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpScheme->CalculateLHSContribution(rCondition, LHS_Contribution, rEquationIdVector, rCurrentProcessInfo);
    }

#ifdef KRATOS_NONSQUARE_SUPPORT
    void CalculateLHSContribution(
        ConditionType& rCondition,
        LocalSystemMatrixType& LHS_Contribution,
        typename ConditionType::EquationIdVectorType& rRowEquationIdVector,
        typename ConditionType::EquationIdVectorType& rColEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpScheme->CalculateLHSContribution(rCondition, LHS_Contribution, rRowEquationIdVector, rColEquationIdVector, rCurrentProcessInfo);
    }
#endif

    ///@}
    ///@name Access
    ///@{

    /// Set the element weight for a specific element ID.
    void SetElementWeight(IndexType ElementId, WeightType Weight)
    {
        mElementWeights[ElementId] = Weight;
    }

    /// Set the element weight map, where the key is the element ID and the value is the weight.
    void SetElementWeights(const std::map<IndexType, WeightType>& rElementWeights)
    {
        mElementWeights = rElementWeights;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        if (mpScheme == nullptr)
            return "ElementWeightingScheme<nullptr>";
        else
            return "ElementWeightingScheme<" + mpScheme->Info() + ">";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        if (mpScheme != nullptr)
            mpScheme->PrintData(rOStream);
    }

    ///@}

private:

    BaseType::Pointer mpScheme = nullptr;
    std::map<IndexType, WeightType> mElementWeights;

}; /* Class ElementWeightingScheme */

}  /* namespace Kratos.*/

#endif /* KRATOS_ERSATZ_ANWENDUNG_ELEMENT_WEIGHTING_SCHEME_H_INCLUDED  defined */
