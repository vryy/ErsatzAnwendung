//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 28, 2026 $
//
//


#if !defined(KRATOS_PROJECTION_BASED_POD_BUILDER_AND_SOLVER_H_INCLUDED )
#define  KRATOS_PROJECTION_BASED_POD_BUILDER_AND_SOLVER_H_INCLUDED


/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/**
 * This builder and solver devotes exclusively to construct the reduced linear system
 * employing the projection method. It leverages the scheme to perform the local projection,
 * then the assembly can be performed normally by the base builder and solver
 */
template<class TBuilderAndSolverType>
class ProjectionBasedPodBuilderAndSolver : public TBuilderAndSolverType
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( ProjectionBasedPodBuilderAndSolver );

    typedef TBuilderAndSolverType BaseType;

    typedef typename BaseType::ModelPartType ModelPartType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;
    typedef typename BaseType::TDenseSpaceType TDenseSpaceType;

    typedef typename BaseType::TLinearSolverType TLinearSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::ValueType ValueType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::ElementType ElementType;
    typedef typename BaseType::ConditionType ConditionType;

    typedef typename BaseType::NodesContainerType NodesContainerType;
    typedef typename BaseType::ElementsContainerType ElementsContainerType;
    typedef typename BaseType::ConditionsContainerType ConditionsContainerType;

    typedef typename MatrixVectorTypeSelector<TDataType>::ZeroVectorType ZeroVectorType;

    typedef SkylineLUFactorizationSolver<TSparseSpaceType, TDenseSpaceType, ModelPartType> SkylineLUFactorizationSolverType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Default constructor.
     */
    ProjectionBasedPodBuilderAndSolver()
    : BaseType(typename TLinearSolverType::Pointer(new SkylineLUFactorizationSolverType()))
    {
    }

    /** Copy constructor.
     */
    ProjectionBasedPodBuilderAndSolver(const BaseType& rOther)
    : BaseType(rOther)
    {
    }

    /** Copy constructor.
     */
    ProjectionBasedPodBuilderAndSolver(const ProjectionBasedPodBuilderAndSolver& rOther)
    : BaseType(rOther), mpPhi(rOther.mpPhi)
    {
    }

    /** Destructor.
    */
    ~ProjectionBasedPodBuilderAndSolver() override
    {}

    /*@} */
    /**@name Operators
    */
    /*@{ */


    /*@} */
    /**@name Operations */
    /*@{ */

    void SetProjectionOperator(const LocalSystemMatrixType& Phi)
    {
        mpPhi = &Phi;
    }

    void ResizeAndInitializeVectors(
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ElementsContainerType& rElements,
        ConditionsContainerType& rConditions,
        const ProcessInfo& CurrentProcessInfo
    ) override
    {
        if (mpPhi == nullptr)
            KRATOS_ERROR << "The projection operator is not set";

        BaseType::ResizeAndInitializeVectors(pA, pDx, pb, rElements, rConditions, CurrentProcessInfo);

        // resize the system accordingly
        SizeType reduced_system_size = TDenseSpaceType::Size2(*mpPhi);

        if (pA != nullptr)
        {
            TSystemMatrixType& A = *pA;
            if (TSparseSpaceType::Size1(A) != reduced_system_size || TSparseSpaceType::Size2(A) != reduced_system_size)
                TSparseSpaceType::Resize(A, reduced_system_size, reduced_system_size);
        }
        if (pb != nullptr)
        {
            TSystemVectorType& b = *pb;
            if (TSparseSpaceType::Size(b) != reduced_system_size)
                TSparseSpaceType::Resize(b, reduced_system_size);
        }
    }

    void SystemSolve(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
    ) override
    {
        if (mpPhi == nullptr)
            KRATOS_ERROR << "The projection operator is not set";

        SizeType reduced_system_size = TDenseSpaceType::Size2(*mpPhi);

        TSystemVectorType Dx;
        TSparseSpaceType::Resize(Dx, reduced_system_size);

        BaseType::SystemSolve(rA, Dx, rb);
        TDenseSpaceType::Mult(*mpPhi, Dx, rDx);
    }

    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& r_model_part,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
    ) override
    {
        if (mpPhi == nullptr)
            KRATOS_ERROR << "The projection operator is not set";

        SizeType full_system_size = TDenseSpaceType::Size1(*mpPhi);

        TSystemVectorType b;
        if (TSparseSpaceType::Size(b) != full_system_size)
            TSparseSpaceType::Resize(b, full_system_size);
        TSparseSpaceType::SetToZero(b);

        BaseType::CalculateReactions(pScheme, r_model_part, rA, rDx, b);
    }

    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */


    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    const LocalSystemMatrixType* mpPhi;

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class ProjectionBasedPodBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_PROJECTION_BASED_POD_BUILDER_AND_SOLVER_H_INCLUDED  defined */
