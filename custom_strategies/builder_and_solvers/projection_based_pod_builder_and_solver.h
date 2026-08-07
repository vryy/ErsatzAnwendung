//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 07, 2026 $
//
//


#if !defined(KRATOS_PROJECTION_BASED_POD_BUILDER_AND_SOLVER_H_INCLUDED )
#define  KRATOS_PROJECTION_BASED_POD_BUILDER_AND_SOLVER_H_INCLUDED


/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "custom_utilities/pod_utils.h"


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
 * This builder and solver plays as a wrapper of other builder and solver to enable POD support.
 * The default mode is to construct the reduced linear system by means of projection technique.
 * Nevertheless, it can collect the snapshot and construct the reduced basis on demand.
 * Various ways to construct the reduces basis are supported.
 */
template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
         class TModelPartType
         >
class ProjectionBasedPODBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver, TModelPartType >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( ProjectionBasedPODBuilderAndSolver );

    typedef BuilderAndSolver<TSparseSpace,TDenseSpace, TLinearSolver, TModelPartType> BaseType;

    typedef typename BaseType::ModelPartType ModelPartType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TSparseSpaceType TSparseSpaceType;

    typedef typename BaseType::TLinearSolverType TLinearSolverType;

    typedef typename BaseType::TDataType TDataType;

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

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    ProjectionBasedPODBuilderAndSolver(typename BaseType::Pointer pBuilderAndSolver)
    : BaseType()
    {}

    /** Destructor.
    */
    ~ProjectionBasedPODBuilderAndSolver() override
    {}

    /*@} */
    /**@name Operators
    */
    /*@{ */


    /*@} */
    /**@name Operations */
    /*@{ */

    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPartType& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
    ) override
    {
        KRATOS_TRY

        const std::size_t rsize = mPhi.size2(); // size of reduced system

        // construct the reduced linear system

        LocalSystemMatrixType Ared;
        POD_Utils::VtKV(Ared, mPhi, rA);

        LocalSystemVectorType bred(rsize);
        noalias(bred) = prod(trans(mPhi), rb);

        // solve the reduced linear system

        LocalSystemVectorType xred(rsize);
        POD_Utils::Solve(Ared, xred, bred);

        // project back the solution to full space
        noalias(rDx) = prod(mPhi, xred);

        KRATOS_CATCH("")
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

    typename BaseType::Pointer mpBuilderAndSolver; // the computing BuiderAndSolver

    LocalSystemMatrixType mPhi; // the projection matrix

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

}; /* Class ProjectionBasedPODBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_PROJECTION_BASED_POD_BUILDER_AND_SOLVER_H_INCLUDED  defined */
