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
#include "custom_processes/pod_process.h"
#include "custom_processes/snapshot_collecting_process.h"


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

    typedef PodProcess<TSparseSpaceType, TDenseSpaceType, ModelPartType> PodProcessType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Default constructor.
     */
    ProjectionBasedPodBuilderAndSolver(typename TLinearSolverType::Pointer pLinearSystemSolver)
    : BaseType(pLinearSystemSolver)
    {
        /// The default PodProcess is the one collecting snapshots only. On the solve, it relies on FOM linear solver.
        this->SetPodProcess(typename PodProcessType::Pointer(
            new SnapshotCollectingProcess<TSparseSpaceType, TDenseSpaceType, ModelPartType>(pLinearSystemSolver)));
    }

    /** Copy constructor.
     */
    ProjectionBasedPodBuilderAndSolver(const BaseType& rOther)
    : BaseType(rOther)
    {
        /// The default PodProcess is the one collecting snapshots only. On the solve, it relies on FOM linear solver.
        this->SetPodProcess(typename PodProcessType::Pointer(
            new SnapshotCollectingProcess<TSparseSpaceType, TDenseSpaceType, ModelPartType>(rOther.GetLinearSystemSolver())));
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

    void SetUpSystem(ModelPartType& rModelPart) override
    {
        BaseType::SetUpSystem(rModelPart);

        mpPodProcess->SetModelPart(rModelPart);
        mpPodProcess->SetDofSet(BaseType::GetDofSet());
        mpPodProcess->ExecuteInitialize();
    }

    void InitializeSolutionStep(
        ModelPartType& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
    ) override
    {
        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);
        mpPodProcess->ExecuteInitializeSolutionStep();
    }

    void FinalizeSolutionStep(
        ModelPartType& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
    ) override
    {
        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);
        mpPodProcess->ExecuteFinalizeSolutionStep();
    }

    void Clear() override
    {
        BaseType::Clear();
    }

    /*@} */
    /**@name Access */
    /*@{ */

    void SetPodProcess(typename PodProcessType::Pointer pPodProcess)
    {
        mpPodProcess = pPodProcess;
        BaseType::SetLinearSystemSolver(pPodProcess);
    }

    typename PodProcessType::Pointer pGetPodProcess() const
    {
        return mpPodProcess;
    }

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

    typename PodProcessType::Pointer mpPodProcess;
        // The PodBuilderAndSolver maintains a pointer of PodProcess, also as the linear solver.
        // By this way, the methods BuildAndSolve of the base class can be reused non-intrusively.

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
