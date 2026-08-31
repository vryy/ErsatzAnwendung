//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 07, 2026 $
//
//

#if !defined(KRATOS_POD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_POD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "custom_strategies/builder_and_solvers/pod_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/projection_based_pod_builder_and_solver.h"

namespace Kratos
{

struct PodBuilderAndSolverFactory
{
    KRATOS_CLASS_POINTER_DEFINITION( PodBuilderAndSolverFactory );

    template<class TBuilderAndSolverType>
    static typename PodBuilderAndSolver<TBuilderAndSolverType>::Pointer Create(typename TBuilderAndSolverType::Pointer pBuilderAndSolver)
    {
        if (pBuilderAndSolver == nullptr)
            KRATOS_ERROR << "The builder and solver is nullptr";
        typedef PodBuilderAndSolver<TBuilderAndSolverType> PodBuilderAndSolverType;
        return typename PodBuilderAndSolverType::Pointer(new PodBuilderAndSolverType(*pBuilderAndSolver));
    }
};

struct ProjectionBasedPodBuilderAndSolverFactory
{
    KRATOS_CLASS_POINTER_DEFINITION( ProjectionBasedPodBuilderAndSolverFactory );

    template<class TBuilderAndSolverType>
    static typename ProjectionBasedPodBuilderAndSolver<TBuilderAndSolverType>::Pointer Create(typename TBuilderAndSolverType::Pointer pBuilderAndSolver)
    {
        if (pBuilderAndSolver == nullptr)
            KRATOS_ERROR << "The builder and solver is nullptr";
        typedef ProjectionBasedPodBuilderAndSolver<TBuilderAndSolverType> PodBuilderAndSolverType;
        return typename PodBuilderAndSolverType::Pointer(new PodBuilderAndSolverType(*pBuilderAndSolver));
    }
};

} // end namespace Kratos

#endif // KRATOS_POD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED
