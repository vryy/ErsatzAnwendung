//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 07, 2026 $
//
//


#include "includes/define.h"
#include "includes/kratos_exception.h"
#if !defined(KRATOS_PROJECTION_BASED_POD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_PROJECTION_BASED_POD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED


/* System includes */

/* External includes */

/* Project includes */
#include "custom_strategies/builder_and_solvers/projection_based_pod_builder_and_solver.h"

namespace Kratos
{

struct ProjectionBasedPodBuilderAndSolverFactory
{
	KRATOS_CLASS_POINTER_DEFINITION( ProjectionBasedPodBuilderAndSolverFactory );

	template<class TBuilderAndSolverType>
	static typename ProjectionBasedPodBuilderAndSolver<TBuilderAndSolverType>::Pointer Create(typename TBuilderAndSolverType::Pointer pBuilderAndSolver)
	{
		if (pBuilderAndSolver == nullptr)
			KRATOS_ERROR << "The builder and solver is nullptr";
		typedef ProjectionBasedPodBuilderAndSolver<TBuilderAndSolverType> ProjectionBasedPodBuilderAndSolverType;
		return typename ProjectionBasedPodBuilderAndSolverType::Pointer(new ProjectionBasedPodBuilderAndSolverType(*pBuilderAndSolver));
	}
};

} // end namespace Kratos

#endif // KRATOS_PROJECTION_BASED_POD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED
