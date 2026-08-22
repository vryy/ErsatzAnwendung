//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 07 Aug 2026 $
//
//


// System includes


// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

// builder_and_solvers
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_deactivation.h"
#include "custom_strategies/builder_and_solvers/projection_based_pod_builder_and_solver.h"
#include "custom_strategies/schemes/element_weighting_scheme.h"
#include "custom_utilities/projection_based_pod_builder_and_solver_factory.h"
#include "custom_python/add_custom_strategies_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void ErsatzAnwendung_AddCustomStrategiesToPython()
{
    typedef UblasSpace<KRATOS_DOUBLE_TYPE, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<KRATOS_DOUBLE_TYPE, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType, ModelPart> LinearSolverType;

    typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType, ModelPart> BuilderAndSolverType;

    typedef Scheme<SparseSpaceType, LocalSpaceType, ModelPart> SchemeType;

    typedef ResidualBasedEliminationBuilderAndSolverDeactivation<SparseSpaceType, LocalSpaceType, LinearSolverType, ModelPart> ResidualBasedEliminationBuilderAndSolverDeactivationType;

    typedef ProjectionBasedPodBuilderAndSolver<ResidualBasedEliminationBuilderAndSolverDeactivationType> ProjectionBasedPODResidualBasedEliminationBuilderAndSolverDeactivationType;

    typedef ElementWeightingScheme<SparseSpaceType, LocalSpaceType, ModelPart> ElementWeightingSchemeType;

    //********************************************************************
    //********************************************************************

    class_< ProjectionBasedPODResidualBasedEliminationBuilderAndSolverDeactivationType,
            ProjectionBasedPODResidualBasedEliminationBuilderAndSolverDeactivationType::Pointer,
            bases<ResidualBasedEliminationBuilderAndSolverDeactivationType>,
            boost::noncopyable
          > (
                "ProjectionBasedPODResidualBasedEliminationBuilderAndSolverDeactivation",
                init<typename LinearSolverType::Pointer>()
            )
    .def("SetPodProcess", &ProjectionBasedPODResidualBasedEliminationBuilderAndSolverDeactivationType::SetPodProcess)
    .def("GetPodProcess", &ProjectionBasedPODResidualBasedEliminationBuilderAndSolverDeactivationType::pGetPodProcess)
    ;

    //********************************************************************
    //********************************************************************

    class_<ProjectionBasedPodBuilderAndSolverFactory, ProjectionBasedPodBuilderAndSolverFactory::Pointer, boost::noncopyable>
    ("ProjectionBasedPodBuilderAndSolverFactory", init<>())
    .def("Create", &ProjectionBasedPodBuilderAndSolverFactory::Create<ResidualBasedEliminationBuilderAndSolverDeactivationType>)
    .staticmethod("Create");
    ;

    //********************************************************************
    //********************************************************************

    class_<ElementWeightingSchemeType, ElementWeightingSchemeType::Pointer, bases<SchemeType>, boost::noncopyable>
    ("ElementWeightingScheme", init<typename SchemeType::Pointer>())
    .def("SetElementWeight", &ElementWeightingSchemeType::SetElementWeight)
    ;
}

}  // namespace Python.

} // Namespace Kratos
