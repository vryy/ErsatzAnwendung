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
#include "custom_strategies/builder_and_solvers/pod_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/projection_based_pod_builder_and_solver.h"
#include "custom_strategies/schemes/element_weighting_scheme.h"
#include "custom_strategies/schemes/rayleigh_ritz_projection_scheme.h"
#include "custom_strategies/convergencecriterias/multiphaseflow_pod_criteria.h"
#include "custom_utilities/pod_builder_and_solver_factory.h"
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

    typedef Scheme<SparseSpaceType, LocalSpaceType, ModelPart> SchemeType;

    typedef ElementWeightingScheme<SparseSpaceType, LocalSpaceType, ModelPart> ElementWeightingSchemeType;
    typedef RayleighRitzProjectionScheme<SparseSpaceType, LocalSpaceType, ModelPart> RayleighRitzProjectionSchemeType;

    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType, ModelPart> ConvergenceCriteriaBaseType;
    typedef MultiPhaseFlowPodCriteria< SparseSpaceType, LocalSpaceType, ModelPart> MultiPhaseFlowPodCriteriaType;

    typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType, ModelPart> BuilderAndSolverType;

    typedef ResidualBasedEliminationBuilderAndSolverDeactivation<SparseSpaceType, LocalSpaceType, LinearSolverType, ModelPart> ResidualBasedEliminationBuilderAndSolverDeactivationType;

    typedef PodBuilderAndSolver<ResidualBasedEliminationBuilderAndSolverDeactivationType> PodResidualBasedEliminationBuilderAndSolverDeactivationType;
    typedef ProjectionBasedPodBuilderAndSolver<ResidualBasedEliminationBuilderAndSolverDeactivationType> ProjectionBasedPodResidualBasedEliminationBuilderAndSolverDeactivationType;

    //********************************************************************
    //********************************************************************

    class_< MultiPhaseFlowPodCriteriaType,
            bases<ConvergenceCriteriaBaseType>,
            boost::noncopyable
          > (
                "MultiPhaseFlowPodCriteria",
                init<double, double>()
            )
    .def("SetType", &MultiPhaseFlowPodCriteriaType::SetType)
    ;

    //********************************************************************
    //********************************************************************

    class_< PodResidualBasedEliminationBuilderAndSolverDeactivationType,
            PodResidualBasedEliminationBuilderAndSolverDeactivationType::Pointer,
            bases<ResidualBasedEliminationBuilderAndSolverDeactivationType>,
            boost::noncopyable
          > (
                "PodResidualBasedEliminationBuilderAndSolverDeactivation",
                init<typename LinearSolverType::Pointer>()
            )
    .def("SetPodProcess", &PodResidualBasedEliminationBuilderAndSolverDeactivationType::SetPodProcess)
    .def("GetPodProcess", &PodResidualBasedEliminationBuilderAndSolverDeactivationType::pGetPodProcess)
    ;

    class_< ProjectionBasedPodResidualBasedEliminationBuilderAndSolverDeactivationType,
            ProjectionBasedPodResidualBasedEliminationBuilderAndSolverDeactivationType::Pointer,
            bases<ResidualBasedEliminationBuilderAndSolverDeactivationType>,
            boost::noncopyable
          > (
                "ProjectionBasedPodResidualBasedEliminationBuilderAndSolverDeactivation",
                init<>()
            )
    .def("SetProjectionOperator", &ProjectionBasedPodResidualBasedEliminationBuilderAndSolverDeactivationType::SetProjectionOperator)
    ;

    //********************************************************************
    //********************************************************************

    class_<PodBuilderAndSolverFactory, typename PodBuilderAndSolverFactory::Pointer, boost::noncopyable>
    ("PodBuilderAndSolverFactory", init<>())
    .def("Create", &PodBuilderAndSolverFactory::Create<ResidualBasedEliminationBuilderAndSolverDeactivationType>)
    .staticmethod("Create");
    ;

    class_<ProjectionBasedPodBuilderAndSolverFactory, typename ProjectionBasedPodBuilderAndSolverFactory::Pointer, boost::noncopyable>
    ("ProjectionBasedPodBuilderAndSolverFactory", init<>())
    .def("Create", &ProjectionBasedPodBuilderAndSolverFactory::Create<ResidualBasedEliminationBuilderAndSolverDeactivationType>)
    .staticmethod("Create");
    ;

    //********************************************************************
    //********************************************************************

    class_<ElementWeightingSchemeType, ElementWeightingSchemeType::Pointer, bases<SchemeType>, boost::noncopyable>
    ("ElementWeightingScheme", init<typename SchemeType::Pointer>())
    .def("SetElementWeight", &ElementWeightingSchemeType::SetElementWeight)
    .def("SetProjectionOperator", &ElementWeightingSchemeType::SetProjectionOperator)
    ;

    class_<RayleighRitzProjectionSchemeType, RayleighRitzProjectionSchemeType::Pointer, bases<SchemeType>, boost::noncopyable>
    ("RayleighRitzProjectionScheme", init<typename SchemeType::Pointer>())
    .def("SetProjectionOperator", &RayleighRitzProjectionSchemeType::SetProjectionOperator)
    ;
}

} // namespace Python.

} // Namespace Kratos
