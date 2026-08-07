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
#include "custom_strategies/builder_and_solvers/projection_based_pod_builder_and_solver.h"

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

    typedef ProjectionBasedPODBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType, ModelPart> ProjectionBasedPODBuilderAndSolverType;

    //********************************************************************
    //********************************************************************

    class_< ProjectionBasedPODBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable>
    ("ProjectionBasedPODBuilderAndSolver", init<typename BuilderAndSolverType::Pointer>())
    ;
}

}  // namespace Python.

} // Namespace Kratos
