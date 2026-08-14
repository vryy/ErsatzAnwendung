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
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

// processes
#include "custom_processes/pod_process.h"
#include "custom_processes/snapshot_collecting_process.h"
#include "custom_processes/rayleigh_ritz_projection_process.h"
#include "custom_processes/petrov_galerkin_projection_process.h"
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void ErsatzAnwendung_AddCustomProcessesToPython()
{
    typedef UblasSpace<KRATOS_DOUBLE_TYPE, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<KRATOS_DOUBLE_TYPE, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType, ModelPart> LinearSolverType;

    typedef PodProcess<SparseSpaceType, LocalSpaceType, ModelPart> PodProcessType;

    typedef SnapshotCollectingProcess<SparseSpaceType, LocalSpaceType, ModelPart> SnapshotCollectingProcessType;

    typedef RayleighRitzProjectionProcess<SparseSpaceType, LocalSpaceType, ModelPart> RayleighRitzProjectionProcessType;
    typedef PetrovGalerkinProjectionProcess<SparseSpaceType, LocalSpaceType, ModelPart> PetrovGalerkinProjectionProcessType;

    class_<PodProcessType, typename PodProcessType::Pointer, bases<Process, LinearSolverType>, boost::noncopyable>
    ("PodProcess", init<>())
    ;

    class_<SnapshotCollectingProcessType, typename SnapshotCollectingProcessType::Pointer, bases<PodProcessType>, boost::noncopyable>
    ("SnapshotCollectingProcess", init<typename LinearSolverType::Pointer>())
    .def("SavePrincipalComponents", &SnapshotCollectingProcessType::SavePrincipalComponents)
    ;

    class_<RayleighRitzProjectionProcessType, typename RayleighRitzProjectionProcessType::Pointer, bases<PodProcessType>, boost::noncopyable>
    ("RayleighRitzProjectionProcess", init<const std::string&>())
    ;

    class_<PetrovGalerkinProjectionProcessType, typename PetrovGalerkinProjectionProcessType::Pointer, bases<PodProcessType>, boost::noncopyable>
    ("PetrovGalerkinProjectionProcess", init<const std::string&, const std::string&>())
    ;
}

}  // namespace Python.

} // Namespace Kratos
