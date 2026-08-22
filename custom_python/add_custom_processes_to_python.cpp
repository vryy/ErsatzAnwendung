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
#include "custom_processes/ecsw_snapshot_collecting_process.h"
#include "custom_processes/rayleigh_ritz_projection_process.h"
#include "custom_processes/petrov_galerkin_projection_process.h"
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

template<class TProcessType>
boost::python::dict EcswSnapshotCollectingProcess_ConstructSystem(TProcessType& self,
    Matrix& rG, Vector& rb, std::size_t number_of_modes)
{
    typedef typename TProcessType::IndexType IndexType;
    std::map<IndexType, IndexType> element_weight_index;
    self.ConstructSystem(rG, rb, element_weight_index, number_of_modes);
    boost::python::dict result;
    for (auto it = element_weight_index.begin(); it != element_weight_index.end(); ++it)
        result[it->first] = it->second;
    return result;
}

void ErsatzAnwendung_AddCustomProcessesToPython()
{
    typedef UblasSpace<KRATOS_DOUBLE_TYPE, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<KRATOS_DOUBLE_TYPE, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType, ModelPart> LinearSolverType;

    typedef PodProcess<SparseSpaceType, LocalSpaceType, ModelPart> PodProcessType;

    typedef SnapshotCollectingProcess<SparseSpaceType, LocalSpaceType, ModelPart> SnapshotCollectingProcessType;
    typedef EcswSnapshotCollectingProcess<SparseSpaceType, LocalSpaceType, ModelPart> EcswSnapshotCollectingProcessType;

    typedef RayleighRitzProjectionProcess<SparseSpaceType, LocalSpaceType, ModelPart> RayleighRitzProjectionProcessType;
    typedef PetrovGalerkinProjectionProcess<SparseSpaceType, LocalSpaceType, ModelPart> PetrovGalerkinProjectionProcessType;

    class_<PodProcessType, typename PodProcessType::Pointer, bases<Process, LinearSolverType>, boost::noncopyable>
    ("PodProcess", init<>())
    ;

    class_<SnapshotCollectingProcessType, typename SnapshotCollectingProcessType::Pointer, bases<PodProcessType>, boost::noncopyable>
    ("SnapshotCollectingProcess", init<typename LinearSolverType::Pointer>())
    .def("GetPrincipalValues", &SnapshotCollectingProcessType::GetPrincipalValues)
    .def("GetPrincipalComponents", &SnapshotCollectingProcessType::GetPrincipalComponents)
    .def("SavePrincipalComponents", &SnapshotCollectingProcessType::SavePrincipalComponents)
    ;

    class_<EcswSnapshotCollectingProcessType, typename EcswSnapshotCollectingProcessType::Pointer, bases<SnapshotCollectingProcessType>, boost::noncopyable>
    ("EcswSnapshotCollectingProcess", init<typename LinearSolverType::Pointer>())
    .def("ConstructSystem", &EcswSnapshotCollectingProcess_ConstructSystem<EcswSnapshotCollectingProcessType>)
    ;

    class_<RayleighRitzProjectionProcessType, typename RayleighRitzProjectionProcessType::Pointer, bases<PodProcessType>, boost::noncopyable>
    ("RayleighRitzProjectionProcess", init<const Matrix&>())
    ;

    class_<PetrovGalerkinProjectionProcessType, typename PetrovGalerkinProjectionProcessType::Pointer, bases<PodProcessType>, boost::noncopyable>
    ("PetrovGalerkinProjectionProcess", init<const Matrix&, const Matrix&>())
    ;
}

}  // namespace Python.

} // Namespace Kratos
