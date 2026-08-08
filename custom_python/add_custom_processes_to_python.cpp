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

// builder_and_solvers
#include "custom_processes/pod_process.h"
#include "custom_processes/snapshot_collecting_process.h"
#include "custom_processes/pod_mode_reading_process.h"
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void ErsatzAnwendung_AddCustomProcessesToPython()
{
    class_<PodProcess, PodProcess::Pointer, bases<Process>, boost::noncopyable>("PodProcess")
    ;

    class_<SnapshotCollectingProcess, SnapshotCollectingProcess::Pointer, bases<PodProcess>, boost::noncopyable>
    ("SnapshotCollectingProcess", init<>())
    .def("SavePrincipalComponents", &SnapshotCollectingProcess::SavePrincipalComponents)
    ;

    class_<PodModeReadingProcess, PodModeReadingProcess::Pointer, bases<PodProcess>, boost::noncopyable>
    ("PodModeReadingProcess", init<const std::string&>())
    ;
}

}  // namespace Python.

} // Namespace Kratos
