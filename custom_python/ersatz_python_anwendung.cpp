//
//   Project Name:        KratosErsatzAnwendung
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Jan 20, 2017 $
//
//


// System includes


// External includes
#if defined(KRATOS_PYTHON)
#include <boost/python.hpp>


// Project includes
#include "includes/define_python.h"
#include "ersatz_anwendung.h"
#include "ersatz_anwendung_variables.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos
{

namespace Python
{

    using namespace boost::python;
    BOOST_PYTHON_MODULE(KratosErsatzAnwendung)
    {

        class_<KratosErsatzAnwendung, KratosErsatzAnwendung::Pointer,
               bases<KratosApplication>, boost::noncopyable>
               ("KratosErsatzAnwendung")
        .def("Has", &KratosErsatzAnwendung::Has)
            .staticmethod("Has")
        ;

        ErsatzAnwendung_AddCustomProcessesToPython();
        ErsatzAnwendung_AddCustomStrategiesToPython();
        ErsatzAnwendung_AddCustomUtilitiesToPython();

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELEMENT_WEIGHT )

    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
