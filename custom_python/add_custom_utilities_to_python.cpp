//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 10 Aug 2017 $
//
//


// Project includes
#include "includes/element.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/pod_utils.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

boost::python::list ErsatzAnwendung_SVD(POD_Utils& rDummy, const Matrix& rA)
{
    Matrix U;
    Vector S;
    Matrix VT;

    int error = POD_Utils::SVD(rA, U, S, VT);
    if(error != 0)
        KRATOS_ERROR << "SVD failed, error code = " << error;

    boost::python::list Output;
    Output.append(U);
    Output.append(S);
    Output.append(VT);

    return Output;
}

void ErsatzAnwendung_AddCustomUtilitiesToPython()
{
    class_<POD_Utils, POD_Utils::Pointer, boost::noncopyable>
    ("POD_Utils", init<>())
    .def("SVD", ErsatzAnwendung_SVD)
    ;
}

}  // namespace Python.

}  // namespace Kratos.
