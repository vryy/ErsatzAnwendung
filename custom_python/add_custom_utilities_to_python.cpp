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

typedef MatrixVectorTypeSelector<int>::VectorType IntegerVector;

boost::python::list POD_Utils_SVD(POD_Utils& rDummy, const Matrix& rA)
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

void POD_Utils_WriteMat(POD_Utils& rDummy, const std::string& filename, const std::string& variable_name, const Matrix& rA, bool append)
{
#ifdef ERSATZ_APP_USE_MATIO
    POD_Utils::WriteMat(filename, variable_name, rA, append);
#else
    KRATOS_ERROR << "MATIO support is not enabled";
#endif
}

void POD_Utils_WriteVec(POD_Utils& rDummy, const std::string& filename, const std::string& variable_name, const Vector& rb, bool append)
{
#ifdef ERSATZ_APP_USE_MATIO
    POD_Utils::WriteVec(filename, variable_name, rb, append);
#else
    KRATOS_ERROR << "MATIO support is not enabled";
#endif
}

void POD_Utils_WriteIntVec(POD_Utils& rDummy, const std::string& filename, const std::string& variable_name,
        const IntegerVector& rb, bool append)
{
#ifdef ERSATZ_APP_USE_MATIO
    POD_Utils::WriteIntegerVec(filename, variable_name, rb, append);
#else
    KRATOS_ERROR << "MATIO support is not enabled";
#endif
}

Matrix POD_Utils_ReadMat(POD_Utils& rDummy, const std::string& filename, const std::string& variable_name)
{
#ifdef ERSATZ_APP_USE_MATIO
    return POD_Utils::ReadMat<Matrix>(filename, variable_name);
#else
    KRATOS_ERROR << "MATIO support is not enabled";
    return Matrix();
#endif
}

Vector POD_Utils_ReadVec(POD_Utils& rDummy, const std::string& filename, const std::string& variable_name)
{
#ifdef ERSATZ_APP_USE_MATIO
    return POD_Utils::ReadVec<Vector>(filename, variable_name);
#else
    KRATOS_ERROR << "MATIO support is not enabled";
    return Vector();
#endif
}

IntegerVector POD_Utils_ReadIntVec(POD_Utils& rDummy, const std::string& filename, const std::string& variable_name)
{
#ifdef ERSATZ_APP_USE_MATIO
    return POD_Utils::ReadVec<IntegerVector>(filename, variable_name);
#else
    KRATOS_ERROR << "MATIO support is not enabled";
    return Vector();
#endif
}

void ErsatzAnwendung_AddCustomUtilitiesToPython()
{
    class_<POD_Utils, POD_Utils::Pointer, boost::noncopyable>
    ("POD_Utils", init<>())
    .def("SVD", POD_Utils_SVD)
#ifdef ERSATZ_APP_USE_MATIO
    .def("ListVariables", POD_Utils::ListVariables)
    .staticmethod("ListVariables")
#endif
    .def("WriteMat", POD_Utils_WriteMat)
    .def("WriteVec", POD_Utils_WriteVec)
    .def("WriteIntVec", POD_Utils_WriteIntVec)
    .def("ReadMat", POD_Utils_ReadMat)
    .def("ReadVec", POD_Utils_ReadVec)
    .def("ReadIntVec", POD_Utils_ReadIntVec)
    ;
}

}  // namespace Python.

}  // namespace Kratos.
