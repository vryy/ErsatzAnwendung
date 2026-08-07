//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 10 Aug 2017 $
//
//


#if !defined(KRATOS_ERSATZ_ANWENDUNG_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ERSATZ_ANWENDUNG_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED


// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"

namespace Kratos
{

namespace Python
{

void  ErsatzAnwendung_AddCustomUtilitiesToPython();

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ERSATZ_ANWENDUNG_ADD_CUSTOM_UTILITIES_TO_PYTHON_H_INCLUDED  defined
