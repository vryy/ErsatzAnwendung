//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Jan 21, 2016$
//
//


// System includes


// External includes


// Project includes
#include "ersatz_anwendung.h"
#include "ersatz_anwendung_variables.h"

namespace Kratos
{
    // create the application variables here

    // constructor
    KratosErsatzAnwendung::KratosErsatzAnwendung()
        : KratosApplication("ErsatzAnwendung")
    {
        KRATOS_REGISTER_VARIABLE( ELEMENT_WEIGHT )
    }

    // register the application to the Kratos kernel
    void KratosErsatzAnwendung::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosErsatzAnwendung..." << std::endl;
    }

    bool KratosErsatzAnwendung::Has(const std::string& feature)
    {
#ifdef ERSATZ_APP_USE_MATIO
        if (feature == "Matio")
            return true;
#endif
        return false;
    }

} // namespace Kratos
