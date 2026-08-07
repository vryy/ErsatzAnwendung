//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Jan 20, 2017 $
//
//

#if !defined(KRATOS_ERSATZ_ANWENDUNG_H_INCLUDED)
#define KRATOS_ERSATZ_ANWENDUNG_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

namespace Kratos
{

    ///@name Kratos Globals
    ///@{

    // Variables definition

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Enum's
    ///@{

    ///@}
    ///@name Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// Short class definition.
    /** Detail class definition.
    */
    class KRATOS_API(ERSATZ_ANWENDUNG) KratosErsatzAnwendung : public KratosApplication
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosMultiphaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(KratosErsatzAnwendung);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        KratosErsatzAnwendung();

        /// Destructor.
        ~KratosErsatzAnwendung() override {}

        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        void Register() override;

        ///@}
        ///@name Access
        ///@{


        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const override
        {
            return "Application for POD & surrogate modelling";
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << Info();
            PrintData(rOStream);
        }

        ///// Print object's data.
        void PrintData(std::ostream& rOStream) const override
        {
            rOStream << "in KratosErsatzAnwendung:";
            KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
            rOStream << "Variables:" << std::endl;
            KratosComponents<VariableData>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Elements:" << std::endl;
            KratosComponents<Element>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Conditions:" << std::endl;
            KratosComponents<Condition>().PrintData(rOStream);
        }

        ///@}
        ///@name Friends
        ///@{


        ///@}

    protected:
        ///@name Protected static Member Variables
        ///@{


        ///@}
        ///@name Protected member Variables
        ///@{


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


        ///@}
        ///@name Protected  Access
        ///@{


        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
        ///@{


        ///@}

    private:
        ///@name Static Member Variables
        ///@{


        ///@}
        ///@name Member Variables
        ///@{


        ///@}
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{


        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.
        KratosErsatzAnwendung& operator=(KratosErsatzAnwendung const& rOther);

        /// Copy constructor.
        KratosErsatzAnwendung(KratosErsatzAnwendung const& rOther);

        ///@}

    }; // Class KratosErsatzAnwendung

    ///@}


    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}

} // namespace Kratos

#endif // KRATOS_ERSATZ_ANWENDUNG_H_INCLUDED defined
