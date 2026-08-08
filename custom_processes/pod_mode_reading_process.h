//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 07, 2026 $
//
//


#if !defined(KRATOS_ERSATZ_ANWENDUNG_POD_MODE_READING_PROCESS_H_INCLUDED )
#define  KRATOS_ERSATZ_ANWENDUNG_POD_MODE_READING_PROCESS_H_INCLUDED


/* System includes */

/* External includes */

/* Project includes */
#include "custom_processes/pod_process.h"
#include "custom_utilities/pod_utils.h"


namespace Kratos
{

/**
 * This process reads the POD mode from file
 */
class PodModeReadingProcess : public PodProcess
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( PodModeReadingProcess );

    typedef PodProcess BaseType;

    PodModeReadingProcess(const std::string& filename) : BaseType()
    {
        this->ReadPrincipalComponents(mPhi, filename);
    }

    void InitializeProjectionMatrix(Matrix& rPhi) override
    {
        if (rPhi.size1() != mPhi.size1() || rPhi.size2() != mPhi.size2())
            rPhi.resize(mPhi.size1(), mPhi.size2(), false);
        noalias(rPhi) = mPhi;
    }

private:

    void ReadPrincipalComponents(Matrix& Phi, const std::string& filename) const
    {
        std::ifstream file(filename, std::ios::binary);

        std::size_t number_of_modes, m;
        file.read(reinterpret_cast<char*>(&number_of_modes), sizeof(number_of_modes));
        file.read(reinterpret_cast<char*>(&m), sizeof(m));

        Phi.resize(m, number_of_modes, false);

        for (std::size_t i = 0; i < number_of_modes; ++i)
        {
            Vector T(m);
            file.read(reinterpret_cast<char*>(&T[0]), T.size() * sizeof(double));

            noalias(column(Phi, i)) = T;
        }

        file.close();
    }

    Matrix mPhi;

}; /* Class PodProcess */

}  /* namespace Kratos.*/

#endif /* KRATOS_ERSATZ_ANWENDUNG_POD_MODE_READING_PROCESS_H_INCLUDED  defined */
