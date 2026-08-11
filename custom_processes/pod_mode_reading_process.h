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


namespace Kratos
{

/**
 * This process reads the POD mode from file
 */
template<class TSparseSpace, class TModelPart>
class PodModeReadingProcess : public PodProcess<TSparseSpace, TModelPart>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( PodModeReadingProcess );

    typedef PodProcess<TSparseSpace, TModelPart> BaseType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    PodModeReadingProcess(const std::string& filename) : BaseType()
    {
        this->ReadPrincipalComponents(mPhi, filename);
    }

    void ApplyProjection(TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb) override
    {
        BaseType::ApplyProjection(mPhi, rA, rDx, rb);
    }

private:

    /**
     * Read in the principal vectors stored in the data file
     */
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

    Matrix mPhi; // projection matrix

}; /* Class PodModeReadingProcess */

}  /* namespace Kratos.*/

#endif /* KRATOS_ERSATZ_ANWENDUNG_POD_MODE_READING_PROCESS_H_INCLUDED  defined */
