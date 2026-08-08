//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 07, 2026 $
//
//


#if !defined(KRATOS_ERSATZ_ANWENDUNG_SNAPSHOT_COLLECTING_PROCESS_H_INCLUDED )
#define  KRATOS_ERSATZ_ANWENDUNG_SNAPSHOT_COLLECTING_PROCESS_H_INCLUDED


/* System includes */

/* External includes */

/* Project includes */
#include "custom_processes/pod_process.h"
#include "custom_utilities/pod_utils.h"


namespace Kratos
{

/**
 * This process only collects the snapshot and computing the POD basis after the simulation ends
 */
class SnapshotCollectingProcess : public PodProcess
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( SnapshotCollectingProcess );

    typedef PodProcess BaseType;

    SnapshotCollectingProcess() : BaseType()
    {
    }

    void ExecuteInitialize() override
    {
        // record the initial condition of the system
        mSnapshot.push_back(BaseType::TakeSnapshot());
    }

    void InitializeProjectionMatrix(Matrix& rPhi) override
    {
        const std::size_t EquationSystemSize = this->GetDofSet().size();
        if (rPhi.size1() != EquationSystemSize || rPhi.size2() != EquationSystemSize)
            rPhi.resize(EquationSystemSize, EquationSystemSize, false);
        noalias(rPhi) = IdentityMatrix(EquationSystemSize, EquationSystemSize);
    }

    void ExecuteFinalizeSolutionStep() override
    {
        mSnapshot.push_back(BaseType::TakeSnapshot());
    }

    void SavePrincipalComponents(const std::string& filename, const std::size_t number_of_modes) const
    {
        const std::size_t n = mSnapshot.size();
        if (n > 0)
        {
            const std::size_t m = mSnapshot[0].size();
            Matrix Q(m, n);

            for (std::size_t i = 0; i < n; ++i)
            {
                noalias(column(Q, i)) = mSnapshot[i];
            }

            Matrix U, VT;
            Vector S;
            POD_Utils::SVD(Q, U, S, VT);
            KRATOS_WATCH(S)

            //

            std::ofstream file(filename, std::ios::binary);
            file.write(reinterpret_cast<const char*>(&number_of_modes), sizeof(number_of_modes));
            file.write(reinterpret_cast<const char*>(&m), sizeof(m));

            for (std::size_t i = 0; i < number_of_modes; ++i)
            {
                Vector T(m);
                noalias(T) = column(U, i);
                file.write(reinterpret_cast<const char*>(&T[0]), T.size() * sizeof(double));
            }

            file.close();
        }
    }

private:

    std::vector<Vector> mSnapshot;

}; /* Class PodProcess */

}  /* namespace Kratos.*/

#endif /* KRATOS_ERSATZ_ANWENDUNG_SNAPSHOT_COLLECTING_PROCESS_H_INCLUDED  defined */
