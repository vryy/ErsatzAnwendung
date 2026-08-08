//
//   Project Name:        KratosErsatzAnwendung
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Aug 07, 2026 $
//
//


#if !defined(KRATOS_POD_PROCESS_H_INCLUDED )
#define  KRATOS_POD_PROCESS_H_INCLUDED


/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"


namespace Kratos
{

/**
 * This base class provides abstract functions for POD operations
 */
class PodProcess : public Process
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( PodProcess );

    typedef Process BaseType;
    typedef ModelPart::DofsArrayType DofsArrayType;

    PodProcess() : BaseType()
    {
    }

    /**
     * Assign the dof set
     */
    void SetDofSet(const DofsArrayType& rDofSet)
    {
        mpDofSet = &rDofSet;
    }

    /**
     * Get the dof set
     */
    const DofsArrayType& GetDofSet() const
    {
        if (mpDofSet == nullptr)
            KRATOS_ERROR << "Dof set is not assigned";
        return *mpDofSet;
    }

    /**
     * Initialize the projection matrix
     */
    virtual void InitializeProjectionMatrix(Matrix& rPhi)
    {}

    /**
     * Compute the projection matrix
     */
    virtual void UpdateProjectionMatrix(Matrix& rPhi)
    {}

protected:

    Vector TakeSnapshot() const
    {
        const auto& DofSet = this->GetDofSet();
        const std::size_t SystemSize = DofSet.size();
        Vector dof_vector(SystemSize);

        std::size_t i = 0;
        for (auto it = DofSet.begin(); it != DofSet.end(); ++it, ++i)
        {
            dof_vector[i] = it->GetSolutionStepValue();
        }

        return dof_vector;
    }

private:

    const DofsArrayType* mpDofSet = nullptr;

}; /* Class PodProcess */

}  /* namespace Kratos.*/

#endif /* KRATOS_POD_PROCESS_H_INCLUDED  defined */
