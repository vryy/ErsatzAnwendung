/* ***********************************************************
 *
 *   Created by:          $Author: hbui $
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: Sep 01, 2026 $
 *
 * ***********************************************************/

#if !defined(KRATOS_MULTIPHASEFLOW_POD_CRITERIA )
#define  KRATOS_MULTIPHASEFLOW_POD_CRITERIA

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */


/*@} */
/**@name Kratos Classes */
/*@{ */

/**
 * Define a convergence criteria for multiphase problem, including water, air and displacement.
 * This criteria is used for POD-based reduced order model, where the residual vector is considered
 * in reduced form. Hence, the energy norm of each phase is not checked individually.
 */
template<class TSparseSpace,
         class TDenseSpace,
         class TModelPartType
         >
class MultiPhaseFlowPodCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace, TModelPartType >
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( MultiPhaseFlowPodCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace, TModelPartType > BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::ValueType ValueType;

    typedef typename BaseType::ModelPartType ModelPartType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    MultiPhaseFlowPodCriteria(ValueType RelativeTolerance, ValueType AbsoluteTolerance)
        : BaseType()
    {
        mRelativeTolerance = RelativeTolerance;
        mAbsoluteTolerance = AbsoluteTolerance;

        mCheckType = 4;
    }

    /** Destructor.
     */
    ~MultiPhaseFlowPodCriteria() override
    {}

    /// Set the type of convergence criteria to be checked for each phase
    //  Type = 1 - only check the absolute criteria
    //  Type = 2 - only check the relative criteria
    //  Type = 3 - check both absolute and relative criteria
    //  Type = 4 - check either absolute or relative criteria
    void SetType(const int Type)
    {
        mCheckType = Type;
    }

    /*@} */
    /**@name Operators
     */
    /*@{ */

    /*Criterias that need to be called before getting the solution */
    bool PreCriteria(
        ModelPartType& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) override
    {
        return true;
    }

    /*Criterias that need to be called after getting the solution */
    bool PostCriteria(
        ModelPartType& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) override
    {
        if (Dx.size() != 0) // if we are solving for something
        {
            double norm_Dx = 0.0;
            double norm_b = TSparseSpace::TwoNorm(b);
            double norm_x = 0.0;
            double norm_Dx_WATER = 0.0;
            double norm_x_WATER = 0.0;
            double norm_Dx_AIR = 0.0;
            double norm_x_AIR = 0.0;

            bool HasDisplacement = false;
            bool HasWaterPres = false;
            bool HasAirPres = false;

            for (auto i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof)
            {
                if (i_dof->IsFree())
                {
                    const auto eq_id = i_dof->EquationId();

                    if (i_dof->GetVariable() == DISPLACEMENT_X)
                    {
                        HasDisplacement = true;
                        double aux = i_dof->GetSolutionStepValue(DISPLACEMENT_X);

                        norm_Dx += Dx[eq_id] * Dx[eq_id];
                        norm_x += aux * aux;
                    }
                    if (i_dof->GetVariable() == DISPLACEMENT_Y)
                    {
                        HasDisplacement = true;
                        double aux = i_dof->GetSolutionStepValue(DISPLACEMENT_Y);

                        norm_Dx += Dx[eq_id] * Dx[eq_id];
                        norm_x += aux * aux;
                    }
                    if (i_dof->GetVariable() == DISPLACEMENT_Z)
                    {
                        HasDisplacement = true;
                        double aux = i_dof->GetSolutionStepValue(DISPLACEMENT_Z);

                        norm_Dx += Dx[eq_id] * Dx[eq_id];
                        norm_x += aux * aux;
                    }
                    if (i_dof->GetVariable() == WATER_PRESSURE)
                    {
                        HasWaterPres = true;
                        double aux = i_dof->GetSolutionStepValue(WATER_PRESSURE);

                        norm_Dx_WATER += Dx[eq_id] * Dx[eq_id];
                        norm_x_WATER += aux * aux;
                    }
                    if (i_dof->GetVariable() == AIR_PRESSURE)
                    {
                        HasAirPres = true;
                        double aux = i_dof->GetSolutionStepValue(AIR_PRESSURE);

                        norm_Dx_AIR += Dx[eq_id] * Dx[eq_id];
                        norm_x_AIR += aux * aux;
                    }
                }
            }

            norm_x = std::sqrt(norm_x);
            norm_Dx = std::sqrt(norm_Dx);

            if (HasWaterPres)
            {
                norm_x_WATER = std::sqrt(norm_x_WATER);
                norm_Dx_WATER = std::sqrt(norm_Dx_WATER);
            }

            if (HasAirPres)
            {
                norm_x_AIR = std::sqrt(norm_x_AIR);
                norm_Dx_AIR = std::sqrt(norm_Dx_AIR);
            }

            double ratioDisp = 1.0;

            double ratioWater = 1.0;

            double ratioAir = 1.0;

            if (norm_x > 0)
                ratioDisp = norm_Dx / norm_x;
            if (norm_Dx == 0.0)
                ratioDisp = 0.0;

            if (norm_x_WATER > 0)
                ratioWater = norm_Dx_WATER / norm_x_WATER;
            if (norm_Dx_WATER == 0.0)
                ratioWater = 0.0;

            if (norm_x_AIR > 0)
                ratioAir = norm_Dx_AIR / norm_x_AIR;
            if (norm_Dx_AIR == 0.0)
                ratioAir = 0.0;

            std::cout << "***********************************************"
                      << "CONVERGENCE CRITERIA FOR MULTIPHASE PROBLEMS"
                      << "***********************************************"
                      << std::endl;
            std::cout.precision(6);
            std::cout.setf(std::ios::scientific);
            if (HasDisplacement)
            {
                std::cout << "** expected values: \t\t\t\t\t\tabs_tol = " << mAbsoluteTolerance
                          << "\trel_tol = " << mRelativeTolerance << "\t\t\t\t**" << std::endl;
                std::cout << "** obtained values displacement:\tratio = " << ratioDisp
                          << "\t||Dx|| = " << norm_Dx
                          << "\t||x|| = " << norm_x
                          << "\t\t\t\t**" << std::endl;
            }
            if (HasWaterPres)
            {
                std::cout << "** obtained values water pressure:\tratio = " << ratioWater
                          << "\t||Dx|| = " << norm_Dx_WATER
                          << "\t||x|| = " << norm_x_WATER
                          << "\t\t\t\t**" << std::endl;
                if (HasAirPres)
                {
                    std::cout << "** obtained values air pressure:\tratio = " << ratioAir
                              << "\t||Dx|| = " << norm_Dx_AIR
                              << "\t||x|| = " << norm_x_AIR
                              << "\t\t\t\t**" << std::endl;

                    std::cout << "** obtained values total:\t\tratio = " << ratioAir + ratioWater + ratioDisp
                              << "\tchange = " << norm_Dx + norm_Dx_WATER + norm_Dx_AIR
                              << "\tabsolute = " << norm_x_AIR + norm_x_WATER + norm_x
                              << "\tenergy = " << norm_b
                              << "\t**" << std::endl;
                }
                else
                {
                    std::cout << "** obtained values total:\t\tratio = " << ratioWater + ratioDisp
                              << "\tchange = " << norm_Dx + norm_Dx_WATER
                              << "\tabsolute = " << norm_x_WATER + norm_x
                              << "\tenergy = " << norm_b
                              << "\t**" << std::endl;
                }
            }
            else
            {
                if (HasAirPres)
                {
                    std::cout << "** obtained values air pressure:\tratio = " << ratioAir
                              << "\t||Dx|| = " << norm_Dx_AIR
                              << "\t||x|| = " << norm_x_AIR
                              << "\t\t\t\t**" << std::endl;

                    std::cout << "** obtained values total:\t\tratio = " << ratioAir + ratioDisp
                              << "\tchange = " << norm_Dx + norm_Dx_AIR
                              << "\tabsolute = " << norm_x_AIR + norm_x
                              << "\tenergy = " << norm_b
                              << "\t**" << std::endl;
                }
                else
                {
                    std::cout << "** obtained values total:\t\tratio = " << ratioWater + ratioDisp
                              << "\tchange = " << norm_Dx
                              << "\tabsolute = " << norm_x
                              << "\tenergy = " << norm_b
                              << "\t**" << std::endl;
                }
            }
            std::cout << "*****************************************************"
                      << "*****************************************************"
                      << "********************************" << std::endl;

            bool all_converged = (norm_b <= mAbsoluteTolerance);

            bool disp_reason_2 = false;
            int disp_reason_2_case = 0;
            if (ratioDisp <= mRelativeTolerance)
            {
                disp_reason_2 = true;
                disp_reason_2_case = 1;
            }
            else
            {
                if (norm_x < mAbsoluteTolerance) // case that zero state is solved
                {
                    if (norm_Dx < mAbsoluteTolerance)
                    {
                        disp_reason_2 = true;
                        disp_reason_2_case = 2;
                    }
                }
            }
            bool disp_converged = false;

            if (mCheckType == 1)
                disp_converged = all_converged;
            else if (mCheckType == 2)
                disp_converged = disp_reason_2;
            else if (mCheckType == 3)
                disp_converged = all_converged && disp_reason_2;
            else if (mCheckType == 4)
                disp_converged = all_converged || disp_reason_2;

            bool water_converged = false;
            bool water_reason_2 = false;
            if(HasWaterPres)
            {
                water_reason_2 = (ratioWater <= mRelativeTolerance);
                if (mCheckType == 1)
                    water_converged = all_converged;
                else if (mCheckType == 2)
                    water_converged = water_reason_2;
                else if (mCheckType == 3)
                    water_converged = all_converged && water_reason_2;
                else if (mCheckType == 4)
                    water_converged = all_converged || water_reason_2;
            }
            else
                water_converged = true;

            bool air_converged = false;
            bool air_reason_2 = false;
            if(HasAirPres)
            {
                air_reason_2 = (ratioAir <= mRelativeTolerance);
                if (mCheckType == 1)
                    air_converged = all_converged;
                else if (mCheckType == 2)
                    air_converged = air_reason_2;
                else if (mCheckType == 3)
                    air_converged = all_converged && air_reason_2;
                else if (mCheckType == 4)
                    air_converged = all_converged || air_reason_2;
            }
            else
                air_converged = true;

            if(disp_converged && water_converged && air_converged)
            {
                std::cout << "Congratulations the solution strategy is converged." << std::endl;
                std::cout << "Reason for converged displacement:";
                if (mCheckType == 1)
                {
                    std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                }
                else if (mCheckType == 2)
                {
                    if (disp_reason_2_case == 1)
                        std::cout << " {(||Dx||/||x|| = " << ratioDisp << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    else if (disp_reason_2_case == 2)
                        std::cout << " {(||Dx|| = " << norm_Dx << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                  << " and {(||x|| = " << norm_x << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                  << std::endl;
                }
                else if (mCheckType == 3)
                {
                    std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                    if (disp_reason_2 == 1)
                        std::cout << "and {(||Dx||/||x|| = " << ratioDisp << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    else if (disp_reason_2_case == 2)
                        std::cout << "and {(||Dx|| = " << norm_Dx << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                  << " and {(||x|| = " << norm_x << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                  << std::endl;
                }
                else if (mCheckType == 4)
                {
                    if(all_converged) std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                    if (disp_reason_2)
                    {
                        if (disp_reason_2_case == 1)
                            std::cout << " {(||Dx||/||x|| = " << ratioDisp << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                        else if (disp_reason_2_case == 2)
                            std::cout << " {(||Dx|| = " << norm_Dx << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                      << " and {(||x|| = " << norm_x << ") <= (abs_tol = " << mAbsoluteTolerance << ")}"
                                      << std::endl;
                    }
                }

                if(HasWaterPres)
                {
                    std::cout << "Reason for converged water pressure:";
                    if (mCheckType == 1)
                    {
                        std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 2)
                    {
                        std::cout << " {(||Dx||/||x|| (water) = " << ratioWater << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 3)
                    {
                        std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mRelativeTolerance << ")}" << std::endl;
                        std::cout << "and {(||Dx||/||x|| (water) = " << ratioWater << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 4)
                    {
                        if(all_converged) std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                        if(water_reason_2) std::cout << " {(||Dx||/||x|| (water) = " << ratioWater << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                }

                if(HasAirPres)
                {
                    std::cout << "Reason for converged air pressure:";
                    if (mCheckType == 4)
                    {
                        std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 2)
                    {
                        std::cout << " {(||Dx||/||x|| (air) = " << ratioAir << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 3)
                    {
                        std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                        std::cout << "and {(||Dx||/||x|| = " << ratioAir << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                    else if (mCheckType == 4)
                    {
                        if(all_converged) std::cout << " {(||b|| = " << norm_b << ") <= (expected ||b|| = " << mAbsoluteTolerance << ")}" << std::endl;
                        if(air_reason_2) std::cout << " {(||Dx||/||x|| (air) = " << ratioAir << ") <= (expected ||Dx||/||x|| = " << mRelativeTolerance << ")}" << std::endl;
                    }
                }

                return true;
            }
            else
                return false;
        }
        else // in this case all the displacements are imposed!
        {
            return true;
        }
    }

    void Initialize(
        ModelPartType& r_model_part
    ) override
    {
    }

    void InitializeSolutionStep(
        ModelPartType& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) override
    {
    }

    void FinalizeSolutionStep(
        ModelPartType& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    ) override
    {
    }

    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "MultiPhaseFlowPodCriteria";
    }

    ///@}

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */


    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */

    /*@{ */

    ValueType mRelativeTolerance;
    ValueType mAbsoluteTolerance;
    int mCheckType;

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_MULTIPHASEFLOW_POD_CRITERIA  defined */
