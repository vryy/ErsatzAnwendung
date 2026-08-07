/* *********************************************************
*
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 10 Aug 17 $
*
* ***********************************************************/


#if !defined(KRATOS_POD_UTILS )
#define  KRATOS_POD_UTILS


/* System includes */

/* External includes */
#include <boost/numeric/bindings/traits/sparse_traits.hpp>
#include <boost/numeric/bindings/traits/matrix_traits.hpp>
#include <boost/numeric/bindings/traits/vector_traits.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>

/* Project includes */
#include "includes/ublas_interface.h"


// template for LAPACK function
extern "C" void dgesvd_(char* JOBU,
             char* JOBVT,
             int* M,
             int* N,
             double* A,
             int* LDA,
             double* S,
             double* U,
             int* LDU,
             double* VT,
             int* LDVT,
             double* WORK,
             int* LWORK,
             int* INFO);

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
 * Mathematical functions for POD operations
 */
class POD_Utils
{
public:
    /**@name Type Definitions */
    KRATOS_CLASS_POINTER_DEFINITION(POD_Utils);

    /*@{ */

    /*@} */
    /**@name Life Cycle
             */
    /*@{ */

    /** Constructor.
     */
    /// Default constructor.
    POD_Utils() {}

    /** Destructor.
     */
    /// Destructor.
    virtual ~POD_Utils() {}

    /*@} */
    /**@name Operators
         */
    /*@{ */

    /*@} */
    /**@name Operations */
    /*@{ */

    /// perform the SVD decomposition on matrix A
    /// On return: A = U * diag(S) * VT
    /// U and VT are orthogonal; S is descending
    static int SVD(const Matrix& rA, Matrix& rU, Vector& rS, Matrix& rVT)
    {
        char JOBU = 'A';
        char JOBVT = 'A';

        int M = rA.size1();
        int N = rA.size2();

        std::vector<double> A(M*N);
        for(int j = 0; j < N; ++j)
            for(int i = 0; i < M; ++i)
                A[j*M + i] = rA(i, j);

        int LDA = M;
        int LDS = std::min(M, N);
        std::vector<double> S(LDS);

        int LDU = M;
        std::vector<double> U(LDU*M);

        int LDVT = N;
        std::vector<double> VT(LDVT*N);

        // determine the size of WORK
        int LWORK = -1;
        double WORKS = 0;
        int INFO;

        dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, VT.data(), &LDVT, &WORKS, &LWORK, &INFO);

        LWORK = static_cast<int>(WORKS);
        std::cout << "Required buffer for SVD operations: " << LWORK << std::endl;

        // perform SVD
        std::vector<double> WORK(LWORK);
        dgesvd_(&JOBU, &JOBVT, &M, &N, A.data(), &LDA, S.data(), U.data(), &LDU, VT.data(), &LDVT, WORK.data(), &LWORK, &INFO);

        // export the data
        rU.resize(LDU, M, false);
        for(int j = 0; j < M; ++j)
            for(int i = 0; i < LDU; ++i)
                rU(i, j) = U[j*LDU + i];

        rS.resize(LDS, false);
        for(int i = 0; i < LDS; ++i)
            rS(i) = S[i];

        rVT.resize(LDVT, N, false);
        for(int j = 0; j < N; ++j)
            for(int i = 0; i < LDVT; ++i)
                rVT(i, j) = VT[j*LDVT + i];

        return INFO;
    }

    /// Perform triple matrix multiplication for POD-FEM simulation
    /// A = V^t K V
    /// V shall be a dense matrix coming from SVD and K is the unreduced stiffness matrix of the system
    static void VtKV(Matrix& rA, const Matrix& rV, const CompressedMatrix& rK)
    {
        const std::size_t fsize = rV.size1(); // full size
        const std::size_t rsize = rV.size2(); // reduced size

        if (rK.size1() != fsize)
            KRATOS_ERROR << "Incompatible size between stiffness matrix and the reduced matrix size";

        if (rA.size1() != rsize || rA.size2() != rsize)
            rA.resize(rsize, rsize, false);
        noalias(rA) = ZeroMatrix(rsize, rsize);

        typedef boost::numeric::bindings::traits::sparse_matrix_traits<const CompressedMatrix> matraits;
        const double* values = matraits::value_storage(rK);

        // ABC = a_ij x b_jk x c_kl
        std::size_t k, nz;
        for(std::size_t l = 0; l < rsize; ++l)
        {
            for(std::size_t j = 0; j < fsize; ++j)
            {
                nz = rK.index1_data()[j+1] - rK.index1_data()[j]; // do we need to make check here?
                for(std::size_t inz = 0; inz < nz; ++inz)
                {
                    k = rK.index2_data()[rK.index1_data()[j] + inz];
                    const double val = values[rK.index1_data()[j] + inz];
                    for(std::size_t i = 0; i < rsize; ++i)
                    {
                        rA(i, l) += rV(j, i) * val * rV(k, l);
                    }
                }
            }
        }
    }

    /**
    * Solve a (small) linear system
    * @param A the given lhs
    * @param x the solution
    * @param b the given rhs
    */
    static int Solve( const Matrix& A, Vector& x, const Vector& b )
    {
        using namespace boost::numeric::ublas;
        typedef permutation_matrix<std::size_t> pmatrix;
        Matrix Acopy(A);
        const std::size_t size = A.size1();
        pmatrix pm(size);
        const int singular = lu_factorize(Acopy, pm);
        Matrix inverse;
        inverse.assign(IdentityMatrix(size));
        lu_substitute(Acopy, pm, inverse);
        noalias(x) = prod(inverse, b);
        return singular;
    }

    /*@} */
    /**@name Acces */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */


    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Acces */
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

}  /* namespace Kratos.*/

#endif /* KRATOS_POD_UTILS  defined */
