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

#ifdef ERSATZ_APP_USE_MATIO
#include <matio.h>
#endif

/* Project includes */
#include "includes/ublas_interface.h"

// template for LAPACK function
#ifdef KRATOS_USE_BLAS_LAPACK
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
#endif

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

#ifdef ERSATZ_APP_USE_MATIO

template<typename TDataType>
struct MatIoTypeHelper {};

template<>
struct MatIoTypeHelper<double>
{
    static const auto MatlabClassType = MAT_C_DOUBLE;
    static const auto MatlabDataType = MAT_T_DOUBLE;
};

template<>
struct MatIoTypeHelper<int>
{
    static const auto MatlabClassType = MAT_C_INT32;
    static const auto MatlabDataType = MAT_T_INT32;
};

template<>
struct MatIoTypeHelper<long>
{
    static const auto MatlabClassType = MAT_C_INT64;
    static const auto MatlabDataType = MAT_T_INT64;
};

#endif

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
#ifdef KRATOS_USE_BLAS_LAPACK
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
#else
        KRATOS_ERROR << "Lapack is needed for SVD";
#endif
    }

    /// Perform triple matrix multiplication for POD-FEM simulation
    /// A = W^t K V
    /// W and V shall be dense matrices coming from SVD and K is the unreduced stiffness matrix of the system
    static void WtKV(Matrix& rA, const Matrix& rW, const CompressedMatrix& rK, const Matrix& rV)
    {
        const std::size_t fsize = rW.size1(); // full size
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
        for (std::size_t l = 0; l < rsize; ++l)
        {
            for (std::size_t j = 0; j < fsize; ++j)
            {
                nz = rK.index1_data()[j + 1] - rK.index1_data()[j]; // do we need to make check here?
                for (std::size_t inz = 0; inz < nz; ++inz)
                {
                    k = rK.index2_data()[rK.index1_data()[j] + inz];
                    const double val = values[rK.index1_data()[j] + inz];
                    for (std::size_t i = 0; i < rsize; ++i)
                    {
                        rA(i, l) += rW(j, i) * val * rV(k, l);
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
    static int Solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using namespace boost::numeric::ublas;
        typedef permutation_matrix<std::size_t> pmatrix;
        Matrix Acopy(A);
        const std::size_t size = A.size1();
        pmatrix pm(size);
        const int singular = lu_factorize(Acopy, pm);
        if (!singular)
        {
            noalias(x) = b;
            lu_substitute(Acopy, pm, x);
        }
        return singular;
    }

#ifdef ERSATZ_APP_USE_MATIO
    ///  List all the variables in the mat file
    static void ListVariables(const std::string& filename)
    {
        // 1. Open MAT file in Read-Only mode
        mat_t* matfp = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
        if (!matfp) {
            KRATOS_ERROR << "Error opening file: " << filename;
        }

        // Force pointer back to first variable
        Mat_Rewind(matfp);

        std::cout << "\nListing variables in '" << filename << "':\n";
        std::cout << "--------------------------------------------------\n";
        std::cout << "Name\t\tType\t\tDimensions\n";
        std::cout << "--------------------------------------------------\n";

        matvar_t* matvar = nullptr;

        // 2. Loop through variable headers using Mat_VarReadNextInfo
        while ((matvar = Mat_VarReadNextInfo(matfp)) != nullptr) {
            std::cout << matvar->name << "\t\t";
            std::cout << get_matlab_class_string(matvar->class_type) << "\t\t";

            // Print dimensions (e.g. 10x20)
            for (int i = 0; i < matvar->rank; ++i) {
                std::cout << matvar->dims[i] << (i < matvar->rank - 1 ? "x" : "");
            }
            std::cout << "\n";

            // 3. Free header metadata memory before reading next
            Mat_VarFree(matvar);
        }

        // 4. Close file handle
        Mat_Close(matfp);
    }

    /// Read Matlab matrix stored in .mat file
    /// It is noted that the variable must be saved using v7.3 format, .e.g,
    /// > save('ew.mat','x',"-v7.3")
    template<typename TMatrixType>
    static TMatrixType ReadMat(const std::string& filename,
                               const std::string& var_name)
    {
        typedef typename TMatrixType::value_type DataType;

        // 1. Open MAT file in Read-Only mode
        mat_t* matfp = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
        if (!matfp) {
            KRATOS_ERROR << "Error opening file: " << filename;
        }

        // 2. Read variable metadata + data from file
        matvar_t* matvar = Mat_VarRead(matfp, var_name.c_str());
        if (!matvar) {
            Mat_Close(matfp);
            KRATOS_ERROR << "Variable '" << var_name << "' not found in file.";
        }

        // 3. Verify variable shape and type
        if (matvar->rank != 2) {
            Mat_VarFree(matvar);
            Mat_Close(matfp);
            KRATOS_ERROR << "Variable is not a 2D matrix.";
        }

        if (matvar->class_type != MatIoTypeHelper<DataType>::MatlabClassType) {
            Mat_VarFree(matvar);
            Mat_Close(matfp);
            KRATOS_ERROR << "Variable class type and matrix type is not compatible.";
        }

        // Dimensions in MATLAB (Rows x Cols)
        std::size_t rows = matvar->dims[0];
        std::size_t cols = matvar->dims[1];

        // 4. Access raw memory buffer
        // matio stores double arrays as a contiguous column-major double*
        const DataType* data_ptr = static_cast<const DataType*>(matvar->data);

        // 5. Convert Column-Major raw buffer -> Row-Major Boost uBLAS Matrix
        TMatrixType mat(rows, cols);
        std::size_t idx = 0;
        for (std::size_t c = 0; c < cols; ++c) {
            for (std::size_t r = 0; r < rows; ++r) {
                mat(r, c) = data_ptr[idx++];
            }
        }

        // 6. Clean up pointers
        Mat_VarFree(matvar);
        Mat_Close(matfp);

        return mat;
    }

    /// Read Matlab vector stored in .mat file
    template<typename TVectorType>
    static TVectorType ReadVec(const std::string& filename,
                               const std::string& var_name)
    {
        typedef typename TVectorType::value_type DataType;
        typedef typename MatrixVectorTypeSelector<DataType>::MatrixType MatrixType;

        MatrixType mat = ReadMat<MatrixType>(filename, var_name);
        if (mat.size1() != 1 && mat.size2() != 1) {
            std::cout << "Variable '" << var_name << "' is not a vector.";
            return Vector();
        }

        if (mat.size1() == 1) {
            // Row vector
            TVectorType vec(mat.size2());
            for (std::size_t j = 0; j < mat.size2(); ++j) {
                vec(j) = mat(0, j);
            }
            return vec;
        }
        else {
            // Column vector
            TVectorType vec(mat.size1());
            for (std::size_t i = 0; i < mat.size1(); ++i) {
                vec(i) = mat(i, 0);
            }
            return vec;
        }
    }

    /// Write a matrix to Matlab's mat file. The output file can be loaded
    /// directly in Matlab using load(filename)
    static void WriteMat(const std::string& filename,
                         const std::string& var_name,
                         const Matrix& A,
                         bool append)
    {
        const std::size_t m = A.size1();
        const std::size_t n = A.size2();

        std::vector<double> col_major_buffer;
        col_major_buffer.reserve(m * n);

        for (size_t col = 0; col < n; ++col) {
            for (size_t row = 0; row < m; ++row) {
                col_major_buffer.push_back(A(row, col));
            }
        }

        WriteMat(filename, var_name, col_major_buffer, m, n, append);
    }

    /// Write a vector to Matlab's mat file. The output file can be loaded
    /// directly in Matlab using load(filename)
    static void WriteVec(const std::string& filename,
        const std::string& var_name,
        const Vector& b,
        bool append)
    {
        const std::size_t m = b.size();

        std::vector<double> col_major_buffer;
        col_major_buffer.reserve(m);

        for (size_t row = 0; row < m; ++row) {
            col_major_buffer.push_back(b(row));
        }

        WriteMat(filename, var_name, col_major_buffer, m, 1, append);
    }

    /// Write an integer vector to Matlab's mat file. The output file can be loaded
    /// directly in Matlab using load(filename)
    static void WriteIntegerVec(const std::string& filename,
        const std::string& var_name,
        const MatrixVectorTypeSelector<int>::VectorType& b,
        bool append)
    {
        const std::size_t m = b.size();

        std::vector<int> col_major_buffer;
        col_major_buffer.reserve(m);

        for (size_t row = 0; row < m; ++row) {
            col_major_buffer.push_back(b(row));
        }

        WriteMat(filename, var_name, col_major_buffer, m, 1, append);
    }

    /// Write a matrix to Matlab's mat file. The output file can be loaded
    /// directly in Matlab using load(filename)
    template<typename TDataType>
    static void WriteMat(const std::string& filename,
                         const std::string& var_name,
                         const std::vector<TDataType>& col_major_data, // Flattened matrix
                         std::size_t rows,
                         std::size_t cols,
                         bool append = false)
    {
        // 1. Create MAT file (MAT_FT_MAT73 for v7.3 / HDF5, or MAT_FT_MAT5 for v5)
        mat_t* matfp = nullptr;
        if (append)
        {
            std::ifstream check_file(filename.c_str());
            bool file_exists = check_file.good();
            check_file.close();

            if (file_exists)
            {
                // if the file exists, open it in read-write mode to append
                matfp = Mat_Open(filename.c_str(), MAT_ACC_RDWR);
            }
            else
            {
                // if the file does not exist, create a new MAT file
                matfp = Mat_CreateVer(filename.c_str(), NULL, MAT_FT_MAT73);
            }
        }
        else
        {
            // create a new MAT file (overwrites if exists)
            matfp = Mat_CreateVer(filename.c_str(), NULL, MAT_FT_MAT73);
        }
        if (!matfp) {
            KRATOS_ERROR << "Error creating MAT file\n";
            return;
        }

        // 2. Define dimensions: [rows, cols]
        std::vector<std::size_t> dims = {rows, cols};

        // 3. Create MATLAB variable (Column-major array)
        matvar_t *matvar = Mat_VarCreate(
            var_name.c_str(),
            MatIoTypeHelper<TDataType>::MatlabClassType,
            MatIoTypeHelper<TDataType>::MatlabDataType,
            2,
            dims.data(),
            col_major_data.data(),
            0 // flags
        );

        // 4. Write variable and close file
        Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
        Mat_VarFree(matvar);
        Mat_Close(matfp);

        std::cout << "Successfully written " << var_name << " to Matlab's MAT-file"
                  << " " << filename << std::endl;
    }
#endif

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
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */

#ifdef ERSATZ_APP_USE_MATIO
    // Helper to convert matio enum types to readable strings
    static std::string get_matlab_class_string(enum matio_classes class_type) {
        switch (class_type) {
        case MAT_C_DOUBLE:   return "double";
        case MAT_C_SINGLE:   return "single";
        case MAT_C_INT8:     return "int8";
        case MAT_C_UINT8:    return "uint8";
        case MAT_C_INT16:    return "int16";
        case MAT_C_UINT16:   return "uint16";
        case MAT_C_INT32:    return "int32";
        case MAT_C_UINT32:   return "uint32";
        case MAT_C_INT64:    return "int64";
        case MAT_C_UINT64:   return "uint64";
        case MAT_C_CHAR:     return "char / string";
        case MAT_C_STRUCT:   return "struct";
        case MAT_C_CELL:     return "cell";
        default:             return "unknown / complex";
        }
    }
#endif

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
