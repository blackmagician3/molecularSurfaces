#ifndef SOLVER_NEWTON_HPP
#define SOLVER_NEWTON_HPP

#include <stdio.h>
#include <helper_math.h>
#include <math.h>

#ifndef ZERO
#define ZERO 1e-7
#endif
////////////////////////////////////////////////////////////////////////////////
// custom matrix struct from float3 rows
////////////////////////////////////////////////////////////////////////////////
struct MatrixDim3
{
    float3 rows[3];

    // multiplication from right A*b
    inline __host__ __device__ MatrixDim3 operator*(float b) const
    {
        MatrixDim3 result;
        for (int i = 0; i < 3; ++i)
        {
            result.rows[i] = make_float3(rows[i].x * b, rows[i].y * b, rows[i].z * b);
        }
        return result;
    }
    // multiplication from left b*A
    friend inline __host__ __device__ MatrixDim3 operator*(float b, const MatrixDim3 &a)
    {
        MatrixDim3 result;
        for (int i = 0; i < 3; ++i)
        {
            result.rows[i] = make_float3(b * a.rows[i].x, b * a.rows[i].y, b * a.rows[i].z);
        }
        return result;
    }

    // multiplication in place A*=b
    inline __host__ __device__ MatrixDim3 &operator*=(float b)
    {
        for (int i = 0; i < 3; ++i)
        {
            rows[i].x *= b;
            rows[i].y *= b;
            rows[i].z *= b;
        }
        return *this;
    }
    // addition A+B
    inline __host__ __device__ MatrixDim3 operator+(const MatrixDim3 &other) const
    {
        MatrixDim3 result;
        for (int i = 0; i < 3; ++i)
        {
            result.rows[i] = make_float3(rows[i].x + other.rows[i].x, rows[i].y + other.rows[i].y, rows[i].z + other.rows[i].z);
        }
        return result;
    }
    // addition in place A+=B
    inline __host__ __device__ MatrixDim3 &operator+=(const MatrixDim3 &other)
    {
        for (int i = 0; i < 3; ++i)
        {
            rows[i].x += other.rows[i].x;
            rows[i].y += other.rows[i].y;
            rows[i].z += other.rows[i].z;
        }
        return *this;
    }

    // Matrix-vector product
    inline __host__ __device__ float3 operator*(const float3 &vec) const
    {
        float3 result;
        result.x = rows[0].x * vec.x + rows[0].y * vec.y + rows[0].z * vec.z;
        result.y = rows[1].x * vec.x + rows[1].y * vec.y + rows[1].z * vec.z;
        result.z = rows[2].x * vec.x + rows[2].y * vec.y + rows[2].z * vec.z;
        return result;
    }
};

__device__ MatrixDim3 outer_product(float3 a, float3 b)
{
    MatrixDim3 result;
    result.rows[0].x = a.x * b.x;
    result.rows[0].y = a.x * b.y;
    result.rows[0].z = a.x * b.z;
    result.rows[1].x = a.y * b.x;
    result.rows[1].y = a.y * b.y;
    result.rows[1].z = a.y * b.z;
    result.rows[2].x = a.z * b.x;
    result.rows[2].y = a.z * b.y;
    result.rows[2].z = a.z * b.z;

    return result;
}

/**
 * @brief  calculating the inverse of a matrix using the adjugate method
 *         A_inverse = (1/det(A))*adj(A), where adj(A) is the adjugate of A (the transposed cofactor matrix of A)
 *
 * @param A
 * @return __device__
 */
__device__ MatrixDim3 matrix_inverse1(MatrixDim3 A, bool getDet = false, float *det = nullptr)
{
    MatrixDim3 result;

    // Calculate the determinant
    float m_det = A.rows[0].x * (A.rows[1].y * A.rows[2].z - A.rows[1].z * A.rows[2].y) -
                  A.rows[0].y * (A.rows[1].x * A.rows[2].z - A.rows[1].z * A.rows[2].x) +
                  A.rows[0].z * (A.rows[1].x * A.rows[2].y - A.rows[1].y * A.rows[2].x);

    // Check if the determinant is near zero
    if (getDet)
    {
        *det = fabs(m_det);
    }

    if (fabs(m_det) < ZERO)
    {
        // Handle the singular matrix case by returning the identity matrix (or any error indication you prefer)
        // return A; // or return an identity matrix or any other indicator of failure
        m_det = ZERO;
    }

    // Calculate the adjugate matrix (cofactor matrix transposed)
    result.rows[0].x = A.rows[1].y * A.rows[2].z - A.rows[1].z * A.rows[2].y;
    result.rows[0].y = A.rows[0].z * A.rows[2].y - A.rows[0].y * A.rows[2].z;
    result.rows[0].z = A.rows[0].y * A.rows[1].z - A.rows[0].z * A.rows[1].y;

    result.rows[1].x = A.rows[1].z * A.rows[2].x - A.rows[1].x * A.rows[2].z;
    result.rows[1].y = A.rows[0].x * A.rows[2].z - A.rows[0].z * A.rows[2].x;
    result.rows[1].z = A.rows[0].z * A.rows[1].x - A.rows[0].x * A.rows[1].z;

    result.rows[2].x = A.rows[1].x * A.rows[2].y - A.rows[1].y * A.rows[2].x;
    result.rows[2].y = A.rows[0].y * A.rows[2].x - A.rows[0].x * A.rows[2].y;
    result.rows[2].z = A.rows[0].x * A.rows[1].y - A.rows[0].y * A.rows[1].x;

    // Divide each element of the adjugate matrix by the determinant
    result = result * (1 / m_det);

    return result;
}

__device__ MatrixDim3 matrix_inverse2(MatrixDim3 A)
{
    MatrixDim3 result;

    float m_det = A.rows[0].x * (A.rows[1].y * A.rows[2].z - A.rows[1].z * A.rows[2].y) - A.rows[0].y * (A.rows[1].x * A.rows[2].z - A.rows[1].z * A.rows[2].x) - A.rows[0].z * (A.rows[1].x * A.rows[2].y - A.rows[1].y * A.rows[2].x);
    if (abs(m_det) < ZERO)
        m_det = ZERO;
    result.rows[0] = make_float3(A.rows[1].y * A.rows[2].z - A.rows[1].z * A.rows[2].y, A.rows[0].z * A.rows[2].y - A.rows[0].y * A.rows[2].z, A.rows[0].y * A.rows[1].z - A.rows[0].z * A.rows[1].y);
    result.rows[1] = make_float3(A.rows[1].z * A.rows[2].x - A.rows[1].x * A.rows[2].z, A.rows[0].x * A.rows[2].z - A.rows[0].z * A.rows[2].x, A.rows[0].z * A.rows[1].x - A.rows[0].x * A.rows[1].z);
    result.rows[2] = make_float3(A.rows[1].x * A.rows[2].y - A.rows[1].y * A.rows[2].x, A.rows[0].y * A.rows[2].x - A.rows[0].x * A.rows[2].y, A.rows[0].x * A.rows[1].y - A.rows[0].y * A.rows[1].x);
    result *= (1 / m_det);
    return result;
}

__device__ void swap_rows(float3 *a, float3 *b)
{
    float3 temp = *a;
    *a = *b;
    *b = temp;
}

/**
 * @brief
 *
 * @param A pointer to input matrix
 * @param P
 */
__device__ void lu_factorization(float3 *A, int *P)
{
    for (int i = 0; i < 3; i++)
    {
        P[i] = i;
    }

    for (int k = 0; k < 3; k++)
    {
        // Partial pivoting
        int max = k;
        if (fabs(A[k].x) > fabs(A[max].x))
            max = k;
        if (fabs(A[k].y) > fabs(A[max].y))
            max = k;
        if (fabs(A[k].z) > fabs(A[max].z))
            max = k;

        if (max != k)
        {
            // Swap rows in permutation matrix
            int temp = P[k];
            P[k] = P[max];
            P[max] = temp;

            // Swap rows in A
            swap_rows(&A[k], &A[max]);
        }

        // LU decomposition
        for (int i = k + 1; i < 3; i++)
        {
            if (k == 0)
                A[i].x /= A[k].x;
            if (k == 1)
                A[i].y /= A[k].y;
            if (k == 2)
                A[i].z /= A[k].z;

            if (k == 0)
            {
                A[i].y -= A[i].x * A[k].y;
                A[i].z -= A[i].x * A[k].z;
            }
            else if (k == 1)
            {
                A[i].z -= A[i].y * A[k].z;
            }
        }
    }
}

__device__ void solve_lu(float3 *LU, int *P, float *b, float *x)
{
    float y[3];

    // Forward substitution
    for (int i = 0; i < 3; i++)
    {
        y[i] = b[P[i]];
        for (int j = 0; j < i; j++)
        {
            y[i] -= LU[i].x * y[j];
        }
    }

    // Backward substitution
    for (int i = 2; i >= 0; i--)
    {
        x[i] = y[i];
        for (int j = i + 1; j < 3; j++)
        {
            if (i == 0)
                x[i] -= LU[i].x * x[j];
            if (i == 1)
                x[i] -= LU[i].y * x[j];
            if (i == 2)
                x[i] -= LU[i].z * x[j];
        }
        if (i == 0)
            x[i] /= LU[i].x;
        if (i == 1)
            x[i] /= LU[i].y;
        if (i == 2)
            x[i] /= LU[i].z;
    }
}

/**
 * @brief calculate the inverse of a matrix using LU
 *
 * @param LU
 * @param P       permutation vector
 * @param inv_out pointer to output matrix
 */
__device__ void invert_matrix(float3 *LU, int *P, float3 *inv_out)
{
    float b[3];
    float x[3];

    for (int i = 0; i < 3; i++)
    {
        b[0] = b[1] = b[2] = 0.0f;
        b[i] = 1.0f;

        solve_lu(LU, P, b, x);

        inv_out[0].x = x[0];
        inv_out[1].x = x[1];
        inv_out[2].x = x[2];
    }
}

__device__ float backtrack_length(float3 point, const MatrixDim3 &atoms, const float r[3])
{
    float3 grad = make_float3(.0f);

    for (int i = 0; i < 3; i++)
    {

        float3 v = point - atoms.rows[i];
        float l = length(v);
        grad += 2 * v * (1 - (r[i] / l));
    }
    return length(grad);
}

#endif