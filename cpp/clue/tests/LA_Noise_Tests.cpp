#include "Linalg.hpp"

void transpose_empty_cvec_error_catch()
{
    dd::CVec V;

    try
    {
        auto V_tranposed = transpose(V);
    }
    catch (const std::logic_error &e)
    {
        return;
    }
    throw std::runtime_error("No error was caught when trying to transpose empty vector");
}

void transpose_empty_cmat_error_catch()
{
    dd::CMat M;

    try
    {
        auto M_tranposed = transpose(M);
    }
    catch (const std::logic_error &e)
    {
        return;
    }
    throw std::runtime_error("No error was caught when trying to transpose empty matrix");
}

void transpose_cvec_correct_size()
{
    dd::CVec V = dd::CVec({CC(1), CC(2), CC(3)});
    auto V_T = transpose(V);

    int V_T_rows = V_T.size();
    int V_T_cols = V_T[0].size();

    if (V_T_rows != V.size())
        throw std::runtime_error("The number of rows in the transposed vector doesn't match the expected number. Got " + std::to_string(V_T_rows) + ", expected " + std::to_string(V.size()));
    if (V_T_cols != 1)
        throw std::runtime_error("The number of columns in the transposed vector doesn't match the expected number. Got " + std::to_string(V_T_cols) + ", expected 1");
}

void transpose_cvec_correct_values()
{
    dd::CVec V = dd::CVec({CC(1), CC(2), CC(3)});
    auto V_T = transpose(V);

    for (int i = 0; i < 3; i++)
    {
        if (V[i] != V_T[i][0])
            throw std::runtime_error("The values in the transposed and original vec does not match");
    }
}

void transpose_cmat_correct_size()
{
    dd::CMat M = {{CC(1), CC(2), CC(3), CC(4)}, {CC(5), CC(6), CC(7), CC(8)}, {CC(9), CC(10), CC(11), CC(12)}};
    auto M_T = transpose(M);

    int M_T_rows = M_T.size();
    int M_T_cols = M_T[0].size();

    if (M_T_rows != M[0].size())
        throw std::runtime_error("The number of rows in the transposed matrix doesn't match the expected number. Got " + std::to_string(M_T_rows) + ", expected " + std::to_string(M[0].size()));
    if (M_T_cols != M.size())
        throw std::runtime_error("The number of columns in the transposed matrix doesn't match the expected number. Got " + std::to_string(M_T_cols) + ", expected 1");
}

void transpose_cmat_correct_values()
{
    dd::CMat M = {{CC(1), CC(2), CC(3), CC(4)}, {CC(5), CC(6), CC(7), CC(8)}, {CC(9), CC(10), CC(11), CC(12)}};
    auto M_T = transpose(M);

    int M_T_rows = M_T.size();
    int M_T_cols = M_T[0].size();

    for (int i = 0; i < M_T_rows; i++)
    {
        for (int j = 0; j < M_T_cols; j++)
        {
            if (M[j][i] != M_T[i][j])
                throw std::runtime_error("The values in the transposed matrix does not match the value in the original matrix");
        }
    }
}

void density_matrix_from_empty_cvec_error_catch()
{
    dd::CVec empty_V;

    try
    {
        auto result = get_density_matrix(empty_V);
    }
    catch (const std::logic_error &e)
    {
        return;
    }
    throw std::runtime_error("Error was not thrown when trying to get the density matrix from empty vector.");
}

void density_matrix_from_empty_cmat_error_catch()
{
    dd::CMat empty_M;

    try
    {
        auto result = get_density_matrix(empty_M);
    }
    catch (const std::logic_error &e)
    {
        return;
    }
    throw std::runtime_error("Error was not thrown when trying to get the density matrix from empty matrix.");
}

void density_matrix_from_non_cvec_matrix_error_catch()
{
    dd::CMat M = {{CC(1), CC(2), CC(3), CC(4)}, {CC(5), CC(6), CC(7), CC(8)}, {CC(9), CC(10), CC(11), CC(12)}};

    try
    {
        auto result = get_density_matrix(M);
    }
    catch (const std::logic_error &e)
    {
        return;
    }
    throw std::runtime_error("Error was not thrown when trying to get density matrix from non-vector matrix.");
}

void density_matrix_from_cvec_incorrect_size()
{
    dd::CVec V = {CC(1), CC(0), CC(1)};

    auto result = get_density_matrix(V);

    if (result.size() != V.size())
        throw std::runtime_error("The density matrix does not have the expected rows. Got " + to_string(result.size()) + "expected " + to_string(V.size()));
    else if (result[0].size() != V.size())
        throw std::runtime_error("The density matrix does not have the expected columns. Got " + to_string(result[0].size()) + "expected " + to_string(V.size()));
}

void density_matrix_from_cmat_incorrect_size()
{
    dd::CMat M = {{CC(1)}, {CC(1)}, {CC(0)}}; // this is a 3x1 "vector"

    auto result = get_density_matrix(M);

    if (result.size() != M.size())
        throw std::runtime_error("The density matrix does not have the expected rows. Got " + to_string(result.size()) + "expected " + to_string(M.size()));
    else if (result[0].size() != M.size())
        throw std::runtime_error("The density matrix does not have the expected columns. Got " + to_string(result[0].size()) + "expected " + to_string(M.size()));
}

void density_matrix_from_cvec_correct()
{
    dd::CVec V = {CC(1), CC(0), CC(1)};
    dd::CMat expected = {{CC(1), CC(0), CC(1)}, {CC(0), CC(0), CC(0)}, {CC(1), CC(0), CC(1)}};

    auto result = get_density_matrix(V);

    if (expected == result)
    {
        // This is the case we want, so do nothing.
    }
    else
        throw std::runtime_error("The density matrix does not have the expected values in it.");
}

int main()
{
    transpose_empty_cvec_error_catch();
    transpose_empty_cmat_error_catch();
    transpose_cvec_correct_size();
    transpose_cvec_correct_values();
    transpose_cmat_correct_size();
    transpose_cmat_correct_values();
    density_matrix_from_empty_cvec_error_catch();
    density_matrix_from_empty_cmat_error_catch();
    density_matrix_from_non_cvec_matrix_error_catch();
    density_matrix_from_cvec_incorrect_size();
    density_matrix_from_cmat_incorrect_size();
    density_matrix_from_cvec_correct();
    // Need to implement the above for CMat also.
    return 0;
}