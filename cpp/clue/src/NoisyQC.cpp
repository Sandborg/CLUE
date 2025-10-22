#include "NoisyQC.hpp"
#include <random>
#include "Types.hpp"

NoiseModel::NoiseModel(double eP1, double eP2, double eP3)
{

    if (eP1 > 1 or eP1 < 0)
        throw std::logic_error("The probability for depolarization should not be higher than 1 nor lower than 0, was given " + std::to_string(eP1));
    else
        this->p_depolarization = eP1;

    if (eP2 > 1 or eP2 < 0)
        throw std::logic_error("The probability for depolarization should not be higher than 1 nor lower than 0, was given " + std::to_string(eP2));
    else
        this->p_phaseflip = eP2;

    if (eP3 > 1 or eP3 < 0)
        throw std::logic_error("The probability for depolarization should not be higher than 1 nor lower than 0, was given " + std::to_string(eP3));
    else
        this->p_amplitude_damp = eP3;
}

/*  We apply all three error types but in sequence. Specifically, assume we have a one-qubit gate U and let us fix the probabilities p, p_1, p_2.
    Then, the noisy version of U is modelled via the following random experiment:

    Steps:
    1. Compute z1=U*z0, where z0 is the input (meaning we apply the correct gate)
    2. Set z2=z1 with probability 1-3p (identity), set z2=X*z1 with probability p, set z2=Y*z1 with probability p, set z2=Z*z1 with probability p.
    3. Set z3=E_0*z2 with probability 1-p1, set z3=E_1*z2 with probability p1 (T1 decoherence, Eq. 6 in Wille's paper).
    4. Set z4=z3 with probability 1-p2, set z4=Z*z3 with probability p2 (T2 decoherence, Eq. 7 in Wille's paper --- the one with the typo).

    Step 1 where the gate is applied, step 2 is depolarization, step 3 is T1 decoherence (amplitude damping) and step 4 is T2 decoherence (phase flip)
*/
qc::QuantumComputation *NoiseModel::build_noisy_qc(qc::QuantumComputation &qc)
{
    auto noisy_qc = new qc::QuantumComputation(qc.getNqubits());

    double I = 1 - ((3 * this->p_depolarization) / 4);
    double X = this->p_depolarization / 4;
    double Y = this->p_depolarization / 4;
    double Z = this->p_depolarization / 4;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> depolarization_dist({I, X, Y, Z});
    std::discrete_distribution<> T1_dist({1 - this->p_amplitude_damp, this->p_amplitude_damp});
    std::discrete_distribution<> T2_dist({1 - this->p_phaseflip, this->p_phaseflip});

    for (const auto &op : qc)
    {

        // Step 1: Compute z1=U*z0, where z0 is the input (meaning we apply the correct gate)
        noisy_qc->emplace_back(op->clone());

        // step 2: Depolarization
        for (const auto &target : op->getTargets())
        {
            switch (depolarization_dist(gen))
            {
            case 0:
                noisy_qc->i(target);
                break;
            case 1:
                noisy_qc->x(target);
                break;
            case 2:
                noisy_qc->y(target);
                break;
            case 3:
                noisy_qc->z(target);
                break;
            default:
                break;
            }

            // step 3: T1 decoherene (amplitude damping) **NOT THE CORRECT GATES**
            switch (T1_dist(gen))
            {
            case 0:
                noisy_qc->i(target);
                break;
            case 1:
                noisy_qc->z(target);
                break;
            default:
                break;
            }

            // step 4: T2 decoherence (phase flip)
            switch (T2_dist(gen))
            {
            case 0:
                noisy_qc->i(target);
                break;
            case 1:
                noisy_qc->z(target);
                break;
            default:
                break;
            }
        }
    }

    return noisy_qc;
}

dd::CMat collect_eta_gamma(const std::vector<std::vector<dd::fp>> &inner_products)
{
    luint d = inner_products.size() - 1; // size - 1 because inner product matrix is 1 bigger than d.
    dd::CMat A_hat(d, dd::CVec(d, clue::CC(0)));

    for (luint i = 0; i < d; i++)
    {
        A_hat[0][i] = inner_products[0][i + 1];
    }

    for (int k = 1; k < d; k++)
    {
        std::cerr << "k = " << k + 1 << std::endl;
        // Calculate eta_k
        std::complex<double> coeff_sum = 0.0;
        for (int i = 0; i <= k - 1; i++)
        {
            auto temp = std::pow(std::abs(A_hat[i][k - 1]), 2);
            std::cerr << "Ahat[" << i << "][" << k - 1 << "] = " << A_hat[i][k - 1] << " result: " << temp << std::endl;
            // std::cerr << "coeff_sum k = " << k << " and i = " << i << std::endl;
            coeff_sum += temp;
        }
        std::cerr << "coeff_sum = " << coeff_sum << std::endl;

        std::cerr << "<A^" << k << ", A^" << k << "> = " << inner_products[k][k] << std::endl;
        auto eta_squared = inner_products[k][k] - coeff_sum;

        if (eta_squared.real() < 0)
            eta_squared *= -1;

        A_hat[k][k - 1] = sqrt(eta_squared); // eta_k
        // std::cerr << "eta_" << k + 1 << " = A_hat[" << k << "][" << k - 1 << "] = " << A_hat[k][k - 1] << "\n\n";
        std::cerr << "eta_" << k + 1 << " = " << inner_products[k][k] << " - " << coeff_sum << " = " << eta_squared << " = A_hat[" << k << "][" << k - 1 << "] = " << A_hat[k][k - 1] << std::endl;
        std::cerr << "Squareroot of eta_" << k + 1 << " = " << A_hat[k][k - 1] << "\n\n";

        for (int l = k; l < d; l++)
        {
            // Calculate gamma_{l,k}
            std::complex<double> gamma_products = 0.0;

            for (int i = 0; i <= k - 1; i++)
            {
                auto temp = std::conj(A_hat[i][k - 1]) * A_hat[i][l];
                gamma_products += temp;
                std::cerr << "Ahat[" << i << "][" << k - 1 << "] = " << A_hat[i][k - 1] << ", Ahat[" << i << "][" << l << "] = " << A_hat[i][l] << ", multipled = " << temp << std::endl;
            }

            A_hat[k][l] = (inner_products[k][l + 1] / A_hat[k][k - 1]) - (gamma_products / A_hat[k][k - 1]); // gamma_{l,k}
            if (A_hat[k][l].real() < 0)
                A_hat[k][l] *= -1;

            std::cerr << "gamma_sums = " << gamma_products << std::endl;
            std::cerr << "<A^" << l + 1 << ", A^" << k << "> = " << inner_products[k][l + 1] << std::endl;
            std::cerr << "gamma_{" << l + 1 << "," << k + 1 << "} = A_hat[" << k << "][" << l << "] = " << (inner_products[k][l + 1] / A_hat[k][k - 1]) << " - " << (gamma_products / A_hat[k][k - 1]) << " = " << A_hat[k][l] << "\n\n";
        }
    }

    return A_hat;
}
/*
Main loops as presented by Max in overleaf text.
In the text when seeing gamme_{l,k}, it it means we are looking at the k'th row and l'th column.

For the scalar products in the text, < A^l, A^k>, where d >= l >= k, are calculated for all combinations,
stored in the inner_products vector.

So in the end, say we want <A²,A¹>, it's stored at inner_products[1][2].

This part is confusingly made atm: The code starts from the second iteration.
this means i don't use k - 1 in <A^k-1,A^k-1> for example, since k = 1 in this loop actually is k = 2 in the pseuodo code on overleaf
On line 211, i use l + 1, because that vector is made so that indexes tell the number of times a circuit have been applied to a state,
thus, inner_products[1][2] is <A²,A¹>, for this reason we have to use l + 1.
*/
dd::CMat get_A_hat(const std::vector<std::vector<dd::fp>> &inner_products)
{

    luint d = inner_products.size() - 1; // size - 1 because inner product matrix is 1 bigger than d.
    dd::CMat A_hat(d, dd::CVec(d, clue::CC(0)));

    for (int i = 1; i <= d; i++)
    {
        A_hat[0][i - 1] = inner_products[0][i];
    }

    for (luint k = 1; k < d; k++)
    {
        // Calculate eta_k
        std::complex<double> coeff_sum = 0.0;
        for (luint i = 0; i < k; i++)
        {
            coeff_sum += std::pow(std::abs(A_hat[i][k - 1]), 2);
        }
        A_hat[k][k - 1] = sqrt(inner_products[k][k] - coeff_sum); // eta_k
        // std::cerr << "eta_" << k + 1 << "= A_hat[" << k << "][" << k - 1 << "] = " << A_hat[k][k - 1] << "\n";
        for (luint l = k; l < d; l++)
        {

            // Calculate gamma_{l,k}
            std::complex<double> gamma_products = 0.0;

            for (luint i = 0; i < k; i++)
            {
                gamma_products += std::conj(A_hat[i][k - 1]) * A_hat[i][l];
            }

            A_hat[k][l] = (inner_products[k][l + 1] / A_hat[k][k - 1]) - (gamma_products / A_hat[k - 1][k]); // gamma_{l,k}
            // std::cerr << "gamma_{" << l + 1 << "," << k + 1 << "} = A_hat[" << k << "][" << l << "] = " << A_hat[k][l] << "\n\n";
        }
    }

    return A_hat;
}

/*
We want to remove the last column from A_hat and the place e1 as the first column in I_gscb.
*/
dd::CMat get_I_gscb(const dd::CMat &A_hat)
{

    luint d = A_hat.size();
    dd::CMat I_gscb(d, dd::CVec(d, clue::CC(0)));
    I_gscb[0][0] = 1;

    for (int i = 0; i < d; i++)
    {
        for (int j = 1; j < d; j++)
        {
            I_gscb[i][j] = A_hat[i][j - 1];
        }
    }

    return I_gscb;
}