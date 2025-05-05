#include "NoisyQC.hpp"
#include <random>

NoisyQuantumComputation::NoisyQuantumComputation(luint _nQubits, double eP1, double eP2, double eP3)
{
    this->nQubits = _nQubits;

    if (eP1 > 1 or eP1 < 0)
        throw std::logic_error("The probability for depolarization should not be higher than 1 nor lower than 0, was given" + std::to_string(eP1));
    else
        this->p_depolarization = eP1;

    if (eP2 > 1 or eP2 < 0)
        throw std::logic_error("The probability for depolarization should not be higher than 1 nor lower than 0, was given" + std::to_string(eP2));
    else
        this->p_phaseflip = eP2;

    if (eP3 > 1 or eP3 < 0)
        throw std::logic_error("The probability for depolarization should not be higher than 1 nor lower than 0, was given" + std::to_string(eP3));
    else
        this->p_amplitude_damp = eP3;
}

/*  We apply all three error typos but in sequence. Specifically, assume we have a one-qubit gate U and let us fix the probabilities p, p_1, p_2.
    Then, the noisy version of U is modelled via the following random experiment:

    Steps:
    1. Compute z1=U*z0, where z0 is the input (meaning we apply the correct gate)
    2. Set z2=z1 with probability 1-3p (identity), set z2=X*z1 with probability p, set z2=Y*z1 with probability p, set z2=Z*z1 with probability p.
    3. Set z3=E_0*z2 with probability 1-p1, set z3=E_1*z2 with probability p1 (T1 decoherence, Eq. 6 in Wille's paper).
    4. Set z4=z3 with probability 1-p2, set z4=Z*z3 with probability p2 (T2 decoherence, Eq. 7 in Wille's paper --- the one with the typo).

    Step 1 where the gate is applied, step 2 is depolarization, step 3 is T1 decoherence (amplitude damping) and step 4 is T2 decoherence (phase flip)
*/
qc::QuantumComputation *NoisyQuantumComputation::build_noisy_qc()
{
    auto qc = new qc::QuantumComputation(this->nQubits);

    double I = 1 - ((3 * this->p_depolarization) / 4);
    double X = this->p_depolarization / 4;
    double Y = this->p_depolarization / 4;
    double Z = this->p_depolarization / 4;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> depolarization_dist({I, X, Y, Z});
    std::discrete_distribution<> T1_dist({1 - this->p_amplitude_damp, this->p_amplitude_damp});
    std::discrete_distribution<> T2_dist({1 - this->p_phaseflip, this->p_phaseflip});

    for (const auto &layer : this->layers)
    {
        // Step 1: Compute z1=U*z0, where z0 is the input (meaning we apply the correct gate)
        qc->emplace_back(layer.front()->clone());

        for (auto target : layer.front()->getTargets())
        {
            // int gate_idx = depolarization_dist(gen);

            // step 2: Depolarization
            switch (depolarization_dist(gen))
            {
            case 0: // the intended gate to be applied
                qc->i(target);
                break;
            case 1: // The X gate is applied.
                qc->x(target);
                break;
            case 2: // The Y (not gate) is applied.
                qc->y(target);
                break;
            case 3: // The Z gate is applied.
                qc->z(target);
                break;
            default:
                break;
            }

            // step 3: T1 decoherene (amplitude damping) **NOT THE CORRECT GATES**
            switch (T1_dist(gen))
            {
            case 0:
                qc->i(target);
                break;
            case 1:
                qc->z(target);
                break;
            default:
                break;
            }

            // step 4: T2 decoherence (phase flip)
            switch (T2_dist(gen))
            {
            case 0:
                qc->i(target);
                break;
            case 1:
                qc->z(target);
                break;
            default:
                break;
            }
        }
    }

    return qc;
}

qc::QuantumComputation *NoisyQuantumComputation::build_non_noisy_qc()
{
    auto qc = new qc::QuantumComputation(this->nQubits);

    for (const auto &layer : this->layers)
    {
        qc->emplace_back(layer.front()->clone());
    }

    return qc;
}