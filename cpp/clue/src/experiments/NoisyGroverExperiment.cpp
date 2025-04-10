#include "experiments/NoisyGroverExperiment.hpp"

#include <cstdlib>
#include "dd/Package.hpp"

NoisyQuantumSearch::NoisyQuantumSearch(luint nQbits, vector<luint> success, luint eIterations, ExperimentType eType, dd::Package<> *ePackage) : Experiment("Grover", "H", eIterations, eType, ePackage)
{
    this->qbits = nQbits;
    luint bound = static_cast<luint>(pow(2UL, nQbits - 1));
    for (luint el : success)
    {
        if (el >= bound)
        {
            throw domain_error("The value for success is out of bound");
        }
        else
        {
            this->success_set.insert(el);
        }
    }
}
/*method to instead of having random succes values we search for, we want to use a trivial case of n nQbits - 1 ones*/
/*static*/ NoisyQuantumSearch *NoisyQuantumSearch::ones_string(luint nQbits, ExperimentType eType, dd::Package<> *ePackage)
{
    luint value = static_cast<luint>(pow(2UL, nQbits - 1));
    luint iterations = static_cast<luint>(ceil(pow(2., static_cast<double>(nQbits - 1) / 2.)));

    auto success_set = vector<luint>();
    success_set.push_back(value - 1UL);

    return new NoisyQuantumSearch(nQbits, success_set, iterations, eType, ePackage);
}
/*static*/ NoisyQuantumSearch *NoisyQuantumSearch::random(luint nQbits, ExperimentType eType, dd::Package<> *ePackage)
{
    luint half_size = static_cast<luint>(pow(2UL, nQbits - 1));
    luint iterations = static_cast<luint>(ceil(pow(2., static_cast<double>(nQbits - 1) / 2.)));

    // RANDOM WITH SEVERAL SUCCESS VALUES
    luint number_of_successes = 1U; //(static_cast<luint>(rand())%(nQbits-1))+1;
    vector<luint> success_set = vector<luint>();
    for (luint i = 0; i < number_of_successes; i++)
    {
        success_set.push_back(static_cast<luint>(rand()) % half_size);
    }
    return new NoisyQuantumSearch(nQbits, success_set, iterations, eType, ePackage);
}

bool NoisyQuantumSearch::oracle(boost::dynamic_bitset<> bitchain)
{
    luint value = 0UL, to_add = 1UL;
    for (luint i = 0; i < bitchain.size(); i++)
    {
        if (bitchain[i])
        {
            value += to_add;
        }
        to_add *= 2UL;
    }
    return this->oracle(value);
}

/**
 * Since succes_set does not allow for duplicate values, this should work as intended,
 * returning 0 when the value is not in the set and 1 when it is.
 */
bool NoisyQuantumSearch::oracle(luint value)
{
    return this->success_set.count(value);
    // return this->success_set.contains(value);
}

/**a
 * @brief Applies the Quantum Oracle associated with `this`.
 *
 * This method applies to the given ``circuit`` the Quantum Oracle associated with the
 * success set defined in `this`. For getting a description of the circuit itself, look
 * into the notes in https://cnot.io/quantum_algorithms/grover/grovers_algorithm.html,
 * where the oracle is described.
 *
 * As a summary, for each `element` in the success set, we add a controlled `X` gate
 * over the ancillary qubit, where the controls on the other qubits are as the bit chain
 * representing `element`. More precisely, for `3` in 4 bits, we have `3 = 0101`, so we would
 * apply a controlled `X` to the 5-th qubit with controls `C(-)C(+)C(-)C(+)`.
 */
void NoisyQuantumSearch::quantum_oracle(NoisyQuantumComputation &circuit, double epsilon)
{
    vector<boost::dynamic_bitset<>> success_bitchains = vector<boost::dynamic_bitset<>>(this->success_set.size());
    luint j = 0;
    for (luint success : this->success_set)
    {
        success_bitchains[j] = boost::dynamic_bitset<>(this->size() - 1, success);
        j++;
    }
    qc::Controls controls{};
    for (boost::dynamic_bitset<> element : success_bitchains)
    {
        controls.clear();
        for (luint i = 0; i < this->size() - 1; ++i)
        {
            controls.emplace(static_cast<qc::Qubit>(i), (element[i] ? qc::Control::Type::Pos : qc::Control::Type::Neg));
        }
        auto qc = qc::QuantumComputation(this->size());
        qc.mcx(controls, static_cast<qc::Qubit>(this->size() - 1));
        circuit.push_back(qc, epsilon);
    }
}
/**
 * @brief Applies the Quantum Diffusion operator for Grover's algorithm.
 *
 * This method applies to the given ``circuit`` the Quantum Diffusion operator. For getting
 * a description of the circuit itself, look into the notes in
 * https://cnot.io/quantum_algorithms/grover/grovers_algorithm.html,
 * where the diffusion operator is described.
 *
 * As a summary, we apply the hadamard gate to each qubit and then a controlled `X` gate
 * over the ancillary qubit with negative controls all over other qubits.
 *
 */
void NoisyQuantumSearch::quantum_diffusion(NoisyQuantumComputation &circuit, double epsilon)
{
    // Code taken from mqt-core/algorithms/Grover.cpp
    for (luint i = 1; i < this->size() - 1; ++i)
    {
        auto qc = qc::QuantumComputation(this->size());
        qc.h(static_cast<qc::Qubit>(i));
        circuit.push_back(qc, epsilon);
    }

    qc::Controls controls{};
    for (qc::Qubit j = 1; j < this->size() - 1; ++j)
    {
        controls.emplace(j, qc::Control::Type::Neg);
    }

    auto qc1 = qc::QuantumComputation(this->size());
    auto qc2 = qc::QuantumComputation(this->size());
    auto qc3 = qc::QuantumComputation(this->size());
    qc1.z(0); // X-H-X
    qc2.mcx(controls, 0);
    qc3.z(0); // X-H-X

    circuit.push_back(qc1, epsilon);
    circuit.push_back(qc2, epsilon);
    circuit.push_back(qc3, epsilon);

    for (luint i = this->size() - 2; i > 0; --i)
    {
        auto qc = qc::QuantumComputation(this->size());
        qc.h(static_cast<qc::Qubit>(i));
        circuit.push_back(qc, epsilon);
    }
}

/* Overriden methods from Experiment */
CCSparseVector NoisyQuantumSearch::clue_observable()
{
    luint full_size = static_cast<luint>(pow(2UL, this->size())), half = full_size / 2UL;
    CCSparseVector result = CCSparseVector(full_size);
    CC coeff = CC(sqrt(1. / static_cast<double>(half)));
    for (luint i = half; i < full_size; i++)
    {
        result.set_value(i, coeff);
    }

    return result;
}
/**
 * @brief Computes the initial state for the Grover algorithm.
 *
 * As described in https://cnot.io/quantum_algorithms/grover/grovers_algorithm.html,
 * the initial state for Grover's Algorithm is an entangled state in the non-ancillary qubits
 * with a |-> state for the ancillary qubit.
 *
 */
dd::vEdge NoisyQuantumSearch::dd_observable()
{
    vector<dd::BasisStates> states;
    for (luint i = 0; i < this->size() - 1; i++)
    {
        states.push_back(dd::BasisStates::plus);
    }
    states.push_back(dd::BasisStates::minus);

    return this->package->makeBasisState(this->size(), states);
}
/* Virtual methods from Experiment */
array<dd::CMat, 2U> NoisyQuantumSearch::direct()
{
    return {}; // TODO: This is incomplete
}
vector<CCSparseVector> NoisyQuantumSearch::matrix()
{
    luint full_size = static_cast<luint>(pow(2UL, this->size()));
    luint half_size = full_size / 2UL;
    vector<CCSparseVector> result = vector<CCSparseVector>(full_size, full_size);

    CC coeff = -CC(1 / pow(2UL, this->size() - 2)), coeff_one = coeff + CC(1.);
    for (luint i = 0; i < half_size; i++)
    {
        for (luint j = 0; j < half_size; j++)
        {
            if (i == j)
            {
                result[i].set_value(j, coeff_one);
                result[i + half_size].set_value(j + half_size, (this->success_set.count(j)) ? coeff_one : -coeff_one);
            }
            else
            {
                result[i].set_value(j, coeff);
                result[i + half_size].set_value(j + half_size, (this->success_set.count(j)) ? coeff : -coeff);
            }
        }
    }
    return result; // TODO: This is incomplete
}

dd::CMat NoisyQuantumSearch::matrix_B(dd::CMat &U)
{
    return identity_matrix(U.size()); // There is no begin hamiltonian: we use the identity
}
qc::QuantumComputation *NoisyQuantumSearch::quantum(double)
{
    double epsilon = 0.3;
    NoisyQuantumComputation circuit = NoisyQuantumComputation(this->size());

    this->quantum_oracle(circuit, epsilon);
    this->quantum_diffusion(circuit, epsilon);

    return circuit.build_noisy_qc();
}
qc::QuantumComputation *NoisyQuantumSearch::quantum_B(double)
{
    qc::QuantumComputation *circuit = new qc::QuantumComputation(this->size());
    return circuit; // Identity circuit
}
NoisyQuantumSearch *NoisyQuantumSearch::change_exec_type(ExperimentType new_type)
{
    vector<luint> to_copy;
    to_copy.reserve(this->success_set.size());
    for (std::unordered_set<luint>::iterator it = this->success_set.begin(); it != this->success_set.end(); it++)
    {
        to_copy.push_back(*it);
    }

    return new NoisyQuantumSearch(this->size() - 1, to_copy, this->iterations, new_type, this->package);
}

string NoisyQuantumSearch::to_string()
{
    stringstream stream;
    stream << "\"Grover Search Algorithm, with noisy circuit, of " << this->size() << " q-bits (1 is a flag) of [";
    std::unordered_set<luint>::iterator it = this->success_set.begin();
    if (it != this->success_set.end())
    {
        stream << "(" << *it << ", " << boost::dynamic_bitset<>(this->size() - 1UL, *it) << ")";
        it++;
    }
    while (it != this->success_set.end())
    {
        stream << ", " << "(" << *it << ", " << boost::dynamic_bitset<>(this->size() - 1UL, *it) << ")";
        ;
        it++;
    }
    stream << "]\"";
    return stream.str();
}
