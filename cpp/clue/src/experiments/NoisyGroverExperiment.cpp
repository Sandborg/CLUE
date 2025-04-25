#include "experiments/NoisyGroverExperiment.hpp"

#include <cstdlib>
#include "dd/Simulation.hpp"

NoisyQuantumSearch::NoisyQuantumSearch(luint nQbits, vector<luint> success, luint eIterations, ExperimentType eType, dd::Package<> *ePackage, unordered_map<string, vector<double>> P) : Experiment("Grover", "H", eIterations, eType, ePackage)
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

    // TODO: Check if all distributions in P equal 1.
}

/*method to instead of having random succes values we search for, we want to use a trivial case of n nQbits - 1 ones*/
/*static*/ NoisyQuantumSearch *NoisyQuantumSearch::ones_string(luint nQbits, ExperimentType eType, dd::Package<> *ePackage, unordered_map<string, vector<double>> P)
{
    luint value = static_cast<luint>(pow(2UL, nQbits - 1));
    luint iterations = static_cast<luint>(ceil(pow(2., static_cast<double>(nQbits - 1) / 2.))) - 1;

    auto success_set = vector<luint>();
    success_set.push_back(value - 1UL);

    return new NoisyQuantumSearch(nQbits, success_set, iterations, eType, ePackage, P);
}
/*static*/ NoisyQuantumSearch *NoisyQuantumSearch::random(luint nQbits, ExperimentType eType, dd::Package<> *ePackage, unordered_map<string, vector<double>> P)
{
    luint half_size = static_cast<luint>(pow(2UL, nQbits - 1));
    luint iterations = static_cast<luint>(ceil(pow(2., static_cast<double>(nQbits - 1) / 2.))) - 1;

    // RANDOM WITH SEVERAL SUCCESS VALUES
    luint number_of_successes = 1U; //(static_cast<luint>(rand())%(nQbits-1))+1;
    vector<luint> success_set = vector<luint>();
    for (luint i = 0; i < number_of_successes; i++)
    {
        success_set.push_back(static_cast<luint>(rand()) % half_size);
    }
    return new NoisyQuantumSearch(nQbits, success_set, iterations, eType, ePackage, P);
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
void NoisyQuantumSearch::quantum_oracle(NoisyQuantumComputation &circuit)
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
        circuit.push_back(qc);
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
void NoisyQuantumSearch::quantum_diffusion(NoisyQuantumComputation &circuit)
{
    // Code taken from mqt-core/algorithms/Grover.cpp
    for (luint i = 1; i < this->size() - 1; ++i)
    {
        auto qc = qc::QuantumComputation(this->size());
        qc.h(static_cast<qc::Qubit>(i));
        circuit.push_back(qc);
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

    circuit.push_back(qc1);
    circuit.push_back(qc2);
    circuit.push_back(qc3);

    for (luint i = this->size() - 2; i > 0; --i)
    {
        auto qc = qc::QuantumComputation(this->size());
        qc.h(static_cast<qc::Qubit>(i));
        circuit.push_back(qc);
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
    NoisyQuantumComputation circuit = NoisyQuantumComputation(this->size());

    this->quantum_oracle(circuit);
    this->quantum_diffusion(circuit);

    return circuit.build_noisy_qc(this->P);
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

    return new NoisyQuantumSearch(this->size() - 1, to_copy, this->iterations, new_type, this->package, this->P);
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

void NoisyQuantumSearch::convert_succes_set_qstate()
{
    for (const auto &succes : success_set)
    {
        auto bitchain = boost::dynamic_bitset<>(this->size() - 1UL, succes);
        vector<dd::BasisStates> states;

        for (int i = 0; i < bitchain.size(); i++)
        {
            if (bitchain[i])
                states.push_back(dd::BasisStates::one);
            else
                states.push_back(dd::BasisStates::zero);
        }

        states.push_back(dd::BasisStates::minus);
        succes_states.push_back(this->package->makeBasisState(this->size(), states));
    }
}

dd::fp fid_test(dd::Package<> *package, dd::vEdge dd)
{

    auto fidelity = package->fidelity(dd, dd);

    return fidelity;
}

/* Method that simulates a quantum circuit without reduction (only used when this->type == DDSIM_ALONE) */
void NoisyQuantumSearch::run_ddsim_alone()
{
    cerr << "+++ [ddsim-only @ " << this->name << "] Computing DDSIM ONLY execution for " << this->name << endl;
    clock_t begin = clock();
    cerr << "+++ [ddsim-only @ " << this->name << "] Setting up observable (" << this->observable << ") and system..." << endl;
    dd::vEdge obs = this->dd_observable();
    this->convert_succes_set_qstate();
    double par_value = 1. / (pow(2., static_cast<double>(this->size())) * static_cast<double>(10 * this->iterations));
    qc::QuantumComputation *U_P = this->quantum(par_value);
    qc::QuantumComputation *U_B = this->quantum_B(par_value);

    cerr << "Circuit Created:" << endl;
    cerr << *U_P << endl;

    cerr << "+++ [ddsim-only @ " << this->name << "] Computing the iteration (U_P*U_B)^iterations..." << endl;
    clock_t b_iteration = clock();
    dd::vEdge current = obs; // We create a new vector for the current

    for (luint i = 0; i < this->iterations; i++)
    {
        current = dd::simulate<>(U_P, current, *package);
        current = dd::simulate<>(U_B, current, *package);
    }

    /*cerr << "state vector after simulation" << endl;
    current.printVector();
    cerr << "the state vector that we are looking for" << endl;
    this->succes_states[0].printVector();

    // auto fid = fid_test(this->package, current);
    // cerr << "The fidelity between the state after simulation and itself: " << fid_test << endl;
    */
    this->fidelity = this->package->fidelity(current, this->succes_states[0]);
    cerr << "The fidelity between the expected state and the result from the simulation: " << this->fidelity << endl;
    clock_t a_iteration = clock();
    clock_t end = clock();

    // We store the data
    this->red_time = 0.0;
    this->it_time = time_to_double(b_iteration, a_iteration);
    this->tot_time = time_to_double(begin, end);

    delete U_P;
    delete U_B;

    return;
}

/*
    The implementation uses discrete_distribution which techincally does not need to be given values that sum to 1,
    because it will be normalised so that it does.
    Checking for it just in case.

    Can always be removed without causing errors.

    The probabilities in the distribution will be understood as follows:
        {intended gate, I, X, Y, Z}, the first element is the probability for the intended gate to be applied,
        the second is the probability for the I gate to be applied, third for X gate and so on.
*/
void NoisyQuantumSearch::add_distribution(string type, vector<double> probabilities)
{
    auto sum = 0.0;

    for (const auto &prob : probabilities)
        sum += prob;

    if (sum != 1.0)
        throw std::runtime_error("The probabilities given does not sum to 1.");

    this->P[type] = probabilities;
}

string NoisyQuantumSearch::to_csv(char delimiter)
{
    stringstream stream;
    stream << this->size()
           << delimiter
           << this->bound_size()
           << delimiter
           << this->name
           << delimiter
           << this->observable
           << delimiter
           << this->red_time
           << delimiter
           << this->red_ratio
           << delimiter
           << this->iterations
           << delimiter
           << this->it_time
           << delimiter
           << this->tot_time
           << delimiter
           << 0 // It's the epsilon value in the csv file, change this to something appropriate
           << delimiter
           << fidelity
           << delimiter
           << this->to_string();

    return stream.str();
}