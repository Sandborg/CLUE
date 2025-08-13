#ifndef CLUE_EX_SEARCH
#define CLUE_EX_SEARCH

#include "NoiseExperiment.hpp"
#include "dd/Package.hpp"
#include "boost/dynamic_bitset.hpp"
#include "NoisyQC.hpp"

using namespace std;

/**
 * Class for clauses in a formula.
 *
 * A clause is a a formula with the following shape: "x1 v x2 v x3...", where `vi` is a variable (possibly negated).
 * The variables can only appear once (since appearing together with different value implies the clause is always true).
 */
class NoisyQuantumSearch : public NoiseExperiment
{
protected:
    luint qbits;
    unordered_set<luint> success_set;
    vector<dd::vEdge> succes_states;
    dd::fp fidelity = 0;
    double epsilon;

    /* Method that serves as an oracle for the search function */
    bool oracle(boost::dynamic_bitset<>);
    bool oracle(luint);
    void quantum_oracle(qc::QuantumComputation &);    // Adds the oracle part to the circuit given as input
    void quantum_diffusion(qc::QuantumComputation &); // Adds the diffusion operator to the circuit given as input

    /* Overriden methods from Experiment */
    CCSparseVector clue_observable();
    dd::vEdge dd_observable();
    /* Virtual methods from Experiment */
    luint size() { return this->qbits; }
    luint correct_size() { return 2UL; }
    luint bound_size() { return 2UL; }
    array<dd::CMat, 2U> direct();
    vector<CCSparseVector> matrix();
    dd::CMat matrix_B(dd::CMat &);
    qc::QuantumComputation *quantum(double);
    qc::QuantumComputation *quantum_B(double);
    NoisyQuantumSearch *change_exec_type(ExperimentType);
    dd::fp calc_fidelity(dd::vEdge);
    dd::fp calc_fidelity(dd::vEdge &, dd::vEdge &);

public:
    NoisyQuantumSearch(luint, vector<luint>, luint, ExperimentType, dd::Package<> *, NoiseModel *, luint);

    /* Methods to create a succes state, namely the number we are looking for*/
    static NoisyQuantumSearch *random(luint, ExperimentType, dd::Package<> *, NoiseModel *, luint);
    static NoisyQuantumSearch *ones_string(luint, ExperimentType, dd::Package<> *, NoiseModel *, luint);
    void convert_succes_set_qstate(); // Convert the succes values to a quantum state

    /* Method to get the string out of an experiment */
    string to_string();
};

#endif