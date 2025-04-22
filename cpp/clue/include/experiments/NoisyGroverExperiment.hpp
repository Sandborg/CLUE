#ifndef CLUE_EX_SEARCH
#define CLUE_EX_SEARCH

#include "Experiment.hpp"
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
class NoisyQuantumSearch : public Experiment
{
protected:
    luint qbits;
    double epsilon;
    unordered_set<luint> success_set;
    vector<dd::vEdge> succes_states;
    dd::fp fidelity = 0;

    /* Method that serves as an oracle for the search function */
    bool oracle(boost::dynamic_bitset<>);
    bool oracle(luint);
    void quantum_oracle(NoisyQuantumComputation &, double);    // Adds the oracle part to the circuit given as input
    void quantum_diffusion(NoisyQuantumComputation &, double); // Adds the diffusion operator to the circuit given as input

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
    void run_ddsim_alone() override;

public:
    NoisyQuantumSearch(luint, vector<luint>, luint, ExperimentType, dd::Package<> *, double);

    static NoisyQuantumSearch *random(luint, ExperimentType, dd::Package<> *, double);
    static NoisyQuantumSearch *ones_string(luint, ExperimentType, dd::Package<> *, double);
    void convert_succes_set_qstate(); // Convert the succes values to a quantum state

    string to_csv(char = ',') override;
    /* Method to get the string out of an experiment */
    string to_string();
};

#endif