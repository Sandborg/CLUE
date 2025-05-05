#pragma once

#include "QuantumComputation.hpp"

typedef long unsigned int luint;

/*************************************************************************/

/* Class for noisy quantum circuits
 */
class NoisyQuantumComputation
{
protected:
    luint nQubits;
    double p_depolarization = 0.0;
    double p_phaseflip = 0.0;
    double p_amplitude_damp = 0.0;
    std::vector<qc::QuantumComputation> layers;

public:
    /* Constructor*/
    NoisyQuantumComputation(luint, double, double, double);

    /* Helper functions*/
    int size() { return this->layers.size(); }
    luint getNqubits() { return this->nQubits; }  // Number of qubits in the circuit
    qc::QuantumComputation *build_noisy_qc();     // Build a noisy quantum circuit based on the epsilon values. Have to return a pointer to match with experiment class.
    qc::QuantumComputation *build_non_noisy_qc(); // Get the quantum circuit as if no noise is present.

    /* Add a layer to the Noisy Quantum Computation
       Primary use of the noisy quantum computation. With this implementation we can each gate a specific epsilon value.

       This is a bit more cumbersome of an implementation, but it allows us to build it in similar fashion to the python implementation.
    */
    void push_back(const qc::QuantumComputation qc)
    {
        if (qc.size() == 0)
            throw std::runtime_error("The quantum circuit is empty.");
        if (qc.getNqubits() != this->getNqubits())
            throw std::runtime_error("The qubits of the circuits does not match the expected number.");

        layers.push_back(qc);
    }
};