#pragma once

#include "QuantumComputation.hpp"

typedef long unsigned int luint;

/*************************************************************************/

/* Class for noisy quantum circuits
 */
class NoiseModel
{
protected:
    luint nQubits;
    double p_depolarization = 0.0;
    double p_phaseflip = 0.0;
    double p_amplitude_damp = 0.0;

public:
    /* Constructor*/
    NoiseModel(luint, double, double, double);

    virtual ~NoiseModel() = default;

    /* Helper functions*/
    double getDepolarization() { return p_depolarization; }
    double getPhaseFlip() { return p_phaseflip; }
    double getAmpltitudeDamping() { return p_amplitude_damp; }
    luint getNqubits() { return this->nQubits; }                      // Number of qubits in the circuit
    qc::QuantumComputation *build_noisy_qc(qc::QuantumComputation &); // Build a noisy quantum circuit based on the epsilon values. Have to return a pointer to match with experiment class.
};