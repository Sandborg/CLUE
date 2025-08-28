#pragma once

#include "QuantumComputation.hpp"
#include "dd/Package.hpp"

typedef long unsigned int luint;

/*************************************************************************/

/* Class for noisy quantum circuits
 */
class NoiseModel
{
protected:
    double p_depolarization = 0.0;
    double p_phaseflip = 0.0;
    double p_amplitude_damp = 0.0;

public:
    /* Constructor*/
    NoiseModel(double, double, double);

    virtual ~NoiseModel() = default;

    /* Helper functions*/
    double getDepolarization() { return p_depolarization; }
    double getPhaseFlip() { return p_phaseflip; }
    double getAmpltitudeDamping() { return p_amplitude_damp; }
    qc::QuantumComputation *build_noisy_qc(qc::QuantumComputation &); // Build a noisy quantum circuit based on the epsilon values. Have to return a pointer to match with experiment class.
};

/* Methods for bisimulation of noisy matrices*/

dd::CMat get_A_hat(const std::vector<std::vector<dd::fp>> &);
dd::CMat get_I_gscb(const dd::CMat &);