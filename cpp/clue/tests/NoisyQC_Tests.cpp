#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "NoisyQC.hpp"

// Test if a layer is added to the noisy_qc currectly with the push_back function
void add_layer_succes(int n, luint qbits)
{
    double p1 = 0.0;
    double p2 = 0.0;
    double p3 = 0.0;
    auto nqc = NoiseModel(qbits, p1, p2, p3);

    // Add layers to the Noisy QuantumComputation
    for (int i = 0; i < n; i++)
    {
        auto layer = qc::QuantumComputation(qbits);
        layer.h(0);
        nqc.push_back(layer);
    }

    // Throw error if we do not have the expected number of layers
    if (nqc.size() != n)
    {
        throw std::runtime_error("The noisy quantum computation layers did not get added currectly. Expected " + std::to_string(n) + " got " + std::to_string(nqc.size()));
    }
}

void add_layer_failure_incorrect_size()
{
    luint nqc_size = 3;
    luint qc_size = 4;
    double p1 = 0.0;
    double p2 = 0.0;
    double p3 = 0.0;

    auto nqc = NoiseModel(nqc_size, p1, p2, p3);
    auto invalid_qc = qc::QuantumComputation(qc_size);

    try
    {
        nqc.push_back(invalid_qc);
    }
    catch (std::runtime_error)
    {
        return;
    }
    throw std::runtime_error("Expected runtime_error was not thrown for mismatched sizes.");
}

void invalid_depolarization_probabilities()
{
    double invalid_probability_pos = 1.1;
    double invalid_probability_neg = -0.5;

    try
    {
        auto nqc1 = NoiseModel(3, invalid_probability_pos, 0.0, 0.0);
        auto nqc2 = NoiseModel(3, invalid_probability_neg, 0.0, 0.0);
    }
    catch (const std::logic_error)
    {
        return;
    }
    throw std::runtime_error("A wrong depolarization probability was not caught");
}

void invalid_amplitude_damping_probabilities()
{
    double invalid_probability_pos = 1.1;
    double invalid_probability_neg = -0.5;

    try
    {
        auto nqc1 = NoiseModel(3, 0.0, invalid_probability_pos, 0.0);
        auto nqc2 = NoiseModel(3, 0.0, invalid_probability_pos, 0.0);
    }
    catch (const std::logic_error)
    {
        return;
    }
    throw std::runtime_error("A wrong amplitude damping probability was not caught");
}

void invalid_phaseflip_probabilities()
{
    double invalid_probability_pos = 1.1;
    double invalid_probability_neg = -0.5;

    try
    {
        auto nqc1 = NoiseModel(3, 0.0, 0.0, invalid_probability_pos);
        auto nqc2 = NoiseModel(3, 0.0, 0.0, invalid_probability_neg);
    }
    catch (const std::logic_error)
    {
        return;
    }
    throw std::runtime_error("A wrong phase flip probability was not caught");
}

void add_layer_failure_invalid_epsilon()
{
    luint nqc_size = 3;
    double invalid_epsilon = 1.5; // Invalid epsilon value (greater than 1.0)   throw std::runtime_error("Incorrect size of Noisy QC and QC didn't throw error");

    auto nqc = NoiseModel(nqc_size, 0.0, 0.0, 0.0);
    auto qc = qc::QuantumComputation(nqc_size);

    try
    {
        nqc.push_back(qc);
    }
    catch (std::runtime_error)
    {
        return;
    }
    throw std::runtime_error("Expected runtime_error was not thrown for invalid epsilon.");
}

void add_layer_failure_empty_qc()
{
    luint nqc_size = 3;
    double epsilon = 0.01;

    auto nqc = NoiseModel(nqc_size, 0.0, 0.0, 0.0);
    auto empty_qc = qc::QuantumComputation(nqc_size);

    try
    {
        nqc.push_back(empty_qc);
    }
    catch (std::runtime_error)
    {
        return;
    }
    throw std::runtime_error("Expected runtime_error was not thrown for empty quantum circuit.");
}

void build_non_noisy_circuit(luint qubits)
{
    double epsilon = 0.01;

    auto nqc = NoiseModel(qubits, 0.0, 0.0, 0.0);

    // Add layers to the Noisy QuantumComputation
    for (int i = 0; i < qubits; i++)
    {
        auto qc = qc::QuantumComputation(qubits);
        qc.h(i);
        nqc.push_back(qc);
    }

    // Build the non-noisy circuit
    auto non_noisy_qc = nqc.build_non_noisy_qc();

    // Check if the size of the non-noisy circuit is equal to the number of layers
    if (non_noisy_qc->size() != nqc.size())
    {
        throw std::runtime_error("The size of the non-noisy circuit does not match the number of layers.");
    }
}

/*
    Need to figure out how to test this.
    void build_noisy_circuit(luint qubits)
    {
        double epsilon = 0.01;

        auto nqc = NoisyQuantumComputation(qubits);

        // Add layers to the Noisy QuantumComputation
        for (int i = 0; i < qubits; i++)
        {
            auto qc = qc::QuantumComputation(qubits);
            qc.h(i);
            nqc.push_back(qc);
        }

        std::map<std::string, std::vector<double>> P;

        P["x"] = {0.88, 0.05, 0.02, 0.03, 0.02};
        P["h"] = {0.92, 0.01, 0.02, 0.03, 0.02};
        P["z"] = {0.96, 0.01, 0.01, 0.01, 0.01};
        // Build the noisy circuit
        auto noisy_qc = nqc.build_noisy_qc(P);

        // Check if the size of the noisy circuit is equal to the number of layers
        if (noisy_qc->size() != nqc.size())
        {
            throw std::runtime_error("The size of the noisy circuit does not match the number of layers.");
        }
    }
    */

int main()
{
    add_layer_succes(3, 3);
    add_layer_failure_incorrect_size();
    add_layer_failure_invalid_epsilon();
    add_layer_failure_empty_qc();
    build_non_noisy_circuit(3);
    invalid_depolarization_probabilities();
    invalid_amplitude_damping_probabilities();
    invalid_phaseflip_probabilities();
    // build_noisy_circuit(3);

    return 0;
}