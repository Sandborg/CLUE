#include "NoisyQC.hpp"
#include <random>

NoisyQuantumComputation::NoisyQuantumComputation(luint _nQubits)
{
    this->nQubits = _nQubits;
}

/*  When we build the noisy qc, we simply add the operation with probability 1-epsilon or else we add the identity gate to the intended target.
    With this implementation, we are also requiring the qc to have only one operation in the layer.
    This is a bit more cumbersome of an implementation, but it allows us to build it in similar fashion to the python implementation.?*/
qc::QuantumComputation *NoisyQuantumComputation::build_noisy_qc(const std::unordered_map<std::string, std::vector<double>> &P)
{
    auto qc = new qc::QuantumComputation(this->nQubits);

    // Generate random number between [0,1)
    std::random_device rd;
    std::mt19937 gen(rd());

    for (const auto &layer : this->layers)
    {
        auto op_name = layer.front()->getName();
        std::discrete_distribution<> d(P[op_name]);

        int gate_idx = 0;
        // std::cerr << layer.front()->getName() << " " << qc::toString(layer.front()->getType()) << std::endl;

        switch (gate_idx)
        {
        case 0: // the intended gate to be applied
            qc->emplace_back(layer.front()->clone());
            break;
        case 1: // The I gate is applied.
            qc->i(layer.front()->getTargets()[0]);
            break;
        case 2: // The X (not gate) is applied.
            qc->x(layer.front()->getTargets()[0]);
            break;
        case 3: // The Y gate is applied.
            qc->y(layer.front()->getTargets()[0]);
            break;
        case 4: // The Z gate is applied.
            qc->z(layer.front()->getTargets()[0]);
            break;
        default:
            break;
        }
    }

    /*
    std::cerr << "Built the following circuit based on the given epsilon values: \n"
              << *qc
              << "The circuit without noise: \n"
              << this->build_non_noisy_qc()
              << std::endl;
    */
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