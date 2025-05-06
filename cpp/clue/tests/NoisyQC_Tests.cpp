#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "NoisyQC.hpp"

void failure_noise_model_build_invalid_depolarization()
{
    luint size = 3;
    double invalid_depolarization_above = 1.1;  // Invalid epsilon value (greater than 1.0)   throw std::runtime_error("Incorrect size of Noisy QC and QC didn't throw error");
    double invalid_depolarization_below = -0.1; // Invalid epsilon value (greater than 1.0)   throw std::runtime_error("Incorrect size of Noisy QC and QC didn't throw error");

    try
    {
        auto nqc_above = NoiseModel(invalid_depolarization_above, 0.0, 0.0);
    }
    catch (const std::logic_error &e)
    {
        try
        {
            auto nqc_below = NoiseModel(invalid_depolarization_below, 0.0, 0.0);
        }
        catch (const std::logic_error &e)
        {
            return;
        }
        throw std::runtime_error("Expected logic_error was not thrown for invalid invalid depolarization value below 0.");
    }
    throw std::runtime_error("Expected logic_error was not thrown for invalid invalid depolarization value above 1.");
}

void failure_noise_model_build_invalid_T1()
{
    luint size = 3;
    double invalid_T1_above = 1.1;  // Invalid epsilon value (greater than 1.0)   throw std::runtime_error("Incorrect size of Noisy QC and QC didn't throw error");
    double invalid_T1_below = -0.1; // Invalid epsilon value (greater than 1.0)   throw std::runtime_error("Incorrect size of Noisy QC and QC didn't throw error");

    try
    {
        auto nqc_above = NoiseModel(0.0, invalid_T1_above, 0.0);
    }
    catch (const std::logic_error &e)
    {
        try
        {
            auto nqc_below = NoiseModel(0.0, invalid_T1_below, 0.0);
        }
        catch (const std::logic_error &e)
        {
            return;
        }
        throw std::runtime_error("Expected logic_error was not thrown for invalid invalid T1 value below 0.");
    }
    throw std::runtime_error("Expected logic_error was not thrown for invalid invalid T1 value above 1.");
}

void failure_noise_model_build_invalid_T2()
{
    luint size = 3;
    double invalid_T2_above = 1.1;  // Invalid epsilon value (greater than 1.0)   throw std::runtime_error("Incorrect size of Noisy QC and QC didn't throw error");
    double invalid_T2_below = -0.1; // Invalid epsilon value (greater than 1.0)   throw std::runtime_error("Incorrect size of Noisy QC and QC didn't throw error");

    try
    {
        auto nqc_above = NoiseModel(0.0, 0.0, invalid_T2_above);
    }
    catch (const std::logic_error &e)
    {
        try
        {
            auto nqc_below = NoiseModel(0.0, 0.0, invalid_T2_below);
        }
        catch (const std::logic_error &e)
        {
            return;
        }
        throw std::runtime_error("Expected logic_error was not thrown for invalid invalid T2 value below 0.");
    }
    throw std::runtime_error("Expected logic_error was not thrown for invalid invalid T2 value above 1.");
}

void build_circuit()
{
    luint size = 3;

    // build circuit to make noisy
    auto qc = qc::QuantumComputation(size);
    qc.h(0);
    qc.x(1);
    qc.y(2);
    qc.h(0);
    qc.swap(0, 2);

    // build noise model
    double depolarization = 0.1;
    double T1 = 0.2;
    double T2 = 0.1;
    auto noise_model = NoiseModel(depolarization, T1, T2);

    // check current size, should be number of targets in qc ops * 4
    auto nqc = noise_model.build_noisy_qc(qc);

    luint n_targets = 0;
    for (auto &op : qc)
        n_targets += op->getNtargets();

    luint circuit_size = n_targets * 4;
    luint expected_size = 6 * 4; // we have 6 target, 1 for h,x,y, and 2 for swap

    if (expected_size != circuit_size)
    {
        delete nqc;
        throw std::runtime_error("The built circuit does not have the currect size, expected: " + std::to_string(expected_size) + ", built: " + std::to_string(circuit_size) + " with " + std::to_string(n_targets) + " targets");
    }

    delete nqc;
    return;
}

int main()
{
    failure_noise_model_build_invalid_depolarization();
    failure_noise_model_build_invalid_T1();
    failure_noise_model_build_invalid_T2();
    build_circuit();

    return 0;
}