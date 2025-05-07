#include "experiments/NoiseExperiment.hpp"

#include <boost/algorithm/string.hpp>
#include "dd/FunctionalityConstruction.hpp"
#include "dd/Simulation.hpp"

// IMPLEMENTATION OF GENERIC METHODS OF CLASS Experiment

// PROTECTED METHODS
/* Method to get the observable for use with CLUE */
CCSparseVector
NoiseExperiment::clue_observable()
{
    CCSparseVector result = CCSparseVector(static_cast<luint>(pow(2, this->size())));
    if (this->observable == "H")
    {
        CC c = CC(1 / sqrt(result.dimension()));
        for (luint i = 0; i < result.dimension(); i++)
        {
            result.set_value(i, c);
        }
    }
    else
    {
        result.set_value(stoul(this->observable), CC(1));
    }
    return result;
}
/* Method to get the observable for use with DD */
dd::vEdge NoiseExperiment::dd_observable()
{
    dd::vEdge obs;
    if (this->observable == "H")
    {
        obs = this->package->makeBasisState(this->size(), vector<dd::BasisStates>(this->size(), dd::BasisStates::plus), 0);
    }
    else
    {
        int val = stoi(this->observable);
        vector<bool> binary = vector<bool>(this->size(), false);
        for (luint i = this->size(); i > 0 && val > 0; i--)
        {
            binary[i - 1] = (val % 2 == 1);
            val /= 2;
        }
        obs = this->package->makeBasisState(this->size(), binary, 0);
    }
    return obs;
}

// PRIVATE METHODS
/* Method that runs the CLUE reduction (only used when this->type == DDSIM_ALONE) */
void NoiseExperiment::run_ddsim_alone()
{
    cerr << "+++ [ddsim-only @ " << this->name << "] Computing DDSIM ONLY execution for " << this->name << endl;
    clock_t begin = clock();
    cerr << "+++ [ddsim-only @ " << this->name << "] Setting up observable (" << this->observable << ") and system..." << endl;
    dd::vEdge obs = this->dd_observable();
    double par_value = 1. / (pow(2., static_cast<double>(this->size())) * static_cast<double>(10 * this->iterations));

    cerr << "+++ [ddsim-only @ " << this->name << "] Computing the iteration (U_P*U_B)^iterations..." << endl;
    clock_t b_iteration = clock();
    dd::vEdge current = obs; // We create a new vector for the current

    for (luint i = 0; i < this->iterations; i++)
    {
        qc::QuantumComputation *U_P = this->quantum(par_value);
        qc::QuantumComputation *U_B = this->quantum_B(par_value);
        // cerr << "Circuit Created:" << endl;
        // cerr << *U_P << endl;
        current = dd::simulate<>(U_P, current, *package);
        current = dd::simulate<>(U_B, current, *package);
        delete U_P;
        delete U_B;
    }

    clock_t a_iteration = clock();
    this->fidelity = this->calc_fidelity(current);
    cerr << "The fidelity between the expected state and the result from the simulation: " << this->fidelity << endl;
    clock_t end = clock();

    // We store the data
    this->red_time = 0.0;
    this->it_time = time_to_double(b_iteration, a_iteration);
    this->tot_time = time_to_double(begin, end);

    return;
}

// PUBLIC METHODS
/* Method that runs the experiment */
void NoiseExperiment::run()
{
    if (!this->executed)
    {
        cerr << "+++ Executing experiment " << this->name << " of type " << ExperimentType_toString(this->type) << endl;
        switch (this->type)
        {
        case ExperimentType::DDSIM_ALONE:
            run_ddsim_alone();
            break;
        default:
            throw logic_error("Unexpected value for type of Experiment");
        }
        this->executed = true; // We avoid repeating computations
    }
    return;
}
/* Method that generate the CSV row for this experiment */
string NoiseExperiment::to_csv(char delimiter)
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
           << this->noise_model->getDepolarization()
           << delimiter
           << this->noise_model->getAmpltitudeDamping()
           << delimiter
           << this->noise_model->getPhaseFlip()
           << delimiter
           << fidelity
           << delimiter
           << this->to_string();

    return stream.str();
}