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
/* Method that runs simulation of a circuit without any reduction. Collects fidelity between success state and outcome
(only used when this->type == DDSIM_ALONE) */
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
        // Noisy case: We create the circuit anew for each iteration
        qc::QuantumComputation *U_P = this->quantum(par_value);
        qc::QuantumComputation *U_B = this->quantum_B(par_value);
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

/*Method used to apply a circuit n times to a state*/
void NoiseExperiment::sim_n_iterations(dd::vEdge &state, luint n)
{
    double par_value = 1. / (pow(2., static_cast<double>(this->size())) * static_cast<double>(10 * this->iterations));

    for (luint i = 0; i < n; i++)
    {
        // We need to create the circuit anew for each application of the circuit
        qc::QuantumComputation *U_P = this->quantum(par_value);
        qc::QuantumComputation *U_B = this->quantum_B(par_value);
        state = dd::simulate<>(U_P, state, *package);
        state = dd::simulate<>(U_B, state, *package);
        delete U_P;
        delete U_B;
    }

    return;
}

void NoiseExperiment::run_ddsim_noise()
{
    cerr << "+++ [ddsim-noise @ " << this->name << "] Computing DDSIM NOISY  execution for " << this->name << endl;
    clock_t begin = clock();
    cerr << "+++ [ddsim-noise @ " << this->name << "] Setting up observable (" << this->observable << ") and system..." << endl;
    dd::vEdge obs = this->dd_observable();
    double par_value = 1. / (pow(2., static_cast<double>(this->size())) * static_cast<double>(10 * this->iterations));

    cerr << "+++ [ddsim-only @ " << this->name << "] Computing the iteration (U_P*U_B)^iterations..." << endl;
    clock_t b_iteration = clock();

    /*Setup for bisimulation: create a matrix with the inner product results.*/
    luint d = this->iterations; // In the beginning we can just guess what the best reduction is.
    luint M = 5000;             // The samples needed to make an accurate guess
    std::vector<std::vector<dd::fp>> inner_products(d + 1, std::vector<dd::fp>(d + 1, 0));

    for (luint k = 0; k <= d; k++)
    {
        for (luint l = k; l <= d; l++)
        {
            for (luint m = 0; m < M; m++)
            {
                cerr << "Currently working on k = " << k << ", l = " << l << ", m = " << m << "\n";
                dd::vEdge w = obs;
                dd::vEdge v = obs;

                this->sim_n_iterations(w, k);
                this->sim_n_iterations(v, l);

                inner_products[k][l] += this->calc_fidelity(v, w);
            }
            inner_products[k][l] /= M;
        }
    }

    for (const auto &r : inner_products)
    {
        for (const auto &c : r)
        {
            std::cerr << c << ", ";
        }
        std::cerr << "\n";
    }

    /*
    Main loops as presented by Max in overleaf text.
    In the text when seeing gamme_{l,k}, it it means we are looking at the k'th row and l'th column.

    For the scalar products in the text, < A^l, A^k>, where d >= l >= k, are calculated for all combinations,
    stored in the inner_products vector.

    So in the, say we want <A²,A¹>, it's stored at inner_products[1][2].
    */
    std::vector<std::vector<dd::fp>> coeff_matrix(d + 1, std::vector<dd::fp>(d + 1, 0)); // d + 1 because we want to check up to and including d.

    for (luint k = 1; k <= d; k++)
    {
        for (luint l = k; l <= d; l++)
        {
            // Calculate eta_k
            dd::fp coeff_sum = 0.0;
            for (luint i = 1; i <= k - 1; i++)
            {
                coeff_sum += std::pow(std::fabs(coeff_matrix[i][k - 1]), 2);
            }
            dd::fp eta = inner_products[k - 1][k - 1] - coeff_sum;
            coeff_matrix[k][k] = eta;
            cerr << "eta_k = " << eta << "\n";

            // Calculate gamma_{l,k}
            dd::fp eta_sqrt = std::sqrt(coeff_matrix[k][k]);

            dd::fp gamma_products = 0.0;

            for (luint i = 1; i <= k - 1; i++)
            {
                gamma_products += coeff_matrix[i][k - 1] * coeff_matrix[i][l];
            }

            coeff_matrix[k][l] = (inner_products[k - 1][l] / eta_sqrt) - (gamma_products / eta_sqrt);
        }
    }

    cerr << "\n\n";
    for (const auto &r : coeff_matrix)
    {
        for (const auto &c : r)
        {
            std::cerr << c << ", ";
        }
        std::cerr << "\n";
    }

    /*This is just placeholder atm.*/
    clock_t a_iteration = clock();
    this->fidelity = this->calc_fidelity(obs); // Should be the fidelity between the result from the reduced system compared to what we are looking for.
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
        case ExperimentType::DDSIM_NOISE:
            run_ddsim_noise();
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