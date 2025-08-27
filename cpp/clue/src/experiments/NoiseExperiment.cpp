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

    cerr << "+++ [ddsim-noise @ " << this->name << "] Computing the iteration (U_P*U_B)^iterations..." << endl;
    clock_t b_iteration = clock();

    /*Setup for bisimulation: create a matrix with the inner product results.
      In the beginning we can just guess what the best reduction is.
    */
    luint d = 0;
    if (this->drg == 0)
    {
        d = this->iterations; // If no guess is made, we just use the expected number of iterations to apply a circuit
    }
    else
    {
        d = this->drg; // A guess is made, try to reduce to this dimension.
    }
    luint M = 5000; // The samples needed to make an accurate guess
    std::vector<std::vector<dd::fp>> inner_products(d + 1, std::vector<dd::fp>(d + 1, 0));

    cerr << "+++ Beginning calculation of inner products with d = " << d << "...\n";
    for (luint k = 0; k <= d; k++)
    {
        for (luint l = k; l <= d; l++)
        {
            cerr << "Currently working on k = " << k << ", l = " << l << "\n";
            for (luint m = 0; m < M; m++)
            {
                dd::vEdge w = obs;
                dd::vEdge v = obs;

                this->sim_n_iterations(w, k);
                this->sim_n_iterations(v, l);

                inner_products[k][l] += this->calc_fidelity(v, w);
            }
            inner_products[k][l] /= M;
        }
    }
    cerr << "Finished calculation of inner products, moving onto reduction... \n\n";

    cerr << "Calculated inner products for <A^l,A^k>, with 0 <= k <= l <= " << d << "\n";
    for (const auto &r : inner_products)
    {
        for (const auto &c : r)
        {
            std::cerr << c << ", ";
        }
        std::cerr << "\n";
    }
    cerr << "\n";

    dd::CMat A_hat(d, dd::CVec(d, 0));

    for (int i = 1; i <= d; i++)
    {
        A_hat[0][i - 1] = inner_products[0][i];
    }
    cerr << "\n";

    /*
    Main loops as presented by Max in overleaf text.
    In the text when seeing gamme_{l,k}, it it means we are looking at the k'th row and l'th column.

    For the scalar products in the text, < A^l, A^k>, where d >= l >= k, are calculated for all combinations,
    stored in the inner_products vector.

    So in the, say we want <A²,A¹>, it's stored at inner_products[1][2].
    */

    // This part is confusingly made atm: The code starts from the second iteration.
    // this means i don't use k - 1 in <A^k-1,A^k-1> for example, since k = 1 in this loop actually is k = 2 in the pseuodo code on overleaf
    // On line 211, i use l + 1, because that vector is made so that indexes tell the number of times a circuit have been applied to a state,
    // thus, inner_products[1][2] is <A²,A¹>, for this reason we have to use l + 1.
    for (luint k = 1; k < d; k++)
    {
        for (luint l = k; l < d; l++)
        {
            // Calculate eta_k
            complex<double> coeff_sum = 0.0;
            for (luint i = 0; i < k; i++)
            {
                coeff_sum += std::pow(std::abs(A_hat[i][k - 1]), 2);
            }
            A_hat[k][k - 1] = sqrt(inner_products[k][k] - coeff_sum); // eta_k
            cerr << "eta_" << k + 1 << "= A_hat[" << k << "][" << k - 1 << "] = " << A_hat[k][k - 1] << "\n";

            // Calculate gamma_{l,k}
            complex<double> gamma_products = 0.0;

            for (luint i = 0; i < k; i++)
            {
                gamma_products += conj(A_hat[i][k - 1]) * A_hat[i][l];
            }

            A_hat[k][l] = (inner_products[k][l + 1] / A_hat[k][k - 1]) - (gamma_products / A_hat[k - 1][k]); // gamma_{l,k}
            cerr << "gamma_{" << l + 1 << "," << k + 1 << "} = A_hat[" << k << "][" << l << "] = " << A_hat[k][l] << "\n\n";
        }
    }

    /* The C_hat matrix from overleaf.
        We want to remove the last column from A_hat and the place e1 as the first column in C_hat.
    */
    dd::CMat I_gscb(d, dd::CVec(d, 0));
    I_gscb[0][0] = 1;

    for (int i = 0; i < d; i++)
    {
        for (int j = 1; j < d; j++)
        {
            I_gscb[i][j] = A_hat[i][j - 1];
        }
    }

    // Not finished, need to calc the inverse of I_gscb

    dd::CMat C_hat = matmul(A_hat, I_gscb);
    cerr << "A_hat  = " << matrix_to_string(A_hat) << endl;
    cerr << "I_gscb = " << matrix_to_string(I_gscb) << endl;
    cerr << "C_hat  = " << matrix_to_string(C_hat) << endl;

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