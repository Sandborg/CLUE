#include "experiments/NoiseExperiment.hpp"

#include <boost/algorithm/string.hpp>
#include "dd/FunctionalityConstruction.hpp"
#include "dd/Simulation.hpp"

// IMPLEMENTATION OF GENERIC METHODS OF CLASS Experiment

// PROTECTED METHODS
/* Method to get the observable for use with CLUE */
CCSparseVector NoiseExperiment::clue_observable()
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
void NoiseExperiment::sim_n_iterations(dd::vEdge &state, luint iterations)
{
    double par_value = 1. / (pow(2., static_cast<double>(this->size())) * static_cast<double>(10 * this->iterations));

    for (luint i = 0; i < iterations; i++)
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

std::vector<std::vector<dd::fp>> NoiseExperiment::collect_inner_products(const dd::vEdge &obs, luint dim, luint samples)
{
    std::vector<std::vector<dd::fp>> inner_products(dim + 1, std::vector<dd::fp>(dim + 1, 0));

    double par_value = 1. / (pow(2., static_cast<double>(this->size())) * static_cast<double>(10 * this->iterations));

    for (luint m = 0; m < samples; m++)
    {
        if (m % 250 == 0 || m == 0)
            cerr << "Currently working on m = " << m + 1 << "\n";
        auto curr_w = obs;
        for (luint k = 0; k <= dim; k++)
        {
            qc::QuantumComputation *U_P_w = this->quantum(par_value);
            qc::QuantumComputation *U_B_w = this->quantum_B(par_value);
            if (k != 0)
            {
                curr_w = dd::simulate<>(U_P_w, curr_w, *package);
                curr_w = dd::simulate<>(U_B_w, curr_w, *package);
            }
            delete U_P_w;
            delete U_B_w;

            auto curr_v = obs;
            for (luint l = k; l <= dim; l++)
            {
                qc::QuantumComputation *U_P_v = this->quantum(par_value);
                qc::QuantumComputation *U_B_v = this->quantum_B(par_value);
                if (l != 0)
                {
                    curr_v = dd::simulate<>(U_P_v, curr_v, *package);
                    curr_v = dd::simulate<>(U_B_v, curr_v, *package);
                }
                delete U_P_v;
                delete U_B_v;

                inner_products[k][l] += this->calc_fidelity(curr_v, curr_w);
            }
        }
    }
    for (auto &row : inner_products)
        for (auto &x : row)
            x /= M;

    return inner_products;
}

std::vector<dd::fp> NoiseExperiment::collect_expected_inner_products(const dd::vEdge &obs, luint dim, luint samples)
{
    std::vector<dd::fp> inner_products(dim, 0);
    double par_value = 1. / (pow(2., static_cast<double>(this->size())) * static_cast<double>(10 * this->iterations));

    for (luint m = 0; m < samples; m++)
    {
        if (m % 250 == 0 || m == 0)
            cerr << "Currently working on m = " << m + 1 << "\n";
        auto curr = obs;
        for (luint k = 0; k < dim; k++)
        {
            qc::QuantumComputation *U_P = this->quantum(par_value);
            qc::QuantumComputation *U_B = this->quantum_B(par_value);
            if (k != 0)
            {
                curr = dd::simulate<>(U_P, curr, *package);
                curr = dd::simulate<>(U_B, curr, *package);
            }

            delete U_P;
            delete U_B;
            inner_products[k] += this->calc_fidelity(curr);
        }
    }

    for (auto &x : inner_products)
        x /= M;

    return inner_products;
}

void NoiseExperiment::run_ddsim_noise()
{
    cerr << "+++ [ddsim-noise @ " << this->name << "] Computing DDSIM NOISY  execution for " << this->name << endl;
    clock_t begin = clock();
    cerr << "+++ [ddsim-noise @ " << this->name << "] Setting up observable (" << this->observable << ") and system..." << endl;
    dd::vEdge obs = this->dd_observable();
    double par_value = 1. / (pow(2., static_cast<double>(this->size())) * static_cast<double>(10 * this->iterations));

    luint d = 0;
    if (this->drg == 0)
    {
        d = this->iterations; // No guess is made, use the expected number of iterations to apply a circuit
    }
    else
    {
        d = this->drg; // A guess is made, try to reduce to this dimension.
    }
    luint M = this->M; // The samples needed to make an accurate guess

    cerr << "+++ [ddsim-noise @ " << this->name << "] Computing the reduced system with dimension reduction guess " << d << ", using " << this->M << " Samples..." << endl;
    clock_t b_iteration = clock();

    /*
    Setup for reduction: create a matrix with the inner product results.
    */

    clock_t r_begin = clock();
    cerr << "+++ [ddsim-noise @ " << this->name << "] Beginning calculation of inner products with d = " << d << "...\n";
    auto inner_products = collect_inner_products(obs, d, M);
    cerr << "+++ [ddsim-noise @ " << this->name << "] Finished calculation of inner products, moving onto expected inner products... \n\n";
    cerr << "+++ [ddsim-noise @ " << this->name << "] Beginning calculation of expected inner products with d = " << d << "...\n";
    auto avg_expected_fids = collect_expected_inner_products(obs, d, M);
    cerr << "+++ [ddsim-noise @ " << this->name << "] Finished calculation of expected inner products, moving onto reduction... \n\n";

    cerr << "IP1: " << endl;
    for (const auto &r : inner_products)
    {
        for (const auto &x : r)
        {

            cerr << x << ", ";
        }
        cerr << endl;
    }

    cerr << "IP2: " << endl;
    for (const auto &x : avg_expected_fids)
        cerr << x << ", ";

    cerr << "\n";

    dd::CMat Ahat = this->noise_model->collect_eta_gamma(inner_products);
    // dd::CMat Ahat = get_A_hat(inner_products);
    dd::CMat I_gscb = this->noise_model->get_I_gscb(Ahat);
    dd::CMat I_gscb_inverse = get_inverse(I_gscb);
    dd::CMat Chat = matmul(Ahat, I_gscb_inverse);
    dd::CMat Chat_k = matrix_power(Chat, this->iterations);
    dd::CVec Chat_e1(d);

    // We collect the first column in Chat.
    for (luint i = 0; i < d; i++)
    {
        Chat_e1[i] = Chat_k[i][0];
    }
    cerr << "Ahat: " << matrix_to_string(Ahat) << endl;
    cerr << "I_gscb: " << matrix_to_string(I_gscb) << endl;
    cerr << "I_gscb^-1: " << matrix_to_string(I_gscb_inverse) << endl;
    cerr << "Chat: " << matrix_to_string(Chat) << endl;
    cerr << "Chat^k: " << matrix_to_string(Chat_k) << endl;
    cerr << "Chat_e1: " << vector_to_string(Chat_e1) << endl;

    dd::CVec results(d + 1);
    CC result = CC(0);

    cerr << "+++ [ddsim-noise @ " << this->name << "] Beginning calculation of fidelity..." << endl;
    for (luint i = 1; i <= d; i++)
    {
        auto eta = I_gscb[i - 1][i - 1]; // We use I_gscb here since that have all the eta in a convinient placement

        CC t1 = (CC(1) / eta) * avg_expected_fids[i - 1];

        CC t2_sum = 0;
        for (luint k = 0; k < i - 1; k++)
        {
            t2_sum += Ahat[k][i - 1] * results[k]; // We use Ahat here for the gamma's, since they are placed nicely here.
                                                   //   cerr << "Ahat[k][j-1] = Ahat[" << k + 1 << "][" << i << "] = " << Ahat[k][i - 1] << endl;
        }

        CC t2 = (CC(1) / eta) * t2_sum;

        // cerr << "t1 = " << t1 << ", t2 = " << t2 << endl;
        results[i - 1] = t1 - t2;
        result += Chat_e1[i - 1] * (t1 - t2);
        // cerr << "Results[" << i - 1 << "] = " << results[i - 1] << endl;
    }

    /*This is just placeholder atm.*/
    clock_t a_iteration = clock();
    clock_t r_end = clock();
    this->fidelity = result.real(); // Should be the fidelity between the result from the reduced system compared to what we are looking for.
    // this->fidelity = 0; // Should be the fidelity between the result from the reduced system compared to what we are looking for.
    cerr << "The fidelity between the expected state and the result from the simulation: " << this->fidelity << endl;
    clock_t end = clock();

    // We store the data
    this->red_time = time_to_double(r_begin, r_end);
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