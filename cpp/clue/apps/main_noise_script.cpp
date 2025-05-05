#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <bits/stdc++.h>
#include <time.h>
#include <cstdlib>
#include <string>

#include <experiments/NoiseExperiment.hpp>
#include "experiments/NoisyGroverExperiment.hpp"
/*
#include "experiments/SATExperiment.hpp"
#include "experiments/CUTExperiment.hpp"
#include "experiments/BenchmarkExperiment.hpp"
*/

using namespace std;

NoiseExperiment *generate_example(string name, luint size, ExperimentType type, string observable, dd::Package<> *package, NoiseModel *noise_model)
{
    string upper = boost::to_upper_copy<std::string>(name);
    /*
    if (upper == "SAT")
     {
         return SATFormula::random(size, static_cast<luint>(rand()) % (2 * size) + size, true, 3UL, type, package);
     }
     else if (upper == "MAXCUT")
     {
         return UndirectedGraph::random(size, 1. / 3., type, package);
     }
     else if (upper == "SEARCH")
     {
         auto distribution = generate_random_grover_distribution(epsilon);
         auto search = NoisyQuantumSearch::ones_string(size, type, package, distribution, epsilon);
         search->convert_succes_set_qstate();
         return search;
     }
     else
     {
         try
         {
             return new BenchmarkExperiment(size, name, observable, type, package);
         }
         catch (const domain_error &err)
         {
             throw logic_error("The given class of experiments is not recognized.");
         }
     }
    */
    auto search = NoisyQuantumSearch::ones_string(size, type, package, noise_model);
    search->convert_succes_set_qstate();
    return search;
}

vector<string> generate_observables(string observable, luint size)
{
    vector<string> result = vector<string>();
    if (observable == "all")
    {
        result.push_back("H");
        for (luint i = 0; i < static_cast<luint>(pow(2, size)); i++)
        {
            result.push_back(std::to_string(i));
        }
    }
    else
    {
        result.push_back(observable);
    }
    return result;
}

int main_script(string name, ExperimentType type, luint m, luint M, luint repeats, string observable)
{
    double total_time = 0.;
    ofstream out;
    // PROCESSING THE OUTPUT FILE
    string upper = boost::to_upper_copy<std::string>(name);
    string lower = boost::to_lower_copy<std::string>(name);
    string filename;
    if (upper != "SAT" && upper != "MAXCUT" && upper != "SEARCH")
    {
        filename = "benchmark";
    }
    else
    {
        filename = lower;
    }

    filesystem::path out_path = filesystem::path("../../../../tests/quantum/noise_results/[result-cpp]noisy_q_" + filename + "_" + ExperimentType_toString(type) + ".csv");
    out.open(out_path, std::ios::app);
    cout << "##################################################################################" << endl;
    cout << "### EXECUTION ON " << boost::to_upper_copy<std::string>(name) << "[m=" << m << ", M=" << M << ", repeats=" << repeats << ", method=" << type << "]" << endl;
    cout << "##################################################################################" << endl;

    for (luint size = m; size <= M; size++)
    { // We repeat for each size.
      // For each size we also go from starting to finishing epsilon values.

        double starting = 0.001;
        double finishing = 1.0;
        double depolarization_noise = starting;
        double phaseflip_noise = 0.1;
        double amp_damping_noise = 0.2;

        while (depolarization_noise <= finishing)
        {

            for (string obs : generate_observables(observable, size))
            {
                for (luint execution = 1; execution <= repeats; execution++)
                { // We repeat "repeats" times
                    try
                    {
                        dd::Package<> *package = new dd::Package<>(size);
                        auto noise_model = new NoiseModel(size, depolarization_noise, phaseflip_noise, amp_damping_noise);
                        NoiseExperiment *experiment = generate_example(name, size, type, obs, package, noise_model);
                        cout << "##################################################################################" << endl;
                        cout << "Generated example\n\t" << experiment->to_string() << endl;
                        cout << "Current Errors: Depolarization: "
                             << depolarization_noise
                             << ", Amplitude Damping: "
                             << amp_damping_noise
                             << ", Phase flip: "
                             << phaseflip_noise
                             << endl;

                        experiment->run();

                        cout << "### -- Finished execution " << execution << "/" << repeats << "(size=" << size << "): took " << experiment->total_time() << "s." << endl;

                        total_time += experiment->total_time();
                        out << experiment->to_csv() << endl;
                        delete experiment;
                        delete package;
                        delete noise_model;
                    }
                    catch (qc::QFRException &e)
                    {
                        cout << "### -- Error in execution " << execution << "/" << repeats << "(size=" << size << "): " << e.what() << endl;
                    }
                }
            }
            depolarization_noise = depolarization_noise + pow(10, floor(log10(depolarization_noise)));
        }
    }
    double average_time = total_time / static_cast<double>((M - m + 1) * repeats);
    cout << "### Average execution time: " << average_time << endl;
    cout << "##################################################################################" << endl;
    return 0;
}

enum ArgumentValues
{
    type,
    min,
    max,
    repeats,
    observable
};

std::map<std::string, ArgumentValues> create_argument_map()
{
    std::map<std::string, ArgumentValues> m;
    m["-t"] = ArgumentValues::type;
    m["-m"] = ArgumentValues::min;
    m["-M"] = ArgumentValues::max;
    m["-repeats"] = ArgumentValues::repeats;
    m["-r"] = ArgumentValues::repeats;
    m["-obs"] = ArgumentValues::observable;
    return m;
}
static std::map<std::string, ArgumentValues> s_mapArgumentValues = create_argument_map();

int main(int argc, char **argv)
{
    srand(static_cast<unsigned>(time(NULL)));
    string test = "search";
    ExperimentType type = ExperimentType::DDSIM_ALONE;
    luint m = 5, M = 5, repeats = 1;
    string observable = "H";

    if (argc > 1)
    {
        test = argv[1];
        int i = 2;
        while (i < argc)
        {
            switch (s_mapArgumentValues[argv[i]])
            {
            case ArgumentValues::type:
                type = ExperimentType_fromString(string(argv[i + 1]));
                i += 2;
                break;
            case ArgumentValues::min:
                m = stoul(argv[i + 1]);
                i += 2;
                break;
            case ArgumentValues::max:
                M = stoul(argv[i + 1]);
                i += 2;
                break;
            case ArgumentValues::repeats:
                repeats = stoul(argv[i + 1]);
                i += 2;
                break;
            case ArgumentValues::observable:
                observable = argv[i + 1];
                i += 2;
                break;
            default:
                cout << "Error in arguments: found " << argv[i];
                return -1;
            }
        }
    }

    return main_script(test, type, m, M, repeats, observable);
}