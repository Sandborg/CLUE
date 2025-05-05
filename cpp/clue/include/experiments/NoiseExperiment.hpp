#ifndef CLUE_NOISE_EX_EX
#define CLUE_NOISE_EX_EX

#include "Linalg.hpp"
#include "dd/Package.hpp"
#include "Experiment.hpp"

using namespace std;

/**
 * Class that define an interface for noise experiments.
 *
 * This class is abstract and the methods required for an experiment to work are the pure virtual methods.
 *
 * Taken from Experiments.hpp and added noise parameters.
 */
class NoiseExperiment
{
protected:
    /** ABSTRACT METHODS FOR AN EXPERIMENT **/
    /* Method to get the number of q-bits for the experiment*/
    virtual luint size() = 0;
    /* Method to get the correct size of the lumping for this experiment */
    virtual luint correct_size() = 0;
    /* Method to get a quivk bound for the lumping size */
    virtual luint bound_size() = 0;
    /* Method to compute the direct lumping of the experiment (usually faster) */
    virtual array<dd::CMat, 2U> direct() = 0;
    /**
     *  Method to obtain the matrix to apply a lumping.
     *
     * This matrix (in the case of QAOA) will be the problem matrix.
     */
    virtual vector<CCSparseVector> matrix() = 0;
    /* Method to obtain the REDUCED begin matrix (if necessary) */
    virtual dd::CMat matrix_B(dd::CMat &) = 0;
    /**
     * Method to obtain a quantum computation to apply a lumping.
     *
     * This method obtains the "problem" circuit in the case of QAOA.
     */
    virtual qc::QuantumComputation *quantum(double) = 0;
    /* Method to obtain the BEGIN circuit (if necessary) */
    virtual qc::QuantumComputation *quantum_B(double) = 0;
    /* Method to change the type of experiment */
    virtual NoiseExperiment *change_exec_type(ExperimentType) = 0;

    // Protected attributes
    string name;            // Name of the experiment
    string observable;      // String representing the observable for the lumping
    luint iterations;       // Number of iterations to perform in an experiment.
    ExperimentType type;    // Type of the experiment. Depending on the type, different methods will be run
    dd::Package<> *package; // dd::Package with the cache information for the size for  whole execution.

    /* Execution attributes */
    bool executed = false;  // Flag indicating if the experiment has been executed or not
    double red_ratio = 1.0; // The reduction ratio obtained in the execution. If no reduction, the value is 1.
    double red_time = -1.0; // Execution time of the reduction.
    double it_time = -1.0;  // Execution time of the iteration.
    double tot_time = -1.0; // Execution time of the iteration.
    double mem_used = -1.0; // Memory usage of the execution

    /* The three probabilities that noise is applied
            Example: P1 = depolarization probability
                     P2 = amplitude damping probability
                     P3 = phase flip probability

            Depolarization    = [I: 1-(3*P1), X: P1, Y: P1, Z: P1]
            Amplitude Damping = [E_0: 1-P2, E_1: P2] E_0 and E_1 is T1 decoherence eq. 6 in Wille's paper
            Phase Flip        = [I (or nothing): 1-P3, Z: P3] This is the T2 decoherence, eq. 7 in Wille's paper
    */
    double p_depolarization = 0.0;
    double p_phaseflip = 0.0;
    double p_amplitude_damp = 0.0;
    dd::fp fidelity = 0;

    // Protected methods
    /* Method to get the observable for use with CLUE */
    virtual CCSparseVector clue_observable();
    /* Method to get the observable for use with DD */
    virtual dd::vEdge dd_observable(); // TODO Currently not working

    /* Method that runs the CLUE reduction (only used when this->type == CLUE) */
    virtual void run_clue();
    /* Method that runs the DDSIM reduction (only used when this->type == DDSIM) */
    virtual void run_ddsim();
    /* Method that runs the DIRECT reduction (only used when this->type == DIRECT) */
    virtual void run_direct();
    /* Method that runs the CLUE reduction (only used when this->type == DDSIM_ALONE) */
    virtual void run_ddsim_alone();

public:
    /** CONSTRUCTORS **/
    NoiseExperiment(string eName, string eObservable, luint eIterations, ExperimentType eType, dd::Package<> *ePackage, double _p_depolarization, double _p_phaseflip, double _p_amplitude_damp)
    {
        this->name = eName;
        this->observable = eObservable;
        this->iterations = eIterations;
        this->type = eType;
        this->package = ePackage;

        if (_p_depolarization > 1)
            throw std::logic_error("The probability for depolarization should not be higher than 1, was given" + std::to_string(_p_depolarization));
        else
            this->p_depolarization = _p_depolarization;

        if (_p_phaseflip > 1)
            throw std::logic_error("The probability for depolarization should not be higher than 1, was given" + std::to_string(_p_phaseflip));
        else
            this->p_phaseflip = _p_phaseflip;

        if (_p_amplitude_damp > 1)
            throw std::logic_error("The probability for depolarization should not be higher than 1, was given" + std::to_string(_p_amplitude_damp));
        else
            this->p_amplitude_damp = _p_amplitude_damp;
    }

    virtual ~NoiseExperiment() = default;
    /* Method that runs the experiment */
    void run();
    /* Method to clean the execution run */
    void clean_exec() { this->executed = false; }
    /* Method to get the string out of an experiment */
    virtual string to_string() = 0;
    /* Method that generate the CSV row for this experiment */
    virtual string to_csv(char = ',');
    /* Method to get the fidelity between expected result and succes state (might only be useful for Grover?)*/
    virtual dd::fp calc_fidelity(dd::vEdge) = 0;

    /* Method to get the total execution time */
    double total_time() { return this->tot_time; }
    /* Method to get the total execution time */
    double reduction_ratio() { return this->red_ratio; }
    /* Method to get the total execution time */
    double reduction_time() { return this->red_time; }
    /* Method to get the total execution time */
    double iteration_time() { return this->it_time; }
    /* Auxiliar method to convert time clocks into double time */
    double time_to_double(clock_t &init, clock_t &end)
    {
        return (double(end - init) / double(CLOCKS_PER_SEC));
    }
};

#endif