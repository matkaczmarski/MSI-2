package javabbob;

import java.util.Random;
import java.util.Calendar;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.Set;

/**
 * Wrapper class running an entire BBOB experiment. It illustrates the
 * benchmarking of MY_OPTIMIZER on the noise-free testbed or the noisy testbed
 * (change the ifun loop in this case as given below). This class reimplements
 * the content of exampleexperiment.c from the original C version of the BBOB
 * code.
 */
public class ExampleExperiment {

    final static double F = 0.5;
    final static double CR = 0.5;
    final static double w = 0.64; //0.75;
    final static double c1 = 1.4; //0.4;//0.2;
    final static double c2 = 1.4; //0.4;//0.5;
    final static double THRESHOLD = 0.5;
    final static int BOUND = 5;
    final static int populationSize = 50;
    final static int NO_IMPROVEMENT_ITERATIONS_THRESHOLD = 5;

    /**
     * This optimiser takes as argument an instance of JNIfgeneric which have
     * all the information on the problem to solve. Only the methods
     * getFtarget() and evaluate(x) of the class JNIfgeneric are used.<p>
     * This method also takes as argument an instance of Random since one might
     * want to set the seed of the random search.<p>
     * The optimiser generates random vectors evaluated on fgeneric until the
     * number of function evaluations is greater than maxfunevals or a function
     * value smaller than the target given by fgeneric.getFtarget() is attained.
     * The parameter maxfunevals to avoid problem when comparing it to
     * 1000000000*dim where dim is the dimension of the problem.
     *
     * @param fgeneric an instance JNIfgeneric object
     * @param dim an integer giving the dimension of the problem
     * @param maxfunevals the maximum number of function evaluations
     * @param rand an instance of Random
     */
    public static void PSO_DE(JNIfgeneric fgeneric, int dim, double maxfunevals, Random rand) {

        double[] globalBest = {0, 1};
        double globalBestEvaluation = 0;
        double ftarget = fgeneric.getFtarget();

        ArrayList<double[]> generation = new ArrayList<>();
        ArrayList<double[]> velocities = new ArrayList<>();
        ArrayList<double[]> bests = new ArrayList<>();
        ArrayList<Double> evaluations = new ArrayList<>();
        ArrayList<Double> bestsEvaluations = new ArrayList<>();

        for (int i = 0; i < populationSize; i++) {
            double[] individual = new double[dim];
            for (int j = 0; j < dim; j++) {
                individual[j] = rand.nextDouble() * 2 * BOUND - BOUND;
            }
            generation.add(individual);
            bests.add(individual);

            double[] velocity = new double[dim];
            for (int j = 0; j < dim; j++) {
                velocity[j] = rand.nextDouble() * 2 * BOUND - BOUND;
            }
            velocities.add(velocity);
        }

        for (int i = 0; i < generation.size(); i++) {
            double evaluation = fgeneric.evaluate(generation.get(i));
            evaluations.add(evaluation);
            bestsEvaluations.add(evaluation);

            if (i == 0) {
                globalBest = generation.get(i);
                globalBestEvaluation = fgeneric.evaluate(globalBest);
            } else {
                if (evaluations.get(i) < globalBestEvaluation) {
                    globalBest = generation.get(i);
                    globalBestEvaluation = evaluations.get(i);
                }
            }
        }

        int iteration = 0;
        while (true) {
            if (iteration >= maxfunevals) {
                break;
            }

            iteration += generation.size();
            for (int i = 0; i < generation.size(); i++) {
                Set<Integer> randomSet = new LinkedHashSet<>();
                randomSet.add(i);
                while (randomSet.size() < 4) {
                    randomSet.add(rand.nextInt(populationSize));
                }

                Object[] randoms = randomSet.toArray();
                int r1 = (int) randoms[1];
                int r2 = (int) randoms[2];
                int r3 = (int) randoms[3];

                double[] m = new double[dim];
                for (int j = 0; j < dim; j++) {
                    m[j] = generation.get(r1)[j] + F * (generation.get(r2)[j] - generation.get(r3)[j]);
                }

                double[] u = new double[dim];
                for (int j = 0; j < dim; j++) {
                    int jRand = Math.abs(rand.nextInt()) % dim + 1;
                    if (rand.nextDouble() < CR || j == jRand) {
                        u[j] = m[j];
                    } else {
                        u[j] = generation.get(i)[j];
                    }
                }
                double uEvaluation = fgeneric.evaluate(u);
                if (uEvaluation < evaluations.get(i)) {
                    generation.remove(i);
                    generation.add(i, u);
                    evaluations.remove(i);
                    evaluations.add(i, uEvaluation);
                } else {
                    double[] TX = new double[dim];
                    double[] velocity = velocities.get(i);

                    double R1 = rand.nextDouble();
                    double R2 = rand.nextDouble();

                    for (int k = 0; k < dim; k++) {
                        velocity[k] = w * velocity[k] + c1 * R1 * (bests.get(i)[k] - generation.get(i)[k]) + c2 * R2 * (globalBest[k] - generation.get(i)[k]);
                        TX[k] = generation.get(i)[k] + velocity[k];
                    }

                    velocities.remove(i);
                    velocities.add(i, velocity);

                    double TXEvaluation = fgeneric.evaluate(TX);
                    if (TXEvaluation < evaluations.get(i)) {
                        generation.remove(i);
                        generation.add(i, TX);
                        evaluations.remove(i);
                        evaluations.add(i, TXEvaluation);
                    }
                }

                if (evaluations.get(i) < bestsEvaluations.get(i)) {
                    bests.remove(i);
                    bests.add(i, generation.get(i));
                    bestsEvaluations.remove(i);
                    bestsEvaluations.add(i, evaluations.get(i));
                }

                if (evaluations.get(i) < globalBestEvaluation) {
                    globalBest = generation.get(i);
                    globalBestEvaluation = evaluations.get(i);
                }
            }

            if (globalBestEvaluation < ftarget) {
                break;
            }
        }
    }

    public static void PSO_DE_withResets(JNIfgeneric fgeneric, int dim, double maxfunevals, Random rand) {

        double[] globalBest = {0, 1};
        double globalBestEvaluation = 0;
        double ftarget = fgeneric.getFtarget();

        ArrayList<double[]> generation = new ArrayList<>();
        ArrayList<double[]> velocities = new ArrayList<>();
        ArrayList<double[]> bests = new ArrayList<>();
        ArrayList<Double> evaluations = new ArrayList<>();
        ArrayList<Double> bestsEvaluations = new ArrayList<>();

        for (int i = 0; i < populationSize; i++) {
            double[] individual = new double[dim];
            for (int j = 0; j < dim; j++) {
                individual[j] = rand.nextDouble() * 2 * BOUND - BOUND;
            }
            generation.add(individual);
            bests.add(individual);

            double[] velocity = new double[dim];
            for (int j = 0; j < dim; j++) {
                velocity[j] = rand.nextDouble() * 2 * BOUND - BOUND;
            }
            velocities.add(velocity);
        }

        for (int i = 0; i < generation.size(); i++) {
            double evaluation = fgeneric.evaluate(generation.get(i));
            evaluations.add(evaluation);
            bestsEvaluations.add(evaluation);

            if (i == 0) {
                globalBest = generation.get(i);
                globalBestEvaluation = fgeneric.evaluate(globalBest);
            } else {
                if (evaluations.get(i) < globalBestEvaluation) {
                    globalBest = generation.get(i);
                    globalBestEvaluation = evaluations.get(i);
                }
            }
        }

        int iteration = 0;
        boolean improved;
        int iterationsWithoutImprovement = 0;
        while (iteration < maxfunevals) {

            improved = false;
            iteration += generation.size();
            for (int i = 0; i < generation.size(); i++) {

                Set<Integer> randomSet = new LinkedHashSet<>();
                randomSet.add(i);
                while (randomSet.size() < 4) {
                    randomSet.add(rand.nextInt(populationSize));
                }

                Object[] randoms = randomSet.toArray();
                int r1 = (int) randoms[1];
                int r2 = (int) randoms[2];
                int r3 = (int) randoms[3];

                double[] m = new double[dim];
                for (int j = 0; j < dim; j++) {
                    m[j] = generation.get(r1)[j] + F * (generation.get(r2)[j] - generation.get(r3)[j]);
                }

                double[] u = new double[dim];
                for (int j = 0; j < dim; j++) {
                    int jRand = Math.abs(rand.nextInt()) % dim + 1;
                    if (rand.nextDouble() < CR || j == jRand) {
                        u[j] = m[j];
                    } else {
                        u[j] = generation.get(i)[j];
                    }
                }
                double uEvaluation = fgeneric.evaluate(u);
                if (uEvaluation < evaluations.get(i)) {
                    generation.remove(i);
                    generation.add(i, u);
                    evaluations.remove(i);
                    evaluations.add(i, uEvaluation);
                } else {
                    double[] TX = new double[dim];
                    double[] velocity = velocities.get(i);

                    double R1 = rand.nextDouble();
                    double R2 = rand.nextDouble();

                    for (int k = 0; k < dim; k++) {
                        velocity[k] = w * velocity[k] + c1 * R1 * (bests.get(i)[k] - generation.get(i)[k]) + c2 * R2 * (globalBest[k] - generation.get(i)[k]);
                        TX[k] = generation.get(i)[k] + velocity[k];
                    }

                    velocities.remove(i);
                    velocities.add(i, velocity);

                    double TXEvaluation = fgeneric.evaluate(TX);
                    if (TXEvaluation < evaluations.get(i)) {
                        generation.remove(i);
                        generation.add(i, TX);
                        evaluations.remove(i);
                        evaluations.add(i, TXEvaluation);
                    }
                }

                if (evaluations.get(i) < bestsEvaluations.get(i)) {
                    bests.remove(i);
                    bests.add(i, generation.get(i));
                    bestsEvaluations.remove(i);
                    bestsEvaluations.add(i, evaluations.get(i));
                }

                if (evaluations.get(i) < globalBestEvaluation) {
                    globalBest = generation.get(i);
                    globalBestEvaluation = evaluations.get(i);
                    improved = true;
                }
            }

            if (!improved) {
                if (++iterationsWithoutImprovement >= NO_IMPROVEMENT_ITERATIONS_THRESHOLD) {
                    //liczenie rozproszenia itede???????????????????
                    return;
                }
            } else {
                iterationsWithoutImprovement = 0;
            }

            if (globalBestEvaluation < ftarget) {
                break;
            }
        }
    }

    public static void PSO_DE_modified(JNIfgeneric fgeneric, int dim, double maxfunevals, Random rand) {

        double[] globalBest = {0, 1};
        double globalBestEvaluation = 0;
        double ftarget = fgeneric.getFtarget();

        ArrayList<double[]> generation = new ArrayList<>();
        ArrayList<double[]> velocities = new ArrayList<>();
        ArrayList<double[]> bests = new ArrayList<>();
        ArrayList<Double> evaluations = new ArrayList<>();
        ArrayList<Double> bestsEvaluations = new ArrayList<>();

        for (int i = 0; i < populationSize; i++) {
            double[] individual = new double[dim];
            for (int j = 0; j < dim; j++) {
                individual[j] = rand.nextDouble() * 2 * BOUND - BOUND;
            }
            generation.add(individual);
            bests.add(individual);

            double[] velocity = new double[dim];
            for (int j = 0; j < dim; j++) {
                velocity[j] = rand.nextDouble() * 2 * BOUND - BOUND;
            }
            velocities.add(velocity);
        }

        for (int i = 0; i < generation.size(); i++) {
            double evaluation = fgeneric.evaluate(generation.get(i));
            evaluations.add(evaluation);
            bestsEvaluations.add(evaluation);

            if (i == 0) {
                globalBest = generation.get(i);
                globalBestEvaluation = fgeneric.evaluate(globalBest);
            } else {
                if (evaluations.get(i) < globalBestEvaluation) {
                    globalBest = generation.get(i);
                    globalBestEvaluation = evaluations.get(i);
                }
            }
        }

        int iteration = 0;
        while (true) {
            if (iteration >= maxfunevals) {
                break;
            }

            iteration += generation.size();
            for (int i = 0; i < generation.size(); i++) {
                Set<Integer> randomSet = new LinkedHashSet<>();
                randomSet.add(i);
                while (randomSet.size() < 4) {
                    randomSet.add(rand.nextInt(populationSize));
                }

                Object[] randoms = randomSet.toArray();
                int r1 = (int) randoms[1];
                int r2 = (int) randoms[2];
                int r3 = (int) randoms[3];

                double[] m = new double[dim];
                for (int j = 0; j < dim; j++) {
                    m[j] = generation.get(r1)[j] + F * (generation.get(r2)[j] - generation.get(r3)[j]);
                }

                double[] u = new double[dim];
                for (int j = 0; j < dim; j++) {
                    int jRand = Math.abs(rand.nextInt()) % dim + 1;
                    if (rand.nextDouble() < CR || j == jRand) {
                        u[j] = m[j];
                    } else {
                        u[j] = generation.get(i)[j];
                    }
                }
                double uEvaluation = fgeneric.evaluate(u);

                if (uEvaluation < evaluations.get(i)) {

                    double[] velocity = u.clone();

                    for (int k = 0; k < dim; k++) {
                        velocity[k] -= generation.get(i)[k];
                    }

                    velocities.remove(i);
                    velocities.add(i, velocity);

                    generation.remove(i);
                    generation.add(i, u);
                    evaluations.remove(i);
                    evaluations.add(i, uEvaluation);

                } else {
                    double[] TX = new double[dim];
                    double[] velocity = velocities.get(i);

                    double R1 = rand.nextDouble();
                    double R2 = rand.nextDouble();

                    for (int k = 0; k < dim; k++) {
                        velocity[k] = w * velocity[k] + c1 * R1 * (bests.get(i)[k] - generation.get(i)[k]) + c2 * R2 * (globalBest[k] - generation.get(i)[k]);
                        TX[k] = generation.get(i)[k] + velocity[k];
                    }

                    velocities.remove(i);
                    velocities.add(i, velocity);

                    double TXEvaluation = fgeneric.evaluate(TX);
                    if (TXEvaluation < evaluations.get(i)) {
                        generation.remove(i);
                        generation.add(i, TX);
                        evaluations.remove(i);
                        evaluations.add(i, TXEvaluation);
                    }
                }

                if (evaluations.get(i) < bestsEvaluations.get(i)) {
                    bests.remove(i);
                    bests.add(i, generation.get(i));
                    bestsEvaluations.remove(i);
                    bestsEvaluations.add(i, evaluations.get(i));
                }

                if (evaluations.get(i) < globalBestEvaluation) {
                    globalBest = generation.get(i);
                    globalBestEvaluation = evaluations.get(i);
                }

            }
            if (globalBestEvaluation < ftarget) {
                break;
            }
        }
    }

    public static void DE(JNIfgeneric fgeneric, int dim, double maxfunevals, Random rand) {

        double[] globalBest = {0, 1};
        double globalBestEvaluation = 0;
        double ftarget = fgeneric.getFtarget();

        ArrayList<double[]> generation = new ArrayList<>();
        ArrayList<double[]> bests = new ArrayList<>();
        ArrayList<Double> evaluations = new ArrayList<>();
        ArrayList<Double> bestsEvaluations = new ArrayList<>();

        for (int i = 0; i < populationSize; i++) {
            double[] individual = new double[dim];
            for (int j = 0; j < dim; j++) {
                individual[j] = rand.nextDouble() * 2 * BOUND - BOUND;
            }
            generation.add(individual);
            bests.add(individual);

        }

        for (int i = 0; i < generation.size(); i++) {
            double evaluation = fgeneric.evaluate(generation.get(i));
            evaluations.add(evaluation);
            bestsEvaluations.add(evaluation);

            if (i == 0) {
                globalBest = generation.get(i);
                globalBestEvaluation = fgeneric.evaluate(globalBest);
            } else {
                if (evaluations.get(i) < globalBestEvaluation) {
                    globalBest = generation.get(i);
                    globalBestEvaluation = evaluations.get(i);
                }
            }
        }

        int iteration = 0;
        while (true) {
            if (iteration >= maxfunevals) {
                break;
            }

            iteration += generation.size();
            for (int i = 0; i < generation.size(); i++) {
                Set<Integer> randomSet = new LinkedHashSet<>();
                randomSet.add(i);
                while (randomSet.size() < 4) {
                    randomSet.add(rand.nextInt(populationSize));
                }

                Object[] randoms = randomSet.toArray();
                int r1 = (int) randoms[1];
                int r2 = (int) randoms[2];
                int r3 = (int) randoms[3];

                double[] m = new double[dim];
                for (int j = 0; j < dim; j++) {
                    m[j] = generation.get(r1)[j] + F * (generation.get(r2)[j] - generation.get(r3)[j]);
                }

                double[] u = new double[dim];
                for (int j = 0; j < dim; j++) {
                    int jRand = Math.abs(rand.nextInt()) % dim + 1;
                    if (rand.nextDouble() < CR || j == jRand) {
                        u[j] = m[j];
                    } else {
                        u[j] = generation.get(i)[j];
                    }
                }
                double uEvaluation = fgeneric.evaluate(u);
                if (uEvaluation < evaluations.get(i)) {
                    generation.remove(i);
                    generation.add(i, u);
                    evaluations.remove(i);
                    evaluations.add(i, uEvaluation);
                }

                if (evaluations.get(i) < bestsEvaluations.get(i)) {
                    bests.remove(i);
                    bests.add(i, generation.get(i));
                    bestsEvaluations.remove(i);
                    bestsEvaluations.add(i, evaluations.get(i));
                }

                if (evaluations.get(i) < globalBestEvaluation) {
                    globalBest = generation.get(i);
                    globalBestEvaluation = evaluations.get(i);
                }
            }

            if (globalBestEvaluation < ftarget) {
                break;
            }
        }
    }

    /**
     * Main method for running the whole BBOB experiment. Executing this method
     * runs the experiment. The first command-line input argument is
     * interpreted: if given, it denotes the data directory to write in (in
     * which case it overrides the one assigned in the preamble of the method).
     */
    public static void main(String[] args) {

        /* run variables for the function/dimension/instances loops */
        final int dim[] = {2, 3, 5, 10, 20};// {5, 20};//{2, 3, 5, 10, 20, 40};
        final int instances[] = {1, 2, 3, 4, 5, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50};
        int idx_dim, ifun, idx_instances, independent_restarts;
        double maxfunevals;
        String outputPath;

        JNIfgeneric fgeneric = new JNIfgeneric();
        /* The line above loads the library cjavabbob at the core of
         * JNIfgeneric. It will throw runtime errors if the library is not
         * found.
         */

        /**
         * ************************************************
         * BBOB Mandatory initialization *
         * ***********************************************
         */
        JNIfgeneric.Params params = new JNIfgeneric.Params();
        /* Modify the following parameters, choosing a different setting
         * for each new experiment */
        params.algName = "PSO-DE Hybrid";
        params.comments = "Particle Swarm Optimization hybridized with Differential Evolution algorithm.";
        outputPath = "PSO_DE_withResets";

        if (args.length > 0) {
            outputPath = args[0]; // Warning: might override the assignment above.
        }

        /* Creates the folders for storing the experimental data. */
        if (JNIfgeneric.makeBBOBdirs(outputPath, false)) {
            System.out.println("BBOB data directories at " + outputPath
                    + " created.");
        } else {
            System.out.println("Error! BBOB data directories at " + outputPath
                    + " was NOT created, stopping.");
            return;
        };


        /* External initialization of MY_OPTIMIZER */
        long seed = System.currentTimeMillis();
        Random rand = new Random(seed);
        System.out.println("PSO_DE seed: " + seed);

        /* record starting time (also useful as random number generation seed) */
        long t0 = System.currentTimeMillis();

        /* To make the noise deterministic, uncomment the following block. */
        /* int noiseseed = 30; // or (int)t0
         * fgeneric.setNoiseSeed(noiseseed);
         * System.out.println("seed for the noise set to: "+noiseseed); */

        /* now the main loop */
        for (idx_dim = 0; idx_dim < dim.length; idx_dim++)//idx_dim < 6; idx_dim++) 
        {
            /* Function indices are from 1 to 24 (noiseless) or from 101 to 130 (noisy) */
            /* for (ifun = 101; ifun <= 130; ifun++) { // Noisy testbed */
            for (ifun = 1; ifun <= 24; ifun++) { //Noiseless testbed
                for (idx_instances = 0; idx_instances < 15; idx_instances++) {

                    /* Initialize the objective function in fgeneric. */
                    fgeneric.initBBOB(ifun, instances[idx_instances],
                            dim[idx_dim], outputPath, params);
                    /* Call to the optimizer with fgeneric as input */
                    maxfunevals = 1000. * dim[idx_dim]; /* PUT APPROPRIATE MAX. FEVALS */
                    /* 5. * dim is fine to just check everything */

                    independent_restarts = -1;
                    while (fgeneric.getEvaluations() < maxfunevals) {
                        independent_restarts++;
                        PSO_DE_withResets(fgeneric, dim[idx_dim], maxfunevals - fgeneric.getEvaluations(), rand);
                        if (fgeneric.getBest() < fgeneric.getFtarget()) {
                            break;
                        }
                    }

                    System.out.printf("  f%d in %d-D, instance %d: FEs=%.0f with %d restarts,", ifun, dim[idx_dim],
                            instances[idx_instances], fgeneric.getEvaluations(), independent_restarts);
                    System.out.printf(" fbest-ftarget=%.4e, elapsed time [h]: %.2f\n", fgeneric.getBest() - fgeneric.getFtarget(),
                            (double) (System.currentTimeMillis() - t0) / 3600000.);

                    /* call the BBOB closing function to wrap things up neatly */
                    fgeneric.exitBBOB();
                }

                System.out.println("\ndate and time: " + (new SimpleDateFormat("dd-MM-yyyy HH:mm:ss")).format(
                        (Calendar.getInstance()).getTime()));

            }
            System.out.println("---- dimension " + dim[idx_dim] + "-D done ----");
        }
    }
}
