package javabbob;

import java.util.Random;
import java.util.Calendar;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.Set;

/** Wrapper class running an entire BBOB timing experiment.
 * Measure the CPU-time of an experiment using an algorithm (in the example
 * MY_OPTIMIZER) on the function f8 of the BBOB testbeds.
 *
 * Information provided by running the timing experiment using your algorithm
 * are required in your BBOB workshop paper.
 */
public class ExampleTiming {

    /** optimiser.
     * This optimiser takes as argument an instance of JNIfgeneric
     * which have all the information on the problem to solve.
     * Only the methods getFtarget() and evaluate(x) of the class JNIfgeneric
     * are used.<p>
     *
     * This method also takes as argument an instance of Random since one
     * might want to set the seed of the random search.<p>
     *
     * The optimiser generates random vectors evaluated on fgeneric until
     * the number of function evaluations is greater than maxfunevals or
     * a function value smaller than the target given by fgeneric.getFtarget()
     * is attained.
     * @param fgeneric an instance JNIfgeneric object
     * @param dim an integer giving the dimension of the problem
     * @param maxfunevals an integer giving the maximum number of function evaluations
     * @param rand an instance of Random
     */
    
     final static double F = 0.5;
    final static double CR = 0.5;
    final static double w = 0.64; //0.75;
    final static double c1 = 1.4; //0.4;//0.2;
    final static double c2 = 1.4; //0.4;//0.5;
    final static double THRESHOLD = 0.5;
    final static int BOUND = 5;
    final static int populationSize = 50;

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

    /** Main method for running the whole BBOB timing experiment.
     *  It will run an optimizer (in this case MY_OPTIMIZER) on the function
     *  8 of the BBOB testbeds in dimensions from 2 to 40. The run will last
     *  at least 30 seconds for each dimension (if the algorithm stopped,
     *  it is restarted). The method will display (in the standard output)
     *  information such as the number of runs and the time per function
     *  evaluations in seconds for each dimension.
     */
    public static void main(String[] args) {

        /* your own internal stuff */
        /* External initialization of MY_OPTIMIZER */
        Random rand = new Random(System.currentTimeMillis());

        final int dim[] = {2, 3, 5, 10, 20};//, 40};
        double timings[] = new double[dim.length];
        int runs[] = new int[dim.length];
        int dims[] = new int[dim.length];
        long t0;
        int idx_dim, nbrun;
        JNIfgeneric fgeneric = new JNIfgeneric();

        JNIfgeneric.Params params = new JNIfgeneric.Params();
        params.algName = "PSO-DE Hybrid";
        params.comments = "Particle Swarm Optimization hybridized with Differential Evolution algorithm.";
        String outputPath = "PSO_DE timing";

        if ( JNIfgeneric.makeBBOBdirs(outputPath, false) ) {
            System.out.println("BBOB data directories at " + outputPath
                    + " created.");
        } else {
            System.out.println("Error! BBOB data directories at " + outputPath
                    + " NOT created, stopping.");
            return; 
        };

        for (idx_dim = 0; idx_dim < dim.length; idx_dim++) {
            nbrun = 0;
            fgeneric.initBBOB(8, 1, dim[idx_dim], outputPath, params);
            t0 = System.currentTimeMillis();
            while (System.currentTimeMillis() - t0 < 30000) {
                PSO_DE(fgeneric, dim[idx_dim], dim[idx_dim]*100000, rand); /*adjust maxfunevals*/
                nbrun ++;
            }

            timings[idx_dim] = (double)(System.currentTimeMillis() - t0)/ 3600000. / (double)fgeneric.getEvaluations();
            dims[idx_dim] = dim[idx_dim];
            runs[idx_dim] = nbrun;
            fgeneric.exitBBOB();
            System.out.print("Dimensions:");
            for (int i = 0; i <= idx_dim; i++) {
                System.out.printf(" %11d ", dims[i]);
            }
            System.out.print("\n      runs:");
            for (int i = 0; i <= idx_dim; i++) {
                System.out.printf(" %11d ", runs[i]);
            }
            System.out.print("\n times [s]:");
            for (int i = 0; i <= idx_dim; i++) {
                System.out.printf(" %11.1e ", timings[i]);
            }
            System.out.println();
        }
    }
}
