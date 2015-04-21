/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 *
 * @author Kuba
 */
public class algorytm 
{
    private static double[] globalBest;
    private static double globalBestEvaluation;
    private static int populationSize = 50;
    
    private static double evaluate(double[] individual)
    {
        return individual[0] * individual[0] + individual[1] * individual[1];
    }
    
    public static void PSO_DE(/*JNIfgeneric fgeneric,*/ int dim, double maxfunevals, Random rand) 
    {
    	double F = 0.5;
        double CR = 0.5;
        double w = 0.75;
        double c1 = 0.2;
        double c2 = 0.5;
        
        //int maxIteration = 50000;
        
        ArrayList<double[]> generation = new ArrayList<double[]>();
        ArrayList<double[]> velocities = new ArrayList<double[]>();
        ArrayList<double[]> bests = new ArrayList<double[]>();
        ArrayList<Double> evaluations = new ArrayList<Double>();
        ArrayList<Double> bestsEvaluations = new ArrayList<Double>();
        
        for (int i = 0; i < populationSize; i++)
        {
            double[] individual = new double[dim];
            for (int j = 0; j < dim; j++)
                individual[j] = rand.nextDouble() * rand.nextInt();
            generation.add(individual);
            bests.add(individual);
            
            double[] velocity = new double[dim];
            for (int j = 0; j < dim; j++)
                velocity[j] = rand.nextDouble() * rand.nextInt();
            velocities.add(velocity);
        }
        
        for (int i = 0; i < generation.size(); i++)
        {
            double evaluation = evaluate(generation.get(i));
            evaluations.add(evaluation);
            bestsEvaluations.add(evaluation);
            
            if (i == 0)
            {
                globalBest = generation.get(i);
		globalBestEvaluation = evaluate(globalBest);
            }
            else
            {
                if (evaluations.get(i) < globalBestEvaluation)
                {
                    globalBest = generation.get(i);
                    globalBestEvaluation = evaluations.get(i);
                }
            }
        }
        
        int iteration = 0;
        boolean stop = false;
        while (true)
        {
            if (iteration == maxfunevals)
                break;
            System.out.println(iteration + "");
            iteration++;
            for (int i = 0; i < generation.size(); i++)
            {
	             Set<Integer> randomSet = new LinkedHashSet<Integer>();
	             while (randomSet.size() < 3)
	                 randomSet.add(rand.nextInt(populationSize));
	             
	             Object[] randoms = randomSet.toArray();
	             int r1 = (int)randoms[0];
	             int r2 = (int)randoms[1];
	             int r3 = (int)randoms[2];
                
                double[] m = new double[dim];
                for (int j = 0; j < dim; j++)
                    m[j] = generation.get(r1)[j] + F * (generation.get(r2)[j] - generation.get(r3)[j]);
                
                double[] u = new double[dim];
                for (int j = 0; j < dim; j++)
                {
                    int jRand = Math.abs(rand.nextInt()) % dim + 1;
                    if (rand.nextDouble() < CR || j == jRand)
                        u[j] = m[j];
                    else
                        u[j] = generation.get(i)[j];
                    
                    double uEvaluation = evaluate(u);
                    if (uEvaluation < evaluations.get(i))
                    {
                        generation.remove(i);
                        generation.add(i, u);
                        evaluations.remove(i);
                        evaluations.add(i, uEvaluation);
                    }
                    else
                    {
                        double[] TX = new double[dim];
                        double[] velocity = velocities.get(i);
                        
                        double R1 = rand.nextDouble();
                        double R2 = rand.nextDouble();
                        
                        for (int k = 0; k < dim; k++)
                        {
                            velocity[k] = w * velocity[k] + c1 * R1 * (bests.get(i)[k] - generation.get(i)[k]) + c2 * R2 * (globalBest[k] - generation.get(i)[k]);
                            TX[k] = generation.get(i)[k] + velocity[k];
                        }
                        
                        velocities.remove(i);
                        velocities.add(i, velocity);
                        
                        double TXEvaluation = evaluate(TX);
                        if (TXEvaluation < evaluations.get(i))
                        {
                            generation.remove(i);
                            generation.add(i, TX);
                            evaluations.remove(i);
                            evaluations.add(i, TXEvaluation);
                        }
                    }
                    
                    if (evaluations.get(i) < bestsEvaluations.get(i))
                    {
                        bests.remove(i);
                        bests.add(i, generation.get(i));
                        bestsEvaluations.remove(i);
                        bestsEvaluations.add(i, evaluations.get(i));
                    }
                    
                    if (evaluations.get(i) < globalBestEvaluation)
                    {
                        globalBest = generation.get(i);
                        globalBestEvaluation = evaluations.get(i);
                    }
                }
            }
        }
        
        for (int i = 0; i < dim; i++)
        {
            System.out.print(globalBest[i] + " ");
        }
    }
    
    public static void main(String[] args)
    {
        Random rand = new Random();
        int dim = 2;
        
        PSO_DE(dim, 50000, rand);
        
    }
}

