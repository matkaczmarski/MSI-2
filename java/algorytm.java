/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package javabbob;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 *
 * @author Kuba
 */
public class algorytm 
{
    private static double[] globalBest;
    private static int dim;
    private static int populationSize = 50;
    
    private static double evaluate(double[] individual)
    {
        return individual[0] * individual[0] + individual[1] * individual[1];
    }
    
    public static void main(String[] args)
    {
        Random rand = new Random();
        dim = 2;
        
        double F = 0.5;
        double CR = 0.5;
        double w = 0.75;
        double c1 = 0.2;
        double c2 = 0.5;
        
        int maxIteration = 50000;
        
        ArrayList<double[]> generation = new ArrayList<double[]>();
        ArrayList<double[]> velocities = new ArrayList<double[]>();
        ArrayList<double[]> bests = new ArrayList<double[]>();
        
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
            if (i == 0)
                globalBest = generation.get(i);
            else if (evaluate(generation.get(i)) < evaluate(globalBest))
                globalBest = generation.get(i);
        }
        
        int iteration = 0;
        boolean stop = false;
        while (true)
        {
            if (iteration == maxIteration)
                break;
            System.out.println(iteration + "");
            iteration++;
            for (int i = 0; i < generation.size(); i++)
            {
                int r1, r2, r3;
                while (true)
                {
                    r1 = Math.abs(rand.nextInt()) % populationSize;
                    if (r1 == i)
                        continue;
                    r2 = Math.abs(rand.nextInt()) % populationSize;
                    if ((r2 == i) || (r2 == r1))
                        continue;
                    r3 = Math.abs(rand.nextInt()) % populationSize;
                    if ((r3 == i) || (r3 == r1) || (r3 == r2))
                        continue;
                    break;
                }
                
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
                    
                    if (evaluate(u) < evaluate(generation.get(i)))
                    {
                        generation.remove(i);
                        generation.add(i, u);
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
                        
                        if (evaluate(TX) < evaluate(generation.get(i)))
                        {
                            generation.remove(i);
                            generation.add(i, TX);
                        }
                    }
                    
                    if (evaluate(generation.get(i)) < evaluate(bests.get(i)))
                    {
                        bests.remove(i);
                        bests.add(i, generation.get(i));
                    }
                    
                    if (evaluate(generation.get(i)) < evaluate(globalBest))
                        globalBest = generation.get(i);
                }
            }
        }
        
        for (int i = 0; i < dim; i++)
        {
            System.out.print(globalBest[i] + " ");
        }
    }
}

