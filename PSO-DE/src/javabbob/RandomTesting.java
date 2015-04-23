/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package javabbob;

import java.util.LinkedHashSet;
import java.util.Random;
import java.util.Set;

/**
 *
 * @author Mateusz
 */
public class RandomTesting {
    
    public static void main(String[] args)
    {
        long seed = System.currentTimeMillis();
        Random rand = new Random(seed);
        Set<Integer> randomSet = new LinkedHashSet<Integer>();
        randomSet.add(4);
	             while (randomSet.size() < 4)
	                 randomSet.add(rand.nextInt(50));
	             
	             Object[] randoms = randomSet.toArray();
	             int r1 = (int)randoms[0];
	             int r2 = (int)randoms[1];
	             int r3 = (int)randoms[2];
                     int r4 = (int)randoms[3];
                
                     System.out.println("" + r1);
                     System.out.println("" + r2);
                     System.out.println("" + r3);
                     System.out.println("" + r4);
    }
}
