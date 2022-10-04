package uk.ac.ebi;
import uk.ac.ebi.reactionblast.fingerprints.RandomNumber;

import java.util.*;

public class testRandomNumber {

    public static void main(String[] args) throws Exception{
        RandomNumber randomNumberGen = new RandomNumber();
        Map<Integer, List<Long>> repetitions = new HashMap<Integer, List<Long>>();

        for(long i=0; i<1000000; i++){
            int randomNumber = randomNumberGen.generateMersenneTwisterRandomNumber(1024, i);

            if(!repetitions.containsKey(randomNumber)) repetitions.put(randomNumber, new ArrayList<Long>());
            repetitions.get(randomNumber).add(i);
        }

        System.out.println("Size of features: " + repetitions.size());
    }
}
