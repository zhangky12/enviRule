package eawag.envirule;

import org.apache.commons.io.FileUtils;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.Clustering_rxn_class;
import uk.ac.ebi.sim_test;

import java.io.*;
import java.util.*;

public class rxn_adder {

    public static final SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

    public static Set<String> adding(String new_rxn_file, String database, String new_database) throws Exception {

        File directory = new File(new_database);
        if(!directory.exists()){
            directory.mkdir();
        }

        // Load the centerMaps (fingerprints of each group) and map_file (map from group to file) of old reactions
        Map<String, Map<Integer, Double>> centerMaps = null;

        FileInputStream fileIn = new FileInputStream(database + "centers.dat");
        ObjectInputStream objectIn = new ObjectInputStream(fileIn);
        centerMaps = (Map<String, Map<Integer, Double>>) objectIn.readObject();
        objectIn.close();

        Map<String, String> map_file = null;

        fileIn = new FileInputStream(database + "map_file.dat");
        objectIn = new ObjectInputStream(fileIn);
        map_file = (Map<String, String>) objectIn.readObject();
        objectIn.close();

        // Clustering for new reactions
        Map<String, Map<Integer, Double>> reactionFPs = new HashMap<String, Map<Integer, Double>>();
        List<String> reactions = Clustering_rxn_class.parseReactions(new_rxn_file, sp, reactionFPs);

        Map<String, Set<String>> clusters = Clustering_rxn_class.getClusters(reactions, sp, reactionFPs, new_database, true);

        // Load the centerMaps from new clustering
        Map<String, Map<Integer, Double>> new_centerMaps = null;

        fileIn = new FileInputStream(new_database + "centers.dat");
        objectIn = new ObjectInputStream(fileIn);
        new_centerMaps = (Map<String, Map<Integer, Double>>) objectIn.readObject();
        objectIn.close();

        // Compare the new centerMaps with the old centerMaps
        Map<String, String> rxn_combine = new HashMap<>();
        Map<String, String> center_match = new HashMap<>(); // old center to new center
        for (String center: new_centerMaps.keySet()){
            Map<Integer, Double> fingerprint = new_centerMaps.get(center);

            for (String old_center: centerMaps.keySet()){
                if (sim_test.compareDict(fingerprint, centerMaps.get(old_center))==1.0F){
                    String old_rxn_file = map_file.get(old_center);
                    rxn_combine.put(center, old_rxn_file);
                    center_match.put(old_center, center);
                }
            }
        }

        // Combine two centerMaps
        Map<String, Map<Integer, Double>> combineCenterMaps = combineCenterMaps(new_centerMaps, centerMaps, center_match);
        // Write combined centerMaps
        String centerMaps_file = new_database + "centers.dat";
        FileOutputStream fos = new FileOutputStream(centerMaps_file);
        ObjectOutputStream oos = new ObjectOutputStream(fos);
        oos.writeObject(combineCenterMaps);
        oos.close();

        // Copy all the reaction file to current reaction folder, and also record the largest index
        File[] files = new File(database).listFiles();
        int max_index = 0;

        for (File file: files) {
            if (file.getName().contains("-") && file.getName().contains(".txt")){

                int index = Integer.valueOf(file.getName().split("-")[0]);
                if(index > max_index) max_index = index;

                File source = new File(database + file.getName());
                File dest = new File(new_database + file.getName());
                FileUtils.copyFile(source, dest);
            }
        }

        // Combine the reactions
        Set<String> updated_rxn_files = new HashSet<>();
        Map<String, String> new_map_file = new HashMap<>();

        for (String new_center: rxn_combine.keySet()){
            String rxn_file = rxn_combine.get(new_center);

            String[] parse_file_name = rxn_file.split("/");
            String old_file_name = parse_file_name[parse_file_name.length-1];
            String index = old_file_name.split("-")[0];

            Set<String> update_reactions = parseReactions(rxn_file);
            update_reactions.addAll(clusters.get(new_center));

            String new_file_name = index + "-" + update_reactions.size() + ".txt";
            updated_rxn_files.add(new_database + new_file_name);
            // Delete old rxn file first
            File old_file = new File(new_database+old_file_name);
            old_file.delete();
            writeReactions(update_reactions, new_database+new_file_name);
            new_map_file.put(new_center, new_database+new_file_name);
        }

        // For the rest reactions, write them to new rxn files
        max_index ++;

        for (String new_center: clusters.keySet()){
            if (!rxn_combine.containsKey(new_center)){
                String new_file_name = max_index + "-" + clusters.get(new_center).size() + ".txt";
                writeReactions(clusters.get(new_center), new_database+new_file_name);
                new_map_file.put(new_center, new_database+new_file_name);
                updated_rxn_files.add(new_database+new_file_name);
                max_index ++;
            }
        }

        // Write map_file as Object
        String map_File = new_database + "map_file.dat";
        fos = new FileOutputStream(map_File);
        oos = new ObjectOutputStream(fos);
        oos.writeObject(new_map_file);
        oos.close();

        System.out.println("Done");

        return updated_rxn_files;

    }

    public static Set<String> parseReactions(String file) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(file));
        Set<String> res = new HashSet<>();

        String line;
        while ((line = br.readLine()) != null) {
            res.add(line);
        }

        br.close();

        return res;
    }

    public static void writeReactions(Set<String> reactions, String file) throws IOException {

        FileWriter fw = new FileWriter(file);

        for(String reaction: reactions){
            fw.write(reaction + "\n");
        }
        fw.close();
    }

    public static Map<String, Map<Integer, Double>> combineCenterMaps(Map<String, Map<Integer, Double>> centerMaps_source, Map<String, Map<Integer, Double>> centerMaps_target, Map<String, String> center_match){

        /**
         * centerMaps_source: new_centerMaps
         * centerMaps_target: old_centerMaps
         * rxn_combine: the key in rxn_combine is the same with centerMaps_source
         */

        Map<String, Map<Integer, Double>> combined_centerMaps = new HashMap<>();

        for (String new_center: centerMaps_source.keySet()){
            combined_centerMaps.put(new_center, centerMaps_source.get(new_center));
        }

        for (String old_center: centerMaps_target.keySet()){
            if(!center_match.containsKey(old_center)) {
                combined_centerMaps.put(old_center, centerMaps_target.get(old_center));
            }
        }

        return combined_centerMaps;
    }
}
