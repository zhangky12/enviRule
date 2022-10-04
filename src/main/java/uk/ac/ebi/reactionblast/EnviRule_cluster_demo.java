package uk.ac.ebi.reactionblast;

import uk.ac.ebi.reactionblast.tools.Adding_rxn_class;

import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.Map;
import java.util.Set;

public class EnviRule_cluster_demo {

    public static void main(String[] args) throws Exception {
        String rxn_file = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/autoRule_cluster_exp/reactions.txt";
        String out_dir = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/autoRule_cluster_exp/reactions4rules/";

        String new_rxn_file = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/autoRule_cluster_exp/reactions2.txt";
        String new_database = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/autoRule_cluster_exp/tmp1/";

//        Clustering_rxn_class.clustering(rxn_file, out_dir);
//
//        Map<String, Map<Integer, Double>> centerMaps = null;
//
//        FileInputStream fileIn = new FileInputStream(out_dir + "centers.dat");
//        ObjectInputStream objectIn = new ObjectInputStream(fileIn);
//
//        centerMaps = (Map<String, Map<Integer, Double>>) objectIn.readObject();
//        objectIn.close();
//
//        System.out.println(centerMaps.size());
//
//        Map<String, String> map_file = null;
//
//        fileIn = new FileInputStream(out_dir + "map_file.dat");
//        objectIn = new ObjectInputStream(fileIn);
//        map_file = (Map<String, String>) objectIn.readObject();
//        objectIn.close();
//
//        System.out.println(map_file);
        Set<String> updated_rxn_files = Adding_rxn_class.adding(new_rxn_file, out_dir, new_database);
        System.out.println(updated_rxn_files);




    }


}
