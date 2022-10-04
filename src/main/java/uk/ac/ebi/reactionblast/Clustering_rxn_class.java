package uk.ac.ebi.reactionblast;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS;
import uk.ac.ebi.sim_test;

import java.io.*;
import java.util.*;

import static uk.ac.ebi.reactionblast.mechanism.helper.Utility.getCanonicalisedBondChangePattern;

public class Clustering_rxn_class {

    public static final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());


//    public static void main(String[] args) throws Exception {
//
//        SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
//
//        String reaction1 = "CC(C)(C)OC=O>>CC(C)(C)O.C(=O)[O-]";
//
//        IReaction Reaction1 = smilesParser.parseReactionSmiles(reaction1);
//        IReaction performAtomAtomMapping1 = performAtomAtomMapping(Reaction1, "null");
//        BondChangeCalculator bcc1 = new BondChangeCalculator(performAtomAtomMapping1);
//        bcc1.computeBondChanges(true, false);
//        System.out.println(sg.create(performAtomAtomMapping1));
//
//        String reaction2 = "CCOC(=O)C>>CCO.CC(=O)[O-]";
//
//        IReaction Reaction2 = smilesParser.parseReactionSmiles(reaction2);
//        IReaction performAtomAtomMapping2 = performAtomAtomMapping(Reaction2, "null");
//        BondChangeCalculator bcc2 = new BondChangeCalculator(performAtomAtomMapping2);
//        bcc2.computeBondChanges(true, false);
//        System.out.println(sg.create(performAtomAtomMapping2));
//
//        Map<Integer, Double> fp1 = getChangedBonds(performAtomAtomMapping1, bcc1);
//        Map<Integer, Double> fp2 = getChangedBonds(performAtomAtomMapping2, bcc2);
//
//        System.out.println(fp1);
//        System.out.println(fp2);
//
//    }

    public static void clustering(String rxn_file, String out_dir) throws Exception {

        File directory = new File(out_dir);
        if(!directory.exists()){
            directory.mkdir();
        }

        Map<String, Map<Integer, Double>> reactionFPs = new HashMap<String, Map<Integer, Double>>();
        List<String> reactions = parseReactions(rxn_file, smilesParser, reactionFPs);
        Map<String, Set<String>> clusters = getClusters(reactions, smilesParser, reactionFPs, out_dir, true);
        Map<String, String> map_file = new HashMap<>();

        System.out.println("Writing results to files...");
        int file_index = 1;
        for (Map.Entry<String, Set<String>> entry: clusters.entrySet()){

            if (entry.getValue().size() >= 1){
                String file_name = out_dir + file_index + "-" + entry.getValue().size() + ".txt";
                map_file.put(entry.getKey(), file_name);
                FileWriter fw = new FileWriter(file_name);
                for(String react: entry.getValue()){
                    fw.write(react + "\n");
                }
                fw.close();
                file_index ++;
            }
        }

        // Write map_file as Object
        String map_File = out_dir + "map_file.dat";
        FileOutputStream fos = new FileOutputStream(map_File);
        ObjectOutputStream oos = new ObjectOutputStream(fos);
        oos.writeObject(map_file);
        oos.close();

        System.out.println("Done");
    }

    public static List<String> parseReactions(String file, SmilesParser smilesParser, Map<String, Map<Integer, Double>> reactionFPs) throws Exception{

        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        String[] line_array;
        String reaction;
        List<String> reactions = new ArrayList<String>();
        int count = 0;

        while ((line = br.readLine()) != null){
            count += 1;
            if(count % 50 == 0){
                System.out.println("--------------------------------------------------------------------------------");
                System.out.println("Processing " + count + "th reaction");
                System.out.println("--------------------------------------------------------------------------------");
            }
            // Sanity check for ">>"
            if(line.indexOf(">>") == -1){
                // If it's not a valid reaction SMIRKS, then skip
                continue;
            }
            line_array = line.split("\n");
            reaction = line_array[0];

            // Skip all the reactions with CO2 as product.
            if(reaction.split(">>").length < 2) continue;
            String product = reaction.split(">>")[1];
            if (product.compareTo("C(=O)=O") == 0) continue;

            reactions.add(reaction);

            IReaction Reaction1 = smilesParser.parseReactionSmiles(reaction);
            IReaction performAtomAtomMapping1 = performAtomAtomMapping(Reaction1, "null");
            BondChangeCalculator bcc1 = new BondChangeCalculator(performAtomAtomMapping1);
            bcc1.computeBondChanges(true, false);

            if(!reactionFPs.containsKey(reaction)){
                reactionFPs.put(reaction, getChangedBonds(performAtomAtomMapping1, bcc1));
            }

        }
        return reactions;
    }

    private static IReaction performAtomAtomMapping(IReaction cdkReaction, String reactionName) throws InvalidSmilesException, AssertionError, Exception {
        cdkReaction.setID(reactionName);
        /*
         RMT for the reaction mapping
         */
        boolean forceMapping = true;//Overrides any mapping present int the reaction
        boolean generate2D = true;//2D perception of the stereo centers
        boolean generate3D = false;//2D perception of the stereo centers
        boolean standardizeReaction = true; //Standardize the reaction
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, forceMapping, generate2D, generate3D, standardizeReaction);
        MappingSolution s = rmt.getSelectedSolution();//Fetch the AAM Solution
        IReaction reaction = s.getReaction();//Fetch Mapped Reaction
        return reaction;
    }

    public static Map<String, Set<String>> getClusters(List<String> reactions, SmilesParser smilesParser, Map<String, Map<Integer, Double>> reactionFPs, String out_dir, boolean write_centers) throws Exception{

        Map<String, Set<String>> centers = new HashMap<String, Set<String>>();
        Set<Integer> toRemove = new HashSet<Integer>();
        float sim;

        Map<String, Map<Integer, Double>> centerMaps = new HashMap<>();

        for(int i=0; i<reactions.size(); i++){
            if(toRemove.contains(i)) continue;
            String reaction1 = reactions.get(i);

            centerMaps.put(reaction1, reactionFPs.get(reaction1));

            for(int j=i; j<reactions.size(); j++){

                if(toRemove.contains(j)) continue;
                String reaction2 = reactions.get(j);

                if(i == j) sim = 1.0F;
                else {
                    sim = sim_test.compareDict(reactionFPs.get(reaction1), reactionFPs.get(reaction2));
                }

                if(sim == 1.0F){
                    if(!centers.containsKey(reaction1)) centers.put(reaction1, new HashSet<String>());
                    centers.get(reaction1).add(reaction2);
                    toRemove.add(j);
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////
        // Write centerMaps as Object
        if (write_centers) {
            String centerMaps_file = out_dir + "centers.dat";
            FileOutputStream fos = new FileOutputStream(centerMaps_file);
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(centerMaps);
            oos.close();
        }

        return centers;
    }

    private static Map<Integer, Double> getChangedBonds(IReaction performAtomAtomMapping1, BondChangeCalculator bcc1) throws CDKException {

        IAtomContainerSet reactants = performAtomAtomMapping1.getReactants();
        IAtomContainerSet products = performAtomAtomMapping1.getProducts();

        Set<IBond> formed_cleaved_bonds_reactants = new HashSet<IBond>();
        Set<IBond> formed_cleaved_bonds_products = new HashSet<IBond>();
        Set<IBond> changed_bonds_reactants = new HashSet<IBond>();
        Set<IBond> changed_bonds_products = new HashSet<IBond>();

        // Initialize fingerprints

        IPatternFingerprinter formedCleavedWFingerprint_reactants = new PatternFingerprinter();
        formedCleavedWFingerprint_reactants.setFingerprintID(performAtomAtomMapping1.getID() + ":" + "Bond Cleaved and Formed");
        IPatternFingerprinter formedCleavedWFingerprint_products = new PatternFingerprinter();
        formedCleavedWFingerprint_products.setFingerprintID(performAtomAtomMapping1.getID() + ":" + "Bond Cleaved and Formed");
        IPatternFingerprinter changedWFingerprint_reactants = new PatternFingerprinter();;
        changedWFingerprint_reactants.setFingerprintID(performAtomAtomMapping1.getID() + ":" + "Bond Change");
        IPatternFingerprinter changedWFingerprint_products = new PatternFingerprinter();;
        changedWFingerprint_products.setFingerprintID(performAtomAtomMapping1.getID() + ":" + "Bond Change");

        // Collect MapId of changed atoms
        Set<Integer> changed_atom_tags = new HashSet<>();
        for(IAtom atom: bcc1.getReactionCenterSet()){
            changed_atom_tags.add(atom.getMapIdx());
        }

        // Include functional groups in reactants
        // note: the order of the following two methods cannot be changed.
        addFunctionalGroups(reactants, changed_atom_tags);
        // Include all the atoms that only show up in products
        newAddedGroupsInProducts(performAtomAtomMapping1, changed_atom_tags);

        // Identify atoms of reaction center in reactant
        Map<Integer, Double> reactants_rc_atoms = new HashMap<>();
        for(int i=0; i<reactants.getAtomContainerCount(); i++){
            IAtomContainer reactant = reactants.getAtomContainer(i);
            for(IAtom atom: reactant.atoms()){
                if(changed_atom_tags.contains(atom.getMapIdx())){
                    String symbol = atom.getSymbol();
                    if(symbol == "C" && atom.isAromatic()){
                        symbol = "c";
                    }else if(symbol == "N" && atom.isAromatic()){
                        symbol = "n";
                    }
                    int hash = -symbol.hashCode();
                    if(!reactants_rc_atoms.containsKey(hash)) reactants_rc_atoms.put(hash, 0.0);
                    reactants_rc_atoms.put(hash, reactants_rc_atoms.get(hash) + 1);
                }
            }
        }

        // Identify atoms of reaction center in products
        Map<Integer, Double> products_rc_atoms = new HashMap<>();
        for(int i=0; i<products.getAtomContainerCount(); i++){
            IAtomContainer product = products.getAtomContainer(i);
            for(IAtom atom: product.atoms()){
                if(changed_atom_tags.contains(atom.getMapIdx())){
                    String symbol = atom.getSymbol();
                    if(symbol == "C" && atom.isAromatic()){
                        symbol = "c";
                    }else if(symbol == "N" && atom.isAromatic()){
                        symbol = "n";
                    }
                    int hash = -symbol.hashCode() - 200;
                    if(!products_rc_atoms.containsKey(hash)) products_rc_atoms.put(hash, 0.0);
                    products_rc_atoms.put(hash, products_rc_atoms.get(hash) + 1);
                }
            }
        }

        // Identify changed or formed/cleaved bonds in reactants
        for(int i=0; i<reactants.getAtomContainerCount(); i++){
            IAtomContainer reactant = reactants.getAtomContainer(i);
            for(int j=0; j<reactant.getBondCount(); j++){
                IBond bond = reactant.getBond(j);
                if(bond.getProperties().containsKey(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION)){
                    if(bond.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) == ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED ||
                            bond.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) == ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED ||
                            bond.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) == ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED_OR_CLEAVED){
                        formed_cleaved_bonds_reactants.add(bond);
                    }else{
                        changed_bonds_reactants.add(bond);
                    }
                }
            }
        }
        // Identify changed or formed/cleaved bonds in products
        for(int i=0; i<products.getAtomContainerCount(); i++){
            IAtomContainer product = products.getAtomContainer(i);
            for(int j=0; j<product.getBondCount(); j++){
                IBond bond = product.getBond(j);

                if(bond.getProperties().containsKey(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION)){
                    if(bond.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) == ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED ||
                            bond.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) == ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED ||
                            bond.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) == ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED_OR_CLEAVED){
                        formed_cleaved_bonds_products.add(bond);
                    }else{
                        changed_bonds_products.add(bond);
                    }
                }
            }
        }

        // Add formed or cleaved bonds to fingerprints
        for(IBond bond: formed_cleaved_bonds_reactants){
            formedCleavedWFingerprint_reactants.add(new Feature(getCanonicalisedBondChangePattern(bond), 1.0));
        }

        for(IBond bond: formed_cleaved_bonds_products){
            formedCleavedWFingerprint_products.add(new Feature(getCanonicalisedBondChangePattern(bond), 1.0));
        }

        // Add changed bonds to fingerprints
        for(IBond bond: changed_bonds_reactants){
            changedWFingerprint_reactants.add(new Feature(getCanonicalisedBondChangePattern(bond), 1.0));
        }

        for(IBond bond: changed_bonds_products){
            changedWFingerprint_products.add(new Feature(getCanonicalisedBondChangePattern(bond), 1.0));
        }

        Map<Integer, Double> hashedFP = getHashedFingerPrint(formedCleavedWFingerprint_reactants.getWeightedHashedFingerPrint(),
                formedCleavedWFingerprint_products.getWeightedHashedFingerPrint(),
                changedWFingerprint_reactants.getWeightedHashedFingerPrint(),
                changedWFingerprint_products.getWeightedHashedFingerPrint());

        hashedFP.putAll(reactants_rc_atoms);
        hashedFP.putAll(products_rc_atoms);

        return hashedFP;

    }

    private static void addFunctionalGroups(IAtomContainerSet reactants, Set<Integer> changed_atom_tags){

        Set<Integer> addedFunctionalGroups = new HashSet<>();

        String[] group_templates = {
                "C(=O)Cl", //acid chloride
                "C(=O)O", // carboxylic acid
                "C(=O)[O-]", // carboxylic acid
                "[S](=O)(=O)(Cl)", //sulfonyl chloride
                "[B](O)(O)", //boronic acid
                "[N](=!@C=!@O)", //isocyanate
                "[N]=[N]=[N]", //azide
                "O=C1N(Br)C(=O)CC1", //NBS brominating agent
                "C=O", //carbonyl
                "ClS(Cl)=O", //thionyl chloride
                "[Mg][Br,Cl]", //grinard(non - disassociated)
                "[#6]S(=O)(=O)[O]", //RSO3 leaving group
                "[O]S(=O)(=O)[O]", //SO4 group
                "[N-]=[N+]=[C]", //diazo - alkyl
//                "C=C", // double bond
        };

        for(int i=0; i<reactants.getAtomContainerCount(); i++){
            IAtomContainer mol = reactants.getAtomContainer(i);
            for(String smarts: group_templates){
                Pattern pattern = SmartsPattern.create(smarts);
                Mappings matches = pattern.matchAll(mol);

                if(matches.count() == 0) continue;
                // matches are int[][] of atom index (not the MapId of atom).
                for(int[] match: matches){
                    Set<Integer> group_match = new HashSet<>();
                    for(int m: match) group_match.add(mol.getAtom(m).getMapIdx());

                    for(IAtom atom: mol.atoms()){
                        // Make sure the atom in reactant is in reaction center and also in a functional group.
                        if(changed_atom_tags.contains(atom.getMapIdx()) && group_match.contains(atom.getMapIdx())){
                            addedFunctionalGroups.addAll(group_match);
                        }
                    }
                }
            }
        }
        changed_atom_tags.addAll(addedFunctionalGroups);
    }

    /**
     * For all the atoms added to products, they need to be considered as reaction center as well
     * @param performAtomAtomMapping
     * @param changed_atom_tags
     */
    private static void newAddedGroupsInProducts(IReaction performAtomAtomMapping, Set<Integer> changed_atom_tags){

        IAtomContainerSet reactants = performAtomAtomMapping.getReactants();
        IAtomContainerSet products = performAtomAtomMapping.getProducts();

        Set<Integer> reactants_atom_tags = new HashSet<>();

        for(int i=0; i<reactants.getAtomContainerCount(); i++){
            IAtomContainer mol = reactants.getAtomContainer(i);
            for(IAtom atom: mol.atoms()){
                reactants_atom_tags.add(atom.getMapIdx());
            }
        }

        for(int i=0; i<products.getAtomContainerCount(); i++){
            IAtomContainer mol = products.getAtomContainer(i);
            for(IAtom atom: mol.atoms()){
                if(!reactants_atom_tags.contains(atom.getMapIdx())){
                    changed_atom_tags.add(atom.getMapIdx());
                }
            }
        }
    }

    private static Map<Integer, Double> getHashedFingerPrint(double[] whfp_fc_reactant, double[] whfp_fc_product, double[] whfp_ch_reactant, double[] whfp_ch_product) {

        Map<Integer, Double> countLoc = new HashMap<Integer, Double>();

        BitSet binary = new BitSet(1024 * 4);
        for (int i = 0; i < whfp_fc_reactant.length; i++) {
            if (whfp_fc_reactant[i] > 0.) {
                countLoc.put(i, whfp_fc_reactant[i]);
                binary.set(i, true);
            } else {
                binary.set(i, false);
            }
        }

        for (int i = 1024; i < 1024 + whfp_fc_product.length; i++) {
            if (whfp_fc_product[i-1024] > 0.) {
                countLoc.put(i, whfp_fc_product[i-1024]);
                binary.set(i, true);
            } else {
                binary.set(i, false);
            }
        }

        for (int i = 2048; i < 2048 + whfp_ch_reactant.length; i++) {
            if (whfp_ch_reactant[i-2048] > 0.) {
                countLoc.put(i, whfp_ch_reactant[i-2048]);
                binary.set(i, true);
            } else {
                binary.set(i, false);
            }
        }

        for (int i = 3072; i < 3072 + whfp_ch_product.length; i++){
            if (whfp_ch_product[i-3072] > 0.) {
                countLoc.put(i, whfp_ch_product[i-3072]);
                binary.set(i, true);
            } else {
                binary.set(i, false);
            }
        }

        return countLoc;
    }

}
