package eawag.envirule;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;

import java.io.*;
import java.util.*;

public class rxn_clusterer {

    public static final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

    public static void clustering(String rxn_file, String out_dir) throws Exception {

        File directory = new File(out_dir);
        if(!directory.exists()){
            directory.mkdir();
        }

        Map<String, Map<Integer, Double>> reactionFPs = new HashMap<String, Map<Integer, Double>>();
        List<String> reactions = parseReactions(rxn_file, smilesParser, reactionFPs);

        // clusters: Representative reaction of a group -> all the reactions in this group
        Map<String, Set<String>> clusters = getClusters(reactions, smilesParser, reactionFPs, out_dir, true);
        // map_file: Representative reaction of a group -> file name
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

        // Serialize the mapping results
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

            // Remove "\n" of each line
            line_array = line.split("\n");
            reaction = line_array[0];

            // Skip reactions with only CO2 as product.
            if(reaction.split(">>").length < 2) continue;
            String product = reaction.split(">>")[1];
            if (product.compareTo("C(=O)=O") == 0) continue;

            reactions.add(reaction);

            IReaction Reaction = smilesParser.parseReactionSmiles(reaction);
            IReaction performAtomAtomMapping = performAtomAtomMapping(Reaction, "null");
            BondChangeCalculator bcc1 = new BondChangeCalculator(performAtomAtomMapping);
            bcc1.computeBondChanges(true, false);

            if(!reactionFPs.containsKey(reaction)){
                reactionFPs.put(reaction, getChangedBonds(performAtomAtomMapping, bcc1));
            }
        }
        return reactions;;
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
        return null;
    }

    private static Map<Integer, Double> getChangedBonds(IReaction performAtomAtomMapping, BondChangeCalculator bcc) throws CDKException {

        IAtomContainerSet reactants = performAtomAtomMapping.getReactants();
        IAtomContainerSet products = performAtomAtomMapping.getProducts();

        Set<IBond> formed_cleaved_bonds_reactants = new HashSet<IBond>();
        Set<IBond> formed_cleaved_bonds_products = new HashSet<IBond>();
        Set<IBond> changed_bonds_reactants = new HashSet<IBond>();
        Set<IBond> changed_bonds_products = new HashSet<IBond>();

        // Initialize fingerprints
        IPatternFingerprinter formedCleavedWFingerprint_reactants = new PatternFingerprinter();
        formedCleavedWFingerprint_reactants.setFingerprintID(performAtomAtomMapping.getID() + ":" + "Bond Cleaved and Formed");
        IPatternFingerprinter formedCleavedWFingerprint_products = new PatternFingerprinter();
        formedCleavedWFingerprint_products.setFingerprintID(performAtomAtomMapping.getID() + ":" + "Bond Cleaved and Formed");
        IPatternFingerprinter changedWFingerprint_reactants = new PatternFingerprinter();;
        changedWFingerprint_reactants.setFingerprintID(performAtomAtomMapping.getID() + ":" + "Bond Change");
        IPatternFingerprinter changedWFingerprint_products = new PatternFingerprinter();;
        changedWFingerprint_products.setFingerprintID(performAtomAtomMapping.getID() + ":" + "Bond Change");

        // Collect MapId of changed atoms
        Set<Integer> changed_atom_tags = new HashSet<>();
        for(IAtom atom: bcc.getReactionCenterSet()){
            changed_atom_tags.add(atom.getMapIdx());
        }

        // Include functional groups in reactants
        addFunctionalGroups(reactants, changed_atom_tags);

        // Include all the atoms that only show up in products
        newAddedGroupsInProducts(performAtomAtomMapping, changed_atom_tags);




        return null;

    }

    private static void addFunctionalGroups(IAtomContainerSet reactants, Set<Integer> changed_atom_tags){

        Set<Integer> addedFunctionalGroups = new HashSet<>();

        // The list of functional groups below is adapted from
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


}
