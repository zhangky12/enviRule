package uk.ac.ebi;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.tools.Similarity;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS;

import java.io.*;
import java.util.*;

import static uk.ac.ebi.reactionblast.mechanism.helper.Utility.getCanonicalisedBondChangePattern;

// For now, only use rdt mapping.

public class clustering_sim_test {

    public static void main(String[] args) throws Exception{

//        final PrintStream err = new PrintStream(System.err);
//        System.setErr(new PrintStream("/dev/null"));

        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//        String reaction1 = "CCCO>>CCC=O";
//        String reaction2 = "CC(C)O>>CC(C)=O";
//        String reaction3 = "CCC=O>>CCC(=O)OH";

//        float sim = getSimScore(reaction1, reaction2, smilesParser);
//
//        System.out.println("Similarity of reaction cores is: " + sim);

//        List<String> reactions = new ArrayList<String>();
//        reactions.add(reaction1);
//        reactions.add(reaction2);
//        reactions.add(reaction3);
//
//        System.out.println(getCenters(reactions, smilesParser));

        Map<String, Map<Integer, Double>> reactionFPs = new HashMap<String, Map<Integer, Double>>();
        List<String> reactions = parseReactions2("/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/BBD_reactions/all_reactions.txt", smilesParser, reactionFPs);
//        List<String> reactions = parseReactions2("/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/MultistepTrn_SB_combined_cleaned_010422/missed_reactions.txt", smilesParser, reactionFPs);
//        Map<String, Set<String>> clusters = getClusters(reactions, smilesParser, reactionFPs);
        Map<String, Set<String>> clusters = getClusters3(reactions, smilesParser, reactionFPs);

//        System.setErr(err);

        FileWriter fw = new FileWriter("/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/BBD_reactions/clustered_missed_reactions.txt");
//        FileWriter fw = new FileWriter("/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/MultistepTrn_SB_combined_cleaned_010422/clustered_missed_reactions.txt");
        for(Map.Entry<String, Set<String>> entry: clusters.entrySet()){
            Set<String> members = entry.getValue();
            for(String react: members){
                fw.write(react + "\n");
            }
            fw.write("\n");
        }
        fw.close();
        System.exit(0);

    }

    public static List<String> parseReactions2(String file, SmilesParser smilesParser, Map<String, Map<Integer, Double>> reactionFPs) throws Exception{

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
                reactionFPs.put(reaction, sim_test.getChangedBonds(performAtomAtomMapping1, bcc1));
            }

        }
        return reactions;
    }

    public static List<String> parseReactions(String file, SmilesParser smilesParser, Map<String, BitSet> reactionFPs) throws Exception {

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
            reactions.add(reaction);

            if(!reactionFPs.containsKey(reaction)){
                reactionFPs.put(reaction, getReactionCenterFPs(reaction, smilesParser));
            }

        }

        return reactions;
    }

    public static BitSet getReactionCenterFPs(String reaction1, SmilesParser smilesParser) throws Exception {

        IReaction Reaction1 = smilesParser.parseReactionSmiles(reaction1);
        IReaction performAtomAtomMapping1 = performAtomAtomMapping(Reaction1, "null");
        BondChangeCalculator bcc1 = new BondChangeCalculator(performAtomAtomMapping1);
        bcc1.computeBondChanges(true, false);
        IAtomContainer center1 = getReactionCenterMol(performAtomAtomMapping1, bcc1);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(center1);

        CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(center1.getBuilder());
        hAdder.addImplicitHydrogens(center1);

        SmilesGenerator generator = SmilesGenerator.generic();
        String smiles1 = generator.create(center1);

        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer mol1 = parser.parseSmiles(smiles1);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol1);
        MACCSFingerprinter fingerprinter = new MACCSFingerprinter();
        BitSet bitset1 = fingerprinter.getBitFingerprint(mol1).asBitSet();

        return bitset1;
    }

    public static Map<String, Set<String>> getClusters(List<String> reactions, SmilesParser smilesParser, Map<String, BitSet> reactionFPs) throws Exception {

        Map<String, Set<String>> centers = new HashMap<String, Set<String>>();
        Set<Integer> toRemove = new HashSet<Integer>();
        float sim;

        for(int i=0; i<reactions.size(); i++){
            if(toRemove.contains(i)) continue;
            String reaction1 = reactions.get(i);
            for(int j=i; j<reactions.size(); j++){

                if(toRemove.contains(j)) continue;
                String reaction2 = reactions.get(j);

                if(i == j) sim = 1.0F;
                //else sim = getSimScore(reaction1, reaction2, smilesParser);
                else sim = getSimScore(reactionFPs.get(reaction1), reactionFPs.get(reaction2));

                if(sim == 1.0F){
                    if(!centers.containsKey(reaction1)) centers.put(reaction1, new HashSet<String>());
                    centers.get(reaction1).add(reaction2);
                    toRemove.add(j);
                }
            }
        }
        return centers;
    }

    public static Map<String, Set<String>> getClusters2(List<String> reactions, SmilesParser smilesParser, Map<String, BondsChange> reactionFPs) throws Exception {

        Map<String, Set<String>> centers = new HashMap<String, Set<String>>();
        Set<Integer> toRemove = new HashSet<Integer>();
        float sim;

        for(int i=0; i<reactions.size(); i++){
            if(toRemove.contains(i)) continue;
            String reaction1 = reactions.get(i);
            for(int j=i; j<reactions.size(); j++){

                if(toRemove.contains(j)) continue;
                String reaction2 = reactions.get(j);

                if(i == j) sim = 1.0F;
                else {
                    float sim_change = Similarity.getTanimotoSimilarity(reactionFPs.get(reaction1).getchanged(), reactionFPs.get(reaction2).getchanged());
                    float sim_formCleaved = Similarity.getTanimotoSimilarity(reactionFPs.get(reaction1).getformedCleaved(), reactionFPs.get(reaction2).getformedCleaved());
                    if(sim_change != sim_change && sim_formCleaved == 1.0F) sim = 1.0F;
                    else if(sim_formCleaved != sim_formCleaved && sim_change == 1.0F) sim = 1.0F;
                    else sim = 0.0F;
                }

                if(sim == 1.0F){
                    if(!centers.containsKey(reaction1)) centers.put(reaction1, new HashSet<String>());
                    centers.get(reaction1).add(reaction2);
                    toRemove.add(j);
                }
            }
        }
        return centers;
    }

    public static Map<String, Set<String>> getClusters3(List<String> reactions, SmilesParser smilesParser, Map<String, Map<Integer, Double>> reactionFPs) throws Exception{

        Map<String, Set<String>> centers = new HashMap<String, Set<String>>();
        Set<Integer> toRemove = new HashSet<Integer>();
        float sim;

        for(int i=0; i<reactions.size(); i++){
            if(toRemove.contains(i)) continue;
            String reaction1 = reactions.get(i);
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


        return centers;
    }

    public static float getSimScore(BitSet bitset1, BitSet bitset2) throws Exception {
        float sim = Similarity.getTanimotoSimilarity(bitset1, bitset2);
        return sim;
    }


    public static float getSimScore(String reaction1, String reaction2, SmilesParser smilesParser) throws Exception {

        float sim;

        IReaction Reaction1 = smilesParser.parseReactionSmiles(reaction1);
        IReaction Reaction2 = smilesParser.parseReactionSmiles(reaction2);

        IReaction performAtomAtomMapping1 = performAtomAtomMapping(Reaction1, "null");
        IReaction performAtomAtomMapping2 = performAtomAtomMapping(Reaction2, "null");

        BondChangeCalculator bcc1 = new BondChangeCalculator(performAtomAtomMapping1);
        BondChangeCalculator bcc2 = new BondChangeCalculator(performAtomAtomMapping2);

        bcc1.computeBondChanges(true, false);
        bcc2.computeBondChanges(true, false);

        IAtomContainer center1 = getReactionCenterMol(performAtomAtomMapping1, bcc1);
        IAtomContainer center2 = getReactionCenterMol(performAtomAtomMapping2, bcc2);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(center1);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(center2);

        CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(center1.getBuilder());
        hAdder.addImplicitHydrogens(center1);
        hAdder.addImplicitHydrogens(center2);

        SmilesGenerator generator = SmilesGenerator.generic();
        String smiles1 = generator.create(center1);
        String smiles2 = generator.create(center2);

        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer mol1 = parser.parseSmiles(smiles1);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol1);
        MACCSFingerprinter fingerprinter = new MACCSFingerprinter();
        BitSet bitset1 = fingerprinter.getBitFingerprint(mol1).asBitSet();

        IAtomContainer mol2 = parser.parseSmiles(smiles2);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2);
        BitSet bitset2 = fingerprinter.getBitFingerprint(mol2).asBitSet();

        sim = Similarity.getTanimotoSimilarity(bitset1, bitset2);

        return sim;
    }

    public static IReaction performAtomAtomMapping(IReaction cdkReaction, String reactionName) throws InvalidSmilesException, AssertionError, Exception {
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

    public static IAtomContainer getReactionCenterMol(IReaction performAtomAtomMapping1, BondChangeCalculator bcc1){

        Set<Integer> mapId1 = new HashSet<Integer>();
        IAtomContainer centermol1 = new AtomContainer();

        // IAtom container can contain several separated parts
        for(IAtom center: bcc1.getReactionCenterSet()){
            // First, find center in reactants.
            boolean found = false;
            for(int i=0; i < performAtomAtomMapping1.getReactants().getAtomContainerCount(); i++){
                IAtomContainer mol = performAtomAtomMapping1.getReactants().getAtomContainer(i);
                for(IAtom mol_atom: mol.atoms()){
                    if(mol_atom.getMapIdx() == center.getMapIdx()) {
                        found=true;
                        // Changes can happen in bondchangeannotator. Same MapId but different atoms.
                        center = mol_atom;
                        break;
                    }
                }
                if(found && !mapId1.contains(center.getMapIdx())){
                    //System.out.println("1: " + center.getSymbol() + " " + center.getMapIdx() );
                    mapId1.add(center.getMapIdx());
                    centermol1.addAtom(center);
                    break;
                }
            }
        }

        // Get reaction center bonds from reactants
        for(int i=0; i < performAtomAtomMapping1.getReactants().getAtomContainerCount(); i++){
            IAtomContainer mol = performAtomAtomMapping1.getReactants().getAtomContainer(i);
            for(IBond bond: mol.bonds()){
                if(mapId1.contains(bond.getBegin().getMapIdx()) && mapId1.contains(bond.getEnd().getMapIdx())){
                    centermol1.addBond(bond);
                }
            }
        }

//        System.out.println(mapId1);
//        System.out.println("Reaction center atom counts: " + centermol1.getAtomCount() + "; bond counts: " + centermol1.getBondCount());


        return centermol1;
    }

    public static BondsChange getChangedBonds(IReaction performAtomAtomMapping1, BondChangeCalculator bcc1) throws CDKException {

        IAtomContainerSet reactants = performAtomAtomMapping1.getReactants();
        IAtomContainerSet products = performAtomAtomMapping1.getProducts();
        Set<IBond> changed_bonds = new HashSet<IBond>();
        Set<IBond> formed_cleaved_bonds = new HashSet<IBond>();

        // Initialize fingerprints
        IPatternFingerprinter formedCleavedWFingerprint = new PatternFingerprinter();
        formedCleavedWFingerprint.setFingerprintID(performAtomAtomMapping1.getID() + ":" + "Bond Cleaved and Formed");
        IPatternFingerprinter changedWFingerprint = new PatternFingerprinter();;
        changedWFingerprint.setFingerprintID(performAtomAtomMapping1.getID() + ":" + "Bond Change");

        // Identify changed or formed/cleaved bonds in reactants
        for(int i=0; i<reactants.getAtomContainerCount(); i++){
            IAtomContainer reactant = reactants.getAtomContainer(i);
            for(int j=0; j<reactant.getBondCount(); j++){
                IBond bond = reactant.getBond(j);

                if(bond.getProperties().containsKey(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION)){
                    if(bond.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) == ECBLAST_BOND_CHANGE_FLAGS.BOND_CLEAVED ||
                            bond.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) == ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED ||
                            bond.getProperties().get(ECBLAST_FLAGS.BOND_CHANGE_INFORMATION) == ECBLAST_BOND_CHANGE_FLAGS.BOND_FORMED_OR_CLEAVED){
                        formed_cleaved_bonds.add(bond);
                    }else{
                        changed_bonds.add(bond);
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
                        formed_cleaved_bonds.add(bond);
                    }else{
                        changed_bonds.add(bond);
                    }
                }
            }
        }

        // Add formed or cleaved bonds to fingerprints
        for(IBond bond: formed_cleaved_bonds){
            formedCleavedWFingerprint.add(new Feature(getCanonicalisedBondChangePattern(bond), 1.0));
        }

        // Add changed bonds to fingerprints
        for(IBond bond: changed_bonds){
            changedWFingerprint.add(new Feature(getCanonicalisedBondChangePattern(bond), 1.0));
        }


        return new BondsChange(formedCleavedWFingerprint.getHashedFingerPrint(), changedWFingerprint.getHashedFingerPrint());

    }

}