import org.openscience.cdk.*;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.tools.BondEnergies;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IFeature;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.mapping.CallableAtomMappingTool;
import uk.ac.ebi.reactionblast.mapping.Reactor;
import uk.ac.ebi.reactionblast.mapping.interfaces.IMappingAlgorithm;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;
import uk.ac.ebi.reactionblast.tools.StandardizeReaction;

import java.util.HashMap;
import java.util.Map;

import static java.lang.Math.abs;
import static java.util.logging.Level.SEVERE;
import static org.openscience.cdk.CDKConstants.MAPPED;
import static org.openscience.cdk.interfaces.IBond.Order.*;
import static org.openscience.cdk.interfaces.IBond.Order.QUADRUPLE;
import static org.openscience.smsd.tools.BondEnergies.getInstance;

public class rdt {

    public static void main(String[] args) throws CloneNotSupportedException, CDKException, AssertionError, Exception {
        final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        HashMap<Integer, IAtom> reactants_map = new HashMap<Integer, IAtom>();
        HashMap<Integer, IAtom> products_map = new HashMap<Integer, IAtom>();
        int map_id;
        int reactants_MapCount = 0;
        int products_MapCount = 0;

        String reactionSM = "CCCC(=N)C1=C(CC(CC1=O)C2CCCS(=O)C2)O>>CCCC(=N)C1=C(CC(CC1=O)C2CCCS(=O)(=O)C2)O";
        String reactionName = "Test";

        String mappedReaction = "[CH3:1][CH2:2][CH2:3][C:4](=[NH:5])[C:6]1=[C:7]([OH:8])[CH2:9][CH:10]([CH:11]2[CH2:12][CH2:13][CH2:14][S:15](=[O:16])[CH2:18]2)[CH2:19][C:20]1=[O:17]>>[CH3:1][CH2:2][CH2:3][C:4](=[NH:5])[C:6]1=[C:7]([OH:8])[CH2:9][CH:10]([CH:11]2[CH2:12][CH2:13][CH2:14][S:15](=[O:16])(=[O:17])[CH2:18]2)[CH2:19][C:20]1=[O:21]";
        IReaction cdkReaction = smilesParser.parseReactionSmiles(reactionSM);
        IAtomContainerSet reactants = cdkReaction.getReactants();
        IAtomContainerSet products = cdkReaction.getProducts();

        IReaction re_reaction = new Reaction();
        IAtomContainerSet re_reactants = new AtomContainerSet();
        IAtomContainerSet re_products = new AtomContainerSet();

        for(IAtomContainer mol: reactants.atomContainers()){
            IAtomContainer re_mol = new AtomContainer();
            for(IAtom atom: mol.atoms()){
                IAtom re_atom= atom.clone();
                re_mol.addAtom(re_atom);
            }
            re_reactants.addAtomContainer(re_mol);
            re_reaction.addReactant(re_mol);
        }

        for(IAtomContainer mol: products.atomContainers()){
            IAtomContainer re_mol = new AtomContainer();
            for(IAtom atom: mol.atoms()){
                IAtom re_atom= atom.clone();
                re_mol.addAtom(re_atom);
            }
            re_products.addAtomContainer(re_mol);
            re_reaction.addProduct(re_mol);
        }

        for(IAtomContainer mol: re_reactants.atomContainers()){
            //System.out.println(mol.toString());
            for(IAtom atom: mol.atoms()){
                map_id = atom.getMapIdx();
                if(map_id > reactants_MapCount) reactants_MapCount = map_id;
                atom.setID(String.valueOf(atom.getMapIdx()));
                reactants_map.put(map_id, atom);
                //System.out.println("Symbol: " + atom.getSymbol() + ";Index: " + atom.getIndex() + ";MapIndex: " + atom.getMapIdx() + ";ID: " + atom.getID() + ";Bound count: " + atom.getBondCount());
            }
        }

        System.out.println("-----------------------------------------------");

        for(IAtomContainer mol: re_products.atomContainers()){
            for(IAtom atom: mol.atoms()){
                map_id = atom.getMapIdx();
                if(map_id > products_MapCount) products_MapCount = map_id;
                atom.setID(String.valueOf(atom.getMapIdx()));
                products_map.put(map_id, atom);
                //System.out.println("Symbol: " + atom.getSymbol() + ";Index: " + atom.getIndex() + ";MapIndex: " + atom.getMapIdx() + ";ID: " + atom.getID() + ";Bound count: " + atom.getBondCount());
            }
        }

        for(int i=1; i<=Math.min(reactants_MapCount, products_MapCount); i++){
            IMapping new_map = new Mapping(reactants_map.get(i), products_map.get(i));
            re_reaction.addMapping(new_map);
        }


        IReaction performAtomAtomMapping = performAtomAtomMapping(cdkReaction, reactionName);
        System.out.println("AAM sm: " + sg.create(performAtomAtomMapping));

        int res = compareAAM(performAtomAtomMapping, re_reaction);


        return;

    }

    /**
     *
     * @param cdkReaction
     * @param reactionName
     * @return
     * @throws InvalidSmilesException
     * @throws AssertionError
     * @throws Exception
     */
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

    public static int compareAAM(IReaction mapping1, IReaction mapping2) throws Exception {

        BondChangeCalculator bcc1, bcc2;
        int fragmentDeltaChanges1, fragmentDeltaChanges2;
        boolean generate2D = true;//2D perception of the stereo centers
        boolean generate3D = false;//2D perception of the stereo centers

        bcc1 = new BondChangeCalculator(mapping1);
        bcc2 = new BondChangeCalculator(mapping2);

        bcc1.computeBondChanges(generate2D, generate3D);
        bcc2.computeBondChanges(generate2D, generate3D);

        fragmentDeltaChanges1 = bcc1.getTotalFragmentCount();
        fragmentDeltaChanges2 = bcc2.getTotalFragmentCount();

        int bondChange1 = (int) getTotalBondChange(bcc1.getFormedCleavedWFingerprint());
        int bondChange2 = (int) getTotalBondChange(bcc2.getFormedCleavedWFingerprint());

        bondChange1 += getTotalBondChange(bcc1.getOrderChangesWFingerprint());
        bondChange2 += getTotalBondChange(bcc2.getOrderChangesWFingerprint());

        int stereoChanges1 = (int) getTotalBondChange(bcc1.getStereoChangesWFingerprint());
        int stereoChanges2 = (int) getTotalBondChange(bcc2.getStereoChangesWFingerprint());

        boolean skipHydrogenRealtedBondChanges = true;

        int bondBreakingEnergy1 = getTotalBondChangeEnergy(bcc1.getFormedCleavedWFingerprint(), skipHydrogenRealtedBondChanges);
        int bondBreakingEnergy2 = getTotalBondChangeEnergy(bcc2.getFormedCleavedWFingerprint(), skipHydrogenRealtedBondChanges);

        int totalSmallestFragmentCount1 = bcc1.getTotalSmallestFragmentSize();
        int totalSmallestFragmentCount2 = bcc2.getTotalSmallestFragmentSize();

        int totalCarbonBondChanges1 = getTotalCarbonBondChange(bcc1.getFormedCleavedWFingerprint());
        int totalCarbonBondChanges2 = getTotalCarbonBondChange(bcc2.getFormedCleavedWFingerprint());

        int localScore1 = bondChange1 + fragmentDeltaChanges1;
        int localScore2 = bondChange2 + fragmentDeltaChanges2;

        bcc1.getReaction().setFlag(MAPPED, true);
        bcc2.getReaction().setFlag(MAPPED, true);

        // Compare

        if (bondBreakingEnergy2 == 0.
                && fragmentDeltaChanges2 == 0
                && bondChange2 == 0
                && stereoChanges1 >= stereoChanges2){
            // mapping2 is better
            System.out.println("Less stereo change of mapping 2.");
            return 1;
        }else if(bondChange1 > bondChange2
                && totalCarbonBondChanges1 > 0
                && totalCarbonBondChanges1 > totalCarbonBondChanges2
                && (fragmentDeltaChanges1 > fragmentDeltaChanges2
                || bondBreakingEnergy1 > bondBreakingEnergy2)){
            System.out.println("Less bond change (total/carbon) or fragments change or less bond breaking energy of mapping 2.");
            return 1;
        }else if(totalCarbonBondChanges1 > totalCarbonBondChanges2
                && fragmentDeltaChanges1 > 0
                && fragmentDeltaChanges2 > 0){
            System.out.println("Less carbon bond change.");
            return 1;
        }else if(fragmentDeltaChanges1 >= fragmentDeltaChanges2
                && totalSmallestFragmentCount1 >= totalSmallestFragmentCount2
                && bondBreakingEnergy1 > bondBreakingEnergy2
                && totalCarbonBondChanges1 >= totalCarbonBondChanges2){
            System.out.println("Less fragment change and bond breaking energy and total carbon bond change of mapping 2.");
            return 1;
        }else if(fragmentDeltaChanges1 > fragmentDeltaChanges2
                && totalSmallestFragmentCount1 > totalSmallestFragmentCount2){
            System.out.println("Less fragment change and less smallest fragment counts of mapping 2.");
            return 1;
        }else if(fragmentDeltaChanges1 == fragmentDeltaChanges2
                && totalSmallestFragmentCount1 == totalSmallestFragmentCount2
                && bondBreakingEnergy1 > bondBreakingEnergy2
                && totalCarbonBondChanges1 >= totalCarbonBondChanges2){
            System.out.println("Less bond breaking energy and less carbon bond change of mapping 2.");
            return 1;
        }else if(fragmentDeltaChanges1 > fragmentDeltaChanges2
                && bondBreakingEnergy1 > bondBreakingEnergy2){
            System.out.println("Less fragment change and less bond breaking energy of mapping 2.");
            return 1;
        }else if(totalCarbonBondChanges1 == totalCarbonBondChanges2
                && fragmentDeltaChanges1 > fragmentDeltaChanges2){
            System.out.println("Less fragment change of mapping 2.");
            return 1;
        }else if(fragmentDeltaChanges1 == fragmentDeltaChanges2
                && bondBreakingEnergy1 == bondBreakingEnergy2
                && totalCarbonBondChanges1 > totalCarbonBondChanges2){
            System.out.println("Less total carbon bond change of mapping 2.");
            return 1;
        }else if(bondBreakingEnergy1 == bondBreakingEnergy2
                && bondChange1 == bondChange2
                && stereoChanges1 > stereoChanges2){
            System.out.println("Less stereo change of mapping 2.");
            return 1;
        }else if(bondBreakingEnergy1 > bondBreakingEnergy2
                && totalCarbonBondChanges1 > totalCarbonBondChanges2){
            System.out.println("Less bond breaking energy and less total carbon bond change of mapping 2.");
            return 1;
        }else if(bondChange1 < bondChange2
                && bondBreakingEnergy1 < bondBreakingEnergy2
                && totalCarbonBondChanges1 > 0
                && totalCarbonBondChanges1 >totalCarbonBondChanges2
                && totalSmallestFragmentCount1 > totalSmallestFragmentCount2){
            System.out.println("More bond changes and greater bond breaking energy but less carbon bond changes and less smallest fragments of mapping 2.");
            return 1;
        }else if(bondChange1 > bondChange2
                && totalCarbonBondChanges1 > totalCarbonBondChanges2
                && totalSmallestFragmentCount1 > totalSmallestFragmentCount2){
            System.out.println("Less bond changes and less carbon bond changes and less smallest fragments of mapping 2.");
            return 1;
        }else if(bondChange1 == bondChange2
                && totalCarbonBondChanges1 == totalCarbonBondChanges2
                && bondBreakingEnergy1 > bondBreakingEnergy2){
            System.out.println("Less bond breaking energy of mapping 2.");
            return 1;
        }

//        System.out.println("Total stereo changes: " + stereoChanges1);
//        System.out.println("Bond breaking energy: " + bondBreakingEnergy1);

        return 0;
    }

    public static double getTotalBondChange(IPatternFingerprinter fingerprint) throws CDKException {
        double total = 0;
        total = fingerprint.getFeatures().stream().map((key) -> key.getWeight()).filter((val) -> (val > 0.)).map((val) -> val).reduce(total, (accumulator, _item) -> accumulator + _item); //&& !key.contains("H")
        return total;
    }

    public static int getTotalBondChangeEnergy(IPatternFingerprinter fingerprint, boolean skipHydrogen) {
        int total = 0;
        try {
            BondEnergies be = getInstance();
            for (IFeature feature : fingerprint.getFeatures()) {
                double val = feature.getWeight();
                String key = feature.getPattern();
                if (val > 0) {
//                    System.out.println("BOND BROKEN/FORMED: " + key + " : " + val);
                    if (key.contains("-") || key.contains("%") || key.contains("@")) {
                        String[] temp = null;
                        if (key.contains("-")) {
                            temp = key.split("-");
                        } else if (key.contains("%")) {
                            temp = key.split("%");
                        } else if (key.contains("@")) {
                            temp = key.split("@");
                        }
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        int energy = be.getEnergies(temp[0], temp[1], SINGLE);
                        if (energy > 0) {
                            /*
                             * Ring energy correction factor example:R01081
                             */
                            if (key.contains("%")) {
                                total += val * (energy - 5.0);

                            } else {
                                total += val * energy;
                            }
                        }
                    } else if (key.contains("=")) {
                        String[] temp = key.split("=");
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        int energy = be.getEnergies(temp[0], temp[1], DOUBLE);
                        if (energy > 0) {
                            total += val * energy;
                        }
                    } else if (key.contains("#")) {
                        String[] temp = key.split("#");
                        int energy = be.getEnergies(temp[0], temp[1], TRIPLE);
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        if (energy > 0) {
                            total += val * energy;
                        }
                    } else if (key.contains("$")) {
                        String[] temp = key.split("$");
                        if (skipHydrogen && (temp[0].equals("H") || temp[1].equals("H"))) {
                            continue;
                        }
                        int energy = be.getEnergies(temp[0], temp[1], QUADRUPLE);
                        if (energy > 0) {
                            total += val * energy;
                        }
                    }
                }
            }
        } catch (CDKException ex) {
            System.out.println(ex);
        }
        return abs(total);
    }

    public static int getTotalCarbonBondChange(IPatternFingerprinter fingerprint) throws CDKException {
        double total = 0;
        total = fingerprint.getFeatures().stream().filter((key) -> (key.getPattern().contains("C-C")
                || key.getPattern().contains("C=C")
                || key.getPattern().contains("C#C")
                || key.getPattern().contains("C%C")
                || key.getPattern().contains("C@C"))).map((key) -> key.getWeight()).filter((val) -> (val > 0.)).map((val) -> val).reduce(total, (accumulator, _item) -> accumulator + _item); //&& !key.contains("H")
        return (int) total;
    }

}
