package uk.ac.ebi;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.templates.MoleculeFactory;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import uk.ac.ebi.reactionblast.fingerprints.Feature;
import uk.ac.ebi.reactionblast.fingerprints.PatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.interfaces.IPatternFingerprinter;
import uk.ac.ebi.reactionblast.fingerprints.tools.Similarity;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_BOND_CHANGE_FLAGS;
import uk.ac.ebi.reactionblast.mechanism.interfaces.ECBLAST_FLAGS;

import java.util.*;

import static uk.ac.ebi.map_test.performAtomAtomMapping;
import static uk.ac.ebi.reactionblast.mechanism.helper.Utility.getCanonicalisedBondChangePattern;

public class sim_test {

    public static void main(String[] args) throws Exception {
        //String reaction1 = "CCC=O>>CCC(=O)O";
        //String reaction1 = "CCCO>>CC=C";
        //String reaction1 = "CC=C>>CCCO";
        //String reaction2 = "CC=O>>CC(=O)O"；
        //String reaction2 = "CCC=C>>CCCCO";

//        String reaction1 = "CCCO>>CCC=O";
//        String reaction2 = "CC(C)O>>CC(C)=O";

        String reaction1 = "CC1=C(N=C(NC)N=C1C)OC(=O)N(C)C>>CC1=C(N=C(NC)N=C1C)O";
        String reaction2 = "CC1=C(N=C(N=C1C)N(C)C)OC(=O)N(C)C>>CC1=C(N=C(N=C1C)N(C)C)O";

//        String reaction1 = "CC(C)(C1=CC=C(C=C1)C(=O)OC)C(=O)O>>CC(C)(C1=CC=C(C=C1)C(=O)OC)C(=O)OCCC";
//        String reaction2 = "C1CCC(CC1)N2C(=O)C3=C(CCC3O)NC2=O>>C1CCC(CC1)N2C(=O)C3=C(CCC3OS(=O)(=O)O)NC2=O";

        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
        IReaction Reaction1 = smilesParser.parseReactionSmiles(reaction1);
        IReaction Reaction2 = smilesParser.parseReactionSmiles(reaction2);

        IReaction performAtomAtomMapping1 = performAtomAtomMapping(Reaction1, "null");
        IReaction performAtomAtomMapping2 = performAtomAtomMapping(Reaction2, "null");

        System.out.println("AAM sm1: " + sg.create(performAtomAtomMapping1));
        System.out.println("AAM sm2: " + sg.create(performAtomAtomMapping2));

        BondChangeCalculator bcc1 = new BondChangeCalculator(performAtomAtomMapping1);
        BondChangeCalculator bcc2 = new BondChangeCalculator(performAtomAtomMapping2);

        bcc1.computeBondChanges(true, false);
        bcc2.computeBondChanges(true, false);



        /////////////////////////////////////////////////////////////////
        //
        Map<Integer, Double> result1 = getChangedBonds(performAtomAtomMapping1, bcc1);
        Map<Integer, Double> result2 = getChangedBonds(performAtomAtomMapping2, bcc2);
        //float sim_formCleaved = Similarity.getTanimotoSimilarity(result1, result2);
        float sim_formCleaved = compareDict(result1, result2);
        System.out.println("Similarity: " + sim_formCleaved + " dict1: " + result1 + " dict2: " + result2);
        /////////////////////////////////////////////////////////////////

        IAtomContainer center1 = getReactionCenterMol(performAtomAtomMapping1, bcc1);
        IAtomContainer center2 = getReactionCenterMol(performAtomAtomMapping2, bcc2);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(center1);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(center2);
        //CDKHueckelAromaticityDetector.detectAromaticity(center1);
        CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(center1.getBuilder());
        hAdder.addImplicitHydrogens(center1);
        hAdder.addImplicitHydrogens(center2);

        SmilesGenerator generator = SmilesGenerator.generic();
        String smiles1 = generator.create(center1);
        String smiles2 = generator.create(center2);
        System.out.println(smiles1);
        System.out.println(smiles2);

//        MACCSFingerprinter fingerprinter = new MACCSFingerprinter();
//        BitSet bitset = fingerprinter.getBitFingerprint(center1).asBitSet();
//        System.out.println(bitset);

        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer mol1 = parser.parseSmiles(smiles1);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol1);
        MACCSFingerprinter fingerprinter = new MACCSFingerprinter();
        BitSet bitset1 = fingerprinter.getBitFingerprint(mol1).asBitSet();
        System.out.println(bitset1);

        IAtomContainer mol2 = parser.parseSmiles(smiles2);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2);
        BitSet bitset2 = fingerprinter.getBitFingerprint(mol2).asBitSet();
        System.out.println(bitset2);


//        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(center1);
//        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(center2);


//        IAtomContainer random = MoleculeFactory.makePhenylAmine();
//        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(random);
//        CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(center1.getBuilder());
//        hAdder.addImplicitHydrogens(center1);
//
//        SmilesGenerator generator = SmilesGenerator.generic();
//        System.out.println(generator.createSMILES(center1));


        // Construct mol from RC
        // Calculate FP of RC mol
        // Calculate FP tanomito sim
        // Train/evaluate model multigen



//        SmilesGenerator smiles_gen = SmilesGenerator.generic();
//        System.out.println(smiles_gen.create(centermol1));
//        System.out.println(smiles_gen.create(centermol2));

//        MACCSFingerprinter MACCS_FP = new MACCSFingerprinter();
//        IBitFingerprint bitFingerprint1 = MACCS_FP.getBitFingerprint(centermol1);
//        IBitFingerprint bitFingerprint2 = MACCS_FP.getBitFingerprint(centermol2);
//
//        System.out.println(bitFingerprint1);






//        for(IAtom center: bcc1.getReactionCenterSet()){
//            if(performAtomAtomMapping1.getReactants().getAtomContainer(0).contains(center)){
//                System.out.println("Get the connected bonds count: " + performAtomAtomMapping1.getReactants().getAtomContainer(0).getConnectedBondsList(center).size());
//            }
//            if(!mapId1.contains(center.getMapIdx())) {
//                System.out.println("1: " + center.getSymbol() + " " + center.getMapIdx());
//                mapId1.add(center.getMapIdx());
//                centermol1.addAtom(center);
//
//            }else continue;
//        }
//
//        Set<Integer> mapId2 = new HashSet<Integer>();
//        IAtomContainer centermol2 = new AtomContainer();
//
//        for(IAtom center: bcc2.getReactionCenterSet()){
//            if(!mapId2.contains(center.getMapIdx())) {
//                System.out.println("2: " + center.getSymbol() + " " + center.getMapIdx());
//                mapId2.add(center.getMapIdx());
//                centermol2.addAtom(center);
//            }else continue;
//        }


//        BitSet Molecule1_1 = bcc1.getReactionCenterWFingerprint().getHashedFingerPrint();
//        BitSet Molecule2_1 = bcc2.getReactionCenterWFingerprint().getHashedFingerPrint();

//        float sim = Similarity.getTanimotoSimilarity(Molecule1_1, Molecule2_1);
//        System.out.println("Reaction Center Compare result is: " + sim);
//
//        BitSet Molecule1_2 = bcc1.getFormedCleavedWFingerprint().getHashedFingerPrint();
//        BitSet Molecule2_2 = bcc2.getFormedCleavedWFingerprint().getHashedFingerPrint();
//
//        sim = Similarity.getTanimotoSimilarity(Molecule1_2, Molecule2_2);
//        System.out.println("Formed Cleaved Compare result is: " + sim);
//
//        BitSet Molecule1_3 = bcc1.getStereoChangesWFingerprint().getHashedFingerPrint();
//        BitSet Molecule2_3 = bcc2.getStereoChangesWFingerprint().getHashedFingerPrint();
//
//        sim = Similarity.getTanimotoSimilarity(Molecule1_3, Molecule2_3);
//        System.out.println("Stereo Changes Compare result is: " + sim);
//
//        BitSet Molecule1_4 = bcc1.getOrderChangesWFingerprint().getHashedFingerPrint();
//        BitSet Molecule2_4 = bcc2.getOrderChangesWFingerprint().getHashedFingerPrint();
//
//        sim = Similarity.getTanimotoSimilarity(Molecule1_4, Molecule2_4);
//        System.out.println("Order Changes Compare result is: " + sim);
//
//        System.out.println(bcc1.getReactionCenterFormedCleavedFingerprint().size());
//        System.out.println(bcc2.getReactionCenterFormedCleavedFingerprint().size());
//
//        System.out.println("Sim in RC FCfp:");
//
//        for(Map.Entry<Integer, IPatternFingerprinter> entry: bcc1.getReactionCenterFormedCleavedFingerprint().entrySet()){
//            Integer key = entry.getKey();
//            IPatternFingerprinter fp1 = entry.getValue();
//            IPatternFingerprinter fp2 = bcc2.getReactionCenterFormedCleavedFingerprint().get(key);
//
//            BitSet bs1 = fp1.getHashedFingerPrint();
//            BitSet bs2 = fp2.getHashedFingerPrint();
//
//            sim = Similarity.getTanimotoSimilarity(bs1, bs2);
//
//            System.out.println("Sim: " + sim);
//
//        }
//
//        Set<Integer> mapId = new HashSet<Integer>();
//
//        for(IAtom center: bcc1.getReactionCenterSet()){
//            if(!mapId.contains(center.getMapIdx())) {
//                System.out.println("1: " + center.getSymbol() + " " + center.getMapIdx());
//                mapId.add(center.getMapIdx());
//            }else continue;
//        }
//
//        mapId = new HashSet<Integer>();
//        for(IAtom center: bcc2.getReactionCenterSet()){
//            if(!mapId.contains(center.getMapIdx())) {
//                System.out.println("2: " + center.getSymbol() + " " + center.getMapIdx());
//                mapId.add(center.getMapIdx());
//            }else continue;
//        }



    }

    /**
     * For all the atoms added to products, they need to be considered as reaction center as well
     * @param performAtomAtomMapping
     * @param changed_atom_tags
     */
    public static void newAddedGroupsInProducts(IReaction performAtomAtomMapping, Set<Integer> changed_atom_tags){

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

    public static void addFunctionalGroups(IAtomContainerSet reactants, Set<Integer> changed_atom_tags){

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

    public static Map<Integer, Double> getChangedBonds(IReaction performAtomAtomMapping1, BondChangeCalculator bcc1) throws CDKException {

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

        // Identify atoms of reaction center in reactants
//        Map<Integer, Double> reactants_rc_atoms = new HashMap<Integer, Double>();
//        for(int i=0; i<reactants.getAtomContainerCount(); i++){
//            IAtomContainer reactant = reactants.getAtomContainer(i);
//
//            for(IAtom atom: bcc1.getReactionCenterSet()){
////                System.out.println(atom.getImplicitHydrogenCount() + " " + atom.getSymbol());
//                if(reactant.contains(atom)){
//                    String symbol = atom.getSymbol();
//                    if(symbol == "C" && atom.isAromatic()){
//                        symbol = "c";
//                    }else if(symbol == "N" && atom.isAromatic()){
//                        symbol = "n";
//                    }
//                    int hash = -symbol.hashCode();
//                    if(!reactants_rc_atoms.containsKey(hash)) reactants_rc_atoms.put(hash, 0.0);
//                    reactants_rc_atoms.put(hash, reactants_rc_atoms.get(hash) + 1);
//                }
//            }
//        }
//
//        // Identify atoms of reaction center in products
//        Map<Integer, Double> products_rc_atoms = new HashMap<Integer, Double>();
//        for(int i=0; i<products.getAtomContainerCount(); i++){
//            IAtomContainer product = products.getAtomContainer(i);
//
//            for(IAtom atom: bcc1.getReactionCenterSet()){
//                if(product.contains(atom)){
//                    String symbol = atom.getSymbol();
//                    if(symbol == "C" && atom.isAromatic()){
//                        symbol = "c";
//                    }else if(symbol == "N" && atom.isAromatic()){
//                        symbol = "n";
//                    }
//                    int hash = -symbol.hashCode() - 200;
//                    if(!products_rc_atoms.containsKey(hash)) products_rc_atoms.put(hash, 0.0);
//                    products_rc_atoms.put(hash, products_rc_atoms.get(hash) + 1);
//                }
//            }
//        }

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

    public static Map<Integer, Double> getHashedFingerPrint(double[] whfp_fc_reactant, double[] whfp_fc_product, double[] whfp_ch_reactant, double[] whfp_ch_product) {

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
                    System.out.println("1: " + center.getSymbol() + " " + center.getMapIdx() );
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

        System.out.println(mapId1);
        System.out.println("Reaction center atom counts: " + centermol1.getAtomCount() + "; bond counts: " + centermol1.getBondCount());


        return centermol1;
    }

    public static float compareDict(Map<Integer, Double> map1, Map<Integer, Double> map2){
        if (map1.size() != map2.size()) return 0.0F;

        for(Map.Entry<Integer, Double> entry: map1.entrySet()){
            Integer key = entry.getKey();
            if(!map2.containsKey(key)) return 0.0F;
            if(entry.getValue().doubleValue() != map2.get(key)) return 0.0F;
        }
        return 1.0F;
    }
}

class BondsChange{

    BitSet formedCleavedWFingerprint;
    BitSet changedWFingerprint;

    public BondsChange(BitSet formedCleavedWFingerprint, BitSet changedWFingerprint){
        this.formedCleavedWFingerprint = formedCleavedWFingerprint;
        this.changedWFingerprint = changedWFingerprint;
    }

    public BitSet getformedCleaved(){
        return formedCleavedWFingerprint;
    }

    public BitSet getchanged(){
        return changedWFingerprint;
    }
}