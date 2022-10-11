package eawag.envirule.modules;

import ambit2.core.data.MoleculeTools;
import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import ambit2.smarts.SmartsConst;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.graph.ConnectivityChecker;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class reactor {

    // Parse reaction file and rule file, apply all the rule on each reactant.

    public static List<IAtomContainer[]> apply_singleRule(String reaction, String rule_smirks) throws Exception {

        String reactant = reaction.split(">>")[0];

        IAtomContainer Reactant = preprocessMolecule(reactant);
        System.out.println("Reactant: " + toSmiles(Reactant));

        boolean moleculeValid = validateMolecule(Reactant);

        List<IAtomContainer[]> ProductLists;
        SMIRKSReaction rule = fromSmirks(rule_smirks);
        List<IAtomContainer[]> products = apply(Reactant, true, true, rule);

        return products;

    }

    public static IAtomContainerSet applySmirks(final String smrk,
                                                final IAtomContainer
                                                        preprocessedMolecule,
                                                final boolean singlePos) throws Exception {

        SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder
                .getInstance());
        smrkMan.setFlagSSMode(SmartsConst.SSM_MODE.SSM_NON_IDENTICAL_FIRST);
        smrkMan.setFlagProcessResultStructures(true);
        smrkMan.setFlagClearHybridizationBeforeResultProcess(true);
        smrkMan.setFlagClearImplicitHAtomsBeforeResultProcess(true);
        smrkMan.setFlagClearAromaticityBeforeResultProcess(true);
        smrkMan.setFlagAddImplicitHAtomsOnResultProcess(true);
        smrkMan.setFlagConvertAddedImplicitHToExplicitOnResultProcess(false);
        smrkMan.setFlagConvertExplicitHToImplicitOnResultProcess(true);
        smrkMan.getSmartsParser().mSupportDoubleBondAromaticityNotSpecified = false;

        smrkMan.setFlagApplyStereoTransformation(true);

        ambit2.smarts.SMIRKSReaction reaction = smrkMan.parse(smrk);
        if (!smrkMan.getErrors().equals(""))
            throw new Exception("Invalid SMIRKS: " + smrkMan.getErrors());

        IAtomContainerSet set;
        if (singlePos)
            set = smrkMan.applyTransformationWithSingleCopyForEachPos(
                    preprocessedMolecule, null, reaction,
                    SmartsConst.SSM_MODE.SSM_ALL);
        else
            set = smrkMan.applyTransformationWithCombinedOverlappedPos(
                    preprocessedMolecule, null, reaction);

        if (set != null) {
            IAtomContainerSet postProcessedMols = new AtomContainerSet();
            for (IAtomContainer mol : set.atomContainers()) {
                mol = AtomContainerManipulator.suppressHydrogens(mol);
                postProcessedMols.addAtomContainer(mol);
            }
            set = postProcessedMols;
        }
        return set;
    }

    private static List<IAtomContainer[]> apply(IAtomContainer preprocessedMolecule,
                                                boolean splitUnconnectedMolecules,
                                                boolean applyRuleToSinglePosition,
                                                SMIRKSReaction reaction) throws Exception {

        List<IAtomContainer[]> rProducts = new ArrayList<IAtomContainer[]>();
        HashSet<String> rProductSmiles = new HashSet<String>();
        try {
            String smirks = reaction.reactantsSmarts + ">" + reaction.agentsSmarts + ">" + reaction.productsSmarts;

            IAtomContainerSet products = applySmirks(
                    smirks,
                    preprocessedMolecule,
                    applyRuleToSinglePosition);

            if (products != null) {

                for (IAtomContainer p : products.atomContainers()) {
                    if (splitUnconnectedMolecules) {
                        List<IAtomContainer> splitProducts = new ArrayList<IAtomContainer>();
                        for (IAtomContainer m : ConnectivityChecker.partitionIntoMolecules(p).atomContainers()) {
                            String smiles = toSmiles(m);
                            if (!rProductSmiles.contains(smiles)) {
                                rProductSmiles.add(smiles);
                                splitProducts.add(m);
                            }
                        }
                        rProducts.add(splitProducts.toArray(new IAtomContainer[splitProducts.size()]));
                    } else {
                        String smiles = toSmiles(p);
                        if (!rProductSmiles.contains(smiles)) {
                            rProductSmiles.add(smiles);
                            rProducts.add(new IAtomContainer[]{p});
                        }
                    }
                }
            }
        }catch (Exception e){
            throw new Exception("Error while applying Rule to compound " + preprocessedMolecule.toString());
        }

        return rProducts;
    }

    public final SMIRKSReaction setSmirks(final String newSmirks) throws Exception {

        SMIRKSReaction reaction = fromSmirks(newSmirks);
        if (reaction == null) checkValidSmartsException(newSmirks);

        return reaction;
    }

    /**
     * Converts a SMIRKS String to a SMIRKSReaction.
     *
     * @param smirks a SMIRKS String
     * @return a SMIRKSReaction
     * @throws Exception Something went wrong
     */
    public static SMIRKSReaction fromSmirks(final String smirks)
            throws Exception {

        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        final SMIRKSManager ms = new SMIRKSManager(builder);
        SMIRKSReaction react = null;
        try {
            react = ms.parse(smirks);
        } catch (Exception e) {
            throw new Exception("Error parsing SMIRKS '" + smirks +
                    "': "
                    + e.getMessage(), e);
        }
        if (!ms.getErrors().equals("")) {
            throw new Exception("Error parsing SMIRKS '" + smirks
                    + "': " + ms.getErrors());
        }
        return react;
    }

    protected static void checkValidSmartsException(String smarts)
            throws Exception {
        try {
            QueryAtomContainer q = SMARTSParser
                    .parse(smarts, SilentChemObjectBuilder.getInstance());
            if (q == null) {
                throw new NullPointerException();
            }
        } catch (Throwable e) {
            throw new Exception("not a valid smarts string '"
                    + smarts + "', error: " + e.getMessage());
        }
    }

    public static String toSmiles(final IAtomContainer mol)
            throws Exception {
        try {
            return SmilesGenerator.absolute().create(AtomContainerManipulator
                    .copyAndSuppressedHydrogens(mol));
        } catch (final Exception e) {
            throw new Exception("Molecule: " + mol + " " + e);
        }
    }

    public static IAtomContainer fromSmiles(final String smiles)
            throws Exception {
        try {
            final SmilesParser sp = new SmilesParser(
                    SilentChemObjectBuilder.getInstance());
            final IAtomContainer mol = sp.parseSmiles(smiles);
            return mol;
        } catch (final CDKException e) {
            throw new Exception("SMILES: " + smiles);
        }
    }

    public static List<String> parseReactionFile(String file) throws IOException {
        List<String> reactions = new ArrayList<>();

        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        while ((line = br.readLine()) != null) {

            // Remove catalysts
            String[] s = line.split(">>");
//            if(!((!s[0].contains(".") && !s[1].contains(".")) || (s[0].contains(".") && s[0].contains(".")))){
//                continue;
//            }

            if(s[0].contains(".")){
                String[] reactants = s[0].split("\\.");
                String[] products = s[1].split("\\.");
                Set<String> reactants_set = new HashSet<String>();
                for(String reactant: reactants) reactants_set.add(reactant);
                Set<String> products_set = new HashSet<String>();
                for(String product: products) products_set.add(product);
                if(reactants_set.size() != products_set.size()) continue;
                line = "";
                for(String reactant:reactants){
                    if(products_set.contains(reactant)) continue;
                    line += reactant + ".";
                }
                if(line.length() == 0) continue;
                if(line.charAt(line.length()-1) == '.'){
                    line = line.substring(0, line.length()-1);
                }
                line += ">>";
                for(String product:products){
                    if(reactants_set.contains(product)) continue;
                    line += product + ".";
                }
                if(line.charAt(line.length()-1) == '.'){
                    line = line.substring(0, line.length()-1);
                }
            }

            reactions.add(line);
        }

        br.close();
        return reactions;
    }

    public static IAtomContainer preprocessMolecule(final String smiles) throws Exception {

        try {
            IAtomContainer target = fromSmiles(smiles);
            for (IAtom atom : target.atoms())
                if (atom.getFlag(CDKConstants.ISAROMATIC))
                    atom.setFlag(CDKConstants.ISAROMATIC, false);
            for (IBond bond : target.bonds())
                if (bond.getFlag(CDKConstants.ISAROMATIC))
                    bond.setFlag(CDKConstants.ISAROMATIC, false);
            MoleculeTools.convertImplicitToExplicitHydrogens(target);
            Aromaticity aromaticity = new Aromaticity(
                    ElectronDonation.daylight(),
                    Cycles.or(Cycles.all(), Cycles.edgeShort()));
            aromaticity.apply(target);
            return target;
        } catch (CDKException e) {
            throw new Exception("Cannot preprocess compound " + e);
        }
    }

    public static boolean validateMolecule(final IAtomContainer molecule)
            throws CDKException {
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        for (int i = 0; i < molecule.getAtomCount(); i++) {
            if (null == molecule.getAtom(i).getValency()) {
                return false;
            }
        }
        return true;
    }

}
