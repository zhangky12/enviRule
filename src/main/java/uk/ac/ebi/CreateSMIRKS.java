package uk.ac.ebi;

import org.checkerframework.checker.units.qual.A;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class CreateSMIRKS {

    public static void main(String[] args) throws Exception {

        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

//        String reaction1 = "CC(C(=O)O)OC1=CC=C(C=C1)OC2=CC=C(C=C2Cl)Cl>>CC(C(=O)O)OC1=CC=C(C=C1)OC2=CC(=C(C=C2Cl)Cl)O";
        String reaction1 = "CC1=C(N=C(N=C1C)N(C)C)OC(=O)N(C)C>>CC1=C(N=C(N=C1C)N(C)C)O";
        String file = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/MultistepTrn_SBS_combined_missed_reactions_1116/C-O-reactions.txt";
        List<String> reactions = parsereactionFile(file);
        for(String reaction: reactions){
            getReactionCenters(reaction, smilesParser);
        }

        return;

    }

    public static List<String> parsereactionFile(String file) throws IOException {
        List<String> reactions = new ArrayList<>();

        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        while ((line = br.readLine()) != null) {
            reactions.add(line);
        }

        return reactions;
    }

    public static void getReactionCenters(String reaction1, SmilesParser smilesParser) throws Exception {

        IReaction Reaction = smilesParser.parseReactionSmiles(reaction1);
        IReaction performAtomAtomMapping = performAtomAtomMapping(Reaction, null);
        BondChangeCalculator bcc1 = new BondChangeCalculator(performAtomAtomMapping);
        bcc1.computeBondChanges(true, false);

        Collection<IAtom> centerset = bcc1.getReactionCenterSet();

        IAtomContainerSet reactants = bcc1.getReaction().getReactants();
        IAtomContainerSet products = bcc1.getReaction().getProducts();
        IAtomContainer reactant_center_atoms = new AtomContainer();
        IAtomContainer product_center_atoms = new AtomContainer();

        // Identify atoms in reaction center belonging to reactants & products
        for(IAtom center_atom: centerset){
            for(int i=0; i<reactants.getAtomContainerCount(); i++){
                IAtomContainer reactant_mol = reactants.getAtomContainer(i);
                if(reactant_mol.contains(center_atom)){
                    reactant_center_atoms.addAtom(center_atom);
                    break;
                }
            }

            for(int j=0; j<products.getAtomContainerCount(); j++){
                IAtomContainer product_mol = products.getAtomContainer(j);
                if(product_mol.contains(center_atom)){
                    product_center_atoms.addAtom(center_atom);
                }
            }
        }

        // Find bonds from reactants that are in reaction center
        for(int i=0; i<reactants.getAtomContainerCount(); i++){
            IAtomContainer reactant_mol = reactants.getAtomContainer(i);
            for(IBond bond: reactant_mol.bonds()){
                IAtom begin = bond.getBegin();
                IAtom end = bond.getEnd();
                if(reactant_center_atoms.contains(begin) && reactant_center_atoms.contains(end)){
                    reactant_center_atoms.addBond(bond);
                }
            }
        }

        // Find bonds from products that are in reaction center
        for(int j=0; j<products.getAtomContainerCount(); j++){
            IAtomContainer product_mol = products.getAtomContainer(j);
            for(IBond bond: product_mol.bonds()){
                IAtom begin = bond.getBegin();
                IAtom end = bond.getEnd();
                if(product_center_atoms.contains(begin) && product_center_atoms.contains(end)){
                    product_center_atoms.addBond(bond);
                }
            }
        }


        SmilesGenerator generator = SmilesGenerator.generic();
        String smiles_react = generator.create(reactant_center_atoms);
        String smiles_prod = generator.create(product_center_atoms);
        System.out.println(smiles_react);
        System.out.println(smiles_prod);

        String reaction_smiles = smiles_react + ">>" + smiles_prod;


        IReaction Reaction_final = smilesParser.parseReactionSmiles(reaction_smiles);
        IReaction performAtomAtomMapping_final = performAtomAtomMapping(Reaction_final, null);
//        BondChangeCalculator bcc_final = new BondChangeCalculator(performAtomAtomMapping_final);
//        bcc_final.computeBondChanges(true, false);


        final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
        System.out.println("AAM sm: " + sg.create(performAtomAtomMapping_final));

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

}
