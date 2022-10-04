package uk.ac.ebi;

import org.checkerframework.checker.units.qual.A;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;

import java.util.*;

public class Generate_rules {

    public static void main(String[] args) throws Exception{

        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
        IReaction Reaction;
        IReaction performAtomAtomMapping;
        BondChangeCalculator bcc1;

        // Get the mapped reaction
        String unmapped = "C[C@]1(C2=CC=CC=C2)C(=O)N(C(=N1)SC)NC3=CC=CC=C3>>C[C@]1(C2=CC=CC=C2)C(=O)N(C(=N1)SC)NC3=C(C=CC=C3)N(=O)=O";
        Reaction = smilesParser.parseReactionSmiles(unmapped);
        performAtomAtomMapping = performAtomAtomMapping(Reaction, null);

        System.out.println("AAM sm: " + sg.create(performAtomAtomMapping));

        // Get changed atoms
        bcc1 = new BondChangeCalculator(performAtomAtomMapping);
        bcc1.computeBondChanges(true, false);
        Set<Integer> changed_atom_tags = new HashSet<Integer>();

        for(IAtom atom: bcc1.getReactionCenterSet()){
            changed_atom_tags.add(atom.getMapIdx());
        }

        get_fragments_for_changed_atoms(performAtomAtomMapping.getReactants(), changed_atom_tags, 1, "reactants", new HashSet<Integer>());


    }

    public static void get_fragments_for_changed_atoms(IAtomContainerSet compounds, Set<Integer> changed_atom_tags, Integer radius, String type, Set<Integer> expansion){

        String fragments = "";
        Map<Integer, String> symbol_replacement = new HashMap<Integer, String>();
        List<String> symbols = new ArrayList<String>();

        for(int i=0; i<compounds.getAtomContainerCount(); i++){
            IAtomContainer mol = compounds.getAtomContainer(i);

            // TODO: special groups

            Set<Integer> atoms_to_use = new HashSet<Integer>();
            String symbol = "";
            for(IAtom atom: mol.atoms()){
                // First check whether this atom is in the reaction center
                if(changed_atom_tags.contains(atom.getMapIdx())){
                    atoms_to_use.add(atom.getMapIdx());
                    symbol = "[" + atom.getSymbol() + ";H" + atom.getImplicitHydrogenCount() + ":" + atom.getMapIdx() + "]";
                    symbol_replacement.put(atom.getMapIdx(), symbol);
                }
            }

            for(int j=1; j<=radius; j++){
                expand_atoms_to_use(mol, atoms_to_use, symbol_replacement);
            }

            if(type.compareTo("products") == 0){
                if(expansion.size() != 0){
                    for(IAtom atom: mol.atoms()){
//                        if(!changed_atom_tags.contains(atom.getMapIdx())) continue;
                        int label = atom.getMapIdx();
                        if(expansion.contains(label) && !changed_atom_tags.contains(label)){
                            atoms_to_use.add(label);
                            // make the expansion a wildcard
                            symbol_replacement.put(label, convert_atom_to_wildcard(atom));
                        }

                    }
                }
                // Make sure unmapped atoms are included (from products)
                // TODO: whether to include all the unmapped atoms
            }

            // Define new symbols to replace terminal species with wildcard to avoid templates being too strict
            for(IAtom atom: mol.atoms()){
                if(!atoms_to_use.contains(atom.getMapIdx())) continue;
                String atom_symbol = "[" + atom.getSymbol() + ":" + atom.getMapIdx() + "]";
                if(symbol_replacement.containsKey(atom.getMapIdx())){
                    atom_symbol = symbol_replacement.get(atom.getMapIdx());
                }
                symbols.add(atom_symbol);
            }

            if(atoms_to_use.size() == 0) continue;

            System.out.println(symbols);


        }

    }

    public static Set<Integer> expand_changed_atom_tags(Set<Integer> changed_atom_tags, String reactant_fragments){
        Set<Integer> expansion = new HashSet<Integer>();
        Set<Integer> reactant_atom_tags = new HashSet<Integer>();
        boolean start_atom_number = false;
        String atom_number = "";

        for(int i=0; i<reactant_fragments.length(); i++){
            char c = reactant_fragments.charAt(i);
            if(c == ':'){
                start_atom_number = true;
                continue;
            }
            if(c == ']'){
                start_atom_number = false;
                if(atom_number.length() != 0){
                    reactant_atom_tags.add(Integer.valueOf(atom_number));
                    atom_number = "";
                }
                continue;
            }
            if(start_atom_number){
                atom_number += c;
            }
        }
        for(int atom_tag: reactant_atom_tags){
            if(!changed_atom_tags.contains(atom_tag)){
                expansion.add(atom_tag);
            }
        }
        return expansion;
    }

    public static void expand_atoms_to_use(IAtomContainer mol, Set<Integer> atoms_to_use, Map<Integer, String> symbol_replacement){

        for(IAtom atom: mol.atoms()){
            if(!atoms_to_use.contains(atom.getMapIdx())) continue;
            Set<IAtom> neighbors = getNeighbors(mol, atom);
            for(IAtom neighbor: neighbors){
                expand_atoms_to_use_atom(mol, atoms_to_use, neighbor.getMapIdx(), symbol_replacement);
            }
        }
    }

    public static void expand_atoms_to_use_atom(IAtomContainer mol, Set<Integer> atoms_to_use, Integer atom_idx, Map<Integer, String> symbol_replacement){

        if(atoms_to_use.contains(atom_idx)){
            return;
        }

        // TODO: found in groups

        // Include this atom
        atoms_to_use.add(atom_idx);
        for(IAtom atom: mol.atoms()){
            if(atom.getMapIdx() == atom_idx){
                symbol_replacement.put(atom_idx, convert_atom_to_wildcard(atom));
                break;
            }
        }
    }

    public static String convert_atom_to_wildcard(IAtom atom){
        // By default, not super general
        String symbol = "[";

        if(atom.getAtomicNumber() != 6){
            symbol += "#" + atom.getAtomicNumber();
        }else if(atom.isAromatic()){
            symbol += "c";
        }else{
            symbol += "C";
        }

        // Strip extra semicolon (is it necessary?)
        if(symbol.charAt(symbol.length()-1) == ';'){
            symbol = symbol.substring(0, symbol.length()-1);
        }

        // Close with label.
        // TODO: need to think about how to deal with unmapped atoms (but have MapID)
        symbol += ";" + atom.getMapIdx() + "]";

        return symbol;
    }

    public static Set<IAtom> getNeighbors(IAtomContainer mol, IAtom atom){

        Set<IAtom> neighbors = new HashSet<IAtom>();

        for(int i=0; i<mol.getBondCount(); i++){
            IBond bond = mol.getBond(i);
            if(bond.getBegin().getMapIdx() == atom.getMapIdx()){
                neighbors.add(bond.getEnd());
            }else if(bond.getEnd().getMapIdx() == atom.getMapIdx()){
                neighbors.add(bond.getBegin());
            }
        }

        return neighbors;
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
