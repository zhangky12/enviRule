package uk.ac.ebi.reactionblast;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class Generate_r_v3 {

    public static void main(String[] args) throws Exception{

    }

    /**
     * Parse the reaction file, and save each reaction into a list of string
     * @param file
     * @return reactions
     * @throws IOException
     */
    public static List<String> parseReactionFile(String file) throws IOException {
        List<String> reactions = new ArrayList<>();

        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        while ((line = br.readLine()) != null) {
            reactions.add(line);
        }
        br.close();
        return reactions;
    }

    /**
     * Perform the atom-mapping for reaction. Each atom will have a MapId.
     * @param cdkReaction
     * @param reactionName
     * @return reaction
     * @throws AssertionError
     * @throws Exception
     */
    public static IReaction performAtomAtomMapping(IReaction cdkReaction, String reactionName) throws AssertionError, Exception {
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

    public static String getNewRule(String unmapped, SmilesParser smilesParser, SmilesGenerator sg, int radius, Boolean unmapUnrelated) throws Exception {

        // Convert reaction string to IReaction
        IReaction Reaction = smilesParser.parseReactionSmiles(unmapped);
        // Perform atom-mapping for the reaction
        IReaction performAtomAtomMapping = performAtomAtomMapping(Reaction, null);

        // Get changed atoms
        BondChangeCalculator bcc1 = new BondChangeCalculator(performAtomAtomMapping);
        bcc1.computeBondChanges(true, false);
        Set<Integer> changed_atom_tags = new HashSet<Integer>();

        // Get the set of MapId for the atoms in reaction center (including both reactants centers and products centers).
        for(IAtom atom: bcc1.getReactionCenterSet()){
            changed_atom_tags.add(atom.getMapIdx());
        }

        // expand the neighbors in reactants
        for(int i=0; i<radius; i++){
            expandReactantNeighbors(performAtomAtomMapping.getReactants(), changed_atom_tags);
        }

        // include all the coming groups in products if there are any
        newAddedGroupsInProducts(performAtomAtomMapping, changed_atom_tags);

        IAtomContainer reactant_fragments = getFragmentsMol(performAtomAtomMapping.getReactants(), changed_atom_tags);
        IAtomContainer product_fragments = getFragmentsMol(performAtomAtomMapping.getProducts(), changed_atom_tags);

//        System.out.println("New rule: " + cleanSMIRKS(reactant_fragments, product_fragments, sg, changed_atom_tags));

        return cleanSMIRKS(reactant_fragments, product_fragments, sg, changed_atom_tags, performAtomAtomMapping, unmapUnrelated);
    }

    /**
     * Include atoms in compounds that are right next to changed atoms - atoms with changed_atom_tags
     * @param compounds
     * @param changed_atom_tags
     */
    public static void expandReactantNeighbors(IAtomContainerSet compounds, Set<Integer> changed_atom_tags){

        Set<Integer> add_atom_tags = new HashSet<>();

        for(int i=0; i<compounds.getAtomContainerCount(); i++){
            IAtomContainer mol = compounds.getAtomContainer(i);
            for(IBond bond: mol.bonds()){
                int count = 0;
                if(changed_atom_tags.contains(bond.getEnd().getMapIdx())) count ++;
                if(changed_atom_tags.contains(bond.getBegin().getMapIdx())) count ++;
                if(count == 1){
                    // Only one side of bond is already included
                    add_atom_tags.add(bond.getBegin().getMapIdx());
                    add_atom_tags.add(bond.getEnd().getMapIdx());
                }
            }
        }
        changed_atom_tags.addAll(add_atom_tags);
    }

    /**
     * For all the atoms only show in products (added groups), include them in changed_atom_tags
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

    /**
     * Get the reaction center substructure from compounds. Take it as a new molecule.
     * @param compounds
     * @param changed_atom_tags
     * @return
     */
    public static IAtomContainer getFragmentsMol(IAtomContainerSet compounds, Set<Integer> changed_atom_tags){

        IAtomContainer fragmentsMol = new AtomContainer();

        for(int i=0; i<compounds.getAtomContainerCount(); i++){
            IAtomContainer mol = compounds.getAtomContainer(i);
            for(IAtom atom: mol.atoms()){
                if(changed_atom_tags.contains(atom.getMapIdx())){
                    fragmentsMol.addAtom(atom);
                }
            }

            for(IBond bond: mol.bonds()){
                if(changed_atom_tags.contains(bond.getBegin().getMapIdx()) && changed_atom_tags.contains(bond.getEnd().getMapIdx())){
                    fragmentsMol.addBond(bond);
                }
            }
        }
        return fragmentsMol;
    }

    public static String cleanSMIRKS(IAtomContainer reactant_fragments, IAtomContainer product_fragments, SmilesGenerator sg, Set<Integer> changed_atom_tags, IReaction performAtomAtomMapping, Boolean unmapUnrelated) throws CDKException {

        String reactant_smart = sg.create(reactant_fragments);
        String product_smart = sg.create(product_fragments);


        String new_smirks = correctBonds(performAtomAtomMapping, reactant_smart, product_smart);

        if(!unmapUnrelated){
            return new_smirks;
        }

        reactant_smart = new_smirks.split(">>")[0];
        product_smart = new_smirks.split(">>")[1];

        Set<Integer> unique_atom_tags = new HashSet<>();
        Set<Integer> intersect_atom_tags = getIntersectMapId(reactant_fragments, product_fragments);

        for(int tag: changed_atom_tags){
            if(!intersect_atom_tags.contains(tag)){
                unique_atom_tags.add(tag);
            }
        }

        for(int tag: unique_atom_tags){
            String tag_s = ":" + tag;
            reactant_smart = reactant_smart.replace(tag_s, "");
            product_smart = product_smart.replace(tag_s, "");
        }

        return reactant_smart+">>"+product_smart;
    }

    /**
     * Create reaction SMIRKS with explicit bonds
     * @param performAtomAtomMapping
     * @param reactant_smart
     * @param product_smart
     * @return
     */
    public static String correctBonds(IReaction performAtomAtomMapping, String reactant_smart, String product_smart){

        IAtomContainerSet reactants = performAtomAtomMapping.getReactants();
        IAtomContainerSet products = performAtomAtomMapping.getProducts();

        reactant_smart = replaceWithExplicitBond(reactants, reactant_smart);
        product_smart = replaceWithExplicitBond(products, product_smart);

        String new_smirks = reactant_smart + ">>" + product_smart;
        return new_smirks;
    }

    /**
     * Replace implicit bonds with explicit bonds
     * @param compounds
     * @param smart
     * @return
     */
    public static String replaceWithExplicitBond(IAtomContainerSet compounds, String smart){

        Map<String, String> replace = new HashMap<>();
        List<String> segments = new ArrayList<>();
        List<Integer> segments_mapId = new ArrayList<>();

        Map<IBond.Order, String> bond_symbol= new HashMap<>();
        bond_symbol.put(IBond.Order.SINGLE, "-");
        bond_symbol.put(IBond.Order.DOUBLE, "=");
        bond_symbol.put(IBond.Order.TRIPLE, "#");

        // Find segments "[Atom:MapId]"
        boolean start_atom = false;
        String atom_symbol = "";
        for(int i=0; i<smart.length(); i++){
            char c = smart.charAt(i);
            if(c == '['){
                start_atom = true;
            }
            if(start_atom){
                atom_symbol += c;
            }
            if(c == ']'){
                start_atom = false;
                // Extract the mapId of this atom
                String map_id = atom_symbol.split(":")[1];
                map_id = map_id.substring(0, map_id.length()-1);
                int mapId = Integer.valueOf(map_id);
                segments_mapId.add(mapId);

                segments.add(atom_symbol);
                atom_symbol = "";
            }
        }

        // If there is only one segment, then no need to add explicit bond
        if(segments.size() < 2){
            return smart;
        }

        // In smarts, for each atom, search for the atom before that is connected with this atom.
        for(int i=1; i<segments.size(); i++){
            boolean foundStart = false;
            // Search the atom before
            for(int j=i-1; j>=0; j--){
                Set<Integer> bond_ends = new HashSet<>();
                bond_ends.add(segments_mapId.get(i));
                bond_ends.add(segments_mapId.get(j));
                // Check if the two atoms are connected
                for(int k=0; k<compounds.getAtomContainerCount(); k++){
                    IAtomContainer mol = compounds.getAtomContainer(k);
                    for(int b=0; b<mol.getBondCount(); b++){
                        IBond bond = mol.getBond(b);
                        if(bond_ends.contains(bond.getBegin().getMapIdx()) && bond_ends.contains(bond.getEnd().getMapIdx())){
                            // If yes, then replace the implicit bond with the bond between these two atoms
                            String replace_string;
                            if(bond.isAromatic()){
                                replace_string = ":" + segments.get(i);
                            }else{
                                replace_string = bond_symbol.get(bond.getOrder()) + segments.get(i);
                            }
                            replace.put(segments.get(i), replace_string);
                            foundStart = true;
                            break;
                        }
                    }
                    if(foundStart) break;
                }
                if(foundStart) break;
            }
        }
        // remove previous explicit bonds
        smart = smart.replace("=", "");
        // add explicit bonds
        for(Map.Entry<String, String> entry: replace.entrySet()){
            smart = smart.replace(entry.getKey(), entry.getValue());
        }
        return smart;
    }

    /**
     * Get the intersection of MapIds between reactant fragments and product fragments
     * @param reactant_fragments
     * @param product_fragments
     * @return
     */
    public static Set<Integer> getIntersectMapId(IAtomContainer reactant_fragments, IAtomContainer product_fragments){

        Set<Integer> reactants_atom_tags = new HashSet<>();
        Set<Integer> products_atom_tags = new HashSet<>();

        for(IAtom atom: reactant_fragments.atoms()){
            reactants_atom_tags.add(atom.getMapIdx());
        }

        for(IAtom atom: product_fragments.atoms()){
            products_atom_tags.add(atom.getMapIdx());
        }

        // Now reactants_atom_tags becomes the intersection
        reactants_atom_tags.retainAll(products_atom_tags);

        return reactants_atom_tags;
    }
}
