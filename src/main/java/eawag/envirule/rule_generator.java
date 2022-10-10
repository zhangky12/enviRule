package eawag.envirule;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class rule_generator {

    private boolean generalizeIgnoreHydrogen;
    private boolean includeFunctionalGroups;
    private String file;
    private int radius;

    private final SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
    private final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
    private final Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.edgeShort()));

    public rule_generator(boolean generalizeIgnoreHydrogen, boolean includeFunctionalGroups, String file, int radius){
        this.generalizeIgnoreHydrogen = generalizeIgnoreHydrogen;
        this.includeFunctionalGroups = includeFunctionalGroups;
        this.file = file;
        this.radius = radius;
    }

    public Set<String> generate() throws  Exception {



        List<String> reactions = parseReactionFile(file);
        List<String> baseRules = new ArrayList<>();
        // Generalized form of rules
        List<String> newRules = new ArrayList<>();

        // Mapped reaction
        List<IReaction> atomMappingReactions = new ArrayList<>();

        // Final results
        Set<String> final_rules = new HashSet<>();

        // To cut off COA later
        int coa_count = 0;

        // Get rules for each reaction
        for(int i=0; i<reactions.size(); i++){

            // Unmapped form of each reaction
            String unmapped = reactions.get(i);
            String[] components = unmapped.split(">>"); // components[0]: reactants; components[1]: products
            String reactants = components[0];
            String products = components[1];

            // Remove catalysts from reaction
            // Note: now it can only remove catalysts if they show in both reactants and products.
            Set<String> reactants_set = new HashSet<>();
            Set<String> products_set = new HashSet<>();

            if(reactants.contains(".")){
                reactants_set.addAll(List.of(reactants.split("\\.")));
            }
            if(products.contains(".")){
                products_set.addAll(List.of(products.split("\\.")));
            }

            unmapped = "";

            if(reactants_set.size() != 0){
                for(String reactant: reactants.split("\\.")){
                    if(products_set.contains(reactant)) continue; // If one compound shows in both reactants and products
                    unmapped += reactant + ".";
                }
            }else{
                // It means there is only one reactant
                unmapped = reactants;
            }

            // If all compounds in reactants can also be found in products, then this reaction is problematic and should be skipped
            if(unmapped.length() == 0) continue;

            // Remove the redundant dot for reactants part
            if(unmapped.charAt(unmapped.length()-1) == '.'){
                unmapped = unmapped.substring(0, unmapped.length()-1);
            }

            // Append reaction sign
            unmapped += ">>";

            // Remove catalysts for products
            String unmapped_product = "";

            if(products_set.size() != 0){
                for(String product:products.split("\\.")){
                    if(reactants_set.contains(product)) continue;
                    unmapped_product += product + ".";
                }

                if(unmapped_product.charAt(unmapped_product.length()-1) == '.'){
                    unmapped_product = unmapped_product.substring(0, unmapped_product.length()-1);
                }
            }else{
                // If there is no dot in products, then there is only one single product
                unmapped_product += products;
            }

            // If everything in products can be found in reactants, then this reaction is problematic and should be skipped
            if(unmapped_product.length() == 0) continue;

            unmapped += unmapped_product;

            if(checkCoA(unmapped)) {
                coa_count += 1;
            }

            IReaction Reaction = sp.parseReactionSmiles(unmapped);
            IReaction performAtomAtomMapping = performAtomAtomMapping(Reaction, null);
            atomMappingReactions.add(performAtomAtomMapping);

            String base_rule = getNewRule(performAtomAtomMapping, sg, 0, false);
            baseRules.add(base_rule);
            newRules.add(generalizedForm(base_rule, performAtomAtomMapping, sp, sg, radius, false));









        }


        return null;

    }

    private String generalizedForm(String base_rule, IReaction performAtomAtomMapping, SmilesParser smilesParser, SmilesGenerator sg, int radius, boolean unmapUnrelated) throws Exception {

        String base_reactant = base_rule.split(">>")[0];
        String base_product = base_rule.split(">>")[1];

        String expanded_rule = getNewRule(performAtomAtomMapping, sg, radius, false);
        Set<Integer> unique_id = uniqueMapId(base_reactant, base_product);

        
    }

    /**
     * Get unique mapping ID for reactant and product from SMARTS
     * @param reactant_smart
     * @param product_smart
     * @return
     */
    private Set<Integer> uniqueMapId(String reactant_smart, String product_smart){

        Set<Integer> unique = new HashSet<>();

        Set<Integer> reactant_Id = new HashSet<>();
        Set<Integer> product_Id = new HashSet<>();

        // Get segments and mapId of reactants
        List<String> reactant_segments = new ArrayList<>();
        List<Integer> reactant_segmentsId = new ArrayList<>();
        getSegmentsWithID(reactant_smart, reactant_segments, reactant_segmentsId);
        reactant_Id.addAll(reactant_segmentsId);

        // Get segments and mapId of products
        List<String> product_segments = new ArrayList<>();
        List<Integer> product_segmentsId = new ArrayList<>();
        getSegmentsWithID(product_smart, product_segments, product_segmentsId);
        product_Id.addAll(product_segmentsId);

        for(int id: reactant_Id){
            if(!product_Id.contains(id)) unique.add(id);
        }

        for(int id: product_Id){
            if(!reactant_Id.contains(id)) unique.add(id);
        }

        return unique;
    }

    /**
     *
     * @param performAtomAtomMapping: atom-mapped reaction
     * @param sg
     * @param radius: radius to expand reaction center
     * @param unmapUnrelated: whether to set the mapping number of unmapped atoms to null
     * @return
     * @throws Exception
     */
    private String getNewRule(IReaction performAtomAtomMapping, SmilesGenerator sg, int radius, Boolean unmapUnrelated) throws Exception {

        // Get changed atoms
        Set<Integer> changed_atom_tags = getChangedAtoms(performAtomAtomMapping, radius);

        // Extract reaction center from reactants and products
        IAtomContainer reactant_fragments = getFragmentsMol(performAtomAtomMapping.getReactants(), changed_atom_tags);
        IAtomContainer product_fragments = getFragmentsMol(performAtomAtomMapping.getProducts(), changed_atom_tags);

        return cleanSMIRKS(reactant_fragments, product_fragments, sg, changed_atom_tags, performAtomAtomMapping, unmapUnrelated);
    }

    /**
     * Combine SMARTS of reactant fragments and product fragments into SMIRKS
     * @param reactant_fragments
     * @param product_fragments
     * @param sg
     * @param changed_atom_tags
     * @param performAtomAtomMapping
     * @param unmapUnrelated
     * @return
     * @throws CDKException
     */
    private String cleanSMIRKS(IAtomContainer reactant_fragments, IAtomContainer product_fragments, SmilesGenerator sg, Set<Integer> changed_atom_tags, IReaction performAtomAtomMapping, Boolean unmapUnrelated) throws CDKException{

        String reactant_smart = sg.create(reactant_fragments);
        String product_smart = sg.create(product_fragments);

        // First, replace all bonds with explicit form
        String new_smirks = correctBonds(performAtomAtomMapping, reactant_smart, product_smart);

        if(!unmapUnrelated){
            return new_smirks;
        }

        return unmapUnrelated(reactant_fragments, product_fragments, new_smirks, changed_atom_tags);
    }

    private String unmapUnrelated(IAtomContainer reactant_fragments, IAtomContainer product_fragments, String rule, Set<Integer> changed_atom_tags){
        Set<Integer> unique_atom_tags = new HashSet<>();
        Set<Integer> intersect_atom_tags = getIntersectMapId(reactant_fragments, product_fragments);
        String reactant_smart = rule.split(">>")[0];
        String product_smart = rule.split(">>")[1];

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

    private Set<Integer> getIntersectMapId(IAtomContainer reactant_fragments, IAtomContainer product_fragments){

        Set<Integer> reactants_atom_tags = new HashSet<Integer>();
        Set<Integer> products_atom_tags = new HashSet<Integer>();

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

    private  String correctBonds(IReaction performAtomAtomMapping, String reactant_smart, String product_smart) throws CDKException {

        reactant_smart = replaceWithExplicitBond(performAtomAtomMapping.getReactants(), reactant_smart, false, true);
        product_smart = replaceWithExplicitBond(performAtomAtomMapping.getProducts(), product_smart, false, false);

        String new_smirks = reactant_smart + ">>" + product_smart;

        return new_smirks;

    }

    /**
     * Replace bonds in SMARTs with explicit version
     * @param compounds
     * @param smart
     * @param fragments
     * @param isReactant
     * @return
     * @throws CDKException
     */
    private String replaceWithExplicitBond(IAtomContainerSet compounds, String smart, boolean fragments, boolean isReactant) throws CDKException {

        Map<String, String> replace = new HashMap<>();
        // segments: for example, [C;$(...):1][N;$(...),$(...):2], the segments will be {[C;$(...):1], [N;$(...),$(...):2]}, and the segments_mapId will be {1, 2}.
        List<String> segments = new ArrayList<>();
        List<Integer> segments_mapId = new ArrayList<>();
        Map<IBond.Order, String> bond_symbol = new HashMap<>();
        bond_symbol.put(IBond.Order.SINGLE, "-");
        bond_symbol.put(IBond.Order.DOUBLE, "=");

        getSegmentsWithID(smart, segments, segments_mapId);

        if(segments.size() < 2){
            // If there is only one atom, then there is no need to add explicit bond
            return smart;
        }

        if(!fragments){
            // If it's not only fragments but the whole molecule, then re-assign aromaticity
            for(int k = 0; k<compounds.getAtomContainerCount(); k++){
                IAtomContainer mol = compounds.getAtomContainer(k);
                aromaticity.apply(mol);
            }
            if(isReactant) {
//                wholeMolecule = compounds;
            }
        }

        for(int i=1; i< segments.size(); i++){

            // Find the bond of this segment that can be seen in SMARTs.
            boolean foundStart = false;

            for(int j=i-1; j>=0; j--){
                Set<Integer> bond_ends = new HashSet<>();
                bond_ends.add(segments_mapId.get(i));
                bond_ends.add(segments_mapId.get(j));
                for(int k=0; k<compounds.getAtomContainerCount(); k++){
                    IAtomContainer mol = compounds.getAtomContainer(k);

                    for(int b=0; b<mol.getBondCount(); b++){
                        IBond bond = mol.getBond(b);
                        if(bond_ends.contains(bond.getBegin().getMapIdx()) && bond_ends.contains(bond.getEnd().getMapIdx())){
                            String replace_string;
                            if(bond.isAromatic() && isReactant){
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

        // First, get rid of all the existing explicit bond information
        smart = smart.replace("=", "");
        smart = smart.replace("#", "");

        // Next, replace the bond with explicit version
        for(Map.Entry<String, String> entry: replace.entrySet()){
            smart = smart.replace(entry.getKey(), entry.getValue());
        }

        return smart;

    }

    private void getSegmentsWithID(String smart, List<String> segments, List<Integer> segments_mapId){

        boolean start_atom = false;
        String atom_symbol = "";
        int brackets_count = 0;
        for(int i=0; i<smart.length(); i++){
            char c = smart.charAt(i);
            if(c == '['){
                start_atom = true;
                brackets_count += 1;
            }
            if(start_atom){
                atom_symbol += c;
            }
            if(c == ']'){
                brackets_count -= 1;
                if (brackets_count > 0) continue;

                start_atom = false;
                // Extract the mapId of this atom
                String map_id;
                try {
                    map_id = atom_symbol.split(":(?=[0-9])")[1];
                }catch (Exception e){
                    atom_symbol = "";
                    continue;
                }
                map_id = map_id.substring(0, map_id.length()-1);
                int mapId = Integer.valueOf(map_id);
                segments_mapId.add(mapId);

                segments.add(atom_symbol);
                atom_symbol = "";
            }
        }
    }


    /**
     * Extract changed atoms and bonds from compounds, and put them into a new molecule
     * @param compounds
     * @param changed_atom_tags
     * @return fragmentsMol
     */

    private IAtomContainer getFragmentsMol(IAtomContainerSet compounds, Set<Integer> changed_atom_tags){

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

    private Set<Integer> getChangedAtoms(IReaction performAtomAtomMapping, int radius) throws Exception{
        // Get changed atoms
        BondChangeCalculator bcc = new BondChangeCalculator(performAtomAtomMapping);
        bcc.computeBondChanges(true, false);
        Set<Integer> changed_atom_tags = new HashSet<Integer>();

        for(IAtom atom: bcc.getReactionCenterSet()){
            changed_atom_tags.add(atom.getMapIdx());
        }

        // expand the neighbors in reactants
        for(int i=0; i<radius; i++){
            expandReactantNeighbors(performAtomAtomMapping.getReactants(), changed_atom_tags);
        }

        // include functional groups
        if(includeFunctionalGroups){
            if(radius != 0) addFunctionalGroups(performAtomAtomMapping.getReactants(), changed_atom_tags);
        }

        // include all the added groups in products if there are any
        newAddedGroupsInProducts(performAtomAtomMapping, changed_atom_tags);

        return changed_atom_tags;
    }

    /**
     * For all the atoms only in products, they should be added to the reaction center.
     * @param performAtomAtomMapping
     * @param changed_atom_tags
     */
    private void newAddedGroupsInProducts(IReaction performAtomAtomMapping, Set<Integer> changed_atom_tags){

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

    private void addFunctionalGroups(IAtomContainerSet reactants, Set<Integer> changed_atom_tags){
        // The list of functional groups are from
        // "Prediction of Organic Reaction Outcomes Using Machine Learning"
        // https://github.com/connorcoley/ochem_predict_nn/blob/master/data/generate_reaction_templates.py

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

    private void expandReactantNeighbors(IAtomContainerSet compounds, Set<Integer> changed_atom_tags){

        Set<Integer> add_atom_tags = new HashSet<Integer>();

        for(int i=0; i<compounds.getAtomContainerCount(); i++){
            IAtomContainer mol = compounds.getAtomContainer(i);
            for(IBond bond: mol.bonds()){

                // For a bond, if one end atom (A) is included in reaction center, while the other end atom (B) is not. Then B should be added into reaction center.
                int count = 0;
                if(changed_atom_tags.contains(bond.getEnd().getMapIdx())) count ++;
                if(changed_atom_tags.contains(bond.getBegin().getMapIdx())) count ++;
                if(count == 1){
                    add_atom_tags.add(bond.getBegin().getMapIdx());
                    add_atom_tags.add(bond.getEnd().getMapIdx());
                }
            }
        }

        changed_atom_tags.addAll(add_atom_tags);
    }


    private IReaction performAtomAtomMapping(IReaction cdkReaction, String reactionName) throws InvalidSmilesException, AssertionError, Exception {
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

    private boolean checkCoA(String base_rule) throws CDKException {

        Pattern ptrn1 = SmartsPattern.create("CC(C)(COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(O)C(=O)NCCC(=O)NCCS");
        String reactants = base_rule.split(">>")[0];
        String products = base_rule.split(">>")[1];

        IAtomContainer reactants_mol = sp.parseSmiles(reactants);
        IAtomContainer products_mol = sp.parseSmiles(products);

        aromaticity.apply(reactants_mol);
        aromaticity.apply(products_mol);

        Mappings reactants_res = ptrn1.matchAll(reactants_mol);
        Mappings products_res = ptrn1.matchAll(products_mol);

        if (reactants_res.count() == 0 && products_res.count() != 0){
            return true;
        }

        return false;
    }

    private List<String> parseReactionFile(String file) throws IOException {
        List<String> reactions = new ArrayList<>();

        BufferedReader br = new BufferedReader(new FileReader(file));
        String line;
        while ((line = br.readLine()) != null) {
            reactions.add(line);
        }
        br.close();
        return reactions;
    }


}
