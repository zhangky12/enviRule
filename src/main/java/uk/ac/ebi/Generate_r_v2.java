package uk.ac.ebi;

import com.google.common.collect.Sets;
import org.checkerframework.checker.units.qual.A;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.xmlcml.euclid.Int;
import uk.ac.ebi.beam.Graph;
import uk.ac.ebi.reactionblast.fingerprints.tools.Similarity;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import static com.google.common.base.Preconditions.checkNotNull;

public class Generate_r_v2 {

    public static boolean generalizeIgnoreHydrogen = true;
    public static final Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.edgeShort()));
    public static IAtomContainerSet wholeMolecule= new AtomContainerSet();
    public static boolean isReactant = true;

    @SuppressWarnings("deprecated")
    public static void main(String[] args) throws Exception{

        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);

        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
        int radius = 4;

        // Get the mapped reaction
//        String unmapped = "COC1=NC(=NC(=C1)OC)NC(=O)NS(=O)(=O)CC2=C(C=CC=C2)C(=O)OC>>COC1=CC(=NC(=N1)NC(=O)NS(=O)(=O)CC2=C(C=CC=C2)C(=O)OC)O";
//        String unmapped = "CC1=NN(C)C(=C1C=NOCC2=CC=C(C=C2)C(=O)OC(C)(C)C)OC3=CC=CC=C3>>CC1=NN(C)C(=C1/C=N/OCC2=CC=C(C=C2)C(=O)O)OC3=CC=CC=C3";
//        String unmapped = "CO/N=C(\\C#N)/C(=O)O>>C(#N)/C(=N\\O)/C(=O)O";

//        String file_name = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/MultistepTrn_SB_combined_cleaned_010422/Reactions4Rules/9.txt";
        String file_name = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/BBD_reactions/Reactions4Rules/18-4.txt"; //87-39
        List<String> reactions =  parseReactionFile(file_name);
//        String unmapped = reactions.get(0);

        List<String> newRules = new ArrayList<>();
        List<String> baseRules = new ArrayList<>();

        // Get extended rules and base rules first.
        for(int i=0; i<reactions.size(); i++){
            String unmapped = reactions.get(i);

            // Remove catalysts
            String[] s = unmapped.split(">>");
//            if(!((!s[0].contains(".") && !s[1].contains(".")) || (s[0].contains(".") && s[0].contains(".")))){
//                continue;
//            }
            if(s[0].contains(".")){
                String[] reactants = s[0].split("\\.");
                String[] products = s[1].split("\\.");
                Set<String> reactants_set = new HashSet<String>(List.of(reactants));
                Set<String> products_set = new HashSet<String>(List.of(products));
                if(reactants_set.size() != products_set.size()) continue;
                unmapped = "";
                for(String reactant:reactants){
                    if(products_set.contains(reactant)) continue;
                    unmapped += reactant + ".";
                }
                if(unmapped.length() == 0) continue;
                if(unmapped.charAt(unmapped.length()-1) == '.'){
                    unmapped = unmapped.substring(0, unmapped.length()-1);
                }
                unmapped += ">>";
                for(String product:products){
                    if(reactants_set.contains(product)) continue;
                    unmapped += product + ".";
                }
                if(unmapped.charAt(unmapped.length()-1) == '.'){
                    unmapped = unmapped.substring(0, unmapped.length()-1);
                }
            }
            newRules.add(generalizedForm(unmapped, smilesParser, sg, radius, false));
//            newRules.add(getNewRule(unmapped, smilesParser, sg, radius, true));
            baseRules.add(getNewRule(unmapped, smilesParser, sg, 0, false));
        }

        for(int i=0; i<newRules.size(); i++){
//            System.out.println(newRules.get(i));
            System.out.println(baseRules.get(i));
        }

        if(newRules.isEmpty()) System.exit(1);

        // Get the standardized form of base rule
        String standardized_base = baseStandardize(baseRules.get(0));

        // Correct the mapping of the rules based on standardized base rule
        List<String> standardized_smirks = new ArrayList<>();
        for(int i=0; i<baseRules.size(); i++){
            String base_smirks =  baseRules.get(i);
            Map<String, String> atom_map = standadizedRuleMapping(standardized_base, base_smirks, smilesParser);
            Map<String, String> first_level_map = new HashMap<>();
            Map<String, String> second_level_map = new HashMap<>();

            int start = 0;

            for(Map.Entry<String, String> entry: atom_map.entrySet()){
                start ++;
                first_level_map.put(entry.getKey(), "a"+start);
                second_level_map.put("a"+start, entry.getValue());
            }

            String expanded_smirks = newRules.get(i);

            for(Map.Entry<String, String> entry: first_level_map.entrySet()){
                expanded_smirks = expanded_smirks.replace(entry.getKey(), entry.getValue());
            }

            for(Map.Entry<String, String> entry: second_level_map.entrySet()){
                expanded_smirks = expanded_smirks.replace(entry.getKey(), entry.getValue());
            }
            standardized_smirks.add(expanded_smirks);
        }

        System.out.println("Standardized Form:");

        for(int i=0; i<standardized_smirks.size(); i++){
            System.out.println(standardized_smirks.get(i));
        }

        String starter_rule = cleanGeneralizedRule(standardized_smirks.get(0));
        starter_rule = completeWithHydrogen(starter_rule, smilesParser);
        // Collect the conditions for each atom from each rule
        Map<Integer, Map<String, String>> conditions = new HashMap<>();
        for(int i=0; i<standardized_smirks.size(); i++){
            getConditions(conditions, standardized_smirks.get(i), smilesParser, factory);
        }

        String standardized_base_reactants = starter_rule.split(">>")[0];
        String standardized_base_products = starter_rule.split(">>")[1];
        standardized_base_reactants = standardized_base_reactants.replace(":", ";:");
        for(Map.Entry<Integer, Map<String, String>> entry: conditions.entrySet()){
            for(Map.Entry<String, String> condition_entry: entry.getValue().entrySet()){
                standardized_base_reactants = standardized_base_reactants.replace(":"+entry.getKey(),","+condition_entry.getValue()+":"+entry.getKey());
            }
        }

        standardized_base_reactants = standardized_base_reactants.replace(";:", ":");
        standardized_base_reactants = standardized_base_reactants.replace(";,", ";");

        Set<Integer> unique_id = uniqueMapId(standardized_base_reactants, standardized_base_products);
        for(int tag: unique_id){
            String tag_s = ":" + tag;
            standardized_base_reactants = standardized_base_reactants.replace(tag_s, "");
            standardized_base_products = standardized_base_products.replace(tag_s, "");
        }

        String generalizedRule = standardized_base_reactants + ">>" + standardized_base_products;

        System.out.println("Generalized rule is: \n" + generalizedRule);

    }

    public static String completeWithHydrogen(String smirks, SmilesParser smilesParser) throws InvalidSmilesException {

        String reactants = smirks.split(">>")[0];
        String products = smirks.split(">>")[1];

        String[] reactants_list = reactants.split("\\.");
        String[] products_list = products.split("\\.");

        List<IAtomContainer> reactants_mols = new ArrayList<>();
        List<IAtomContainer> products_mols = new ArrayList<>();

        for(String reactant_mol: reactants_list){
            reactants_mols.add(smilesParser.parseSmiles(reactant_mol));
        }
        for(String product_mol: products_list){
            products_mols.add(smilesParser.parseSmiles(product_mol));
        }

        Map<Integer, Integer> reactants_bonds_count = new HashMap<>();
        Map<Integer, Integer> products_bonds_count = new HashMap<>();

        for(IAtomContainer reactant_mol: reactants_mols){
            for(IAtom atom: reactant_mol.atoms()){
                if(atom.getMapIdx() != 0){
                    int bond_orders = 0;
                    for(IBond bond: reactant_mol.getConnectedBondsList(atom)){
                        if(bond.getOrder() == IBond.Order.SINGLE || bond.isAromatic()){
                            // TODO: need to think whether aromatic bond can be considered as 1
                            bond_orders += 1;
                        }else if(bond.getOrder() == IBond.Order.DOUBLE){
                            bond_orders += 2;
                        }else if(bond.getOrder() == IBond.Order.TRIPLE){
                            bond_orders += 3;
                        }else{
                            // Haven't seen quad bonds yet
                        }

                        // Correct change of valence cases for S and P
                        if((atom.getAtomicNumber() == 16 || atom.getAtomicNumber() == 15)
                                && (bond.getBegin().getAtomicNumber() == 8 || bond.getEnd().getAtomicNumber() == 8)
                                && bond.getOrder() == IBond.Order.DOUBLE){
                            bond_orders -= 2;
                        }
                    }
                    reactants_bonds_count.put(atom.getMapIdx(), bond_orders);
//                    reactants_bonds_count.put(atom.getMapIdx(), atom.getBondCount());
                }
            }
        }

        for(IAtomContainer product_mol: products_mols){
            for(IAtom atom: product_mol.atoms()){
                if(atom.getMapIdx() != 0){
                    int bond_orders = 0;
                    for(IBond bond: product_mol.getConnectedBondsList(atom)){
                        if(bond.getOrder() == IBond.Order.SINGLE || bond.isAromatic()){
                            bond_orders += 1;
                        }else if(bond.getOrder() == IBond.Order.DOUBLE){
                            bond_orders += 2;
                        }else if(bond.getOrder() == IBond.Order.TRIPLE){
                            bond_orders += 3;
                        }else{

                        }
                        // Correct change of valence cases for P and S
                        if((atom.getAtomicNumber() == 16 || atom.getAtomicNumber() == 15)
                                && (bond.getBegin().getAtomicNumber() == 8 || bond.getEnd().getAtomicNumber() == 8)
                                && bond.getOrder() == IBond.Order.DOUBLE){
                            bond_orders -= 2;
                        }
                    }
                    products_bonds_count.put(atom.getMapIdx(), bond_orders);
//                    products_bonds_count.put(atom.getMapIdx(), atom.getBondCount());
                }
            }
        }

        for(Map.Entry<Integer, Integer> entry: reactants_bonds_count.entrySet()){
            if(products_bonds_count.containsKey(entry.getKey())){
                if(products_bonds_count.get(entry.getKey()) - entry.getValue() > 0){
                    int diff = products_bonds_count.get(entry.getKey()) - entry.getValue();
                    String hydrogen = "";
                    for(int j=diff; j>0; j--) hydrogen += "([H])";
                    reactants = reactants.replace(":"+entry.getKey()+"]", ":"+entry.getKey()+"]"+hydrogen);
                }
            }
        }

        return reactants + ">>" + products;

    }

    public static String cleanGeneralizedRule(String rule){

        String reactants = rule.split(">>")[0];
        String products = rule.split(">>")[1];
        Map<String, String> clean_list = new HashMap<>();

        List<String> segments = generalSegments(reactants);
        Map<Integer, String> condition_map = extractConditionFromSegments(segments);

        for(Map.Entry<Integer, String> entry: condition_map.entrySet()){
            String replace = ";" + entry.getValue() + ":" + entry.getKey();
            clean_list.put(replace, ":"+entry.getKey());
        }

        for(Map.Entry<String, String> entry: clean_list.entrySet()){
            reactants = reactants.replace(entry.getKey(), entry.getValue());
        }

        return reactants + ">>" + products;
    }

    public static void getConditions(Map<Integer, Map<String, String>> conditions, String rule, SmilesParser smilesParser, InChIGeneratorFactory factory) throws CDKException {

        String reactants = rule.split(">>")[0];

        List<String> segments = generalSegments(reactants);
        Map<Integer, String> condition_map = extractConditionFromSegments(segments);
        for(Map.Entry<Integer, String> entry: condition_map.entrySet()){
            String condition = entry.getValue();
            if(generalizeIgnoreHydrogen){
                condition = condition.replaceAll("H[0-9]*", "");
            }

            String InchiKey = getInChiKeyOfCondition(condition, smilesParser, factory);
            if(!conditions.containsKey(entry.getKey())) conditions.put(entry.getKey(), new HashMap<>());
            if(!conditions.get(entry.getKey()).containsKey(InchiKey)){
                conditions.get(entry.getKey()).put(InchiKey, condition);
            }
        }
    }

    public static String getInChiKeyOfCondition(String condition, SmilesParser smilesParser, InChIGeneratorFactory factory) throws CDKException {
        String smarts = condition.substring(2);
        smarts = smarts.substring(0, smarts.length()-1);
        IAtomContainer mol = smilesParser.parseSmiles(smarts);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
        InChIGenerator generator = factory.getInChIGenerator(mol);

        ///////////////////////////////////////////////////////////////
        // TODO: this part is not perfect right now
        // It resolves the problem that "[O]-[C](:[C]):[C]" and "[O]-[C](-[C])-[C]" have the same InchiKey
        // Combines the Inchikey with the number of aromatic bonds
        int aromatic_bonds = 0;
        for(IBond bond: mol.bonds()){
            if(bond.isAromatic()) aromatic_bonds ++ ;
        }

        String expanded_inchikey = generator.getInchiKey() + "-" + aromatic_bonds;

        return expanded_inchikey;
    }

    public static Map<Integer, String> extractConditionFromSegments(List<String> segments){

        Map<Integer, String> condition_map = new HashMap<>();

        for(String segment: segments){
            String condition = "";
            String mapId = "";
            boolean start = false;
            List<Character> stack = new ArrayList<>();

            for(int i=0; i<segment.length(); i++){
                char c = segment.charAt(i);
                if(c == '$') {
                    start = true;
                    condition += c;
                    continue;
                }
                if(start) condition += c;

                if(c == '(') stack.add(c);

                if(c == ')') {
                    stack.remove(stack.size()-1);
                    if(stack.size() == 0) {
//                    start = false;
                        mapId = segment.substring(i+1);
                        mapId = mapId.replace(":", "").replace("]", "");
                        break;
                    }
                }
//                if(!start){
//                    if(c!= ')' && c!=':' && c!=']') mapId += c;
//                }
            }

            if(condition.length() < 3) System.exit(1);

            try{
                condition_map.put(Integer.valueOf(mapId), condition);
            }catch (Exception e){
                System.exit(1);
            }
        }

        return condition_map;

    }

    public static List<String> generalSegments(String smarts){

        List<Character> bracket_stack = new ArrayList<>();
        List<String> segments = new ArrayList<>();

        String segment = "";

        for(int i=0; i<smarts.length(); i++){
            char c = smarts.charAt(i);
            if(c == '['){
                bracket_stack.add('[');
                segment += c;
                continue;
            }
            if(bracket_stack.size()!=0) segment += c;

            if(c == ']'){
                bracket_stack.remove(bracket_stack.size()-1);
            }

            if(bracket_stack.size() == 0 && segment.length() != 0){
                segments.add(segment);
                segment = "";
            }
        }
        return segments;
    }

    public static String getNewRule(String unmapped, SmilesParser smilesParser, SmilesGenerator sg, int radius, Boolean unmapUnrelated) throws Exception {

        IReaction Reaction = smilesParser.parseReactionSmiles(unmapped);
        IReaction performAtomAtomMapping = performAtomAtomMapping(Reaction, null);

        // Get changed atoms
        BondChangeCalculator bcc1 = new BondChangeCalculator(performAtomAtomMapping);
        bcc1.computeBondChanges(true, false);
        Set<Integer> changed_atom_tags = new HashSet<Integer>();

        for(IAtom atom: bcc1.getReactionCenterSet()){
            changed_atom_tags.add(atom.getMapIdx());
        }

        // expand the neighbors in reactants
        for(int i=0; i<radius; i++){
            expandReactantNeighbors(performAtomAtomMapping.getReactants(), changed_atom_tags);
        }

        // include functional groups
        if(radius != 0) addFunctionalGroups(performAtomAtomMapping.getReactants(), changed_atom_tags);

        // include all the coming groups in products if there are any
        newAddedGroupsInProducts(performAtomAtomMapping, changed_atom_tags);

        IAtomContainer reactant_fragments = getFragmentsMol(performAtomAtomMapping.getReactants(), changed_atom_tags);
        IAtomContainer product_fragments = getFragmentsMol(performAtomAtomMapping.getProducts(), changed_atom_tags);

//        System.out.println("New rule: " + cleanSMIRKS(reactant_fragments, product_fragments, sg, changed_atom_tags));

        return cleanSMIRKS(reactant_fragments, product_fragments, sg, changed_atom_tags, performAtomAtomMapping, unmapUnrelated);
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
    public static Set<Integer> uniqueMapId(String reactant_smile, String product_smile){
        Set<Integer> unique = new HashSet<>();

        Set<Integer> reactant_Id = new HashSet<>();
        Set<Integer> product_Id = new HashSet<>();

        List<String> reactant_segs = generalSegments(reactant_smile);
        List<String> product_segs = generalSegments(product_smile);

        for(String seg: reactant_segs){
            if(!seg.contains(":")) continue;
            String[] pieces = seg.split(":");
            String last_piece = pieces[pieces.length-1];
            int mapid = -1;
            try{
                mapid = Integer.valueOf(last_piece.substring(0, last_piece.length()-1));
            }catch (Exception e){
                // Do nothing
            }
            if(mapid != -1){
                reactant_Id.add(mapid);
            }
        }

        for(String seg: product_segs){
            if(!seg.contains(":")) continue;
            String[] pieces = seg.split(":");
            String last_piece = pieces[pieces.length-1];
            int mapid = -1;
            try{
                mapid = Integer.valueOf(last_piece.substring(0, last_piece.length()-1));
            }catch (Exception e){
                // Do nothing
            }
            if(mapid != -1){
                product_Id.add(mapid);
            }
        }

        for(int id: reactant_Id){
            if(!product_Id.contains(id)) unique.add(id);
        }

        for(int id: product_Id){
            if(!reactant_Id.contains(id)) unique.add(id);
        }

        return unique;
    }

    public static String generalizedForm(String unmapped, SmilesParser smilesParser, SmilesGenerator sg, int radius, boolean unmapUnrelated) throws Exception {

        String base_rule = getNewRule(unmapped, smilesParser, sg, 0, false);
        String base_reactant = base_rule.split(">>")[0];
        String base_product = base_rule.split(">>")[1];
        String expanded_rule = getNewRule(unmapped, smilesParser, sg, radius, false);
        Set<Integer> unique_id = uniqueMapId(base_reactant, base_product);

        // Reaction center mapIDs.
        Map<Integer, String> segments_map = new HashMap<>();

        boolean start_atom = false;
        String atom_symbol = "";
        for(int i=0; i<base_reactant.length(); i++){
            char c = base_reactant.charAt(i);
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
                segments_map.put(mapId, atom_symbol);
                atom_symbol = "";
            }
        }

        // Convert the expanded reactant to AtomContainer
        String reactant = expanded_rule.split(">>")[0];
        IAtomContainer m = smilesParser.parseSmiles(reactant);

        // Correct aromatic bond based on the whole molecule
        correctAromaticBondInRing(m);

        // Get atom mapId to symbol map
        Map<Integer, String> atom_symbol_map = new HashMap<>();
        for(IAtom atom: m.atoms()){
            atom_symbol_map.put(atom.getMapIdx(), atom.getSymbol());
        }


        // Now, collect all the bonds that connect atoms in reaction center of reactant
        List<IBond> center_bonds = new ArrayList<>();
        for(IBond bond: m.bonds()){
            if(segments_map.containsKey(bond.getBegin().getMapIdx()) && segments_map.containsKey(bond.getEnd().getMapIdx())) center_bonds.add(bond);
        }
        // Remove them from the AtomContainer
        for(IBond bond: center_bonds){
            m.removeBond(bond);
        }

        Map<Integer, String> generalized_map = new HashMap<>();
        for(IAtom atom: m.atoms()){
            if(segments_map.containsKey(atom.getMapIdx())){
                String generalized = getGeneralizedSmilesForAtom(atom, m, sg);
                generalized_map.put(atom.getMapIdx(), generalized);
//                System.out.println(generalized);
            }
        }

        // Replace with the generalized form
        for(Map.Entry<Integer, String> entry: generalized_map.entrySet()){
            base_reactant = base_reactant.replace(segments_map.get(entry.getKey()), entry.getValue());
        }

        if(unmapUnrelated){
            for(int tag: unique_id){
                String tag_s = ":" + tag;
                base_reactant = base_reactant.replace(tag_s, "");
                base_product = base_product.replace(tag_s, "");
            }
        }

        return base_reactant + ">>" + base_product;

    }

    public static void correctAromaticBondInRing(IAtomContainer mol) {

        for(IBond bond: mol.bonds()){

            Set<Integer> bond_ends = new HashSet<>();
            bond_ends.add(bond.getBegin().getMapIdx());
            bond_ends.add(bond.getEnd().getMapIdx());

            for(int i=0; i<wholeMolecule.getAtomContainerCount(); i++){
                IAtomContainer reactant = wholeMolecule.getAtomContainer(i);
                for(IBond reactant_bond: reactant.bonds()){
                    int match = 0;
                    if(bond_ends.contains(reactant_bond.getBegin().getMapIdx())) match++;
                    if(bond_ends.contains(reactant_bond.getEnd().getMapIdx())) match++;

                    if(match == 2){
                        if(bond.getOrder() != reactant_bond.getOrder()){
                            bond.setOrder(reactant_bond.getOrder());
                        }

                        if(bond.isAromatic() != reactant_bond.isAromatic()){
                            bond.setIsAromatic(reactant_bond.isAromatic());
                        }
                    }
                }
            }
        }
        return;
    }

    public static int[] subset(final int[] electrons) {
        int[] vs = new int[electrons.length];
        int n = 0;

        for (int i = 0; i < electrons.length; i++)
            if (electrons[i] >= 0) vs[n++] = i;

        return Arrays.copyOf(vs, n);
    }

    public static String getGeneralizedSmilesForAtom(IAtom center_atom, IAtomContainer expanded_reactant_mol, SmilesGenerator sg) throws CDKException {

        List<IAtom> reachableAtoms = new ArrayList<>();
        Set<Integer> visited = new HashSet<>();
        List<IAtom> atom_queue = new ArrayList<>();
        atom_queue.add(center_atom);
        reachableAtoms.add(center_atom);

        // First, find all the reachable atoms from center_atom
        while(atom_queue.size()!=0){
            IAtom current_atom = atom_queue.remove(0);
            visited.add(current_atom.getMapIdx());
            List<IAtom> neighbors = findNeighborsOfAtom(current_atom, expanded_reactant_mol, visited);
            atom_queue.addAll(neighbors);
            reachableAtoms.addAll(neighbors);
        }

        // Create an ordered new molecule (center_atom is in the first place)
        IAtomContainer ordered_mol = new AtomContainer();
        for(IAtom atom: reachableAtoms){
            ordered_mol.addAtom(atom);
        }

        // Then, add all the bonds from the original molecule that are inside the reachable groups
        // Non-single bond will be first added. In case it is in a ring but no explicitly specified
        // e.g. [C:1]1-C-C-C-[N:2]1. If the bond connecting two ends is actually [N:2]=[C:1], then this SMARTS is wrong
        // Instead, it should be [C:1]1=[N:2]-C-C-C1
        List<IBond> bonds_cache = new ArrayList<>();
        for(IBond bond: expanded_reactant_mol.bonds()){
            if(visited.contains(bond.getBegin().getMapIdx()) && visited.contains(bond.getEnd().getMapIdx())) {
                if(bond.getOrder()!=IBond.Order.SINGLE) ordered_mol.addBond(bond);
                else bonds_cache.add(bond);
            }
        }
        for(IBond bond:bonds_cache) ordered_mol.addBond(bond);

        String fragment_smiles = sg.create(ordered_mol);
        IAtomContainerSet ordered_mol_sets = new AtomContainerSet();
        ordered_mol_sets.addAtomContainer(ordered_mol);

        fragment_smiles = replaceWithExplicitBond(ordered_mol_sets, fragment_smiles, true);

        String generalized = "[" + center_atom.getSymbol() + ";$(" + fragment_smiles.replaceAll(":[0-9]+", "") + ")" + ":" + center_atom.getMapIdx() + "]";

        return generalized;
    }

    public static List<IAtom> findNeighborsOfAtom(IAtom center_atom, IAtomContainer mol, Set<Integer> visited){

        List<IAtom> neighbors = new ArrayList<>();

        for(IBond bond: mol.bonds()){
            if(bond.getBegin().getMapIdx() == center_atom.getMapIdx()){
                if(!visited.contains(bond.getEnd().getMapIdx())){
                    neighbors.add(bond.getEnd());
                }
            }else if(bond.getEnd().getMapIdx() == center_atom.getMapIdx()){
                if(!visited.contains(bond.getBegin().getMapIdx())){
                    neighbors.add(bond.getBegin());
                }
            }
        }

        return neighbors;
    }

    public static IReaction getCenterReaction(IAtomContainer reactant_fragments, IAtomContainer product_fragments) throws CDKException {

        SmilesGenerator generator = SmilesGenerator.generic();
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        String smiles_react = generator.create(reactant_fragments);
        String smiles_prod = generator.create(product_fragments);

        String reaction_smiles = smiles_react + ">>" + smiles_prod;

        IReaction ruleReaction = smilesParser.parseReactionSmiles(reaction_smiles);

        return ruleReaction;

    }

    public static String baseStandardize(String smirks){
        String standard = initializeStandardizedMapping(smirks);
        System.out.println("Initialized standardized string is: " + standard);
        return standard;
    }

    public static String initializeStandardizedMapping(String smirks){

        String reactants1 = smirks.split(">>")[0];
        Map<String, String> level1_map = new HashMap<>();
        Map<String, String> level2_map = new HashMap<>();

        int mapID = 0;
        String tmp_seg = "";
        boolean start = false;


        for(int i=0; i<reactants1.length(); i++){
            char c = reactants1.charAt(i);
            if(c == '[') start = true;
            if(start) tmp_seg += c;
            if(c == ']'){
                start = false;
                int cur_map = -1;
                try{
                    String[] pieces = tmp_seg.split(":");
                    String last_piece = pieces[pieces.length-1];
                    cur_map = Integer.valueOf(last_piece.substring(0, last_piece.length()-1));
                }catch (Exception e){
                    // Do nothing
                }
                if(cur_map != -1){
                    mapID ++;
                    // Two levels, to avoid cases
                    // e.g [1,2] 1->2, 2->1. Want to get [2,1], but will end up with [1,1]
                    level1_map.put(""+cur_map, "a"+mapID);
                    level2_map.put("a"+mapID, ""+mapID);
                }
                tmp_seg = "";
            }
        }

        for(Map.Entry<String, String> entry: level1_map.entrySet()){
            smirks = smirks.replace(":"+entry.getKey()+"]", ":"+entry.getValue()+"]");
        }

        for(Map.Entry<String, String> entry: level2_map.entrySet()){
            smirks = smirks.replace(":"+entry.getKey()+"]", ":"+entry.getValue()+"]");
        }

        return smirks;

    }

    public static boolean CompareAtomsInTwoMols(IAtom atom1, IAtom atom2, IAtomContainer mol1, IAtomContainer mol2) throws Exception {

        if(mol1.getAtomCount() == 1 && mol2.getAtomCount() == 1 && atom1.getAtomicNumber() == atom2.getAtomicNumber()){
            return true;
        }

        MACCSFingerprinter fingerprinter = new MACCSFingerprinter();

        IAtomContainer copied_mol1 = new AtomContainer();
        for(IAtom atom: mol1.atoms()){
            if(atom.getMapIdx() == atom1.getMapIdx()) continue;
            copied_mol1.addAtom(atom);
        }
        for(IBond bond: mol1.bonds()){
            if(bond.getBegin().getMapIdx() == atom1.getMapIdx() || bond.getEnd().getMapIdx() == atom1.getMapIdx()) continue;
            copied_mol1.addBond(bond);
        }

        IAtomContainer copied_mol2 = new AtomContainer();
        for(IAtom atom: mol2.atoms()){
            if(atom.getMapIdx() == atom2.getMapIdx()) continue;
            copied_mol2.addAtom(atom);
        }
        for(IBond bond: mol2.bonds()){
            if(bond.getBegin().getMapIdx() == atom2.getMapIdx() || bond.getEnd().getMapIdx() == atom2.getMapIdx()) continue;
            copied_mol2.addBond(bond);
        }

        // Feb 8 2022
//        IAtomContainer copied_mol1 = new AtomContainer();


        SmilesGenerator g = new SmilesGenerator(SmiFlavor.Canonical);
        SmilesParser sp  = new SmilesParser(SilentChemObjectBuilder.getInstance());
        String smiles1 = g.create(copied_mol1);
        String smiles2 = g.create(copied_mol2);

        copied_mol1 = sp.parseSmiles(smiles1);
        copied_mol2 = sp.parseSmiles(smiles2);




        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(copied_mol1);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(copied_mol2);

        CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(
                DefaultChemObjectBuilder.getInstance()
        );
        adder.addImplicitHydrogens(copied_mol1);
        adder.addImplicitHydrogens(copied_mol2);

        BitSet bitset1 = fingerprinter.getBitFingerprint(copied_mol1).asBitSet();
        BitSet bitSet2 = fingerprinter.getBitFingerprint(copied_mol2).asBitSet();

        float sim_change = Similarity.getTanimotoSimilarity(bitset1, bitSet2);

        if(sim_change != 1.0f) return false;
        else return true;
    }

    public static Map<String, String> standadizedRuleMapping(String smirks1, String smirks2, SmilesParser smilesParser) throws Exception {

        String reactants1 = smirks1.split(">>")[0];
        String reactants2 = smirks2.split(">>")[0];

        String products1 = smirks1.split(">>")[1];
        String products2 = smirks2.split(">>")[1];

        // Remove all the +,-
        reactants1 = reactants1.replaceAll("(?<=[A-Z,a-z])\\+[0-9]*", "").replaceAll("(?<=[A-Z,a-z])-[0-9]*", "");
        reactants2 = reactants2.replaceAll("(?<=[A-Z,a-z])\\+[0-9]*", "").replaceAll("(?<=[A-Z,a-z])-[0-9]*", "");

        products1 = products1.replaceAll("(?<=[A-Z,a-z])\\+[0-9]*", "").replaceAll("(?<=[A-Z,a-z])-[0-9]*", "");
        products2 = products2.replaceAll("(?<=[A-Z,a-z])\\+[0-9]*", "").replaceAll("(?<=[A-Z,a-z])-[0-9]*", "");

        Map<String, String> atom_map = new HashMap<>();

        // Put reactants into molecule
        IAtomContainerSet reactants1_acs = new AtomContainerSet();
        if(reactants1.contains(".")){
            String[] reactants_list1 = reactants1.split("\\.");
            for(String reactant: reactants_list1){
                IAtomContainer mol = smilesParser.parseSmiles(reactant);
                reactants1_acs.addAtomContainer(mol);
            }
        }else{
            IAtomContainer mol = smilesParser.parseSmiles(reactants1);
            reactants1_acs.addAtomContainer(mol);
        }

        IAtomContainerSet reactants2_acs = new AtomContainerSet();
        if(reactants2.contains(".")){
            String[] reactants_list2 = reactants2.split("\\.");
            for(String reactant: reactants_list2){
                IAtomContainer mol = smilesParser.parseSmiles(reactant);
                reactants2_acs.addAtomContainer(mol);
            }
        }else{
            IAtomContainer mol = smilesParser.parseSmiles(reactants2);
            reactants2_acs.addAtomContainer(mol);
        }

        // Put products into molecule
        IAtomContainerSet products1_acs = new AtomContainerSet();
        if(products1.contains(".")){
            String[] products_list1 = products1.split("\\.");
            for(String product: products_list1){
                IAtomContainer mol = smilesParser.parseSmiles(product);
                products1_acs.addAtomContainer(mol);
            }
        }else{
            IAtomContainer mol = smilesParser.parseSmiles(products1);
            products1_acs.addAtomContainer(mol);
        }

        IAtomContainerSet products2_acs = new AtomContainerSet();
        if(products2.contains(".")){
            String[] products_list2 = products2.split("\\.");
            for(String product: products_list2){
                IAtomContainer mol = smilesParser.parseSmiles(product);
                products2_acs.addAtomContainer(mol);
            }
        }else{
            IAtomContainer mol = smilesParser.parseSmiles(products2);
            products2_acs.addAtomContainer(mol);
        }

        SmilesGenerator g = new SmilesGenerator(SmiFlavor.Canonical);

        if(reactants1_acs.getAtomContainerCount() != reactants2_acs.getAtomContainerCount() ||
        products2_acs.getAtomContainerCount() != products2_acs.getAtomContainerCount()) return atom_map;

        // Now, only support at most two reaction centers
        if(reactants1_acs.getAtomContainerCount() == 2){

            Map<String, String> tmp_map1 = compareInAtomContainer(reactants1_acs.getAtomContainer(0), reactants2_acs.getAtomContainer(0));
            Map<String, String> tmp_map2 = compareInAtomContainer(reactants1_acs.getAtomContainer(1), reactants2_acs.getAtomContainer(1));

            Map<String, String> tmp_map3 = compareInAtomContainer(reactants1_acs.getAtomContainer(0), reactants2_acs.getAtomContainer(1));
            Map<String, String> tmp_map4 = compareInAtomContainer(reactants1_acs.getAtomContainer(1), reactants2_acs.getAtomContainer(0));

            if(tmp_map1.size() == 0 || tmp_map2.size() == 0){
                atom_map.putAll(tmp_map3);
                atom_map.putAll(tmp_map4);
            }else if(tmp_map3.size() == 0 || tmp_map4.size() == 0){
                atom_map.putAll(tmp_map1);
                atom_map.putAll(tmp_map2);
            }else{
                if(tmp_map1.size()+tmp_map2.size() >= tmp_map3.size() + tmp_map4.size()){
                    atom_map.putAll(tmp_map1);
                    atom_map.putAll(tmp_map2);
                }else{
                    atom_map.putAll(tmp_map3);
                    atom_map.putAll(tmp_map4);
                }
            }
        }else{

            // First map for products
            Set<String> mapped_atom = new HashSet<>();
            for(int i=0; i<products1_acs.getAtomContainerCount(); i++){
                IAtomContainer mol1 = products1_acs.getAtomContainer(i);
                for(IAtom atom1: mol1.atoms()){
                    if(atom1.getMapIdx() == 0) continue;
                    boolean atom1_found = false;
                    for(int j=0; j<products2_acs.getAtomContainerCount(); j++){
                        IAtomContainer mol2 = products2_acs.getAtomContainer(j);
                        for(IAtom atom2: mol2.atoms()){
                            if(atom2.getMapIdx() == 0) continue;
                            if(CompareAtomsInTwoMols(atom1, atom2, mol1, mol2)){
                                if(atom_map.containsKey(":"+atom2.getMapIdx()+"]")) continue;
                                atom_map.put(":"+atom2.getMapIdx()+"]", ":"+atom1.getMapIdx() +"]");
                                atom1_found = true;
                                mapped_atom.add(":"+atom1.getMapIdx()+"]");
                                break;
                            }
                        }
                        if(atom1_found) break;
                    }
                }
            }

            // Then map for reactants
            for(int i=0; i<reactants1_acs.getAtomContainerCount(); i++){
                IAtomContainer mol1 = reactants1_acs.getAtomContainer(i);
                for(IAtom atom1: mol1.atoms()){
                    if(atom1.getMapIdx() == 0) continue;
                    if(mapped_atom.contains(":"+atom1.getMapIdx()+"]")) continue;
                    boolean atom1_found = false;
                    for(int j=0; j<reactants2_acs.getAtomContainerCount(); j++){
                        IAtomContainer mol2 = reactants2_acs.getAtomContainer(j);
                        for(IAtom atom2: mol2.atoms()){
                            if(atom2.getMapIdx() == 0) continue;
                            if(CompareAtomsInTwoMols(atom1, atom2, mol1, mol2)){
                                if(atom_map.containsKey(":"+atom2.getMapIdx()+"]")) continue;
                                atom_map.put(":"+atom2.getMapIdx()+"]", ":"+atom1.getMapIdx()+"]");
                                atom1_found = true;
                                break;
                            }
                        }
                        if(atom1_found) break;
                    }
                }
            }
        }
        return atom_map;
    }

    public static Map<String, String> compareInAtomContainer(IAtomContainer mol1, IAtomContainer mol2) throws Exception {

        Map<String, String> atom_map = new HashMap<>();

        for(IAtom atom1: mol1.atoms()){
            if(atom1.getMapIdx() == 0) continue;
            for(IAtom atom2: mol2.atoms()){
                if(atom2.getMapIdx() == 0) continue;
                if(CompareAtomsInTwoMols(atom1, atom2, mol1, mol2)){
                    atom_map.put(":"+atom2.getMapIdx()+"]", ":"+atom1.getMapIdx()+"]");
                    break;
                }
            }
        }

        return atom_map;
    }

    public static String completeMapIDs(String smarts){

        int map_id = 0;
        Map<String, String> replace = new HashMap<>();
        Map<String, Integer> segments = getSegments(smarts);

        // Find the largest mapping id
        for(Map.Entry<String, Integer> entry: segments.entrySet()){
            if(entry.getValue() > map_id){
                map_id = entry.getValue();
            }
        }

        for(Map.Entry<String, Integer> entry: segments.entrySet()){
            // Assign mapID to atoms that have no mapID yet.
            if(entry.getValue() == -1){
                map_id ++;
                String replacement = entry.getKey().substring(0, entry.getKey().length()-1) + ":" + map_id + "]";
                replace.put(entry.getKey(), replacement);
            }
        }

        for(Map.Entry<String, String> entry: replace.entrySet()){
            smarts = smarts.replace(entry.getKey(), entry.getValue());
        }

        return smarts;
    }


    public static Map<String, Integer> getSegments(String smarts){

        Map<String, Integer> segments = new HashMap<>();
        boolean start = false;
        String tmp_seg = "";

        for(int i=0; i<smarts.length(); i++){
            char c = smarts.charAt(i);
            if(c == '[') start = true;
            if(start) tmp_seg += c;
            if(c == ']'){
                start = false;
                int map_id = -1;
                try{
                    String[] pieces = tmp_seg.split(":");
                    String last_piece = pieces[pieces.length-1];
                    map_id = Integer.valueOf(last_piece.substring(0, last_piece.length()-1));
                }catch(Exception e){
                    // Do nothing
                }
                segments.put(tmp_seg, map_id);
                tmp_seg = "";
            }
        }

        return segments;
    }

    public static IReaction standardizeRuleMapping(IAtomContainer reactant_fragments, IAtomContainer product_fragments) throws Exception {

        IReaction ruleReaction = getCenterReaction(reactant_fragments, product_fragments);
        IReaction ruleMap = performAtomAtomMapping(ruleReaction, null);

        SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);

        System.out.println("Standardized rule mapping: " + sg.create(ruleMap));

        return ruleMap;
    }

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

    public static String correctBonds(IReaction performAtomAtomMapping, String reactant_smart, String product_smart) throws CDKException {

        IAtomContainerSet reactants = performAtomAtomMapping.getReactants();
        IAtomContainerSet products = performAtomAtomMapping.getProducts();

        isReactant = true;
        reactant_smart = replaceWithExplicitBond(reactants, reactant_smart, false);
        isReactant = false;
        product_smart = replaceWithExplicitBond(products, product_smart, false);

        String new_smirks = reactant_smart + ">>" + product_smart;
        return new_smirks;
    }

    public static String replaceWithExplicitBond(IAtomContainerSet compounds, String smart, boolean fragments) throws CDKException {

        Map<String, String> replace = new HashMap<>();
        List<String> segments = new ArrayList<>();
        List<Integer> segments_mapId = new ArrayList<>();
        Map<IBond.Order, String> bond_symbol= new HashMap<>();
        bond_symbol.put(IBond.Order.SINGLE, "-");
        bond_symbol.put(IBond.Order.DOUBLE, "=");
        bond_symbol.put(IBond.Order.TRIPLE, "#");

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

        if(segments.size() < 2){
            return smart;
        }

        if(!fragments){
            // If it's not only fragments but the whole molecule, then re-assign aromaticity
            for(int k = 0; k<compounds.getAtomContainerCount(); k++){
                IAtomContainer mol = compounds.getAtomContainer(k);
                aromaticity.apply(mol);
            }
            if(isReactant) wholeMolecule = compounds;
        }


        for(int i=1; i<segments.size(); i++){
            boolean foundStart = false;

            for(int j=i-1; j>=0; j--){
                Set<Integer> bond_ends = new HashSet<>();
                bond_ends.add(segments_mapId.get(i));
                bond_ends.add(segments_mapId.get(j));

                for(int k=0; k<compounds.getAtomContainerCount(); k++){
                    IAtomContainer mol = compounds.getAtomContainer(k);

                    /////////////////////////////////////////////////////////////////////////////////////
                    // Determine aromatic bonds with the same method used in enviPath
//                    if(!fragments) aromaticity.apply(mol);
                    /////////////////////////////////////////////////////////////////////////////////////

                    for(int b=0; b<mol.getBondCount(); b++){
                        IBond bond = mol.getBond(b);
                        if(bond_ends.contains(bond.getBegin().getMapIdx()) && bond_ends.contains(bond.getEnd().getMapIdx())){
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

        smart = smart.replace("=", "");
        smart = smart.replace("#", "");
        for(Map.Entry<String, String> entry: replace.entrySet()){
            smart = smart.replace(entry.getKey(), entry.getValue());
        }

        return smart;
    }

    public static void recordBondsInfo(IAtomContainer mol){
        // First MapId, Second MapId, bondtype, isAromatic
        // (FMapId, SMapId) -> [bondtype, isAromatic]
    }

    public static Set<Integer> getIntersectMapId(IAtomContainer reactant_fragments, IAtomContainer product_fragments){

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

    public static void cleanMapIdx(IAtomContainer reactant_fragments, IAtomContainer product_fragments){

        Set<Integer> intersect_atom_tags = getIntersectMapId(reactant_fragments, product_fragments);

        for(IAtom atom: reactant_fragments.atoms()){
            if(!intersect_atom_tags.contains(atom.getMapIdx())){
                atom.setMapIdx(0);
            }
        }

        for(IAtom atom: product_fragments.atoms()){
            if(!intersect_atom_tags.contains(atom.getMapIdx())){
                atom.setMapIdx(0);
            }
        }
    }

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

    public static void expandReactantNeighbors(IAtomContainerSet compounds, Set<Integer> changed_atom_tags){

        Set<Integer> add_atom_tags = new HashSet<Integer>();

        for(int i=0; i<compounds.getAtomContainerCount(); i++){
            IAtomContainer mol = compounds.getAtomContainer(i);
            for(IBond bond: mol.bonds()){
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
