package uk.ac.ebi.reactionblast;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import uk.ac.ebi.reactionblast.fingerprints.tools.Similarity;
import uk.ac.ebi.reactionblast.mechanism.BondChangeCalculator;
import uk.ac.ebi.reactionblast.mechanism.MappingSolution;
import uk.ac.ebi.reactionblast.mechanism.ReactionMechanismTool;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;

public class Generate_r_class {

    private boolean generalizeIgnoreHydrogen;
    private boolean includeFunctionalGroups;
    private String file;
    private int radius;

    public Generate_r_class(boolean generalizeIgnoreHydrogen, boolean includeFunctionalGroups, String file, int radius){
        this.generalizeIgnoreHydrogen = generalizeIgnoreHydrogen;
        this.includeFunctionalGroups = includeFunctionalGroups;
        this.file = file;
        this.radius = radius;
    }

    public Set<String> generate() throws Exception {

        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.AtomAtomMap);

        List<String> reactions =  parseReactionFile(file);
        List<String> baseRules = new ArrayList<>();
        // Generalized form of rules
        List<String> newRules = new ArrayList<>();

        // Mapped reaction SMILES
        List<IReaction> atomMappingReactions = new ArrayList<>();

        // Final results
        Set<String> final_rules = new HashSet<>();

        int coa_count = 0;

        // Get rule for each reaction
        for(int i=0; i<reactions.size(); i++){

            // Unmapped form of each reaction
            String unmapped = reactions.get(i);
            String[] components = unmapped.split(">>"); // components[0]: reactants; components[1]: products
            String reactants = components[0];
            String products = components[1];

            // Remove catalysts from reaction
            // TODO: now it can only remove catalysts if they appear in both sides. Later it should also be able to recognize catalysts in only one side.
            Set<String> reactants_set = new HashSet<>();
            Set<String> products_set = new HashSet<>();

            if(reactants.contains(".")){
                reactants_set.addAll(List.of(reactants.split("\\.")));
            }
            if(products.contains(".")){
                products_set.addAll(List.of(products.split("\\.")));
            }

            // Reconstruct the unmapped reaction without any catalysts
            unmapped = "";

            if(reactants_set.size() != 0){
                for (String reactant: reactants.split("\\.")){
                    if(products_set.contains(reactant)) continue;
                    unmapped += reactant + ".";
                }
            }else{
                // Only one reactant
                unmapped = reactants;
            }

            // If everything in reactants can be found in products, then this reaction is problematic and should be skipped
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

            if (checkCoA(unmapped)) {
//                System.out.println("CoA adding reactions. No rules.");
                coa_count += 1;
            }

            IReaction Reaction = smilesParser.parseReactionSmiles(unmapped);
            IReaction performAtomAtomMapping = performAtomAtomMapping(Reaction, null);
            atomMappingReactions.add(performAtomAtomMapping);

            String base_rule = getNewRule(performAtomAtomMapping, smilesGenerator, 0, false);
            baseRules.add(base_rule);
            newRules.add(generalizedForm(base_rule, performAtomAtomMapping, smilesParser, smilesGenerator, radius, false));
        }

        if (coa_count != 0 && coa_count >= (reactions.size()/2)){
            List<String> newRules_update = new ArrayList<>();
            List<String> baseRules_update = new ArrayList<>();
            for (int i=0; i< newRules.size(); i++){
                String reactants = newRules.get(i).split(">>")[0];
                List<String> reactant_segments = new ArrayList<>();
                List<Integer> reactant_segmentsId = new ArrayList<>();
                getSegmentsWithID(reactants, reactant_segments, reactant_segmentsId);

                for (int j=0; j<reactant_segments.size(); j++){
                    if (reactant_segments.get(j).contains("[C")){
                        String base_rule = baseRules.get(i).split(">>")[0] + ">>" + "CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(N2C=NC3=C2N=CN=C3N)O1)O)OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCS" + "-[C:" + reactant_segmentsId.get(j) + "]";
                        String new_rule = reactants + ">>" + "CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(N2C=NC3=C2N=CN=C3N)O1)O)OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCS" + "-[C:" + reactant_segmentsId.get(j) + "]";
                        baseRules_update.add(base_rule);
                        newRules_update.add(new_rule);
                        break;
                    }
                }
            }
            baseRules = baseRules_update;
            newRules = newRules_update;
        }

        if(newRules.isEmpty()) System.exit(1);

        // Get the standardized form of base rule
        String standardized_base = baseStandardize(baseRules.get(0));

        // Check if it has disconnected centers
        if (standardized_base.split(">>")[0].contains(".")){

            InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
            Set<String> unique_InchiKey = new HashSet<>();

            for (int i=0; i<atomMappingReactions.size(); i++){

                String connected_rule = connectCenters(atomMappingReactions.get(i), smilesGenerator, smilesParser, radius);

                String InchiKey_reactant = getInChiKeyOfCondition(connected_rule.split(">>")[0], smilesParser, factory);
                String InchiKey_product = getInChiKeyOfCondition(connected_rule.split(">>")[1], smilesParser, factory);
                String InchiKey_total = InchiKey_reactant + "+" + InchiKey_product;
                if(unique_InchiKey.contains(InchiKey_total)) continue;
                unique_InchiKey.add(InchiKey_total);

//                connected_rule = completeWithHydrogen(connected_rule, atomMappingReactions.get(i),smilesParser);
//                System.out.println("Complete: " + connected_rule);

                final_rules.add(connected_rule);
//                System.out.println(connected_rule);
            }

//            System.out.println("------------------------------------------");
//
//            for(String final_rule: final_rules){
//                System.out.println("\"" + final_rule + "\",");
//            }

            return final_rules;
        }

        // Correct the mapping of the rules based on standardized base rule
        List<String> standardized_base_smirks = new ArrayList<>();
        List<String> standardized_smirks = new ArrayList<>();
        // Keep record of the changed atom mapping
        List<Map<String, String>> atom_map_list = new ArrayList<>();

        for(int i=0; i<baseRules.size(); i++){
            String base_smirks =  baseRules.get(i);
            Map<String, String> atom_map = standadizedRuleMapping(standardized_base, base_smirks, smilesParser);
            atom_map_list.add(atom_map);
            Map<String, String> first_level_map = new HashMap<>();
            Map<String, String> second_level_map = new HashMap<>();

            int start = 0;

            for(Map.Entry<String, String> entry: atom_map.entrySet()){
                start ++;
                first_level_map.put(entry.getKey(), "a"+start+"]");
                second_level_map.put("a"+start+"]", entry.getValue());
            }

            String base_rule_smirks = baseRules.get(i);
            String expanded_smirks = newRules.get(i);

            for(Map.Entry<String, String> entry: first_level_map.entrySet()){
                base_rule_smirks = base_rule_smirks.replace(entry.getKey(), entry.getValue());
                expanded_smirks = expanded_smirks.replace(entry.getKey(), entry.getValue());
            }

            for(Map.Entry<String, String> entry: second_level_map.entrySet()){
                base_rule_smirks = base_rule_smirks.replace(entry.getKey(), entry.getValue());
                expanded_smirks = expanded_smirks.replace(entry.getKey(), entry.getValue());
            }
            standardized_base_smirks.add(base_rule_smirks);
            standardized_smirks.add(expanded_smirks);
        }

        // Complete rules with hydrogen
        List<String> complete_rules = new ArrayList<>();
        for (int i=0; i< standardized_base_smirks.size(); i++){
            String base_rule = standardized_base_smirks.get(i);
            String new_rule = standardized_smirks.get(i);

            if(generalizeIgnoreHydrogen){
                base_rule = base_rule.replaceAll("H[0-9]*", "");
                new_rule = new_rule.replaceAll("H[0-9]*", "");
            }

            if (coa_count ==0 || coa_count < (reactions.size()/2)){
                String complete_rule = completeWithHydrogen3(base_rule, new_rule, atomMappingReactions.get(i), atom_map_list.get(i), smilesParser);
                if (complete_rule.compareTo("")==0) continue;
                complete_rules.add(complete_rule);
            }else{
                complete_rules.add(new_rule);
            }

//            System.out.println(complete_rule);
        }

        List<Map<Integer, String>> rule_condition = new ArrayList<>();
        // Map from conditions to inchiKey
        Map<String, String> conditions_inchiKey = new HashMap<>();

        for (String complete_rule: complete_rules){
            String reactant = complete_rule.split(">>")[0];

            List<String> reactant_segments = new ArrayList<>();
            List<Integer> reactant_segmentsId = new ArrayList<>();
            getSegmentsWithID(reactant, reactant_segments, reactant_segmentsId);
            rule_condition.add(getConditions(reactant_segments, conditions_inchiKey, smilesParser));
        }

        // Extract the standardized mapping ID
        List<Integer> std_mapId = new ArrayList<>(rule_condition.get(0).keySet());

        // First, for each center atom, find out rules that can be combined
        Map<Integer, Map<String, Set<Integer>>> combine_indexes = new HashMap<>();
        for (int mapId: std_mapId){
            Map<String, Set<Integer>> combine_rule_index = new HashMap<>();

            Map<String, Integer> unique_inchiKey = new HashMap<>();

            for (int i=0; i<rule_condition.size(); i++){
                String condition = rule_condition.get(i).get(mapId);
                if (unique_inchiKey.containsKey(conditions_inchiKey.get(condition))){

                    if(!combine_rule_index.containsKey(conditions_inchiKey.get(condition))) combine_rule_index.put(conditions_inchiKey.get(condition), new HashSet<>());
                    combine_rule_index.get(conditions_inchiKey.get(condition)).add(unique_inchiKey.get(conditions_inchiKey.get(condition)));
                    combine_rule_index.get(conditions_inchiKey.get(condition)).add(i);
                }else{
                    unique_inchiKey.put(conditions_inchiKey.get(condition), i);
                }
            }
            combine_indexes.put(mapId, combine_rule_index);
        }

        // Next, for each center atom, combine rules. Make sure there is no overlapping of applicable domain
        Set<String> applicable_domain = new HashSet<>();
        Set<Integer> included_indexes = new HashSet<>();
        for (int mapId: std_mapId){
            Map<String, Set<Integer>> combine_rule_index = combine_indexes.get(mapId);
            for (Map.Entry<String, Set<Integer>> entry: combine_rule_index.entrySet()){
                int count = 0;
                String start_rule = "";

                Map<Integer, Set<String>> unique_condition = new HashMap<>();
                for (int mapId3: std_mapId){
                    if (mapId3 == mapId) continue;
                    unique_condition.put(mapId3, new HashSet<>());
                }

                for (int index: entry.getValue()){

                    included_indexes.add(index);
                    // Domain is a concatenated inchiKey
                    String domain = "";
                    for (int mapId2: std_mapId) domain += conditions_inchiKey.get(rule_condition.get(index).get(mapId2))+"+";
                    if(applicable_domain.contains(domain)) continue;
                    else applicable_domain.add(domain);

                    if (count == 0){
                        start_rule = complete_rules.get(index);
                        count += 1;
                        continue;
                    }

                    String reactants = start_rule.split(">>")[0];
                    String products = start_rule.split(">>")[1];
                    for (int mapId2: std_mapId) {
                        if(mapId2 == mapId) continue;
                        // Avoid duplicate conditions for each center atom
                        if(unique_condition.get(mapId2).contains(rule_condition.get(index).get(mapId2))) continue;
                        if(!rule_condition.get(index).containsKey(mapId2)) continue;
                        reactants = reactants.replace(":"+mapId2, ",$(" + rule_condition.get(index).get(mapId2)+"):"+mapId2);
                        unique_condition.get(mapId2).add(rule_condition.get(index).get(mapId2));
                    }
                    start_rule = reactants + ">>" + products;

                }
                if (start_rule.compareTo("")==0) continue;
                final_rules.add(unmapUnrelatedAtoms(start_rule));
            }
        }

        for (int i=0; i<complete_rules.size(); i++){
            if(included_indexes.contains(i)) continue;
            final_rules.add(unmapUnrelatedAtoms(complete_rules.get(i)));
        }

//        for(String final_rule: final_rules){
//            System.out.println(final_rule);
//        }
//
//        System.out.println("------------------------------------------");
//
//        for(String final_rule: final_rules){
//            System.out.println("\"" + final_rule + "\",");
//        }
        return final_rules;
    }

    private boolean checkCoA(String base_rule) throws CDKException {

        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.edgeShort()));

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

    private String connectCenters(IReaction performAtomAtomMapping, SmilesGenerator sg, SmilesParser sp, int radius) throws Exception {

        String mappedReaction = sg.create(performAtomAtomMapping);
        Set<Integer> changed_atom_tags = getChangedAtoms(performAtomAtomMapping, radius);
        // Extract reaction center from reactants and products
        IAtomContainer reactant_fragments = getFragmentsMol(performAtomAtomMapping.getReactants(), changed_atom_tags);
        IAtomContainer product_fragments = getFragmentsMol(performAtomAtomMapping.getProducts(), changed_atom_tags);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // TODO: Rearrange the order of bonds. Try to avoid the case that double bonds are ignored when writing rings to SMILES. On May 31 2022
        IAtomContainer ordered_product_fragments = new AtomContainer();
        for (IAtom atom: product_fragments.atoms()){
            ordered_product_fragments.addAtom(atom);
        }
        List<IBond> bonds_cache = new ArrayList<>();
        for (IBond bond: product_fragments.bonds()){
            if (bond.getOrder() == IBond.Order.SINGLE) {
                bonds_cache.add(bond);
            }else{
                ordered_product_fragments.addBond(bond);
            }
        }

        for (IBond bond: bonds_cache) ordered_product_fragments.addBond(bond);
        product_fragments = ordered_product_fragments;


        String rule_original_id = cleanSMIRKS(reactant_fragments, product_fragments, sg, changed_atom_tags, performAtomAtomMapping, false);
        String reactants = rule_original_id.split(">>")[0];

        if(!reactants.contains("."))  {
            String connected_rule = rule_original_id;
            if(generalizeIgnoreHydrogen){
                connected_rule = connected_rule.replaceAll("H[0-9]*", "");
            }
            connected_rule = completeWithHydrogen(connected_rule, performAtomAtomMapping,sp);

            return unmapUnrelated(reactant_fragments, product_fragments, connected_rule, changed_atom_tags);

        };

        // Prepare for the dynamic table. First, get all the MapIDs
        Set<Integer> mapIds = new HashSet<>();
        for (int i=0; i< performAtomAtomMapping.getReactants().getAtomContainerCount(); i++){
            for (IAtom atom: performAtomAtomMapping.getReactants().getAtomContainer(i).atoms()){
                if (atom.getMapIdx() != 0) mapIds.add(atom.getMapIdx());
            }
        }

        List<Integer> mapIds_list = new ArrayList<>(mapIds);
        Map<Integer, Integer> ids2index = new HashMap<>();

        for (int i=0; i<mapIds_list.size(); i++){
            ids2index.put(mapIds_list.get(i), i);
        }

        int[][] reactant_dynamic_table = new int[mapIds_list.size()][mapIds_list.size()];

        // Fill dynamic table with large (unrealistic) value
        for (int i=0; i<mapIds_list.size(); i++){
            for (int j=0; j<mapIds_list.size(); j++){
                reactant_dynamic_table[i][j] = 1000;
            }
        }

        IAtomContainer reactant= sp.parseSmiles(mappedReaction.split(">>")[0]);
        for (IAtom atom: reactant.atoms()){
            Set<Integer> visited = new HashSet<>();
            BFSDistance(reactant_dynamic_table, reactant, ids2index, atom, visited);
        }

        // Get disconnected components from reactants
        String[] components = reactants.split("\\.");
        List<List<Integer>> component_ids = new ArrayList<>();
        for (int i=0; i<components.length; i++){
            List<String> segments = new ArrayList<>();
            List<Integer> segments_id = new ArrayList<>();
            getSegmentsWithID(components[i], segments, segments_id);
            component_ids.add(segments_id);
        }

        // Enumerate all the possible connections among all the disconnected components

        List<List<Integer>> iterations = new ArrayList<>();

        for (int id: component_ids.get(0)){
            List<Integer> tmp_list = new ArrayList<>();
            tmp_list.add(id);
            iterations.add(tmp_list);
        }

        for (int i=1; i<component_ids.size(); i++){
            List<List<Integer>> expanded_lists = new ArrayList<>();

            for (List<Integer> list: iterations){
                for (int id: component_ids.get(i)){
                    List<Integer> tmp_list = new ArrayList<>();
                    tmp_list.addAll(list);
                    tmp_list.add(id);
                    expanded_lists.add(tmp_list);
                }
            }
            iterations = expanded_lists;
        }

        int min_dist = 100;
        List<Integer> shortest_path = new ArrayList<>();

        for (List<Integer> ids: iterations){
            List<List<Integer>> permute = permute(ids);
            for (List<Integer> iter: permute){
                int dist = 0;
                for (int j=0; j<iter.size()-1; j++){
                    dist += reactant_dynamic_table[ids2index.get(iter.get(j))][ids2index.get(iter.get(j+1))];
                }
                if (dist < min_dist){
                    min_dist = dist;
                    shortest_path = iter;
                }
            }
        }

        // Trace back the shortest path
        List<Integer> path = new ArrayList<>();
        for (int i=0; i<shortest_path.size()-1; i++){

            int start = shortest_path.get(i);
            int end = shortest_path.get(i+1);
            int dist = reactant_dynamic_table[ids2index.get(start)][ids2index.get(end)];

            for (int step=1; step<dist; step++){
                for (int inter_id: ids2index.keySet()){
                    if (reactant_dynamic_table[ids2index.get(start)][ids2index.get(inter_id)] == 1 &&
                            reactant_dynamic_table[ids2index.get(inter_id)][ids2index.get(end)] == dist-step){
                        path.add(inter_id);
                        start = inter_id;
                        break;
                    }
                }
            }
        }

        // Add atoms on the path to changed_atom_tags
        changed_atom_tags.addAll(path);
        // Connected reactant fragments
        IAtomContainer new_reactant_fragments = getFragmentsMol(performAtomAtomMapping.getReactants(), changed_atom_tags);
        IAtomContainer new_product_fragments = getFragmentsMol(performAtomAtomMapping.getProducts(), changed_atom_tags);

        String connected_rule = cleanSMIRKS(new_reactant_fragments, new_product_fragments, sg, changed_atom_tags, performAtomAtomMapping, false);

        if(generalizeIgnoreHydrogen){
            connected_rule = connected_rule.replaceAll("H[0-9]*", "");
        }

        connected_rule = completeWithHydrogen(connected_rule, performAtomAtomMapping,sp);

        return unmapUnrelated(new_reactant_fragments, new_product_fragments, connected_rule, changed_atom_tags);

//        return connected_rule;
    }

    // The permutation code is adapted from
    // https://java2blog.com/permutations-array-java/
    private List<List<Integer>> permute(List<Integer> arr) {
        List<List<Integer>> list = new ArrayList<>();
        permuteHelper(list, new ArrayList<>(), arr);
        return list;
    }

    private void permuteHelper(List<List<Integer>> list, List<Integer> resultList, List<Integer> arr){

        // Base case
        if(resultList.size() == arr.size()){
            list.add(new ArrayList<>(resultList));
        }
        else{
            for(int i = 0; i < arr.size(); i++){

                if(resultList.contains(arr.get(i)))
                {
                    // If element already exists in the list then skip
                    continue;
                }
                // Choose element
                resultList.add(arr.get(i));
                // Explore
                permuteHelper(list, resultList, arr);
                // Unchoose element
                resultList.remove(resultList.size() - 1);
            }
        }
    }

    private void BFSDistance(int[][] dynamic_table, IAtomContainer compound, Map<Integer, Integer> ids2index, IAtom atom, Set<Integer> visited){

        List<IAtom> queue = new ArrayList<>();
        int row = ids2index.get(atom.getMapIdx());
        dynamic_table[row][row] = 0;
        queue.add(atom);

        while(!queue.isEmpty()){
            IAtom current_atom = queue.remove(0);
            visited.add(current_atom.getMapIdx());

            for (IBond bond: compound.bonds()){
                if(bond.getBegin().getMapIdx() == current_atom.getMapIdx() && !visited.contains(bond.getEnd().getMapIdx())){
                    queue.add(bond.getEnd());
                    if (dynamic_table[row][ids2index.get(bond.getEnd().getMapIdx())] > dynamic_table[row][ids2index.get(current_atom.getMapIdx())] + 1) {
                        dynamic_table[row][ids2index.get(bond.getEnd().getMapIdx())] = dynamic_table[row][ids2index.get(current_atom.getMapIdx())] + 1;
                    }
                }else if(bond.getEnd().getMapIdx() == current_atom.getMapIdx() && !visited.contains(bond.getBegin().getMapIdx())){
                    queue.add(bond.getBegin());
                    if (dynamic_table[row][ids2index.get(bond.getBegin().getMapIdx())] > dynamic_table[row][ids2index.get(current_atom.getMapIdx())] + 1) {
                        dynamic_table[row][ids2index.get(bond.getBegin().getMapIdx())] = dynamic_table[row][ids2index.get(current_atom.getMapIdx())] + 1;
                    }
                }
            }
        }
    }

    private String unmapUnrelatedAtoms(String complete_rule){

        String reactant = complete_rule.split(">>")[0];
        String product = complete_rule.split(">>")[1];

        Set<Integer> unique_id = uniqueMapId(reactant, product);
        List<Integer> unique_id_list = new ArrayList<>(unique_id);
        Collections.reverse(unique_id_list);


        for (int id: unique_id_list){
            complete_rule = complete_rule.replace(":"+id, "");
        }

        return complete_rule;
    }

    /**
     * Extract the conditions for each center atom
     * @param segments
     * @return
     */
    private Map<Integer, String> getConditions(List<String> segments, Map<String, String> conditions_inchiKey, SmilesParser smilesParser) throws CDKException {

        Map<Integer, String> conditions = new HashMap<>();
        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();

        for (String segment: segments){
            String condition = "";
            int mapId = 0;
            java.util.regex.Pattern condition_pattern = java.util.regex.Pattern.compile("(?<=\\$\\().+(?=\\):[0-9]+])");
            Matcher condition_matcher = condition_pattern.matcher(segment);
            if (condition_matcher.find()) {
                condition = condition_matcher.group(0);
            }

            java.util.regex.Pattern mapId_pattern = java.util.regex.Pattern.compile("(?<=\\):)[0-9]+(?=\\])");
            Matcher mapId_matcher = mapId_pattern.matcher(segment);
            if(mapId_matcher.find()) {
                mapId = Integer.valueOf(mapId_matcher.group(0));
            }

            conditions.put(mapId, condition);
            conditions_inchiKey.put(condition, getInChiKeyOfCondition(condition, smilesParser, factory));
        }
        return conditions;
    }

    private String getInChiKeyOfCondition(String condition, SmilesParser smilesParser, InChIGeneratorFactory factory) throws CDKException {

        IAtomContainer mol = smilesParser.parseSmiles(condition);
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

    private Set<Integer> getChangedAtoms(IReaction performAtomAtomMapping, int radius) throws Exception{
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
        if(includeFunctionalGroups){
            if(radius != 0) addFunctionalGroups(performAtomAtomMapping.getReactants(), changed_atom_tags);
        }

        // include all the coming groups in products if there are any
        newAddedGroupsInProducts(performAtomAtomMapping, changed_atom_tags);

        return changed_atom_tags;
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

    private String generalizedForm(String base_rule, IReaction performAtomAtomMapping, SmilesParser smilesParser, SmilesGenerator sg, int radius, boolean unmapUnrelated) throws Exception {

        String base_reactant = base_rule.split(">>")[0];
        String base_product = base_rule.split(">>")[1];

        String expanded_rule = getNewRule(performAtomAtomMapping, sg, radius, false);
        Set<Integer> unique_id = uniqueMapId(base_reactant, base_product);

        // Get segments and mapId of reactants (reaction center)
        List<String> reactant_segments = new ArrayList<>();
        List<Integer> reactant_segmentsId = new ArrayList<>();
        getSegmentsWithID(base_reactant, reactant_segments, reactant_segmentsId);

        // Store segments of reactants (reaction center) into a map
        // Reaction center mapIDs.
        Map<Integer, String> segments_map = new HashMap<>();
        for (int i=0; i<reactant_segments.size(); i++){
            segments_map.put(reactant_segmentsId.get(i), reactant_segments.get(i));
        }

        // Convert the expanded reaction center into AtomContainer
        String reactant = expanded_rule.split(">>")[0];
        IAtomContainer m = smilesParser.parseSmiles(reactant);

        // Correct aromatic bond in reaction center (fragments) based on the original whole molecule
        IAtomContainerSet wholeMolecule = performAtomAtomMapping.getReactants();
        correctAromaticBondInRing(m, wholeMolecule);

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

        // Remove them from the AtomContainer, so the context information for each center atom is independent.
        for(IBond bond: center_bonds){
            m.removeBond(bond);
        }

        Map<Integer, String> generalized_map = new HashMap<>();
        for(IAtom atom: m.atoms()){
            if(segments_map.containsKey(atom.getMapIdx())){
                String generalized = getGeneralizedSmilesForAtom(atom, m, sg);
                generalized_map.put(atom.getMapIdx(), generalized);
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

    private void expandReactantNeighbors(IAtomContainerSet compounds, Set<Integer> changed_atom_tags){

        Set<Integer> add_atom_tags = new HashSet<Integer>();

        for(int i=0; i<compounds.getAtomContainerCount(); i++){
            IAtomContainer mol = compounds.getAtomContainer(i);
            for(IBond bond: mol.bonds()){

                // For a bond, if one end compound (A) is included in reaction center, while the other end compound (B) is not. Then B should be added into reaction center.
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

    /*
    Extract all the atoms and bonds from a compounds, and put them into a new molecule.
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

    private String cleanSMIRKS(IAtomContainer reactant_fragments, IAtomContainer product_fragments, SmilesGenerator sg, Set<Integer> changed_atom_tags, IReaction performAtomAtomMapping, Boolean unmapUnrelated) throws CDKException{

        String reactant_smart = sg.create(reactant_fragments);
        String product_smart = sg.create(product_fragments);

        // First, replace all the bonds with explicit form
        String new_smirks = correctBonds(performAtomAtomMapping, reactant_smart, product_smart);

        if(!unmapUnrelated){
            return new_smirks;
        }

        return unmapUnrelated(reactant_fragments, product_fragments, new_smirks, changed_atom_tags);
    }

    private String correctBonds(IReaction performAtomAtomMapping, String reactant_smart, String product_smart) throws CDKException {

        reactant_smart = replaceWithExplicitBond(performAtomAtomMapping.getReactants(), reactant_smart, false, true);
        product_smart = replaceWithExplicitBond(performAtomAtomMapping.getProducts(), product_smart, false, false);

        String new_smirks = reactant_smart + ">>" + product_smart;

        return new_smirks;
    }

    private boolean aromacitySanityCheck(IBond aro_bond, IAtomContainer compound){

        // If a carbon has C=O, then it cannot be aromatic
        IAtom begin = aro_bond.getBegin();
        IAtom end = aro_bond.getEnd();

        if(begin.getAtomicNumber() == 6){
            for (IBond bond: compound.bonds()){
                if(bond.getBegin().getMapIdx() == begin.getMapIdx()){
                    if (bond.getEnd().getAtomicNumber() == 8 && bond.getOrder() == IBond.Order.DOUBLE) {
                        return false;
                    }
                }else if(bond.getEnd().getMapIdx() == begin.getMapIdx()){
                    if (bond.getBegin().getAtomicNumber() == 8 && bond.getOrder() == IBond.Order.DOUBLE) {
                        return false;
                    }
                }

            }
        }

        if(end.getAtomicNumber() == 6){
            for (IBond bond: compound.bonds()){
                if(bond.getBegin().getMapIdx() == end.getMapIdx()){
                    if (bond.getEnd().getAtomicNumber() == 8 && bond.getOrder() == IBond.Order.DOUBLE) {
                        return false;
                    }
                }else if(bond.getEnd().getMapIdx() == end.getMapIdx()){
                    if (bond.getBegin().getAtomicNumber() == 8 && bond.getOrder() == IBond.Order.DOUBLE) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    /**
     * Replace bonds in SMARTs with explicit version
     *
     * @param compounds
     * @param smart
     * @param fragments: whether it's a fragment or the whole molecule
     * @param isReactant: whether it belongs to reactants
     * @return
     * @throws CDKException
     */
    private String replaceWithExplicitBond(IAtomContainerSet compounds, String smart, boolean fragments, boolean isReactant) throws CDKException {

        Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.edgeShort()));
        Map<String, String> replace = new HashMap<>();
        // segments: for example, [C;$(...):1][N;$(...),$(...):2], the segments will be {[C;$(...):1], [N;$(...),$(...):2]}, and the segments_mapId will be {1, 2}.
        List<String> segments = new ArrayList<>();
        List<Integer> segments_mapId = new ArrayList<>();
        Map<IBond.Order, String> bond_symbol= new HashMap<>();
        bond_symbol.put(IBond.Order.SINGLE, "-");
        bond_symbol.put(IBond.Order.DOUBLE, "=");
        bond_symbol.put(IBond.Order.TRIPLE, "#");

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

        for(int i=1; i<segments.size(); i++){

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
//                            if(bond.isAromatic()){
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

    private void correctAromaticBondInRing(IAtomContainer mol, IAtomContainerSet wholeMolecule) throws CDKException {

        Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.edgeShort()));

        for(IBond bond: mol.bonds()){

            Set<Integer> bond_ends = new HashSet<>();
            bond_ends.add(bond.getBegin().getMapIdx());
            bond_ends.add(bond.getEnd().getMapIdx());

            for(int i=0; i<wholeMolecule.getAtomContainerCount(); i++){
                IAtomContainer reactant = wholeMolecule.getAtomContainer(i);
                aromaticity.apply(reactant);
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

    private String getGeneralizedSmilesForAtom(IAtom center_atom, IAtomContainer expanded_reactant_mol, SmilesGenerator sg) throws CDKException {

        List<IAtom> reachableAtoms = new ArrayList<>();
        Set<Integer> visited = new HashSet<>();
        List<IAtom> atom_queue = new ArrayList<>();
        atom_queue.add(center_atom);
        reachableAtoms.add(center_atom);

        // BFS search: First, find all the reachable atoms from center_atom
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
        // Non-single bond will be first added. In case it is in a ring but not explicitly specified
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

        fragment_smiles = replaceWithExplicitBond(ordered_mol_sets, fragment_smiles, true, true);
        String generalized = "[" + center_atom.getSymbol() + ";$(" + fragment_smiles.replaceAll(":[0-9]+", "") + ")" + ":" + center_atom.getMapIdx() + "]";

        return generalized;
    }

    private List<IAtom> findNeighborsOfAtom(IAtom center_atom, IAtomContainer mol, Set<Integer> visited){

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

    /**
     * Change the mapping IDs of base rule to make them start from 1
     * @param smirks
     * @return
     */
    private String baseStandardize(String smirks){
        String standard = initializeStandardizedMapping(smirks);
//        System.out.println("Initialized standardized string is: " + standard);
        return standard;
    }

    /**
     * Change the mapping IDs of rule to make them start from 1
     * @param smirks
     * @return
     */
    private String initializeStandardizedMapping(String smirks){

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

    /**
     * Change the atom-mapping ID of smirks2 based on smirks1. The same atom in reaction center should have the same mapping ID in smirks1 and smirks2.
     * @param smirks1
     * @param smirks2
     * @param smilesParser
     * @return
     * @throws Exception
     */
    private Map<String, String> standadizedRuleMapping(String smirks1, String smirks2, SmilesParser smilesParser) throws Exception {

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
            // Has more than one reactant
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

        // If in two smirks, the number of reactants or the number of products are different, then the base rules are different, so there
        // is no need to continue.
        if(reactants1_acs.getAtomContainerCount() != reactants2_acs.getAtomContainerCount() ||
                products2_acs.getAtomContainerCount() != products2_acs.getAtomContainerCount()) return atom_map;

        // Now, only support at most two reaction centers
        if(reactants1_acs.getAtomContainerCount() == 2){

            // The order of two centers might be different in two reactants
            // case 1: orders are the same
            Map<String, String> tmp_map1 = compareInAtomContainer(reactants1_acs.getAtomContainer(0), reactants2_acs.getAtomContainer(0));
            Map<String, String> tmp_map2 = compareInAtomContainer(reactants1_acs.getAtomContainer(1), reactants2_acs.getAtomContainer(1));
            // case 2: orders are different
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

    /**
     * Find the same atom (same symbol & same connectivity) in two atom containers
     * @param mol1
     * @param mol2
     * @return
     * @throws Exception
     */
    private Map<String, String> compareInAtomContainer(IAtomContainer mol1, IAtomContainer mol2) throws Exception {

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

    /**
     * Compare whether two atoms in two atom containers are the same (same symbol & same connectivity)
     * @param atom1
     * @param atom2
     * @param mol1
     * @param mol2
     * @return
     * @throws Exception
     */
    private boolean CompareAtomsInTwoMols(IAtom atom1, IAtom atom2, IAtomContainer mol1, IAtomContainer mol2) throws Exception {

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

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // TODO: from May 31 2022
    // TODO: before the order of atom was calculated based on adjacent bonds that are shown in fragments. Now the adjacent bonds in the whole origianl molecules are counted
    private String completeWithHydrogen(String expanded_smirks, IReaction performAtomAtomMapping, SmilesParser smilesParser) throws CDKException {

        Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.edgeShort()));

        String reactants = expanded_smirks.split(">>")[0];
        String products = expanded_smirks.split(">>")[1];

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

        Map<Integer, Integer> OrigReactants_atoms_orders = new HashMap<>();
        for (int i=0; i<performAtomAtomMapping.getReactants().getAtomContainerCount(); i++){
            IAtomContainer reactant = performAtomAtomMapping.getReactants().getAtomContainer(i);
            aromaticity.apply(reactant);
            for (IAtom atom: reactant.atoms()){
                if (atom.getMapIdx() != 0){
                    int bond_orders = 0;
                    for(IBond bond: reactant.getConnectedBondsList(atom)){
                        if (bond.getOrder() == IBond.Order.DOUBLE){
                            bond_orders += 2;
                        }else if(bond.getOrder() == IBond.Order.SINGLE){
                            bond_orders += 1;
                        }else if(bond.getOrder() == IBond.Order.TRIPLE){
                            bond_orders += 3;
                        }else{
                            // Haven't seen quad bonds
                        }

                        // Correct change of valence cases for S and P
                        if((atom.getAtomicNumber() == 16 || atom.getAtomicNumber() == 15)
                                && (bond.getBegin().getAtomicNumber() == 8 || bond.getEnd().getAtomicNumber() == 8)
                                && bond.getOrder() == IBond.Order.DOUBLE){
                            bond_orders -= 2;
                        }
                    }

                    OrigReactants_atoms_orders.put(atom.getMapIdx(), bond_orders);
                }
            }
        }

        Map<Integer, Integer> OrigProducts_atoms_orders = new HashMap<>();
        for (int i=0; i<performAtomAtomMapping.getProducts().getAtomContainerCount(); i++){
            IAtomContainer product = performAtomAtomMapping.getProducts().getAtomContainer(i);
            aromaticity.apply(product);
            for (IAtom atom: product.atoms()){
                if (atom.getMapIdx() != 0){
                    int bond_orders = 0;
                    for(IBond bond: product.getConnectedBondsList(atom)){
                        if (bond.getOrder() == IBond.Order.DOUBLE){
                            bond_orders += 2;
                        }else if(bond.getOrder() == IBond.Order.SINGLE){
                            bond_orders += 1;
                        }else if(bond.getOrder() == IBond.Order.TRIPLE){
                            bond_orders += 3;
                        }else{
                            // Haven't seen quad bonds
                        }

                        // Correct change of valence cases for S and P
                        if((atom.getAtomicNumber() == 16 || atom.getAtomicNumber() == 15)
                                && (bond.getBegin().getAtomicNumber() == 8 || bond.getEnd().getAtomicNumber() == 8)
                                && bond.getOrder() == IBond.Order.DOUBLE){
                            bond_orders -= 2;
                        }
                    }

                    OrigProducts_atoms_orders.put(atom.getMapIdx(), bond_orders);
                }
            }
        }

        if(expanded_smirks != null) reactants = expanded_smirks.split(">>")[0];

        for(Map.Entry<Integer, Integer> entry: OrigReactants_atoms_orders.entrySet()){
            if(OrigProducts_atoms_orders.containsKey(entry.getKey())){
                if(OrigProducts_atoms_orders.get(entry.getKey()) - entry.getValue() > 0){
                    int diff = OrigProducts_atoms_orders.get(entry.getKey()) - entry.getValue();
                    String hydrogen = "";
                    for(int j=diff; j>0; j--) hydrogen += "([H])";
                    reactants = reactants.replace(":"+entry.getKey()+"]", ":"+entry.getKey()+"]"+hydrogen);
                }
            }
        }

        return reactants + ">>" + products;
    }

    private String completeWithHydrogen3(String smirks, String expanded_smirks, IReaction performAtomAtomMapping, Map<String, String> atom_map, SmilesParser smilesParser) throws CDKException {

        if(atom_map.isEmpty()) return "";

        Map<Integer, Integer> reverse_mapIds = new HashMap<>();
        for (Map.Entry<String, String> entry: atom_map.entrySet()){
            int key = Integer.valueOf(entry.getValue().replace(":", "").replace("]", ""));
            int value = Integer.valueOf(entry.getKey().replace(":", "").replace("]", ""));
            reverse_mapIds.put(key, value);
        }

        // Get implicit hydrogen counts for each atom
        IAtomContainerSet reactants_set = performAtomAtomMapping.getReactants();
        IAtomContainerSet products_set = performAtomAtomMapping.getProducts();
        Map<Integer, Integer> reactants_hydrogen = new HashMap<>();
        Map<Integer, Integer> products_hydrogen = new HashMap<>();

        for (int i=0; i<reactants_set.getAtomContainerCount(); i++){
            IAtomContainer mol = reactants_set.getAtomContainer(i);
            for (IAtom atom: mol.atoms()){
                reactants_hydrogen.put(atom.getMapIdx(), atom.getImplicitHydrogenCount());
            }
        }

        for (int i=0; i<products_set.getAtomContainerCount(); i++){
            IAtomContainer mol = products_set.getAtomContainer(i);
            for (IAtom atom: mol.atoms()){
                products_hydrogen.put(atom.getMapIdx(), atom.getImplicitHydrogenCount());
            }
        }

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

        // implicit hydrogen bond count
        Map<Integer, Integer> reactants_bonds_count = new HashMap<>();
        Map<Integer, Integer> products_bonds_count = new HashMap<>();

        for(IAtomContainer reactant_mol: reactants_mols){
            for(IAtom atom: reactant_mol.atoms()) {
                if(atom.getMapIdx() != 0 && reverse_mapIds.containsKey(atom.getMapIdx())){
                    int hydrogen_count = reactants_hydrogen.get(reverse_mapIds.get(atom.getMapIdx()));
                    reactants_bonds_count.put(atom.getMapIdx(), hydrogen_count);
                }
            }
        }

        for(IAtomContainer product_mol: products_mols){
            for(IAtom atom: product_mol.atoms()) {
                if(atom.getMapIdx() != 0 && reverse_mapIds.containsKey(atom.getMapIdx())){
                    int hydrogen_count = products_hydrogen.get(reverse_mapIds.get(atom.getMapIdx()));
                    products_bonds_count.put(atom.getMapIdx(), hydrogen_count);
                }
            }
        }

        if(expanded_smirks != null) reactants = expanded_smirks.split(">>")[0];

        for(Map.Entry<Integer, Integer> entry: reactants_bonds_count.entrySet()){
            if(products_bonds_count.containsKey(entry.getKey())){
                if(products_bonds_count.get(entry.getKey()) - entry.getValue() < 0){
                    int diff = entry.getValue() - products_bonds_count.get(entry.getKey()) ;
                    String hydrogen = "";
                    for(int j=diff; j>0; j--) hydrogen += "([H])";
                    reactants = reactants.replace(":"+entry.getKey()+"]", ":"+entry.getKey()+"]"+hydrogen);
                }
            }
        }
        return reactants + ">>" + products;
    }

}
