package eawag.envirule;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class rule_generator {

    private boolean generalizeIgnoreHydrogen;
    private boolean includeFunctionalGroups;
    private String file;
    private int radius;

    public rule_generator(boolean generalizeIgnoreHydrogen, boolean includeFunctionalGroups, String file, int radius){
        this.generalizeIgnoreHydrogen = generalizeIgnoreHydrogen;
        this.includeFunctionalGroups = includeFunctionalGroups;
        this.file = file;
        this.radius = radius;
    }

    public Set<String> generate() throws  Exception {

        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        final SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.AtomAtomMap);

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








        }


        return null;

    }

    private boolean checkCoA(String base_rule) throws CDKException {

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
