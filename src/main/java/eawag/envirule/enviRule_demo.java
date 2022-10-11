package eawag.envirule;

import eawag.envirule.modules.rule_generator;

import java.util.Set;

public class enviRule_demo {

    public static void main(String[] args) throws Exception {
        boolean generalizeIgnoreHydrogen = true;
        boolean includeFunctionalGroups = true;
        String file = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/BBD_reactions/Reactions4Rules/99-15.txt";
//        String file = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/BBD_autoRule_exp/run_1/reactions4rules/11-2.txt";
//        String file = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/BBD_SOIL_autoRule_formal/bbd_soil_reactions4rules/699-6.txt";
//        String file = "/Users/kunyang/Downloads/N-ring-break.txt";
        int radius = 2;

        rule_generator generator = new rule_generator(generalizeIgnoreHydrogen, includeFunctionalGroups, file, radius);
        Set<String> rules = generator.generate();

        for (String rule: rules){
            System.out.println(rule);
        }
    }

}
