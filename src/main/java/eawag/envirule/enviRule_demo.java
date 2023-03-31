package eawag.envirule;

import eawag.envirule.modules.rule_generator;

import java.util.Set;

public class enviRule_demo {

    public static void main(String[] args) throws Exception {
        boolean generalizeIgnoreHydrogen = false;
        boolean includeFunctionalGroups = true;
//        String file = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/BBD_reactions/Reactions4Rules/99-15.txt";
//        String file = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/BBD_autoRule_exp/run_1/reactions4rules/11-2.txt";
//        String file = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/BBD_SOIL_autoRule_formal/bbd_soil_reactions4rules/326-3.txt";
        String file = "/Users/kunyang/Documents/Eawag/envipath/envipath-python/enviPath_python/BBD_SOIL_autoRule_formal/reactions4rules/327-5.txt";
//        String file = "/Users/kunyang/Downloads/N-ring-break-2.txt";
//        String file = "/Users/kunyang/Downloads/38-4-0.txt";
        int radius = 1;

        rule_generator generator = new rule_generator(generalizeIgnoreHydrogen, includeFunctionalGroups, file, radius);
        Set<String> rules = generator.generate();

        for (String rule: rules){
            System.out.println(rule);
        }
    }

}
