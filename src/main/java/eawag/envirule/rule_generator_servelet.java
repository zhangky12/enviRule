package eawag.envirule;

import org.apache.commons.cli.*;
import uk.ac.ebi.reactionblast.Generate_r_class;

import java.util.Set;

public class rule_generator_servelet {

    public static Set<String> test(String[] args) throws Exception {

        Options options = new Options();

        Option ignore = new Option("i", "ignoreHydrogen", true, "further generalize rules by ignoring hydrogen");
        ignore.setRequired(true);
        options.addOption(ignore);

        Option fgs = new Option("fg", "functionalGroups", true, "whether to include the whole functional groups if only part of them are already included");
        fgs.setRequired(false);
        options.addOption(fgs);

        Option input = new Option("f", "inputFile", true, "path of input file");
        input.setRequired(true);
        options.addOption(input);

        Option radius = new Option("r", "radius", true, "radius of expanded rules");
        radius.setRequired(true);
        options.addOption(radius);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("utility-name", options);
            throw e;
        }

        boolean generalizeIgnoreHydrogen = Boolean.parseBoolean(cmd.getOptionValue("ignoreHydrogen"));
        boolean includeFunctionalGroups = Boolean.parseBoolean(cmd.getOptionValue("functionalGroups"));
        String file = cmd.getOptionValue("inputFile");
        int r = Integer.valueOf(cmd.getOptionValue("radius"));

        rule_generator generator = new rule_generator(generalizeIgnoreHydrogen, includeFunctionalGroups, file, r);
        Set<String> rules = generator.generate();

        for (String rule: rules){
            System.out.println(rule);
        }

        return generator.generate();
    }

    public Set<String> run(boolean generalizeIgnoreHydrogen, boolean includeFunctionalGroups, String file, int radius) throws Exception {

        Generate_r_class generator = new Generate_r_class(generalizeIgnoreHydrogen, includeFunctionalGroups, file, radius);
        return generator.generate();
    }

}
