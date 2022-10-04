package uk.ac.ebi.reactionblast;

import org.apache.commons.cli.*;
import uk.ac.ebi.reactionblast.tools.Adding_rxn_class;

import java.io.File;
import java.util.Set;

public class RxnAdder_v1 {

    public static Set<String> adder(String[] args) throws Exception{

        Options options = new Options();

        Option reactions = new Option("r", "reactionFile", true, "reaction file to be clustered");
        reactions.setRequired(true);
        options.addOption(reactions);

        Option database = new Option("b", "oldDatabase", true, "directory of old clustered reactions");
        database.setRequired(true);
        options.addOption(database);

        Option newdatabase = new Option("n", "newDatabase", true, "directory for updated reactions");
        newdatabase.setRequired(true);
        options.addOption(newdatabase);

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

        String reactionFile = cmd.getOptionValue("reactionFile");
        String oldDatabase = cmd.getOptionValue("oldDatabase");
        String newDatabase = cmd.getOptionValue("newDatabase");

        File directory = new File(oldDatabase);
        if(!directory.exists()){
            throw new Exception("Directory of clustered reactions is not found!");
        }

        Set<String> changed_rxn_files = Adding_rxn_class.adding(reactionFile, oldDatabase, newDatabase);
        return changed_rxn_files;
    }

}
