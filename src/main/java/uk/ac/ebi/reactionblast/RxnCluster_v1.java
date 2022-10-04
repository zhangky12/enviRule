package uk.ac.ebi.reactionblast;
import org.apache.commons.cli.*;

public class RxnCluster_v1 {

    public static void cluster(String[] args) throws Exception {


        Options options = new Options();

        Option reactions = new Option("r", "reactionFile", true, "reaction file to be clustered");
        reactions.setRequired(true);
        options.addOption(reactions);

        Option directory = new Option("d", "directory", true, "directory of folder to store clustered reactions");
        directory.setRequired(true);
        options.addOption(directory);

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
        String dir = cmd.getOptionValue("directory");

        Clustering_rxn_class.clustering(reactionFile, dir);

    }

}
