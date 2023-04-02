package eawag.envirule;

import eawag.envirule.servelet.rule_generator_servelet;
import eawag.envirule.servelet.rxn_adder_servelet;
import eawag.envirule.servelet.rxn_clusterer_servelet;
import java.util.*;

public class enviRule_toy {

    public static void main(String[] args) throws Exception {

        CommandLineParser clp = new CommandLineParser(args);

        Set<String> operation_set = new HashSet<>();
        operation_set.add("cluster");
        operation_set.add("generator");
        operation_set.add("adder");

        if(!clp.map.containsKey("operation") || clp.getArgumentValue("operation").length != 1 ||
                !operation_set.contains(clp.getArgumentValue("operation")[0])){
            System.out.println("One and only one operation please. (e.g., cluster, generator, adder)");
            System.exit(1);
        }

        if(clp.getArgumentValue("operation")[0].compareTo("cluster")==0){

            // Clustering

            if(!clp.map.containsKey("reactions") || clp.getArgumentValue("reactions").length != 1){
                System.out.println("One and only one reaction file please.");
                System.exit(1);
            }

            if(!clp.map.containsKey("output") || clp.getArgumentValue("output").length != 1){
                System.out.println("One output directory please. End with '/'");
                System.exit(1);
            }

            System.out.println("Start reaction clustering:");

            String cluster_command = "rxnclust -r " + clp.getArgumentValue("reactions")[0] + " -d " + clp.getArgumentValue("output")[0];

            try{
                rxn_clusterer_servelet.cluster(cluster_command.split(" "));
            }catch (Exception e){
                System.out.println(e.getMessage());
            }

            return;

        }else if(clp.getArgumentValue("operation")[0].compareTo("generator")==0){

            // Rule generation

            if(!clp.map.containsKey("group") || clp.getArgumentValue("group").length != 1){
                System.out.println("One and only one reaction group file please.");
                System.exit(1);
            }

            if(!clp.map.containsKey("radius") || clp.getArgumentValue("radius").length != 1){
                System.out.println("One and only one radius for expanding reaction centers please.");
                System.exit(1);
            }

            System.out.println("Start run generation:");

            String generator_command = "autorule -i false -fg true -f " + clp.getArgumentValue("group")[0] + " -r " + clp.getArgumentValue("radius")[0];
            try{
                rule_generator_servelet.test(generator_command.split(" "));
            }catch (Exception e){
                System.out.println(e.getMessage());
            }

            System.out.println("Done");

            return;
        }else{

            // Adding reactions

            if(!clp.map.containsKey("reactions") || clp.getArgumentValue("reactions").length != 1){
                System.out.println("Choose a file with new reactions please.");
                System.exit(1);
            }

            if(!clp.map.containsKey("input") || clp.getArgumentValue("input").length != 1){
                System.out.println("The folder with old reaction groups please. End with '/'");
                System.exit(1);
            }

            if(!clp.map.containsKey("output") || clp.getArgumentValue("output").length != 1){
                System.out.println("The folder with updated reaction groups please. End with '/'");
                System.exit(1);
            }

            System.out.println("Start reaction adding:");

            String add_command = "rxnadder -r " + clp.getArgumentValue("reactions")[0] + " -b " + clp.getArgumentValue("input")[0] + " -n " + clp.getArgumentValue("output")[0];

            Set<String> res = null;
            try{
                res = rxn_adder_servelet.adder(add_command.split(" "));
            }catch (Exception e){
                System.out.println(e.getMessage());
            }
            System.out.println("Updated or created reaction groups:");
            for (String c: res){
                System.out.println(c);
            }

            return;
        }
    }

    public static class CommandLineParser {
        List <String> args;
        HashMap<String, List<String>> map = new HashMap<>();
        Set<String> flags = new HashSet<>();

        CommandLineParser(String arguments[])
        {
            this.args = Arrays.asList(arguments);
            map();
        }

        // Return argument names
        public Set<String> getArgumentNames()
        {
            Set<String> argumentNames = new HashSet<>();
            argumentNames.addAll(flags);
            argumentNames.addAll(map.keySet());
            return argumentNames;
        }

        // Check if flag is given
        public boolean getFlag(String flagName)
        {
            return flags.contains(flagName);
        }

        // Return argument value for particular argument name
        public String[] getArgumentValue(String argumentName)
        {
            if(map.containsKey(argumentName))
                return map.get(argumentName).toArray(new String[0]);
            else
                return null;
        }

        // Map the flags and argument names with the values
        public void map()
        {
            for(String arg: args)
            {
                if(arg.startsWith("-"))
                {
                    if (args.indexOf(arg) == (args.size() - 1))
                    {
                        flags.add(arg.replace("-", ""));
                    }
                    else if (args.get(args.indexOf(arg)+1).startsWith("-"))
                    {
                        flags.add(arg.replace("-", ""));
                    }
                    else
                    {
                        //List of values (can be multiple)
                        List<String> argumentValues = new ArrayList<>();
                        int i = 1;
                        while(args.indexOf(arg)+i != args.size() && !args.get(args.indexOf(arg)+i).startsWith("-"))
                        {
                            argumentValues.add(args.get(args.indexOf(arg)+i));
                            i++;
                        }
                        map.put(arg.replace("-", ""), argumentValues);
                    }
                }
            }
        }
    }
}
