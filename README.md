Introduction
============

`Eawag enviRule`
-----------------------------

`1. Reaction clustering tool`

`2. Automatically extracting rules from biotransformation reactions and `

`3. Automatically updating rules with additional reactions`

Eawag enviRule is built on top of [Reaction Decoder Tool (rdt)] (https://github.com/asad/ReactionDecoder).

Contact
============
Author: Kunyang Zhang
e-mail: kunyang.zhang@eawag.ch

###Note: you can find reactions to play with under the folder "example". You can apply automatically extracted rules on reactions with envipath_reaction (https://github.com/zhangky12/envipath_reaction)


Quick start
======================
Make sure you have java installed on the computer, for example:
```markdown
java version "16.0.1" 2021-04-20
Java(TM) SE Runtime Environment (build 16.0.1+9-24)
Java HotSpot(TM) 64-Bit Server VM (build 16.0.1+9-24, mixed mode, sharing)
```
Please always replace "<...>" with your local address in the following examples. The envirule_toy.jar can be found under output/artifacts/envirule_toy/

1. Reaction clustering
```markdown
java -jar envirule_toy.jar -operation cluster -reactions <path to reaction file reactions.txt> -output <choose a folder (end with '/') to store clustered groups. E.g., reactions4rules/>
```
2. Rule generation
```markdown
java -jar envirule_toy.jar -operation generator -group <path to any clustered group from the last step. E.g., 1-2.txt> -radius 1
```
ps: If you see more than one lines of rules, it suggests the extracted rule is a composite rule.
3. Adding reactions
```markdown
java -jar envirule_toy.jar -operation adder -reactions <path to new reaction file reactions_new.txt> -input <path to existing clustered reaction groups (end with '/'). E.g., reactions4rules/> -output <choose a folder (end with '/') to store updated clustered groups. E.g., updated_reactions4rules/>
```
ps: It will return the reaction groups that are updated or created with your new reactions.

To use enviRule in a more flexible way, you can start enviRule server and send commands to it via a client. The following examples basically replicate what we have just done, but through the server & client interactions. 

Start enviRule server
========================
The envirule.jar can be found under output/artifacts/envirule_server/
```
java -jar envirule.jar
```
You will be able to see the following message, which means the server is ready to accept requests.
```
Welcome to enviRule developed by EAWAG!
Waiting for a client ...
```

Connect to enviRule server from a client
==========================================
```
public static void main(String args[]) throws UnknownHostException, IOException{
		
		// Connect to enviRule server
		Socket socket = new Socket("127.0.0.1", 5000);
		System.out.println("Connected");
		ObjectInputStream objectInput = new ObjectInputStream(socket.getInputStream());
		PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
		socket.close();
}
```

Reaction clustering with enviRule
=====================================
Put reactions into a file (e.g., reactions.txt). After running the script below, clustered reaction groups will be saved in a folder (e.g., reactions4rules)

```
public static void main(String args[]) throws UnknownHostException, IOException{
		
		// Connect to enviRule server
		Socket socket = new Socket("127.0.0.1", 5000);
		System.out.println("Connected");
		ObjectInputStream objectInput = new ObjectInputStream(socket.getInputStream());
		PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
		
		// Replace here with your local directory
		String dir = "...";
		
		// Clustering reactions
		// Reactions to be clustered are stored in "reactions.txt"
		String reaction_file = dir + "reactions.txt";
		
		// Clustered reaction groups will be stored under the folder "reactions4rules"
		String output_dir = dir + "reactions4rules/";
		
		try {
			clusteringReactions(reaction_file, output_dir, objectInput, out);
		}catch(Exception e) {
			System.out.println(e.getMessage());
			socket.close();
			return;
		}
		
		socket.close();
}

public static void clusteringReactions(String reaction_file, String output_dir, 
			ObjectInputStream objectInput, PrintWriter out) throws Exception {
		
		String cluster_command = "rxnclust -r " + reaction_file + " -d " + output_dir;
		out.println(cluster_command);
		String result = (String)objectInput.readObject();
		if(result == null) {
			System.out.println("Null response");
			throw new Exception("Clustering failed");
		}else{
			System.out.println(result);
		}
		
}
```

Automatic rule generation with enviRule
==========================================
```
public static void main(String args[]) throws UnknownHostException, IOException{
		
		// Connect to enviRule server
		Socket socket = new Socket("127.0.0.1", 5000);
		System.out.println("Connected");
		ObjectInputStream objectInput = new ObjectInputStream(socket.getInputStream());
		PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
		
		// Replace here with your local directory
		String dir = "...";
		
		// Generating rules	
		Map<String, List<String>> autoRules = new HashMap<>();
		Set<String> files = new HashSet<>();
		
		// Define file names for the clustered reaction groups
		String file_name = dir + "reactions4rules/" + "1-2.txt";
		files.add(file_name);
		file_name = dir + "reactions4rules/" + "2-3.txt";
		files.add(file_name);
		file_name = dir + "reactions4rules/" + "3-5.txt";
		files.add(file_name);
		
		// Set parameters
		boolean ignoreHydrogen = false;
		boolean functionalGroups = true;
		int radius = 1;
		
		try {
			autoRules = generatingRulesforFiles(files, objectInput, out, 
					ignoreHydrogen, functionalGroups, radius);
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		System.out.println(autoRules);
		
		socket.close();
		
}

public static Map<String, List<String>> generatingRulesforFiles(Set<String> file_names, 
			ObjectInputStream objectInput, PrintWriter out, boolean ignoreHydrogen, 
			boolean functionalGroups, int radius) throws ClassNotFoundException, IOException{
		
		Map<String, List<String>> autoRules = new HashMap<>();
		
		
		String includeHydrogen_arg;
		String functionalGroups_arg;
		
		if(ignoreHydrogen) {
			includeHydrogen_arg = "true";
		}else includeHydrogen_arg = "false";
		
		if(functionalGroups) {
			functionalGroups_arg = "true";
		}else functionalGroups_arg = "false";
		
		
		for (String file_name: file_names) {
			
			List<String> rule_smarts = new ArrayList<>();
			
			String command = "autorule -i " + includeHydrogen_arg + " -fg " + 
			functionalGroups_arg + " -f " + file_name + " -r " + radius;
			out.println(command);
			Set<String> rules = (Set<String>)objectInput.readObject();
			if(rules == null || rules.size()==0) {
				System.out.println("Null response");
				continue;
			}
			
			rule_smarts.addAll(rules);
			
			try {
				String[] file_name_segs = file_name.split("/");
				file_name = file_name_segs[file_name_segs.length-1];
				String rule_name_base = "rule-" + file_name.split("-")[0];
				autoRules.put(rule_name_base, rule_smarts);
				
			}catch (Exception e) {
				System.out.println(e.getMessage());
				continue;
			}	
		}
		
		return autoRules;

}
```
From the client side, you will see the following results. rule-1 and rule-2 are simple rules, while rule-3 is a composite rule, consisting of two simple rules
```
rule-1=[[C;$([CH2]-[C]),$([C](=[CH])-[CH]):1]-[N;$([NH2])]>>[CH2:1]-[OH]],
rule-2=[[N;$([N](:[CH]):[C]):1]-[C;$([CH3])]>>[NH:1]], 
rule-3=[[P;$([P](-[O-])(-[O])-[O]):1]=[S;$([S])]>>[P:1]=[O], [P;$([P](-[S])(-[O])-[O]),$([P](-[O-])(-[O-])-[O]),$([P](-[O])(-[O])-[O]):1]=[S;$([S])]>>[P:1]=[O]]
```
Rule updates with enviRule
=============================
The reaction adder module in enviRule server will backup the old clustered reaction groups (e.g., reactions4rules), and either add new reactions into old reaction groups or into new groups depending on the similarity of their reaction fingerprints. enviRule server will return all the reaction groups that are updated (including both changed and newly created groups).
```
public static void main(String args[]) throws UnknownHostException, IOException{
		
		// Replace here with your local directory
		String dir = "...";
		
		// Updating rules
		String new_rxn_file = dir + "reactions_new.txt";
		String old_database = dir + "reactions4rules/";
		String new_database = dir + "updated_reactions4rules/";
		Set<String> changed_rxn_files = new HashSet<>();
		try {
			// Get updated or created files
			changed_rxn_files = addReactions(new_rxn_file, old_database, new_database);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.out.println(changed_rxn_files);
		
}

public static Set<String> addReactions(String new_rxn_file, String old_database, String new_database) throws Exception {
    	
    	// Connect to enviRule server
		Socket socket = new Socket("127.0.0.1", 5000);
		System.out.println("Connected");
    	
		ObjectInputStream objectInput
		= new ObjectInputStream(socket.getInputStream());
		PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
		
		// Adding reactions
		String add_command = "rxnadder -r " + new_rxn_file + " -b " + old_database + " -n " + new_database;
		out.println(add_command);
		// Get updated or newly created files
		Set<String> result = (Set<String>)objectInput.readObject();
		if(result == null || result.size()==0) {
			System.out.println("Null response");
			throw new Exception("Adding reactions failed");
		}
		
		socket.close();
		
		return result;
		
}
```
Using the example files provided here, you will see the following results, which include a file that has been expanded (i.e., 1-2.txt -> 1-3.txt), and one file that is newly created (i.e., 4-3.txt)
```
1-3.txt
4-3.txt
```
It suggests we only need to run rule generation on these two groups, instead of all the clustered groups, to get updated set of rules