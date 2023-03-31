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

Start enviRule server
========================
```
java -jar enviRule.jar
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
		
		// Connect to enviRule server
		Socket socket = new Socket("127.0.0.1", 5000);
		System.out.println("Connected");
		ObjectInputStream objectInput = new ObjectInputStream(socket.getInputStream());
		PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
		
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
		
		socket.close();
		
}
```
Using the example files provided here, you will see the following results, which include a file that has been expanded (i.e., 1-2.txt -> 1-3.txt), and one file that is newly created (i.e., 4-3.txt)
```
1-3.txt
4-3.txt
```
It suggests we only need to run rule generation on these two groups, instead of all the clustered groups, to get updated set of rules