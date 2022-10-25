Introduction
============

`Eawag enviRule`
-----------------------------

`1. Reaction clustering tool`

`2. Automatic extracting rules from biotransformation reactions and `

`3. Automatic updating rules with additional reactions`

Eawag enviRule is built on top of [Reaction Decoder Tool (rdt)] (https://github.com/asad/ReactionDecoder).

Contact
============
Author: Kunyang Zhang
e-mail: kunyang.zhang@eawag.ch

Atom Atom Mapping using Java API
=================================

View mapped reaction using [CDKDEPICT Tool](http://www.simolecule.com/cdkdepict/depict.html).

```

public static void main(String[] args) throws CloneNotSupportedException, CDKException, AssertionError, Exception {
        final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        String reactionSM = "CC(=O)C=C.CC=CC=C>>CC1CC(CC=C1)C(C)=O";
        String reactionName = "Test";

        IReaction cdkReaction = smilesParser.parseReactionSmiles(reactionSM);

        IReaction performAtomAtomMapping = performAtomAtomMapping(cdkReaction, reactionName);
        System.out.println("AAM sm: " + sg.create(performAtomAtomMapping));
    }

    /**
     *
     * @param cdkReaction
     * @param reactionName
     * @return
     * @throws InvalidSmilesException
     * @throws AssertionError
     * @throws Exception
     */
    public static IReaction performAtomAtomMapping(IReaction cdkReaction, String reactionName) throws InvalidSmilesException, AssertionError, Exception {
        cdkReaction.setID(reactionName);
        /*
         RMT for the reaction mapping
         */
        boolean forceMapping = true;//Overrides any mapping present int the reaction
        boolean generate2D = true;//2D perception of the stereo centers
        boolean generate3D = false;//2D perception of the stereo centers
        StandardizeReaction standardizeReaction = new StandardizeReaction(); //Standardize the reaction
        ReactionMechanismTool rmt = new ReactionMechanismTool(cdkReaction, forceMapping, generate2D, generate3D, standardizeReaction);
        MappingSolution s = rmt.getSelectedSolution();//Fetch the AAM Solution
        IReaction reaction = s.getReaction();//Fetch Mapped Reaction
        return reaction;
    }

```


License
=======

`RDT` is released under the [GNU General Public License version 3](http://www.gnu.org/licenses/gpl.html).

```
Author: Syed Asad Rahman
e-mail: asad@ebi.ac.uk
c/o EMBL-European BioInformatics Institute (EBI)
WTGC, CB10 1SD Hinxton
UK

Note: The copyright of this software belongs to the author
and EMBL-European BioInformatics Institute (EBI).
```

How to Cite RDT?
================

`SA Rahman, G Torrance, L Baldacci, SM Cuesta, F Fenninger, N Gopal, S Choudhary, JW May, GL Holliday, C Steinbeck and JM Thornton: Reaction Decoder Tool (RDT): Extracting Features from Chemical Reactions, Bioinformatics (2016)`

[doi: 10.1093/bioinformatics/btw096](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4920114/)


Sub-commands
===========


`Perform AAM`
-------------

`AAM using SMILES`
  
  ```
  java -jar ReactionDecoder.jar -Q SMI -q "CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O" -g -c -j AAM -f TEXT
  ```

`Perform AAM` for Transporters
-------------

`AAM using SMILES` (accept mapping with no bond changes -b)
  
  ```
  java -jar ReactionDecoder.jar -Q SMI -q "O=C(O)C(N)CC(=O)N.O=C(O)C(N)CS>>C(N)(CC(=O)N)C(=O)O.O=C(O)C(N)CS" -b -g -c -j AAM -f TEXT
  ```
  
`Annotate Reaction using SMILES`
---------------------------------

  ```
  java -jar ReactionDecoder.jar -Q SMI -q "CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O" -g -c -j ANNOTATE -f XML
  ```


`Compare Reactions`
--------------------

`Compare Reactions using SMILES with precomputed AAM mappings`
  
  ```
  java -jar ReactionDecoder.jar -Q RXN -q example/ReactionDecoder_mapped.rxn  -T RXN -t example/ReactionDecoder_mapped.rxn -j COMPARE -f BOTH -u
  ```


`Compare Reactions using RXN files`
  
  ```
  java -jar ReactionDecoder.jar -Q RXN -q example/ReactionDecoder_mapped.rxn  -T RXN -t example/ReactionDecoder_mapped.rxn -j COMPARE -f BOTH
  ```
