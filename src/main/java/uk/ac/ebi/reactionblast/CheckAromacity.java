package uk.ac.ebi.reactionblast;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.*;

public class CheckAromacity {
    public static void main(String[] args) throws Exception{
        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
        String smiles = "[N:1]1[C:2]([NH2:3])[C:4][C:5][N:7]1[CH:8]";

        IAtomContainer mol = smilesParser.parseSmiles(smiles);

        Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.edgeShort()));
        aromaticity.apply(mol);
        System.out.println();

        String s1 = "[O]-[C](:[C]):[C]";
        String s2 = "[O]-[C](-[C])-[C]";

        IAtomContainer mol1 = smilesParser.parseSmiles(s1);
        IAtomContainer mol2 = smilesParser.parseSmiles(s2);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol1);
        InChIGenerator generator1 = factory.getInChIGenerator(mol1);
        System.out.println(generator1.getInchiKey());

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2);
        InChIGenerator generator2 = factory.getInChIGenerator(mol2);
        System.out.println(generator2.getInchiKey());

        int aromatic_bonds1 = 0;
        for(IBond bond: mol1.bonds()){
            if(bond.isAromatic()) aromatic_bonds1 ++ ;
        }

        int aromatic_bonds2 = 0;
        for(IBond bond: mol2.bonds()){
            if(bond.isAromatic()) aromatic_bonds2 ++;
        }

        System.out.println(aromatic_bonds1);
        System.out.println(aromatic_bonds2);

        MACCSFingerprinter fp = new MACCSFingerprinter();
        BitSet b1 = fp.getBitFingerprint(mol1).asBitSet();
        BitSet b2 = fp.getBitFingerprint(mol2).asBitSet();

        System.out.println(b1);
        System.out.println(b2);

        String s3 = "[C:1]1[C:2][C:3][C:4]2[C:5]1[C:6][C:7][C:8][C:9]2";
        rings(s3);


    }

    public static Map<Integer, List<Integer>> rings(String s){

        Map<Integer, List<Integer>> ring_connect_points = new HashMap<>();
        for(int i=1; i<s.length(); i++){
            // Let's assume there is no compounds that have more than 9 rings
            char c = s.charAt(i);
            try{
                int ring_number = Integer.valueOf(c+"");
                // Look back to include the whole []
                if(s.charAt(i-1) != ']') continue;
                int end_mapId = Integer.valueOf(s.charAt(i-2)+"");
                if(!ring_connect_points.containsKey(ring_number)) ring_connect_points.put(ring_number, new ArrayList<>());
                ring_connect_points.get(ring_number).add(end_mapId);
            }catch(Exception e){
                continue;
            }
        }


        return ring_connect_points;
    }

    public static List<String> generalSegments(String smarts){

        List<Character> bracket_stack = new ArrayList<>();
        List<String> segments = new ArrayList<>();

        String segment = "";

        for(int i=0; i<smarts.length(); i++){
            char c = smarts.charAt(i);
            if(c == '['){
                bracket_stack.add('[');
                segment += c;
                continue;
            }
            if(bracket_stack.size()!=0) segment += c;

            if(c == ']'){
                bracket_stack.remove(bracket_stack.size()-1);
            }

            if(bracket_stack.size() == 0 && segment.length() != 0){
                segments.add(segment);
                segment = "";
            }
        }
        return segments;
    }

    public static Map<Integer, String> extractConditionFromSegments(List<String> segments){

        Map<Integer, String> condition_map = new HashMap<>();

        for(String segment: segments){
            String condition = "";
            String mapId = "";
            boolean start = false;
            List<Character> stack = new ArrayList<>();

            for(int i=0; i<segment.length(); i++){
                char c = segment.charAt(i);
                if(c == '$') {
                    start = true;
                    condition += c;
                    continue;
                }
                if(start) condition += c;

                if(c == '(') stack.add(c);

                if(c == ')') {
                    stack.remove(stack.size()-1);
                    if(stack.size() == 0) {
//                    start = false;
                        mapId = segment.substring(i+1);
                        mapId = mapId.replace(":", "").replace("]", "");
                        break;
                    }
                }
//                if(!start){
//                    if(c!= ')' && c!=':' && c!=']') mapId += c;
//                }
            }

            if(condition.length() < 3) System.exit(1);

            try{
                condition_map.put(Integer.valueOf(mapId), condition);
            }catch (Exception e){
                System.exit(1);
            }
        }

        return condition_map;

    }

}
