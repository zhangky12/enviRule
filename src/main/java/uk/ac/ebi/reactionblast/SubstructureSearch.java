package uk.ac.ebi.reactionblast;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Mapping;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmilesParser;

public class SubstructureSearch {
    public static void main(String[] args) throws Exception{
//        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//        IAtomContainer atomContainer = sp.parseSmiles("[C:1][C:2](=[O:3])[O:4][C:5](=[O:6])[C:7]");
//        Pattern ptrn = SmartsPattern.create("O=CO");
//        Mappings res = ptrn.matchAll(atomContainer);
//
//        for(int[] list: res) {
//            System.out.println("-----------------------------");
//            for(int i:list){
//                IAtom atom = atomContainer.getAtom(i);
//                System.out.println(i + "; " + atom.getMapIdx());
//            }
//
//        }
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        //[O:1]=[C:2]([O-:3])[CH:4]=[CH:5][C:6]([Cl:11])=[CH:7][C:8](=[O:9])[O-:10]>>[CH:4]1=[CH:5][C:6]([O:3][C:2]1=[O:1])=[CH:7][C:8]([O-:10])=[O:9]
        IAtomContainer atomContainer = sp.parseSmiles("[O:1]=[C:2]([O-:3])[CH:4]=[CH:5][C:6]([Cl:11])=[CH:7][C:8](=[O:9])[O-:10].[CH:4]1=[CH:5][C:6]([O:3][C:2]1=[O:1])=[CH:7][C:8]([O-:10])=[O:9]");
        //[O;$([O-]-[C](=[O])-[C]=[C]-[C]=[C]-[C](=[O])-[O-]):1]([H]).[C;$([C](=[C]-[C](=[O])-[O-])-[C]=[C]-[C](=[O])-[O-]):2]-[Cl;$([Cl])]>>[O:1]-[C:2]
        Pattern ptrn1 = SmartsPattern.create("[O;$([O-]-[C](=[O])-[C]=[C]-[C]=[C]-[C](=[O])-[O-]):1]");
//        Pattern ptrn1 = SmartsPattern.create("[O;$([O-]-[C](=[O])-[C]=[C](-[Cl])-[C]=[C](-[Cl])-[C](=[O])-[O-]):1]([H])");
        Pattern ptrn2 = SmartsPattern.create("[C;$([C](=[C](-[Cl])-[C](=[O])-[O-])-[C](=[C]-[C](=[O])-[O-])-[Cl]):2]-[Cl;$([Cl])]");

        Mappings res = ptrn1.matchAll(atomContainer);

        for(int[] list: res){
            System.out.println("-----------------------------");
            for(int i:list){
                IAtom atom = atomContainer.getAtom(i);
                System.out.println(i + "; " + atom.getMapIdx());
            }
        }

    }
}
