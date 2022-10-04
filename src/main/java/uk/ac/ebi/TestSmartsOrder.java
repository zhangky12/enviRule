package uk.ac.ebi;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class TestSmartsOrder {

    public static void main(String[] args) throws Exception{

        final SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);

        IAtomContainer mol = new AtomContainer();

        SmilesParser   sp  = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer m   = sp.parseSmiles("[C:1][O:2][N:3][P:4]");

        Map<Integer, IAtom> atom_map = new HashMap<Integer, IAtom>();

        for(int i=0; i<m.getAtomCount(); i++){
            IAtom atom = m.getAtom(i);
            System.out.println(i+": " + atom.getSymbol());
            atom_map.put(i, atom);
            //mol.addAtom(atom);
        }

        mol.addAtom(atom_map.get(1));
        mol.addAtom(atom_map.get(3));
        mol.addAtom(atom_map.get(0));
        mol.addAtom(atom_map.get(2));

//        for(IBond bond: m.bonds()){
//            mol.addBond(bond);
//        }
        mol.addBond(m.getBond(0));
        mol.addBond(m.getBond(2));

        System.out.println("SMILES is: " + sg.create(mol));
    }
}
