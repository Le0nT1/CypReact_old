/**
 * Authors: Yannick Djoumbou Feunang
 * Class Description:
 */


package reactantpredictor.utils;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class MoleculeExplorer {

	
	/**
	 * This function applies some preprocessing operations, such as setting the
	 * flag of atoms from aromatic rings to "ISAROMATIC", and kelulizing
	 * molecules.
	 * 
	 * @param molecule
	 *            : A molecule of interest
	 * @return : A processed molecule (AtomContainer)
	 * @throws CDKException
	 */
	public static IAtomContainer preprocessContainer(IAtomContainer molecule)
			throws CDKException {
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder.getInstance(molecule.getBuilder()).addImplicitHydrogens(molecule);
		 
	    Aromaticity aromaticity = new Aromaticity( ElectronDonation.cdk(), Cycles.cdkAromaticSet());
//		Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.all());
		
		for (IBond bond : molecule.bonds()) {
			if (bond.getFlag(CDKConstants.SINGLE_OR_DOUBLE)) {
				bond.setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(0).setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(1).setFlag(CDKConstants.ISAROMATIC, true);

			} 
		}
//		aromaticity.apply(molecule);

		Kekulization.kekulize(molecule);
		
		
		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(molecule);
		sdg.generateCoordinates();
		
		IAtomContainer layedOutMol = sdg.getMolecule();

//		StringWriter w2 = new StringWriter();
//		MDLWriter mw2 = new MDLWriter(w2);
//		mw2.write(layedOutMol);	
//		System.out.println("After preprocessing\n" + w2.toString() + "\n\n");

		return layedOutMol;
	}

	
}
