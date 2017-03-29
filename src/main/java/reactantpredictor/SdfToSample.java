/** 
 * Authors: Siyang Tian
 * Class Description:
 */


package reactantpredictor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Properties;

import org.apache.commons.lang3.StringUtils;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;


import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;

import reactantpredictor.utils.ChemSearcher;
import reactantpredictor.utils.*;


public class SdfToSample {
	/**
	 * Create IAtomContainerSet from inputFiles(either sdf or smiles)
	 * @param String inputPath
	 * @return IAtomContainerSet moleculeSet        
	 * @throws Exception
	 */
	
	public IAtomContainerSet createIAtomContainerSet(String inputPath) throws Exception{
		IAtomContainerSet moleculeSet = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		
		FileReader fr = new FileReader(inputPath);
		BufferedReader br = new BufferedReader(fr);
		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
		SmilesParser sp = new SmilesParser(builder);
		String oneLine;
		//If the input is SMILEs
		if(inputPath.contains(".csv")){
			while((oneLine = br.readLine())!=null){
				//Create IAtomContainer molecule from smiles
				IAtomContainer mol = sp.parseSmiles(oneLine);
				StructureDiagramGenerator sdg = new StructureDiagramGenerator();
				sdg.setMolecule(mol);
				sdg.generateCoordinates();
				IAtomContainer layedOutMol = sdg.getMolecule();
				moleculeSet.addAtomContainer(layedOutMol);
			}

		}
		//if the input file is a sdf file
		else if(inputPath.contains(".sdf")){
			 moleculeSet = readFile(inputPath);
		}
		return moleculeSet;
	}
	
	/**
	 * Create csv file contains Inchikey, title, SMILEs... Predict Results. It's an output csv File
	 * @param IAtomContainerSet moleculeSet
	 * @return void      
	 * @throws Exception
	 */	
	public void outPutCsv(IAtomContainerSet moleculeSet, String outSdfPath) throws Exception{
		IAtomContainerSet resultMole = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		String outputCsv = outSdfPath;
		//FileReader fr = new FileReader(new File(outputCsv));
		FileWriter fw = new FileWriter(new File(outputCsv));
			
		
		for(int i = 0; i < moleculeSet.getAtomContainerCount(); i++){

			
			IAtomContainer oneMole = moleculeSet.getAtomContainer(i);
			//Create proper list
			Map<Object, Object> CypPre = oneMole.getProperties();
			Iterator it = CypPre.entrySet().iterator();
			String attributes = "";
			while(it.hasNext()){
				Map.Entry entry = (Map.Entry) it.next();
				//Wrtie Attribute names as the first row
				if(i == 0){
					attributes = attributes + entry.getKey().toString() + ",";
		
				}
			}
			fw.write(attributes);
			/*
			String attributes = "";
			//Map<Object, Object> entry = (Map<Object, Object>) it.next();
			if(it.hasNext()){
				attributes = attributes + entry.getValue() + ",";
			}
			else attributes = attributes + entry.getValue() + "\n";
			fw.write(attributes);
			*/
				
		}
			fw.close();
	}	
	
	/**
	 * Given a IAtomContainerSet SS, generate Instances with all features for all molecules in SS
	 * @param IAtomContainerSet set
	 * @return weka Instances with all raw feature values         
	 * @throws Exception
	 */
	public Instances generateAllFeatures(IAtomContainerSet set) throws Exception{
		
		 ArrayList<Attribute> atts = new ArrayList<Attribute>();
		 //IAtomContainerSet set = readFile(pathToInputFile);
		 //Add attribute names
		 //String names = "InChiKey\tPubChemID\tHMDB\tDrugBank\t1A2\t2A6\t2B6\t2C8\t2C9\t2C19\t2D6\t2|E1\t3A4\tName\tIsomericSmiles";
		 String moleculeFeatures = "nHBAcc\tnHBDon\tnaAromAtom\tnAtomP\tnB\tnAromBond\tnRotB\tALogP\tALogp2\tAMR\tXLogP\tMLogP\tapol\tTopoPSA\tMW\tbpol\tATSc1\tATSc2\tATSc3\tATSc4\tATSc5\tATSm1\tATSm2\tATSm3\tATSm4\tATSm5\tnAcid\tnBase\tMOMI-X\tMOMI-Y\tMOMI-Z\tMOMI-XY\tMOMI-XZ\tMOMI-YZ\tMOMI-R\tAllSurfaceArea";
		 LinkedHashMap<String, String> fpatterns = ChemSearcher.getRINFingerprintPatterns();
		 String[] labels = fpatterns.keySet().toArray(new String[fpatterns.size()]);
		 
		

		 String rinFPnames = "\t" + StringUtils.join(labels,"\t");
		
		 String firstNames = moleculeFeatures+rinFPnames;
		 String[] fnames = firstNames.split("\t");
		 
		 for(int j = 0; j<fnames.length; j++){
			 fnames[j] = fnames[j].replace(",", "-");
			 Attribute Attribute = new Attribute(fnames[j]);
			 atts.add(Attribute);
		 }
		 for(int h = 0; h < 881; h++){
				
			 Attribute Attribute = new Attribute(String.format("pubchem_f%03d", h+1));
			 atts.add(Attribute);
		 }

		 for(int h = 0; h < 166; h++){
				
			 Attribute Attribute = new Attribute(String.format("maccs_k%03d", h+1));
			 atts.add(Attribute);
		}
		//Add class attribute as weka request this when predicting
		Attribute classAtt = new Attribute("class"); 
		atts.add(classAtt);
		//Names have been added, now add values 
		Instances userinput = new Instances("Rel", atts, 100000);
		int length = atts.size();
		for(int idx = 0; idx<set.getAtomContainerCount();idx++){
		  String result = generateOneinstance(set.getAtomContainer(idx));
		  String[] temp = result.split("\t");
		  Instance sample = new DenseInstance(length); 
		  for(int vidx = 0; vidx < temp.length; vidx++){
			  Attribute att = atts.get(vidx);
			  //att.type();
			  double vle = Double.parseDouble(temp[vidx]);
			  //FastVector attributes = new FastVector();
			  sample.setValue(att, vle);
		  }
		  sample.setValue(classAtt, 0.0);
		  userinput.add(sample);
		}
		
		return userinput;
		
		
	}
	
	
	/**
	 * Given a sdf file, generate IAtomContainerSet that contains all molecules in the sdf file
	 * @param String pathToInputFile---(sdf file)
	 * @return IAtomContainerSet that contains all molecules in the sdf file        
	 * @throws Exception
	 */
	public IAtomContainerSet readFile(String pathToInputFile)
			throws FileNotFoundException, CDKException {
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(pathToInputFile),
				bldr);
		Properties prop = new Properties();
		prop.setProperty("ForceReadAs3DCoordinates", "true");
		PropertiesListener listener = new PropertiesListener(prop);
		sdfr.addChemObjectIOListener(listener);
		sdfr.customizeJob();
		IAtomContainerSet MOLS = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		while (sdfr.hasNext())
				MOLS.addAtomContainer(sdfr.next());
		return MOLS;

	}
	
	/**
	 * Given an IAtomContainer of a molecule and a smile string, add the smile as property of the molecule
	 * @param IAtomContainer molecule, String smiles
	 * @return IAtomContainer that contains smile property    
	 * @throws Exception
	 */
	
	/*
	public IAtomContainer addSmiles(IAtomContainer mole, String smiles) throws Exception {
		//Try to add properties
		HashMap<Object, Object> properties = new HashMap<Object, Object>();
		HashMap<Object, Object> smileProp = new HashMap<Object, Object>();
		smileProp.put("SMILEs", smiles);
		mole.addProperties(smileProp);
		properties.putAll(mole.getProperties());
		//Finish adding properties
		return mole;
	}
	
	*/
	/**
	 * Given an IAtomContainer of a molecule, generate a string that contains all raw feature values for that molecule
	 * @param IAtomContainer molecule
	 * @return IAtomContainerSet that contains all molecules in the sdf file        
	 * @throws Exception
	 */
	public String generateOneinstance(IAtomContainer mole) throws Exception {
		StringBuffer sb = new StringBuffer();
		ChemSearcher cs = new ChemSearcher();
		PubchemFingerprinter pbf 	= new PubchemFingerprinter(SilentChemObjectBuilder.getInstance());
		MACCSFingerprinter maccs 	=  new MACCSFingerprinter(SilentChemObjectBuilder.getInstance());

		LinkedHashMap<String, String> fpatterns = cs.getRINFingerprintPatterns();
		FeatureGeneration fgen = new FeatureGeneration();
		
		IAtomContainer container = mole;
		
	
		IAtomContainer prepContainer = MoleculeExplorer.preprocessContainer(container);
		String[] gg = fgen.generateExtendedMolecularFeatures(prepContainer).split(",");
			
		String extendedFeatures = StringUtils.join(fgen.generateExtendedMolecularFeatures(prepContainer).split(","), "\t");
			

		ArrayList<Double> bioTransformerFingerprint_bits = cs.generateClassyfireFingerprintAsDouble(prepContainer, fpatterns).getBitValues();
		for(int x = 0; x < bioTransformerFingerprint_bits.size(); x++){
			extendedFeatures =  extendedFeatures + "\t" + String.valueOf(bioTransformerFingerprint_bits.get(x));
				
		}
		
			
		ArrayList<Double> fingerprint_bits = new ArrayList<Double>();
		IBitFingerprint fingerp	= pbf.getBitFingerprint(prepContainer);

		int[] onbits = fingerp.getSetbits();

		for(int kp = 0; kp < 881; kp++){
			fingerprint_bits.add(0.0);
		}
		for(int o = 0; o < onbits.length; o++){
			fingerprint_bits.set(onbits[o], 1.0);
		}
		
		extendedFeatures =  extendedFeatures + "\t" + StringUtils.join(fingerprint_bits,"\t");
			
		ArrayList<Double> maccs_fingerprint_bits = new ArrayList<Double>();
		IBitFingerprint maccs_fingerp		= maccs.getBitFingerprint(prepContainer);
			
		int[] maccs_onbits = maccs_fingerp.getSetbits();
			
		for(int kp = 0; kp < 166; kp++){
			maccs_fingerprint_bits.add(0.0);
		}
		for(int o = 0; o < maccs_onbits.length; o++){
			maccs_fingerprint_bits.set(maccs_onbits[o], 1.0);
		}			
		
		extendedFeatures =  extendedFeatures + "\t" + StringUtils.join(maccs_fingerprint_bits,"\t");
		
		String finalFeatureValues = extendedFeatures;
		String[] temp = extendedFeatures.split("\t");
		return finalFeatureValues;
	}
	
}


