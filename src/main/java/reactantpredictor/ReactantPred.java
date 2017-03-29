/** 
 * Authors: Siyang Tian
 * Class Description:
 */

package reactantpredictor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.smiles.SmilesGenerator;

import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffSaver;

public class ReactantPred {
	public static void main(String[] args) throws Exception{
		Scanner reader = new Scanner(System.in);
		System.out.println("Enter the inputPath: ");
		String input = reader.next();
		System.out.println("Enter the outputPath: ");
		String output = reader.next();
		System.out.println("Choose one from the following 9 CYP450 isoforms: 1A2, 2C9, 2B6, 2C8, 2C9, 2C19, 2D6, 2E1, 3A4");
		String cyp = reader.next();
		//System.out.println("Choose output format, input 1 for SDF, 2 for CSV");
		//String outputMode = reader.next();
		//C:/Users/Tian/Desktop/BioData/testDatas/HMDB_selected_endogenous_metabolites_mix.sdf
		//input = "C:/Users/Tian/Desktop/BioData/testDatas/HMDB_endogenous_mini.csv";
		//output = "C:/Users/Tian/Desktop/BioData/0121NewCYP/9CypsDataset/1A2/yannickTest.csv";
		SdfToSample sf = new SdfToSample();
		ReactantPred test = new ReactantPred();
		IAtomContainerSet inputMolecules = sf.createIAtomContainerSet(input);
		//supposed test input: C:\Users\Tian\Desktop\BioData\testDatas\HMDB_endogenous_mini.sdf

		//supposed test input: C:\Users\Tian\Desktop\BioData\testDatas\HMDB_endogenous_mini.csv
		//supposed output path: C:\Users\Tian\Desktop\BioData\0121NewCYP\9CypsDataset\1A2\yannickTest.sdf
		//supposed output path: C:\Users\Tian\Desktop\BioData\0121NewCYP\9CypsDataset\1A2\yannickTest.csv
		
		//String input = "C:/Users/Tian/Desktop/BioData/0121NewCYP/9CypsDataset/1A2/1A2_Reduced.arff";
		String supportfile = "Determined by the input";
		//The predictedResult is used to store all 9 predicted results;
		ArrayList<HashMap<String,String>> predictedResult = new ArrayList<HashMap<String,String>>();
		
		//Create arrayList predictedResult with size = inputMolecules.getAtomContainerCount(), store results fore every molecule
		predictedResult = test.initPreResults(predictedResult,inputMolecules.getAtomContainerCount());
		
		String model1 = "To be choosen";
		String model2 = "To be choosen";
		if(cyp.contains("1A2")){
			model1 = "supportfiles/CYP1A2/model/1A2_RT.model";
			
			supportfile = "supportfiles/CYP1A2/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("1A2",model1,inputMolecules,supportfile,predictedResult );
			//System.out.println("---------------------Model built------------------------");
		}
		else if(cyp.contains("3A4")){
			model1 = "supportfiles/CYP3A4/model/3A4_RT.model";
			supportfile = "supportfiles/CYP3A4/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("3A4",model1,inputMolecules,supportfile,predictedResult );
			//System.out.println("---------------------Model built------------------------");
		}
		else if(cyp.contains("2C9")){
			model1 = "supportfiles/CYP2C9/model/2C9_RT.model";
			
			supportfile = "supportfiles/CYP2C9/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2C9",model1,inputMolecules,supportfile,predictedResult );
			//System.out.println("---------------------Model built------------------------");
		}
		else if(cyp.contains("2C19")){
			model1 = "supportfiles/CYP2C19/model/2C19_RT.model";
			supportfile = "supportfiles/CYP2C19/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2C19",model1,inputMolecules,supportfile,predictedResult );
			//System.out.println("---------------------Model built------------------------");
		}
		else if(cyp.contains("2E1")){
			model1 = "supportfiles/CYP2E1/model/2E1_RT.model";
			
			supportfile = "supportfiles/CYP2E1/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2E1",model1,inputMolecules,supportfile,predictedResult );
			//System.out.println("---------------------Model built------------------------");
		}
		else if(cyp.contains("2D6")){
			model1 = "supportfiles/CYP2D6/model/2D6_RT.model";

			supportfile = "supportfiles/CYP2D6/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2D6",model1,inputMolecules,supportfile,predictedResult );
			//System.out.println("---------------------Model built------------------------");
		}
		else if(cyp.contains("2A6")){
			model1 = "supportfiles/CYP2A6/model/2A6_RT.model";
			supportfile = "supportfiles/CYP2A6/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2A6",model1,inputMolecules,supportfile,predictedResult );
			//System.out.println("---------------------Model built------------------------");
		}
		else if(cyp.contains("2B6")){
			model1 = "supportfiles/CYP2B6/model/2B6_RT.model";
			supportfile = "supportfiles/CYP2B6/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2B6",model1,inputMolecules,supportfile,predictedResult);
			//System.out.println("---------------------Model built------------------------");
		}
		else if(cyp.contains("2C8")){
			model1 = "supportfiles/CYP2C8/model/2C8_RT.model";
			
			supportfile = "supportfiles/CYP2C8/supportfile.csv";
			//predict and store the results into the predictedResult arrayList
			predictedResult = test.makePrediction("2C8",model1,inputMolecules,supportfile,predictedResult );
			//System.out.println("---------------------Model built------------------------");
		}
		else{
			System.out.println("Be added soon");
			return;
		}
		if(output.contains("sdf")){
			test.outputResultIAtomContainerSet(inputMolecules, output, predictedResult);
		}
		else if(output.contains("csv")){
			test.outPutCsv(inputMolecules, output, predictedResult);
		}
		else System.out.println("Cannot detect the output file format");
		//test.outputResultIAtomContainerSet(inputMolecules, output, predictedResult);
		//test.outPutCsv(inputMolecules, output, predictedResult);
		//String resultPath = "src/main/java/rinpredictor/supportfiles/1A2_predict.csv";
		
		
	}
	
	/**
	 * Predict reactants, inhibitors and non-RIs for the given test molecules in the sdf file
	 * @param String model1(first layer model),String model2(second layer model),String testfile, string supportfileString resultPath(where to store the result)
	 * @return weka Instances with all raw feature values         
	 * @throws Exception
	 */
	public ArrayList<HashMap<String,String>> makePrediction(String cyp, String model1, IAtomContainerSet inputMolecules, String supportfile, ArrayList<HashMap<String,String>> predictedResult) throws Exception{
		
		ArrayList<HashMap<String,String>> resultList = predictedResult;
		SdfToSample sdfTool = new SdfToSample();
		//Instances testSet = sdfTool.generateAllFeatures("C:/Users/Tian/Desktop/BioData/1218Subset/Yannick_New.sdf");
		//Instances testSet = sdfTool.generateAllFeatures("C:/Users/Tian/Desktop/BioData/0117NewCyp/1A2_ALL_N_SMILEs.sdf");
		
		Instances testSet = sdfTool.generateAllFeatures(inputMolecules);
		//Instances testSet = sdfTool.generateAllFeatures("C:/Users/Tian/Desktop/BioData/1218Subset/yannickTest.sdf");
		
		Classifier cls = (Classifier) weka.core.SerializationHelper.read(model1);
		
		
		File supfile = new File(supportfile);
		FileReader spfr = new FileReader(supfile);
		BufferedReader spbr = new BufferedReader(spfr);
		
		//Create attribute arraylist,mean arraylist, max arraylist, min arraylist
		ArrayList<String> attList = new ArrayList<String>();
		ArrayList<String> meanList = new ArrayList<String>();
		ArrayList<String> maxList = new ArrayList<String>();
		ArrayList<String> minList = new ArrayList<String>();
		//Write the output 
		int counter = 0;

		//attribute	mean	max	min. Skip the first line that contains all the titles
		String supLine = spbr.readLine();
		while((supLine = spbr.readLine())!=null){
			String[] elements = supLine.split(",");
			attList.add(elements[0]);
			meanList.add(elements[1]);
			maxList.add(elements[2]);
			minList.add(elements[3]);
		}
		int countR = 0;
		int countT = 0;
		//Test section Below
		/*
		String input = "C:/Users/Tian/Desktop/BioData/0117NewCyp/1A2_N_Test.arff";
		File f = new File(input);
		FileReader fr = new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		Instances wholeData = new Instances(br);
		wholeData.setClassIndex(wholeData.numAttributes()-1);
		*/
		//Test secion above
		Instances matched = matchAttributes(testSet,attList,meanList,maxList,minList);
		//
		//Instances matched = wholeData;
		matched.setClassIndex(matched.numAttributes()-1);
		
		//double result = cls.classifyInstance(matched.get(0));
		for(int i = 0; i<matched.size();i++){
			counter++;
			Instance oneSample = matched.get(i);
			double result = cls.classifyInstance(oneSample);
			String preRIN = "";
		
			if(result==0.0){
				countT++;
				System.out.println("T");
				preRIN = "T";
			}
			else if(result == 1.0){
				countR++;
				System.out.println("R");
				preRIN = "R";
			}
			
	
			
			resultList.get(i).replace(cyp,preRIN);
		}
		//outputWriter.close();
		System.out.println("R: " + countR + " T: " + countT);
		return resultList;
	}
	
	/**
	 * output predicted results for given molecules as a sdf file
	 * @param IAtomContainerSet rawMolecules, String outSdfPath, ArrayList<HashMap<String,String>> predictedResult
	 * @return void      
	 * @throws Exception
	 */	
	public void outputResultIAtomContainerSet(IAtomContainerSet rawMolecules, String outSdfPath, ArrayList<HashMap<String,String>> predictedResult) throws Exception{
		IAtomContainerSet resultMole = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		String sdfOutput = outSdfPath;
		SDFWriter sdfWriter  = new SDFWriter(new FileWriter(sdfOutput)); 
		
		IAtomContainerSet resultIAtomContainerSet = getResultIAtomContainerSet(rawMolecules, predictedResult);
		for(int i = 0; i < resultIAtomContainerSet.getAtomContainerCount(); i++){
			IAtomContainer outMole = resultIAtomContainerSet.getAtomContainer(i);
			sdfWriter.write(outMole);
		}
		sdfWriter.close();
	}
	
	/**
	 * Return IAtomContainerSet that contains all predicted results and molecules properties
	 * @param IAtomContainerSet rawMolecules, ArrayList<HashMap<String,String>> predictedResult
	 * @return IAtomContainerSet resultMolecules       
	 * @throws Exception
	 */	
	public IAtomContainerSet getResultIAtomContainerSet(IAtomContainerSet rawMolecules, ArrayList<HashMap<String,String>> predictedResult) throws Exception{
		IAtomContainerSet resultMole = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		for(int i = 0; i < rawMolecules.getAtomContainerCount(); i++){
			IAtomContainer oneMole = rawMolecules.getAtomContainer(i);
			IAtomContainer outMole = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
			outMole = oneMole.clone();
			outMole.setProperties(null);
			//Create proper list
			Map<Object, Object> CypPre = new LinkedHashMap<Object, Object>();
			//Add inchikey
			// Generate factory - throws CDKException if native code does not load
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			// Get InChIGenerator
			InChIGenerator gen1 = factory.getInChIGenerator(oneMole);
			String inchikey1 = gen1.getInchiKey();
			CypPre.put("InchiKey", inchikey1);
			
			//Add SMILEs
			SmilesGenerator sg = new SmilesGenerator().unique();
			String smile1 = sg.create(oneMole);
			CypPre.put("SMILES", smile1);
			
			//Add Title
			String title1 = oneMole.getProperty(CDKConstants.TITLE); 
			//HashMap<Object, Object> title = new HashMap<Object, Object>();
			/*
			if(i == 5){
				title1 = "";
			}
			if(i == 7){
				title1 = null;
			}
			*/
			/*
			if(i == 1){
				inchikey1 = null;
				CypPre.put("InchiKey", inchikey1);
			}
			*/
			if(title1 != null && !title1.isEmpty()){
				CypPre.put("cdk:Title", title1);
			}
			else {
				title1 = "molecule:" + (i+1);
				CypPre.put("cdk:Title", title1);
			}
				
			//outMole.addProperties(title);
			/*
			 * Comment the section below in order to generate sdf files from smiles
			 */
			
			CypPre.put("1A2", predictedResult.get(i).get("1A2"));
			CypPre.put("2B6", predictedResult.get(i).get("2B6"));
			CypPre.put("2A6", predictedResult.get(i).get("2A6"));
			CypPre.put("2C8", predictedResult.get(i).get("2C8"));
			CypPre.put("2C9", predictedResult.get(i).get("2C9"));
			CypPre.put("2C19",predictedResult.get(i).get("2C19"));
			CypPre.put("2D6", predictedResult.get(i).get("2D6"));
			CypPre.put("2E1", predictedResult.get(i).get("2E1"));
			CypPre.put("3A4", predictedResult.get(i).get("3A4"));
			
			
			//CypPre.put("InchiKey", inchikey1);
			Map<Object, Object> oldProperties = oneMole.getProperties();

			//Filter out the null mappings from the oldProperties
			HashMap<Object,Object> filteredProperties = new HashMap<Object,Object>();
			for(Object key : oldProperties.keySet()){
				//String oldValue = (String) oldProperties.get(key);
				if(oldProperties.get(key)!=null){
					filteredProperties.put(key, oldProperties.get(key));
				}
			}
			oldProperties = filteredProperties;
			// The section below checks whether title, Smiles and InchiKey from the origin sdf file are the same as the one generated by our tool
			if(oldProperties.keySet().contains("cdk:Title")){
				String checkTitle = (String) oldProperties.get("cdk:Title");
				if(!checkTitle.equals(title1)){
					System.out.println("The " + i + "th mole's title is different from the title generated by our tool");
					System.out.println("Origin title: " + checkTitle);
					System.out.println("OurGen title: " + title1);
					//oldProperties.remove("cdk:Title", null);
					oldProperties.put("cdk:Title",title1);
					
				}
			}
			if(oldProperties.keySet().contains("InchiKey")){
				String checkInchikey = (String) oldProperties.get("InchiKey");
				if(!checkInchikey.equals(inchikey1)){
					System.out.println("The " + i + "th mole's InchiKey is different from the InchiKey generated by our tool");
					System.out.println("Origin title: " + checkInchikey);
					System.out.println("OurGen title: " + inchikey1);
					oldProperties.put("InchiKey",inchikey1);
				}
			}
			if(oldProperties.keySet().contains("SMILES")){
				String checkSmiles = (String) oldProperties.get("SMILES");
				if(!checkSmiles.equals(smile1)){
					System.out.println("The " + i + "th mole's InchiKey is different from the InchiKey generated by our tool");
					System.out.println("Origin title: " + checkSmiles);
					System.out.println("OurGen title: " + smile1);
					oldProperties.put("SMILES",smile1);
				}
			}
			
			
			CypPre.putAll(oldProperties);
			outMole.addProperties(CypPre);
			
			////System.out.println(oneMole.getProperty("InchiKey"));
			//Put it into the resultIAtomContainerSet
			resultMole.addAtomContainer(outMole);
			
			//write the molecule into sdf file for test purpose
		}
		return resultMole;
	}
	
	/**
	 * Create csv file contains Inchikey, title, SMILEs... Predict Results. It's an output csv File
	 * @param IAtomContainerSet moleculeSet
	 * @return void      
	 * @throws Exception
	 */	
	public void outPutCsv(IAtomContainerSet rawMolecules, String outPath, ArrayList<HashMap<String,String>> predictedResult) throws Exception{
		IAtomContainerSet resultMole = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		String csvOutput = outPath;
		FileWriter csvWriter  = new FileWriter(new File(csvOutput));
		//Write the first line: Inchiky, SMILEs, Titile/MoleCounter, PredictedResults
		String titleLine = "Inchiky, SMILES, Title, 1A2, 2B6, 2A6, 2C8, 2C9, 2C19, 2D6, 2|E1, 3A4\n";
		csvWriter.write(titleLine);
		for(int i = 0; i < rawMolecules.getAtomContainerCount(); i++){
			
			IAtomContainer oneMole = rawMolecules.getAtomContainer(i);
			IAtomContainer outMole = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
			outMole = oneMole.clone();
			outMole.setProperties(null);
			//Create proper list
			Map<Object, Object> CypPre = new LinkedHashMap<Object, Object>();
			//Add inchikey
			// Generate factory - throws CDKException if native code does not load
			InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
			// Get InChIGenerator
			InChIGenerator gen1 = factory.getInChIGenerator(oneMole);
			String inchikey1 = gen1.getInchiKey();
					
			//Add SMILEs
			SmilesGenerator sg = new SmilesGenerator().unique();
			String smile1 = sg.create(oneMole);
			
			//Add Title
			String title1orCounter = "";
			String title1 = oneMole.getProperty(CDKConstants.TITLE); 
			//HashMap<Object, Object> title = new HashMap<Object, Object>();
			/*
			if(i == 5){
				title1 = "";
			}
			if(i == 7){
				title1 = null;
			}
			*/
			if(title1 != null && !title1.isEmpty()){
				title1orCounter =  title1;
			}
			else title1orCounter = "molecule:" + (i+1);
			
			CypPre.put("1A2", predictedResult.get(i).get("1A2"));
			CypPre.put("2B6", predictedResult.get(i).get("2B6"));
			CypPre.put("2A6", predictedResult.get(i).get("2A6"));
			CypPre.put("2C8", predictedResult.get(i).get("2C8"));
			CypPre.put("2C9", predictedResult.get(i).get("2C9"));
			CypPre.put("2C19",predictedResult.get(i).get("2C19"));
			CypPre.put("2D6", predictedResult.get(i).get("2D6"));
			CypPre.put("2E1", predictedResult.get(i).get("2E1"));
			CypPre.put("3A4", predictedResult.get(i).get("3A4"));
			
			//CypPre.put("InchiKey", inchikey1);
			if(title1orCounter.contains(",")){
				title1orCounter = "\"" + title1orCounter + "\"";
			}
			String oneSampleLine = inchikey1 + "," + smile1 + "," + title1orCounter + 
					"," + predictedResult.get(i).get("1A2") + "," + predictedResult.get(i).get("2B6") + "," + predictedResult.get(i).get("2A6") + 
					"," + predictedResult.get(i).get("2C8") + "," + predictedResult.get(i).get("2C9") + "," + predictedResult.get(i).get("2C19") +
					"," + predictedResult.get(i).get("2D6") + "," + predictedResult.get(i).get("2|E1") + "," + predictedResult.get(i).get("3A4") + "\n"; 
			////System.out.println(oneMole.getProperty("InchiKey"));
			//Put it into the resultIAtomContainerSet
			
			//write the molecule into sdf file for test purpose
			csvWriter.write(oneSampleLine);
			
		}

		csvWriter.close();
	}		
	/**
	 * Apply feature selecion and normalization on the given molecule
	 * @param Instances testSet, ArrayList<String> attList, ArrayList<String> meanList, ArrayList<String> maxList, ArrayList<String> minList
	 * @return weka Instances with all raw feature values         
	 * @throws Exception
	 */
	public Instances matchAttributes(Instances testSet, ArrayList<String> attList, ArrayList<String> meanList, ArrayList<String> maxList, ArrayList<String> minList) throws Exception{
		
		List my_nominal_values = new ArrayList(3); 
		my_nominal_values.add("R"); 
		my_nominal_values.add("T"); 
		Attribute lastAttribute  = new Attribute(attList.get(attList.size()-1), my_nominal_values);
		// Create nominal attribute "position" 
		
		ArrayList<Attribute> atts = new ArrayList<Attribute>();
		for(int i = 0; i<attList.size();i++){
			if(i<attList.size()-1){
				Attribute Attribute = new Attribute(attList.get(i));
				atts.add(Attribute);
			}
			else{
				
				atts.add(lastAttribute);
			}
		}
		
		int length = testSet.size();
		int numAttributes = attList.size();
		Instances matched = new Instances("Rel", atts, length);
		Instance sample = new DenseInstance(numAttributes);
		
		for(int j = 0; j< length; j++){
			ArrayList misFeatures = new ArrayList();
			Instance temp = testSet.get(j);
			for(int vidx = 0; vidx < numAttributes; vidx++){
				  Attribute att = atts.get(vidx);
				  if(vidx<numAttributes-1){
					  
					  String stAtt = att.toString();
					  String[] whatever = stAtt.split(" ");
					  
					  if(whatever.length!=3){
						  String combineSpaces = "";
						  for(int kk = 1; kk<whatever.length-1; kk++){
							  if(kk<whatever.length - 2){
								  combineSpaces = combineSpaces + whatever[kk]+ " ";
							  }
							  else combineSpaces = combineSpaces + whatever[kk];
						  }
						  whatever[1] = combineSpaces;
					  }
					 
					  if(whatever[1].contains("'")){
						  whatever[1] = whatever[1].replace("'", "");
					  }
					  Attribute convAtt = testSet.attribute(whatever[1]);
					  double vle = temp.value(convAtt);
					  if(!Double.isNaN(vle)){
						  double min = Double.parseDouble(minList.get(vidx));
						  double max = Double.parseDouble(maxList.get(vidx));
						  double normedVal = (vle-min)/(max-min);
						  /*
						  if(normedVal > max){
							  throw new Exception("Attribute value error");
						  }
						  */
						  sample.setValue(att, normedVal);
					  }
					  else{
						  //Add the name of the missing feature into the misFeatures ArrayList
						  misFeatures.add(whatever[1]);
						  //System.out.println("Molecule--" + j + "contains missing feature, synthesized");
						  double syn = Double.parseDouble(meanList.get(vidx));
						  sample.setValue(att, syn);

					  }
				  }
				  else{
					  //Set class with random value. Does not matter
					  sample.setValue(lastAttribute, "T");
					  
				  }
			}
			//Output those missing features
			if(misFeatures.size()>0){
				int numMisFeatures = misFeatures.size();
				String misFeatureString = "Molecule" + j + ": " +  numMisFeatures + 
						" features could not be computed from the structure representation, and were thus synthesized: " + misFeatures.get(0);
				for(int i = 1; i < numMisFeatures; i++){
					misFeatureString = misFeatureString + "," + misFeatures.get(i);
				}
				System.out.println(misFeatureString);
			}
			
			
			matched.add(sample);
		}
		
		return matched;
		
	}
	public ArrayList<HashMap<String,String>> initPreResults(ArrayList<HashMap<String,String>> predictedResult, int numOfMoles){
		ArrayList<HashMap<String,String>> afterInit = new ArrayList<HashMap<String,String>>();
		for(int i = 0; i < numOfMoles; i++){
			HashMap<String,String> hm = new HashMap<String,String>();
			hm.put("1A2", "null");
			hm.put("2A6", "null");
			hm.put("2B6", "null");
			hm.put("2C8", "null");
			hm.put("2C9", "null");
			hm.put("2C19", "null");
			hm.put("2D6", "null");
			hm.put("2E1", "null");
			hm.put("3A4", "null");
			afterInit.add(hm);
			
		}
		
		return afterInit;
	}

}
