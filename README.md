Please download the CypReactBundle.zip to run CypReact tool.
The bundle includes: 1. A runnable .jar file; 2. Supporting files(.model files and supporfile.csv files).
To run the CypReact tool, the user should use command in the terminal as:

->java -jar "PathOfCypReactBundle" "PathOfInputMoles" "PathOfOutputFile" "CYP/CYPs"
The user can either input a .sdf file or a .csv. If the user input a .csv file, it should contains the SMILEs of all molecules and split them with ",". 
The user can output a .sdf file or a .csv file.
For example, Assume CypReactBundle folder is extracted onto the desktop, the user can run CypReactBundle as:

C:\Users\Desktop\CypReactBundle>java -jar cypreact.jar C:\Users\Desktop\CypReactBundle\ C:\Users\Desktop\BioData\1A2_Test.sdf C:\Users\Desktop\BioData\1A2_Result.sdf 1A2

The user can also input multiple CYPs. For example, if the user wants to test his/her molecules on CYP1A2,2A6 and 2B6, he/she can simply replace "1A2" with "1A2,2A6,2B6". Please leave no space between "," and a CYP.

The CypReactBundle has been tested on Windows, Linux and macOS.

CypReact is developed jdk 1.8.0 with Eclipse Mars.1 Release (4.5.1).
The user can also download the repository and run it on his/her local machine.

There was a bug when the user run CypReactBundle for CYP2D6 and CYP2C9. 
The bug has been fixed and the current version of CypReactBund is 1.1.
