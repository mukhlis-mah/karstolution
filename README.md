# Karstolution
The First Speleothem δ18O Model Integrating Karst Hydrological and In-Cave Fractionation Processes.  
Coupling of existing KarstFor (Bradley et al., 2010) and ISOLTUION (Deininger et al., 2012) models.   
Manuscript to be submitted to GCA.  
A Windows executable GUI has also been created and will freely be available...   
This python code is to allow furhter customisation of the Karstolution model, improvements and allow more advanced model runs such as customisable plotting and batch runs.  

# Conceptual Figure


# Dependancies
Numpy  
Scipy  
(optional) matplotlib  
(optional) pandas  

# Installation
The current easiest way is to just to download and use the files in the Karstolution folder, executing the model from the __init__ file.  
The configuration and input csv files are described below  

# Configuration File
A csv file containing all the model parameters is currently the way the GUI works (example is provided). To allow mixed-use of Karstolution between the python script and GUI, the same configuration file is used. This has to be modified manually using a text editor.  
The following is the format of the configuration csv file (see above conceptual figure):  
F1,F3,F5,F6,F7,k_diffuse,f8  
i,j,k,m,n  
k_eevap, k_d18o_soil, k_d18o_epi  
Cave temps (monthly: Jan-Dec)  
Drip intervals (monthly: Jan-Dec)  
Soilstore, epicap, ovicap, epikarst, KS1, KS2 (store sizes)  
Lambda, k (weibull)  
Drip pCO2 (monthly: Jan-Dec)  
Cave pCO2 (monthly: Jan-Dec)  
Rel Humidity (monthly: Jan-Dec)  
Ventilation (monthly: Jan-Dec)  
Mixing parameter (phi)  
Store initial values (mm): soil, epikarst, ks1, ks2, diffuse  
δ18O intial values (per mille):  soil, epikarst, ks1, ks2, prev rain, diffuse  

# Input file
The input file is a csv of climatic inputs, a similar format to that of KarstFor (example is provided).  
Note: the model steps are in months and the number of rows represents the number of model steps   
The columns are:  
tt: an id column with numbers from 1 to total number of model steps  
mm: representing the month of that model step (important for seasonality), values 1-12  
evpt: evapotranspiration (mm)  
prp: rainfall amount (mm)  
tempp: surface temperature (degree celcius)  
d18O: the δ18O of rainfall amount  

