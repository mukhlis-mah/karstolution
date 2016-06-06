import csv
import karst_process

#this function unpacks and initialses some of the model parameters, reads the input file
#iterates each step (according to each entry of the input file)
#and writes the output into the defined file at each step

def karstolution(data,input,output):
    #unpacking inital values taken from the configuration file (in list 'data')
    epx18oxp=data[13][1]        #inital d18O in epikarst
    epxstorxp=data[12][1]       #inital water level of epikarst
    soilstorxp=data[12][0]      #inital d18O in soil store
    soil18oxp=data[13][0]       #inital water level of  soil store
    kststor1xp=data[12][2]      #inital water level of  soil store
    kststor2xp=data[12][3]      #inital water level of w soil store
    kststor118oxp=data[13][2]   #inital d18O in soil store
    kststor218oxp=data[13][3]   #inital d18O in soil store
    d18oxp=data[13][4]          #inital d18O in soil store
    prpxp=0                     #initial value for rain from 'previous' step 
    #weibull distribution (diffuse flow) currently set to 12 months, can be increased
    dpdf=[data[12][4]]*12       #inital water quantity in the weibull distribution (diffuse flow)
    epdf=[data[13][5]]*12       #inital d18O in the weibull distribution (diffuse flow)
    
    #the rest of the data which we are not unpacking (but will do so later)
    data_rest=[thing for thing in data[0:12]]
    
    #Do you want to Karstolution to model drip interval from the level of karst stores? 
    #default is yes. If no, it will just use repeat the monthly drip int averages from the config file
    calculate_drip=True
    
    #36 month surface temp list for purposes of coupling surface to cave; currently a dummy list
    tempp=[thing for thing in range(0,36)]
    #finding the average inputted cave value
    avr_cave=sum(data_rest[3])/len(data_rest[3])
    #setting a dummy value that will be overwritten below
    difference=10
    #getting the output file ready for writing
    csvoutput=open(output,'w')
    writer=csv.writer(csvoutput,lineterminator='\n')
    #name of the headings of the output file
    writer.writerow(['tt','mm','f1','f3','f4','f5','f6','f7','soilstor','epxstor',
    'kststor1','kststor2','soil18o','epx18o','kststor118o','kststor218o','dpdf[0]',
    'stal1d18o','stal2d18o','stal3d18o','stal4d18o','stal5d18o','drip_int_stal1',
    'drip_int_stal4','drip_int_stal3','drip_int_stal2','drip_int_stal5','cave_temp'])
    #reading the input file and using each row as one iteration of the model
    with open(input) as csvinput:
        reader=csv.DictReader(csvinput)
        for row in reader:
            #using the headings of the columns
            tt=int(row['tt']) #step number (starts at 1 and ends at # of iterations)
            mm=int(row['mm']) #month: varies from 1-12
            evpt=float(row['evpt']) #value of evapotranspiration (mm)
            prp=float(row['prp'])  #value of precipitation (mm)
            #update the first value in the list with new input value
            tempp[0]=float(row['tempp']) #surface temperature (celcius)
            d18o=float(row['d18o']) #d18O value of rainfall (per mille)
            
            #for the first loop of the program filling tempp with the first input value
            if tt==1:
                tempp=[tempp[0] for thing in range(tt-1,36)]
                #the difference between the surface temp and cave temp
                #this value determines the set diff for rest of model
                difference=tempp[0]-avr_cave
            #seasonlity factor for that month based on GUI cave temp inputs 
            seasonality=data_rest[3][mm-1]-avr_cave
            
            #average surface temp of last 36 months
            avr_surfacet=sum(tempp)/len(tempp)
            #cave temp is the average surface temp - the set surface-cave temp difference from step 1
            #with an adjustment for seasonality
            cave_temp=avr_surfacet-difference+seasonality
            
            #passes each of the paramters through to the karst_process function...
            out=karst_process.karst_process(tt,mm,evpt,prp,prpxp,tempp,d18o,d18oxp,dpdf,epdf,soilstorxp,
            soil18oxp,epxstorxp,epx18oxp,kststor1xp,kststor118oxp,kststor2xp,kststor218oxp,data_rest,
            calculate_drip, cave_temp)
            
            #output data for this timestep is taken and writen to the csv output file    
            #in order: tt,mm,f1,f3,f4,f5,f6,f7,soilstor,epxstor,kststor1,kststor2, soil18o,epx18o,
            #kststor118o,kststor218o,dpdf[0],stal1d18o,stal2d18o,stal3d18o,stal4d18o,stal5d18o,
            #drip_interval_ks2,drip_interval_epi,drip_interval_stal3, drip_interval_stal2,
            #drip_interval_ks1,cave_temp
            writer.writerow([out[i] for i in range(0,28)])

            
            #update model terms for next iteration
            epx18oxp=out[13]
            epxstorxp=out[9]
            soilstorxp=out[8]
            soil18oxp=out[12]
            kststor1xp=out[10]
            kststor2xp=out[11]
            kststor118oxp=out[14]
            kststor218oxp=out[15]
            epdf[1:11]=epdf[0:10]
            epdf[0]=out[13]
            dpdf[1:11]=dpdf[0:10]
            dpdf[0]=out[16]
            d18oxp=d18o
            prpxp=prp
            tempp[1:36]=tempp[0:35]
    #once model is done returns a True so we know its all G
    return True