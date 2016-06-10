from karstolution1_1 import karstolution
import csv
# import pandas as pd
# import matplotlib.pyplot as plt


####################################################
#__init__ file for the Karstolution speleothem d18O module
#
####################################################

#in order to be consistent with the GUI, a configuration file (CSV) is used
#to store all the model parameters. A full description is given in the readme
#name of the configuration csv file
config='config.csv'

#name of the input (csv) file. A full description is given in the readme
input='input.csv'

#output (csv) file name and location. A full description is given in the readme
output='output.csv'

#reading the config file
csvinput=open(config, 'rU')
reader=csv.reader(csvinput)
#all paramters from config file stores in list 'data'
data=[]
for row in reader:
    if row:
        data.append([])
        for thing in row:
            data[-1].append(float(thing))

#Do you want to Karstolution to model drip interval from the level of karst stores? 
#default is yes. If no, it will just use repeat the monthly drip int averages from the config file
calculate_drip=True
            
#running the Karstolution model
PSM=karstolution(data,input,output,calculate_drip)

if PSM:
    print "Finished running Karstolution!"
else:
    print "something happened?!"
    
    
#additional module that will uses pandas and matplot lib to plot the output
# df = pd.read_csv(output)

# time = df.tt
# stal1=df.stal1d18o
# stal2=df.stal2d18o
# stal3=df.stal3d18o    
# stal4=df.stal4d18o
# stal5=df.stal5d18o
# drip1=df.drip_int_stal1
# drip2=df.drip_int_stal2
# drip3=df.drip_int_stal3
# drip4=df.drip_int_stal4
# drip5=df.drip_int_stal5

# fig = plt.figure(figsize=(15,12))
# fig.suptitle('Plots for %s'%(output), fontsize=14, fontweight='bold')
# ax1 = fig.add_subplot(321)
# ax1.set_title("Stal 1")    
# ax1.set_xlabel('Time(months)')
# ax1.set_ylabel('d18O(per mille)')
# ax1.set_ylim([-7,-0])
# ax1.plot(time,stal1, c='r', label='d18O')

# ax2 = fig.add_subplot(322)
# ax2.set_title("Stal 2")    
# ax2.set_xlabel('Time(months)')
# ax2.set_ylabel('d18O(per mille)')
# ax2.set_ylim([-7,-0])
# ax2.plot(time,stal2, c='r', label='d18O')

# ax3 = fig.add_subplot(323)
# ax3.set_title("Stal 3")    
# ax3.set_xlabel('Time(months)')
# ax3.set_ylabel('d18O(per mille)')
# ax3.set_ylim([-7,-0])
# ax3.plot(time,stal3, c='r', label='d18O')

# ax4 = fig.add_subplot(324)
# ax4.set_title("Stal 4")    
# ax4.set_xlabel('Time(months)')
# ax4.set_ylabel('d18O(per mille)')
# ax4.set_ylim([-7,-0])
# ax4.plot(time,stal4, c='r', label='d18O')

# ax5 = fig.add_subplot(325)
# ax5.set_title("Stal 5")    
# ax5.set_xlabel('Time(months)')
# ax5.set_ylabel('d18O(per mille)')
# ax5.set_ylim([-7,-0])
# ax5.plot(time,stal5, c='r', label='d18O')

# ax6 = fig.add_subplot(326)
# ax6.set_title("Drip Interval")    
# ax6.set_xlabel('Time(months)')
# ax6.set_ylabel('Drip Interval (s)')
# ax6.set_ylim([0,1500])
# ax6.plot(time,drip1, c='r', label='Stal 1')
# ax6.plot(time,drip2, c='g', label='Stal 2')
# ax6.plot(time,drip3, c='b', label='Stal 3')
# ax6.plot(time,drip4, c='c', label='Stal 4')
# ax6.plot(time,drip5, c='y', label='Stal 5')
# leg = ax6.legend()

# fig.tight_layout()
# plt.show()
# plt.savefig('%s.png'%output)