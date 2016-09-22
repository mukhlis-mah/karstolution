import wx, csv, os
from functools import partial #used to make some calls to GUI events
from Main import opj #used to import images
import scipy.stats as s
import numpy as np
import matplotlib.pyplot as plt
#importing main karstolution module, ie the actual model
from karstolution1_1 import karstolution 
from advanced import Advanced #importing the advanced layout
from basic import Basic #importing the basic layout
#pandas is used to simplify the plotting of results
#this is optional, so if you can remove the last 'plot_all' module
#if you don't have pandas already installed
import pandas as pd

#**********************************************************************
#GUI created using wxpython
#This module contains the main shell for GUI
#the main features are located in this module, which uses 
#details for the basic view are then in the basic module
#details for the advanced view module are in that module
#**********************************************************************

class simpleapp_wx(wx.Frame):
    
    def __init__(self,parent,id,title):
        wx.Frame.__init__(self,parent,id,title,size=(1000,1000))
        self.parent=parent
        
        self.name_config='config.csv'
        self.csv_reader(self.name_config)
        #this is the configuration file which holds all the parameters
        #which are changable from the GUI. Any changes to the parameters
        #in the GUI will save it automatically to this csv.
        #This file is read by the csv_reader function,
        #which extracts all parameters into the list, self.data
        #The order of the parameters in the csv are:
        #fluxes=data[0] ie. the first row...
        #bypass_comp=data[1]
        #Evaporation_coeffs=data[2]
        #Cave_temps=data[3]
        #drip_interval=data[4]
        #store_sizes=data[5]
        #weibull_param=data[6]
        #drip_pCO2=data[7]
        #cave_pCO2=data[8]
        #humidity=data[9]
        #Wind=data[10]
        #mixing_parameter=data[11]
        #initial_store=data[12]
        #initial_d18o=data[13]
        #the full details are avaialble at the GitHub depository
        # https://github.com/mukhlis-mah/karstolution
        
        #**********************************************************************
        
        #input csv file, this contains the input climatic time seriess
        #and is necessary to run the model but is only accessed from the 
        #karstolution1.1 module (ie not from the GUI)
        self.input_file='input.csv'
        #this is the output file name, also written from karstolution1.1
        self.output_file='output.csv'
        
        #whether to use absolute user-inputted
        self.calculate_drip=True 
        
        #id used for batch_run, determines which parameter has been chosen from
        #the drop-down list. Default set to zero (ie. no parameter)
        self.batch_id=0 
        #list of the iterable parameters for the batch run and their id so they can be
        #called from self.data (eg. f3 is self.data[0][1])
        self.list_options=[['None'],['f1',0,0],['f3',0,1],['f8',0,6],['f5',0,2],['f6',0,3],['f7',0,4],['k_diffuse',0,5],
        ["k_eevap",2,0],['k_d18o_soil',2,1],["k_d18o_epi",2,2],["Soil_Capac",5,0],["Epikarst_Capac",5,3],
        ["KS1_Capac",5,4],["KS2_Capac",5,5],["Epicap",5,1],["Ovicap",5,2],["Init_SoilStore",12,0],["Init_Epikarst",12,1],
        ["Init_KS1",12,2],["Init KS2",12,3],["Init_Diffuse",12,4],["Init_Soil_d8o",13,0],
        ["Init Epikarst_d8o",13,1],["Init KS1 d8o",13,2],["Init KS2 d8o",13,3],["Init_Diffuse_d8o",13,5],
        ["Weibull_Lambda",6,0],["Weibull_k",6,1],['Bypass_i',1,0],['Bypass_m',1,3],["MixingParameter",11,0],["CaveTemp",3,0],
        ["Mindripint",4,0],["Drip_pCO2",7,0],["Cave_pCO2",8,0],["Rel_Humdity",9,0],["Ventilation",10,0]]
        #parameters for iteration [min value, max value, # iterations]
        #eg. the below is iterations starting at 0 till 0.2 in 10 steps
        #ie. 0, 0.02222, 0.04444, 0.06666, 0.08888,...,0.2
        self.batch_p=[0,0.2,10]
        #runs function below, with the wx pythony things
        self.initialise()
        
    def initialise(self):
        
        #title font
        self.font = wx.Font(14, wx.DEFAULT, wx.NORMAL, wx.BOLD)
               
        #use of scroll bar window
        self.scroll = wx.ScrolledWindow(self, -1)
        self.scroll.SetScrollbars(1,1,1000,1000)
        
        #setting up the sizer for use with scroll bar window 
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.scroll.SetSizer(self.sizer)
        
        #the karstolution schematic diagram in the top-left
        png = wx.Image(opj('structure.png'), wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        wx.StaticBitmap(self.scroll, -1, png, (10, 10), (png.GetWidth(), png.GetHeight())) 
        #the run button on the top-left, binded to self.OnButtonClick function
        button=wx.Button(self.scroll,-1,"Run Karstolution",(800,60))
        self.Bind(wx.EVT_BUTTON,partial( self.OnButtonClick,data=self.data,output=self.output_file),button)
        # the rest of the GUI interface takes either from the basic or advanced view
        #which come from their independant modules
        self.basic_panel=Basic(self.scroll,self.data,self.name_config)
        self.advanced_panel=Advanced(self.scroll,self.data,self.name_config)
        #the basic view is hidden first but can be switched into
        self.basic_panel.Hide()
        #adding the two views into the sizer
        self.sizer.Add(self.basic_panel, 1, wx.EXPAND)
        self.sizer.Add(self.advanced_panel, 2, wx.EXPAND)
        
        #creating the status bar at the bottom
        self.CreateStatusBar()
        
        #creating menu at top
        filemenu=wx.Menu()
        menuBar=wx.MenuBar()
        #file is only option at this stage
        menuBar.Append(filemenu,"&File")
        self.SetMenuBar(menuBar)
        
        #switch between the Advanced and Basic view
        toggle=filemenu.Append(wx.ID_ANY,"Switch Basic/Advanced",
        """Switch between the default Basic view or the Advanced view 
        with all possible parameters""")
        #open box with information about Karstolution
        menuAbout=filemenu.Append(wx.ID_ABOUT, "&About",
        "Information about this program")
        #exit the program
        menuExit=filemenu.Append(wx.ID_EXIT,"&Exit","Terminate the program")
        
        #binding the events (functions) to the menu items
        self.Bind(wx.EVT_MENU,self.OnAbout,menuAbout)
        self.Bind(wx.EVT_MENU,self.OnCheck,toggle)
        self.Bind(wx.EVT_MENU,self.OnExit, menuExit)
        
        #hidden features, that you need to expand the window to access
        #the Calculate Drip option (default=yes), allows you to over-ride the
        #drip-rate calculation (based on user-inputted min drip interval & store size)
        #and instead run the model just using the user-inputted values as the
        #absolute drip-rate (used for all time-steps of the model)
        self.cb = wx.CheckBox(self, -1, 'Calculate Drips ?', (1050, 270))
        self.cb.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.GetId, self.cb)
        
        #Allows the model to be run multiple times, iterating over a range
        #of parameter values. 
        batch1=wx.Button(self,-1,"Run Batch Mode",(1050,225))
        self.Bind(wx.EVT_BUTTON,self.OnBatch,batch1)
        self.label=wx.StaticText(self,-1,'Batch Mode: ',(1050,10))
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'Choose parameter to iterate: ',(1050,40))
        lister=wx.ListBox(self,-1,(1050,60),(180,80),[thing[0] for thing in self.list_options],wx.LB_SINGLE)
        lister.SetSelection(0)
        self.Bind(wx.EVT_LISTBOX,self.OnChoose,lister)
        self.min_b=wx.StaticText(self,-1,'Min Value: ',(1050,150))
        self.min_batch=wx.TextCtrl(self,-1,str(self.batch_p[0]),(1130,150))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.min_batch,id1=0),self.min_batch)
        self.max_b=wx.StaticText(self,-1,'Max Value: ',(1050,175))
        self.max_batch=wx.TextCtrl(self,-1,str(self.batch_p[1]),(1130,175))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.max_batch,id1=1),self.max_batch)
        self.it_b=wx.StaticText(self,-1,'# Iterations: ',(1050,200))
        self.it_batch=wx.TextCtrl(self,-1,str(self.batch_p[2]),(1130,200))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.it_batch,id1=2),self.it_batch)
        
        #telling wxpython the layout is all G       
        self.sizer.Layout()
    
    #function used to extract existing values (used in binding events)
    def assign(self,event,name,id1):
        self.batch_p[id1]=name.GetValue()
        
    #reads the config csv file and extracts it to the self.data list       
    def csv_reader(self,config):
        csvinput=open(config, 'rU')
        reader=csv.reader(csvinput)
        self.data=[]
        for row in reader:
            if row:
                self.data.append([])
                for thing in row:
                    self.data[-1].append(float(thing))
        return self.data
        
    #toggling between basic and advanced view
    def OnCheck(self,e):
        if self.basic_panel.IsShown():
            self.scroll.Scroll(1,1)
            self.basic_panel.Hide()
            self.advanced_panel.Show()
        elif self.advanced_panel.IsShown():
            self.scroll.Scroll(1,1)
            self.advanced_panel.Hide()
            self.basic_panel.Show()
        self.sizer.Layout()
    
    #info about Karstolution in the message dialogue
    def OnAbout(self,e):
        dlg=wx.MessageDialog(self,"""The Karstolution integrated karst and in-cave fractionation\
        speleothem d18O PSM was developed by coupling the KarstFor model with the ISOLUTION model.
        For documentation, python scripts, YouTube tutorial, papers and more, visit:\
        www.connectedwaters.unsw.edu.au/karstolution 
        Code written by Mukhlis Mah, UNSW
        Contact: a.baker@unsw.edu.au""",
        "Karstolution Speleothem d18O PSM",wx.OK)
        dlg.ShowModal()
        dlg.Destroy
    
    #exits the program
    def OnExit(self,e):
        self.Close(True)
    
    #allows the choice of a selection from the list (for batch run)
    def OnChoose(self,e):
        self.batch_id=e.GetSelection()
    
    #does all the magic for the batch run
    #nb: you could play around here to iterate across two parameters....
    def OnBatch(self,e):
        #makes sure they choose a parameter
        if self.batch_id==0:
            print "Choose a parameter to run batch mode!"
        choice=self.list_options[self.batch_id] #their choice
        min=float(self.batch_p[0]) #min value
        max=float(self.batch_p[1]) #max value
        iterations=float(self.batch_p[2]) #iterations
        step=(max-min)/(iterations-1) #works out the step based on above
        original=self.output_file #the original output file to revert to after batch
        #gets the current paramater value to revert to after batch
        original_d=self.data[choice[1]][choice[2]]
        #creates a folder to put output files into
        if not os.path.exists('%s' %(choice[0])):
            os.makedirs('%s' %(choice[0]))
        #the loop for batch run
        for number in np.arange(min,max+step,step):
            #names the output file for each iteration as 
            # output_parameter_parametervalue (eg. output_f1_0.05)
            self.output_file="%s/output_%s_%f.csv" %(choice[0],choice[0],number)
            #for the cave parameters, the batch mode adds the value to the existing value
            #this is to preserve seasonality
            if choice[0] in ["CaveTemp","Mindripint","Drip_pCO2","Cave_pCO2",
            "Rel_Humdity","Ventilation"]:
                for i in range(0,12):
                    self.data[choice[1]][i]+=number
            #for bypass i and m, the other bypass parameters are adjusted to everything adds to 1 
            elif choice[0]=='Bypass_i':
                self.data[choice[1]][choice[2]]=number
                self.data[choice[1]][choice[2]+1]=(1-number)/2
                self.data[choice[1]][choice[2]+2]=(1-number)/2
            elif choice[0]=='Bypass_m':
                self.data[choice[1]][choice[2]]=number
                self.data[choice[1]][choice[2]+1]=1-number
            #all other parameters are the absolute value calculated in the batch run
            #ie. independant of what value they have in the configuration file (displayed on the GUI)
            else:
                self.data[choice[1]][choice[2]]=number
            #clicking the to activate batch run!
            self.OnButtonClick(e,self.data,self.output_file)
            #reverting the cave parameter
            if choice[0] in ["CaveTemp","Mindripint","Drip_pCO2","Cave_pCO2","Rel_Humdity",
            "Ventilation"]:
                    for i in range(0,12):
                            self.data[choice[1]][i]+= -number
        #reverting the original output file and other original parameters
        self.output_file=original
        self.data[choice[1]][choice[2]]=original_d
            
    #Do you want to Karstolution to model drip interval from the level of karst stores? 
    #default is yes. If no, it will just use repeat the monthly drip int averages from the config file
    def GetId(self,e):
        self.calculate_drip=self.cb.GetValue()
    
    #runs the actual model!
    def OnButtonClick(self,event,data,output):
        #text in the status bar
        self.SetStatusText('Running Karsolution! Please wait')
        running=karstolution(data,self.input_file,output,self.calculate_drip)
        #text in the status bar
        if running:
            self.SetStatusText('Finished running Karstolution!')
        #plots the output from the csv file(s) created
        self.plot_all()
        self.Refresh()
    #plots the output from the csv file(s) created
    #nb: this was just done as default graphing to get a crude feel for results
    #this part can be edited and added to, such that you can display the appropriate data
    #you interested in!
    def plot_all(self):
        #using pandas here, as it easily pulls out the columns of the output csv
        df = pd.read_csv(self.output_file)
        time = df.tt
        stal1=df.stal1d18o
        stal2=df.stal2d18o
        stal3=df.stal3d18o    
        stal4=df.stal4d18o
        stal5=df.stal5d18o
        drip1=df.drip_int_stal1
        drip2=df.drip_int_stal2
        drip3=df.drip_int_stal3
        drip4=df.drip_int_stal4
        drip5=df.drip_int_stal5
        
        #plotting six figures based on the five output stals + drip intervals
        fig = plt.figure(figsize=(15,12))
        fig.suptitle('Plots for %s'%(self.output_file), fontsize=14, fontweight='bold')
        ax1 = fig.add_subplot(321)
        ax1.set_title("Stal 1")    
        ax1.set_xlabel('Time(months)')
        ax1.set_ylabel('d18O(per mille)')
        ax1.set_ylim([-7,-0])
        ax1.plot(time,stal1, c='r', label='d18O')

        ax2 = fig.add_subplot(322)
        ax2.set_title("Stal 2")    
        ax2.set_xlabel('Time(months)')
        ax2.set_ylabel('d18O(per mille)')
        ax2.set_ylim([-7,-0])
        ax2.plot(time,stal2, c='r', label='d18O')

        ax3 = fig.add_subplot(323)
        ax3.set_title("Stal 3")    
        ax3.set_xlabel('Time(months)')
        ax3.set_ylabel('d18O(per mille)')
        ax3.set_ylim([-7,-0])
        ax3.plot(time,stal3, c='r', label='d18O')

        ax4 = fig.add_subplot(324)
        ax4.set_title("Stal 4")    
        ax4.set_xlabel('Time(months)')
        ax4.set_ylabel('d18O(per mille)')
        ax4.set_ylim([-7,-0])
        ax4.plot(time,stal4, c='r', label='d18O')

        ax5 = fig.add_subplot(325)
        ax5.set_title("Stal 5")    
        ax5.set_xlabel('Time(months)')
        ax5.set_ylabel('d18O(per mille)')
        ax5.set_ylim([-7,-0])
        ax5.plot(time,stal5, c='r', label='d18O')

        ax6 = fig.add_subplot(326)
        ax6.set_title("Drip Interval")    
        ax6.set_xlabel('Time(months)')
        ax6.set_ylabel('Drip Interval (s)')
        for thing in [drip1,drip2,drip3,drip4,drip5]:
            if max(thing)>1000:
                    ax6.set_ylim([0,1000])
        ax6.plot(time,drip1, c='r', label='Stal 1')
        ax6.plot(time,drip2, c='g', label='Stal 2')
        ax6.plot(time,drip3, c='b', label='Stal 3')
        ax6.plot(time,drip4, c='c', label='Stal 4')
        ax6.plot(time,drip5, c='y', label='Stal 5')
        leg = ax6.legend()

        fig.tight_layout()
        plt.show()
        plt.savefig('%s.png'%self.output_file) #saves the plot as well!

app=wx.App(False)
frame=simpleapp_wx(None,-1,'Karstolution')
frame.Show()
app.MainLoop()    