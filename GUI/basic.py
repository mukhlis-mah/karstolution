import wx
import csv
from functools import partial
import scipy.stats as s
import numpy as np
import matplotlib.pyplot as plt

class Basic(wx.Panel):
    def __init__(self,parent,data,name_config):
        self.data1=data
        self.config=name_config
        wx.Panel.__init__(self,parent=parent)
        self.font = wx.Font(14, wx.DEFAULT, wx.NORMAL, wx.BOLD)
        
        #store sizes: all
        #overflow lims: all
        #store init values: except diffuse
        
        text_1="""This is the basic view of Karstolution, with selected parameters and overall mean cave inputs. 
        Change to advanced view for all parameters and cave monthly means"""
        text_2="""Input: "input.csv" in the same directory as Karstolution"""
        text_3="""Output: will create "output.csv" in same directory + plot some values"""
        text_4="""Config: the paramaters in the gui are stored in "config.txt" in the same directory and 
        automatically updates with changes."""
        text_5="""Batch Mode: (full screen) choose a parameter to iterate, # iterations, max and min. 
        Files will be ouputted into a new folder in the same directory"""
        text_6="""Calculate drips: (full screen) ticked=calculate (default); unticked=user inputted drip interval 
        After executing Karstolution, the GUI will freeze, so please be patient."""
        
        posx=350
        #fluxes
        self.flabel=wx.StaticText(self,-1,'Fluxes',(30,posx),(250,-1),wx.ALIGN_CENTER)
        self.flabel.SetFont(self.font)
        self.flabel=wx.StaticText(self,-1,'(1/month)',(110,posx+10),(200,-1))
        self.label=wx.StaticText(self,-1,'Cave Parameters',(220,posx),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'Karstolution Beta 1.1(Basic View)',(500,posx),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        
        posx+=40
        self.F1=wx.StaticText(self,-1,'F1: ',(10,posx))
        self.F1_v=wx.TextCtrl(self,-1,str(self.data1[0][0]),(70,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.F1_v,id1=0,id2=0),self.F1_v)
        self.cave_t=wx.StaticText(self,-1,'Cave Temp (Celcius): ',(190,posx))
        self.cave_temp=wx.TextCtrl(self,-1,str(self.data1[3][0]),(340,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign_cave, name=self.cave_temp,id1=3),self.cave_temp)
        
        
        self.Drip=wx.StaticText(self,-1,text_1,(450,posx))
        
        posx+=25
        self.F3=wx.StaticText(self,-1,'F3: ',(10,posx))
        self.F3_v=wx.TextCtrl(self,-1,str(self.data1[0][1]),(70,posx))
        self.Drip=wx.StaticText(self,-1,'Min Drip interval (s): ',(190,posx))
        self.Drip_interval=wx.TextCtrl(self,-1,str(self.data1[4][0]),(340,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign_cave, name=self.Drip_interval,id1=4),self.Drip_interval)
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.F3_v,id1=0,id2=1),self.F3_v)
        
        
        posx+=25
        self.F8=wx.StaticText(self,-1,'F8: ',(10,posx))
        self.F8_v=wx.TextCtrl(self,-1,str(self.data1[0][6]),(70,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.F8_v,id1=0,id2=6),self.F8_v)
        self.Drip=wx.StaticText(self,-1,'Drip pCO2 (ppmv): ',(190,posx))
        self.Drip_pCO2=wx.TextCtrl(self,-1,str(self.data1[7][0]),(340,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign_cave, name=self.Drip_pCO2,id1=7),self.Drip_pCO2)
        self.Drip=wx.StaticText(self,-1,text_2,(450,posx))
        
        
        posx+=25
        self.F7=wx.StaticText(self,-1,'F7: ',(10,posx))
        self.F7_v=wx.TextCtrl(self,-1,str(self.data1[0][4]),(70,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.F7_v,id1=0,id2=4),self.F7_v)
        self.Cave=wx.StaticText(self,-1,'Cave pCO2 (ppmv): ',(190,posx))
        self.Cave_pCO2=wx.TextCtrl(self,-1,str(self.data1[8][0]),(340,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign_cave, name=self.Cave_pCO2,id1=8),self.Cave_pCO2)
        self.Drip=wx.StaticText(self,-1,text_3,(450,posx))
        
        posx+=25
        self.Wind=wx.StaticText(self,-1,'Wind (m/s): ',(190,posx))
        self.v=wx.TextCtrl(self,-1,str(self.data1[10][0]),(340,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign_cave, name=self.v,id1=10),self.v)
        self.Drip=wx.StaticText(self,-1,text_4,(450,posx))
        posx+=25
        self.Drip=wx.StaticText(self,-1,text_5,(450,posx+10))
        self.Hum=wx.StaticText(self,-1,'Humidity (%; 0<h<=1): ',(190,posx))
        self.Humidity=wx.TextCtrl(self,-1,str(self.data1[9][0]),(340,posx))
        self.Bind(wx.EVT_TEXT,partial( self.assign_cave, name=self.Humidity,id1=9),self.Humidity)
        posx+=45
        self.Drip=wx.StaticText(self,-1,text_6,(450,posx))
        
        #store size
        self.flabel=wx.StaticText(self,-1,'Store Sizes',(490,10),(200,-1),wx.ALIGN_CENTER)
        self.flabel.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(mm)',(630,20))
        self.soils=wx.StaticText(self,-1,'Soil Store: ',(475,40))
        self.soilcap=wx.TextCtrl(self,-1,str(self.data1[5][0]),(575,40))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.soilcap,id1=5,id2=0),self.soilcap)
        self.epics=wx.StaticText(self,-1,'Epikarst: ',(475,65))
        self.epistore=wx.TextCtrl(self,-1,str(self.data1[5][3]),(575,65))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.epistore,id1=5,id2=3),self.epistore)
        self.ks1=wx.StaticText(self,-1,'Karst Store 1: ',(475,90))
        self.ks1_size=wx.TextCtrl(self,-1,str(self.data1[5][4]),(575,90))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.ks1_size,id1=5,id2=4),self.ks1_size)
        self.ks2=wx.StaticText(self,-1,'Karst Store 2: ',(475,115))
        self.ks2_size=wx.TextCtrl(self,-1,str(self.data1[5][5]),(575,115))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.ks2_size,id1=5,id2=5),self.ks2_size)
        
        #Overflow limits
        self.flabel=wx.StaticText(self,-1,'Overflow Limits',(475,140),(200,-1),wx.ALIGN_CENTER)
        self.flabel.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(mm)',(665,150))
        self.epic=wx.StaticText(self,-1,'Epicap: ',(475,175))
        self.epicap=wx.TextCtrl(self,-1,str(self.data1[5][1]),(530,175))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.epicap,id1=5,id2=1),self.epicap)
        self.ovc=wx.StaticText(self,-1,'Ovicap: ',(475,200))
        self.ovicap=wx.TextCtrl(self,-1,str(self.data1[5][2]),(530,200))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.ovicap,id1=5,id2=2),self.ovicap)
        
        pos_init=100
        #Store Initial Values
        self.label=wx.StaticText(self,-1,'Store Initial Values',(700,pos_init),(200,-1),wx.ALIGN_CENTER)
        self.label.SetFont(self.font)
        self.label=wx.StaticText(self,-1,'(mm)',(800,pos_init+30))
        pos_init+=50
        self.ss=wx.StaticText(self,-1,'Soil Store: ',(710,pos_init))
        self.sstore=wx.TextCtrl(self,-1,str(self.data1[12][0]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.sstore,id1=12,id2=0),self.sstore)
        pos_init+=25
        self.es=wx.StaticText(self,-1,'Epikarst: ',(710,pos_init))
        self.estore=wx.TextCtrl(self,-1,str(self.data1[12][1]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.estore,id1=12,id2=1),self.estore)
        pos_init+=25
        self.ks1=wx.StaticText(self,-1,'Karst Store 1: ',(710,pos_init))
        self.kstore1=wx.TextCtrl(self,-1,str(self.data1[12][2]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.kstore1,id1=12,id2=2),self.kstore1)
        pos_init+=25
        self.ks2=wx.StaticText(self,-1,'Karst Store 2: ',(710,pos_init))
        self.kstore2=wx.TextCtrl(self,-1,str(self.data1[12][3]),(810,pos_init))
        self.Bind(wx.EVT_TEXT,partial( self.assign, name=self.kstore2,id1=12,id2=3),self.kstore2)
        
        
        
    def assign(self,event,name,id1,id2):
        self.data1[id1][id2]=float(name.GetValue())
        csvoutput=open(self.config,'w')
        writer=csv.writer(csvoutput)
        for list in self.data1:
            writer.writerow(list)
        
    def assign_cave(self,event,name,id1):
        for counter in range(0,12):
            self.data1[id1][counter]=float(name.GetValue())
        csvoutput=open(self.config,'w')
        writer=csv.writer(csvoutput)
        for list in self.data1:
            writer.writerow(list)
    