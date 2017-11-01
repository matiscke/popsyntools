# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:55:15 2016

@author: bec
"""
import scipy.optimize as opt
import matplotlib 
import PIL as pil
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
import MikeMath as MM
import csv
import glob
import os




def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
    
    
    
class analyseData:
    
    def __init__(self):
        self.filename=None
        self.path=None
        self.rawData=None
        self.data=None
        self.energy=None
        self.wavelength=None
        self.time=None
        self.yValues=None
        self.yValuesFloat=None
        self.yName=None
        self.yUnit=None
        self.FitParameters=None
        self.FitParametersError=None
        self.FitData=None
        self.FitDataList=None
        self.header=[]
        self.image=None
        self.dataList=None
        self.imageList=None
        self.xDistance=None
        self.xImageScaling=None
        self.yImageScaling=None
        self.dataFFT=None
        self.xValues=None
        self.xValuesFloat=None
        self.xName=None
        self.xUnit=None
        self.zValues=None
        self.zValuesFloat=None
        self.zName=None
        self.zUnit=None
        self.format=None
        return

    
    
    def importExportRandomCSV(self,path,plot=True,save=True):
                
        self.path=path
        self.format=".csv"
        for file in os.listdir(self.path):
            if file.endswith(".csv"):
                self.filename=None
                self.path=None
                self.rawData=None
                self.data=None
                self.energy=None
                self.wavelength=None
                self.time=None
                self.yValues=None
                self.yValuesFloat=None
                self.yName=None
                self.yUnit=None
                self.FitParameters=None
                self.FitParametersError=None
                self.FitData=None
                self.header=[]
                print file
                self.importDataCSV(path,file)
                if plot:                
                    self.plotData()
                if save:
                    self.saveDataTXT()
    def importDataTIFF(self,path,filename=None, sortByDate=False, Offset=0, ScalingX=1,ScalingY=1,RescaleFactor=None):
        self.path=path
        self.xImageScaling=ScalingX   
        self.yImageScaling=ScalingY   
        self.dataList=[]
        self.imageList=[]
        self.format=".tif"
        if filename==None: # Load several Single Files                
            print "No filename provided"
            filenames=[]            
            
            if sortByDate:
                filenames=sorted(glob.glob(os.path.join(path, '*.tif')),key=os.path.getmtime)
            else:
                for filename in glob.glob(os.path.join(path, '*.tif')):   
                    filenames.append(filename)
            for i in range(len(filenames)):
                print filenames[i]
                self.filename=filenames[-1].split("\\")[-1]
                self.image=pil.Image.open(filenames[i])
                self.data=pd.DataFrame(np.array(self.image)-float(Offset))
                if RescaleFactor==None:
                    self.image=pil.Image.fromarray(self.data.as_matrix())
                else:
                    self.image=pil.Image.fromarray(self.data.as_matrix()/RescaleFactor)   
                self.dataList.append(self.data)
                self.imageList.append(self.image)
            self.xDistance=np.linspace(0,self.image.size[1]-1,self.image.size[1])*self.xImageScaling
            self.yValues=np.linspace(0,self.image.size[0]-1,self.image.size[0])*self.yImageScaling
            self.set_yValues(yName="y",yUnit="um")
            self.xValues=self.xDistance
            self.set_xValues(self.xDistance,"x","um")
            self.zValues=np.linspace(0,len(filenames)-1,len(filenames))
            self.set_zValues(np.linspace(0,len(filenames)-1,len(filenames)),"Time","arb. units")
        else:
            self.filename=filename
            self.image=pil.Image.open(path+filename)
            self.data=pd.DataFrame(np.array(self.image)-float(Offset))
            self.xDistance=np.linspace(0,self.image.size[1]-1,self.image.size[1])*self.xImageScaling
            self.yValues=np.linspace(0,self.image.size[0]-1,self.image.size[0])*self.yImageScaling
            self.set_yValues(yName="y",yUnit="um")
            self.xValues=self.xDistance
            self.set_xValues(self.xDistance,"x","um")
            self.zValues=np.linspace(0,len(filenames)-1,len(filenames))
            self.set_zValues(np.linspace(0,len(filenames)-1,len(filenames)),"Time","arb. units")
            if RescaleFactor==None:
                self.image=pil.Image.fromarray(self.data.as_matrix())
            else:
                self.image=pil.Image.fromarray(self.data.as_matrix()/RescaleFactor) 
        
            #self.image.show()
    def importDataCSV(self,path,filename=None, sortByDate=False):
        self.path=path
            
        self.format=".csv"
        if filename==None: # Load several Single Files
            print "No filename provided"
            filenames=[]            
            
            if sortByDate:
                filenames=sorted(glob.glob(os.path.join(path, '*.csv')),key=os.path.getmtime)
            else:
                for filename in glob.glob(os.path.join(path, '*.csv')):   
                    filenames.append(filename) 
            
            #Determine Header length
                        
            with open(filenames[0]) as csvfile:
                reader = csv.reader(csvfile)
            
                for row in reader:
                    if row[0]=="Wavelength" or row[0]=="Time":
                        break
                        
                    else:
                        self.header.append(str(row))
    
            self.rawData=pd.read_csv(filenames[0],header=len(self.header))

            for filename in filenames[1:]:
                
                currentfile=pd.read_csv(filename,header=len(self.header))
                self.rawData=pd.concat([self.rawData,currentfile.iloc[:,1]],axis=1)           
                
            self.yValues=self.rawData.columns.values[1:]
            self.set_yValues()  
            self.data=self.rawData.iloc[:,1:]
            if self.rawData.columns.values[0]=='Wavelength':
                self.wavelength=self.rawData.Wavelength
                self.energy=1239.8/self.rawData.Wavelength
                self.xValues=self.energy
                self.set_xValues(self.energy,"Energy","eV")
            elif self.rawData.columns.values[0]=='Time':
                self.time=self.rawData.Time
                self.xValues=self.time                
                self.set_xValues(self.time,"Time","ns")
            else:
                print "Unknown Data Format"    
            self.filename=filenames[-1].split("\\")[-1]
            print str(len(filenames))+" Files loaded"
        
        else: # Load Single File
            self.filename=filename
            
            
            #Determine Header length
            
            with open(path+filename) as csvfile:
                reader = csv.reader(csvfile)
            
                for row in reader:
                    if row[0]=="Wavelength" or row[0]=="Time":
                        break
                        
                    else:
                        self.header.append(str(row))
                        
            
            # Load Data using Pandas
            self.rawData=pd.read_csv(path+filename,header=len(self.header))
            self.data=self.rawData.iloc[:,1:]
            self.yValues=self.rawData.columns.values[1:]
            self.set_yValues()
            
            
            if self.rawData.columns.values[0]=='Wavelength':
                self.wavelength=self.rawData.Wavelength
                self.energy=1239.8/self.rawData.Wavelength
                self.xValues=self.energy
                self.set_xValues(self.energy,"Energy","eV")
            elif self.rawData.columns.values[0]=='Time':
                self.time=self.rawData.Time
                self.xValues=self.time                
                self.set_xValues(self.time,"Time","ns")
            else:
                print "Unknown Data Format" 
                
    def load_yValues(self,pathFilename=None,yName="y Values", yUnit = "arb. Units"):
        if not pathFilename:
            print "Give path + Filename to y-Values Summary"
            return
        yValues=[]
        
        with open(pathFilename) as csvfile:
            reader=csv.reader(csvfile)  
            for row in reader:
                yValues.append(float(row[0]))
        if len(yValues)==len(self.yValues):
                self.yValuesFloat=yValues
                self.yValues=[str(yValues[i]) for i in range(len(self.yValues))]
        else:
            print "Not right amount of y-Values in File"
        print yValues
        self.yName=yName
        self.yUnit=yUnit
        
                
    def set_yValues(self,yValues=[],yName="y Values", yUnit = "arb. Units", filenamepath=None):
        if not yValues:
            if isfloat(self.yValues[0]):
                self.yValuesFloat=self.yValues.astype(np.float)
            else:
                self.yValuesFloat=range(len(self.yValues))
        else:
            if len(yValues)==len(self.yValues):
                self.yValuesFloat=yValues
                self.yValues=[str(yValues[i]) for i in range(len(self.yValues))]
            else:
                print "Wrong amount of entered yValues. "+str(len(self.yValues))+" Values needed!"
        self.yName=yName
        self.yUnit=yUnit
    
    def set_zValues(self,zValues=[],zName="z Values", zUnit = "arb. Units", filenamepath=None):
        
        if len(zValues)==len(self.zValues):
            self.zValuesFloat=zValues
            self.zValues=[str(zValues[i]) for i in range(len(self.zValues))]
        else:
            print "Wrong amount of entered zValues. "+str(len(self.zValues))+" Values needed!"
        self.zName=zName
        self.zUnit=zUnit
    
    def set_xValues(self,xValues=[],xName="x Values", xUnit = "arb. Units", filenamepath=None):
        
        if len(xValues)==len(self.xValues):
            self.xValuesFloat=xValues
            self.xValues=[str(xValues[i]) for i in range(len(self.xValues))]
        else:
            print "Wrong amount of entered xValues. "+str(len(self.xValues))+" Values needed!"
        self.xName=xName
        self.xUnit=xUnit
    
    def shiftToZero(self,Normalize=True):
        if Normalize:
           self.normalizeData()
        for i in range(len(self.data.columns)):
            if i==0:
                maxIndex1=np.argmax(data.data.iloc[:,i])
                xValueMaximum=self.xValuesFloat[np.argmax(data.data.iloc[:,i])]
                
                if self.xValuesFloat is not None:
                    xValues=self.xValuesFloat-xValueMaximum
                    self.xValuesFloat=xValues
                    self.xValues=[str(xValues[i]) for i in range(len(self.xValues))]
                    
            else:
                shift=maxIndex1-np.argmax(data.data.iloc[:,i])
                print shift
                self.data.iloc[:,i]=np.roll(np.array(self.data.iloc[:,i]),shift)
                    
    def setBoundariesData(self,xRange=[-np.inf,np.inf],yRange=[-np.inf,np.inf]):
        if not self.xValuesFloat is None:        
            xminIndex=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-xRange[0]))[0]
            xmaxIndex=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-xRange[1]))[0]
            if xminIndex>xmaxIndex:
                a=xminIndex
                xminIndex=xmaxIndex
                xmaxIndex=a
            self.xValuesFloat=self.xValuesFloat[xminIndex:xmaxIndex]
            self.xValues=self.xValues[xminIndex:xmaxIndex]
            yminIndex=min(enumerate(self.yValuesFloat), key=lambda x: abs(x[1]-yRange[0]))[0]            
            ymaxIndex=min(enumerate(self.yValuesFloat), key=lambda x: abs(x[1]-yRange[1]))[0]
            self.yValues=self.yValues[yminIndex:ymaxIndex]
            self.yValuesFloat=self.yValuesFloat[yminIndex:ymaxIndex]
        
        elif not self.energy is None:        
            xminIndex=min(enumerate(self.energy), key=lambda x: abs(x[1]-xRange[0]))[0]
            xmaxIndex=min(enumerate(self.energy), key=lambda x: abs(x[1]-xRange[1]))[0]
            if xminIndex>xmaxIndex:
                a=xminIndex
                xminIndex=xmaxIndex
                xmaxIndex=a
            self.energy=self.energy[xminIndex:xmaxIndex]
            yminIndex=min(enumerate(self.yValuesFloat), key=lambda x: abs(x[1]-yRange[0]))[0]            
            ymaxIndex=min(enumerate(self.yValuesFloat), key=lambda x: abs(x[1]-yRange[1]))[0]
            self.yValues=self.yValues[yminIndex:ymaxIndex]
            self.yValuesFloat=self.yValuesFloat[yminIndex:ymaxIndex]
        elif not self.time is None:
            xminIndex=min(enumerate(self.time), key=lambda x: abs(x[1]-xRange[0]))[0]
            xmaxIndex=min(enumerate(self.time), key=lambda x: abs(x[1]-xRange[1]))[0]
            if xminIndex>xmaxIndex:
                a=xminIndex
                xminIndex=xmaxIndex
                xmaxIndex=a
            self.time=self.time[xminIndex:xmaxIndex]
            yminIndex=min(enumerate(self.yValuesFloat), key=lambda x: abs(x[1]-yRange[0]))[0]            
            ymaxIndex=min(enumerate(self.yValuesFloat), key=lambda x: abs(x[1]-yRange[1]))[0]
            self.yValues=self.yValues[yminIndex:ymaxIndex]
            self.yValuesFloat=self.yValuesFloat[yminIndex:ymaxIndex]
        
        
        self.data=self.data.iloc[xminIndex:xmaxIndex,yminIndex:ymaxIndex]
        
    def setBoundariesDataList(self,xRange=[-np.inf,np.inf],yRange=[-np.inf,np.inf],zRange=[-np.inf,np.inf]):
        zminIndex=min(enumerate(self.zValuesFloat), key=lambda x: abs(x[1]-zRange[0]))[0]
        zmaxIndex=min(enumerate(self.zValuesFloat), key=lambda x: abs(x[1]-zRange[1]))[0]
        if zminIndex>zmaxIndex:
                a=zminIndex
                zminIndex=zmaxIndex
                zmaxIndex=a
        self.dataList=self.dataList[zminIndex:zmaxIndex] 
        self.zValues=self.zValues[zminIndex:zmaxIndex]
        self.zValuesFloat=self.zValuesFloat[zminIndex:zmaxIndex]
        yminIndex=min(enumerate(self.yValuesFloat), key=lambda x: abs(x[1]-yRange[0]))[0]            
        ymaxIndex=min(enumerate(self.yValuesFloat), key=lambda x: abs(x[1]-yRange[1]))[0]
        xminIndex=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-xRange[0]))[0]
        xmaxIndex=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-xRange[1]))[0]
        if xminIndex>xmaxIndex:
            b=xminIndex
            xminIndex=xmaxIndex
            xmaxIndex=b
        for i in range(len(self.dataList)):
            self.dataList[i]=self.dataList[i].iloc[xminIndex:xmaxIndex,yminIndex:ymaxIndex]
            
        self.data=self.data.iloc[xminIndex:xmaxIndex,yminIndex:ymaxIndex]
        self.xValuesFloat=self.xValuesFloat[xminIndex:xmaxIndex]
        self.xValues=self.xValues[xminIndex:xmaxIndex]   
        self.yValues=self.yValues[yminIndex:ymaxIndex]
        self.yValuesFloat=self.yValuesFloat[yminIndex:ymaxIndex]
        
    def normalizeData(self):
        #for ii in range(len(self.data.iloc[0])):
        self.data = self.data.div(self.data.max())
        #print self.data.iloc[:,ii]
    def plotSimple(self, xValues=None,yValues=None,zValues=None,xName=None,yName=None,zName=None,title=None):
        if len(yValues)==len(xValues) and xValues is not None:
            
                fig, intensity=plt.subplots(1,1,figsize=(12,12))
                
                cmap = plt.cm.jet
                if zValues is not None:
                    sc=intensity.scatter(xValues,yValues, c=zValues, vmin=min(zValues), vmax=max(zValues), s=60, cmap=cmap)
                    cbar=plt.colorbar(sc)
                    cbar.ax.set_ylabel(zName)
                else:
                    intensity.scatter(xValues,yValues, color='b',marker="o",s=60)
                intensity.set_xlabel(xName,fontsize=20)
                intensity.set_ylabel(yName,fontsize=20)
                if title is not None:                
                    intensity.set_title(title,fontsize=20)
                intensity.axis([min(xValues)-0.1*(max(xValues)-min(xValues)),max(xValues)+0.1*(max(xValues)-min(xValues)),min(yValues)-0.1*(max(yValues)-min(yValues)),max(yValues)+0.1*(max(yValues)-min(yValues))])
                
                plt.show()        
                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
                
    def plotData(self, seperateSpectra=False,xRange=[],yRange=[],fits=False):
        if len(yRange)==1 or len(yRange)>2 or len(yRange)==1 or len(yRange)>2:
            print "Enter Window ranges correctly: xRange=[x_min,x_max], yRange=[y_min,y_max]"
        
        if seperateSpectra:        
            for i in range(len(self.data.columns)):
            
                if not self.xValuesFloat is None:
                    fig, intensity=plt.subplots(1,1,figsize=(12,12))
                    intensity.plot(self.xValuesFloat,self.data.iloc[:,i], 'b', lw='2',label=self.yValues[i])
                    if fits:
                        intensity.plot(self.xValuesFloat,self.FitData.iloc[:,i], 'orange', lw='2',label="Fit "+self.yValues[i])
                    intensity.set_xlabel(self.xName+"("+self.xUnit+")")
                    intensity.set_ylabel('Intensity (arb. units)')
                    intensity.set_title('Single spectrum:')
                    if xRange and len(xRange)==2:
                        intensity.set_xlim(xRange)
                    if yRange and len(yRange)==2:
                        intensity.set_ylim(yRange)
                    intensity.legend()
                    plt.show()                
                elif not self.time is None:
                    fig, intensity=plt.subplots(1,1,figsize=(12,12))
                    intensity.plot(self.time,self.data.iloc[:,i], 'b', lw='2',label=self.yValues[i])
                    if fits:
                        intensity.plot(self.time,self.FitData.iloc[:,i], 'orange', lw='2',label="Fit "+self.yValues[i])
                    intensity.set_xlabel('Time (ns)')
                    intensity.set_ylabel('Intensity (arb. units)')
                    intensity.set_title('Single spectrum:')
                    if xRange and len(xRange)==2:
                        intensity.set_xlim(xRange)
                    if yRange and len(yRange)==2:
                        intensity.set_ylim(yRange)
                    intensity.legend()
                    plt.show()
                elif not self.wavelength is None:
                    fig, intensity=plt.subplots(1,1,figsize=(12,12))
                    intensity.plot(self.energy,self.data.iloc[:,i], 'b', lw='2',label=self.yValues[i])
                    if fits:
                        intensity.plot(self.energy,self.FitData.iloc[:,i], 'orange', lw='2',label="Fit "+self.yValues[i])
                    intensity.set_xlabel('Energy (eV)')
                    intensity.set_ylabel('Intensity (arb. units)')
                    intensity.set_title('Single spectrum:')
                    if xRange and len(xRange)==2:
                        intensity.set_xlim(xRange)
                    if yRange and len(yRange)==2:
                        intensity.set_ylim(yRange)
                    intensity.legend()
                    plt.show()
        
        else:
            if not self.xValuesFloat is None:
                cmap = plt.cm.jet
                fig, intensity=plt.subplots(1,1,figsize=(12,12))
                for i in range(len(self.data.columns)):
                    if len(self.data.columns)>1:
                        intensity.plot(self.xValuesFloat,self.data.iloc[:,i], color=cmap(i / float(len(self.data.columns)-1)), lw='2',label=self.yValues[i])
                        if fits:
                            intensity.plot(self.xValuesFloat,self.FitData.iloc[:,i], 'orange', lw='2',label="Fit "+self.yValues[i])
                    else:
                        intensity.plot(self.xValuesFloat,self.data.iloc[:,0], 'b', lw='2',label=self.yValues[0])
                        if fits:
                            intensity.plot(self.xValuesFloat,self.FitData.iloc[:,i], 'orange', lw='2',label="Fit "+self.yValues[i])
                intensity.set_xlabel(self.xName+"("+self.xUnit+")")
                intensity.set_ylabel('Intensity (arb. units)')
                intensity.set_title('Single spectrum:')
                if xRange and len(xRange)==2:
                        intensity.set_xlim(xRange)
                if yRange and len(yRange)==2:
                    intensity.set_ylim(yRange)
                
                plt.show()
            """
            elif not self.time is None:
                cmap = plt.cm.jet
                fig, intensity=plt.subplots(1,1,figsize=(12,12))
                for i in range(len(self.data.columns)):
                    if len(self.data.columns)>1:
                        intensity.plot(self.time,self.data.iloc[:,i], color=cmap(i / float(len(self.data.columns)-1)), lw='2',label=self.yValues[i])
                        if fits:
                            intensity.plot(self.time,self.FitData.iloc[:,i], 'orange', lw='2',label="Fit "+self.yValues[i])
                    else:
                        intensity.plot(self.time,self.data.iloc[:,0], 'b', lw='2',label=self.yValues[0])
                        if fits:
                            intensity.plot(self.time,self.FitData.iloc[:,i], 'orange', lw='2',label="Fit "+self.yValues[i])
                intensity.set_xlabel('Time (ns)')
                intensity.set_ylabel('Intensity (arb. units)')
                intensity.set_title('Single spectrum:')
                if xRange and len(xRange)==2:
                        intensity.set_xlim(xRange)
                if yRange and len(yRange)==2:
                    intensity.set_ylim(yRange)
                intensity.legend()
                plt.show()
            elif not self.wavelength is None:
                cmap = plt.cm.jet
                fig, intensity=plt.subplots(1,1,figsize=(12,12))
                for i in range(len(self.data.columns)):
                    if len(self.data.columns)>1:
                        intensity.plot(self.energy,self.data.iloc[:,i], color=cmap(i / float(len(self.data.columns)-1)), lw='2',label=self.yValues[i])
                        if fits:
                            intensity.plot(self.energy,self.FitData.iloc[:,i], 'orange', lw='2',label=self.yValues[i])
                    else:
                        
                        intensity.plot(self.energy,self.data.iloc[:,0], 'b', lw='2',label=self.yValues[0])
                        if fits:
                            intensity.plot(self.energy,self.FitData.iloc[:,i], 'orange', lw='2',label=self.yValues[i])
                intensity.set_xlabel('Energy (eV)')
                intensity.set_ylabel('Intensity (arb. units)')
                intensity.set_title('Single spectrum:')
                if xRange and len(xRange)==2:
                    intensity.set_xlim(xRange)
                if yRange and len(yRange)==2:
                    intensity.set_ylim(yRange)
                intensity.legend()
                plt.show()
            """
    def plotData2D(self,xRange=[],yRange=[],zRange=[],Colormap="viridis", data=None,fits=False,fitData=None):
        if Colormap=="Hot":        
            cm = plt.cm.afmhot
        elif Colormap=="Jet":
            cm = plt.cm.jet
        elif Colormap=="Seismic":
            cm = plt.cm.seismic
        elif Colormap=="Viridis":
            cm = plt.cm.viridis
        else:
            print " Available Colormaps: \"Jet\" , \"Hot\", \"Seismic\" "
            cm = plt.cm.jet
            
                 
        
        if data is None:
            data=np.array(self.data)
        
        
        fig, Contour = plt.subplots(1,1, figsize=(15,10), sharex=True)
        
        if zRange and len(zRange)==2:
            levels = np.linspace(zRange[0], zRange[1], 256)
        else:
            levels = np.linspace(np.array(data).min(), np.array(data).max(), 256)
        
        if self.xValuesFloat is not None:
            cs = Contour.contourf(self.xValuesFloat, self.yValuesFloat,data.T, levels=levels,cmap= cm)
            Contour.set_xlabel(self.xName+" ("+self.xUnit+")")
        
        if fits:
            if fitData is None:
                fitData=np.array(self.FitData)
            Contour.contour(self.xValuesFloat, self.yValuesFloat,fitData.T,8, colors='w')
        Contour.set_ylabel(self.yName+" ("+self.yUnit+")")
        Contour.set_title('2D Spectrum:')
        if xRange and len(xRange)==2:
            Contour.set_xlim(xRange)
        if yRange and len(yRange)==2:
            Contour.set_ylim(yRange)        
        fig.colorbar(cs, ax=Contour, format="%.2f")
    
    
    
    
    def saveDataTXT(self,data=None,path=None,filename=None,filename_extension="Data"):
        
        if path is None:
            path=self.path
        if filename is None:
            filename=self.filename.replace(self.format,"_"+filename_extension+".txt")
        
        if not self.time is None:        
            Header_line1 = "Time "
            Header_line2 = "ns  "
            Xaxis=self.time
        if not self.energy is None:        
            Header_line1 = "Energy "
            Header_line2 = "eV  "
            Xaxis=self.energy
        Header_line3=self.yName
        for i in range(len(self.yValues)):
            Header_line3=Header_line3+"\t"+str(self.yValuesFloat[i])
            Header_line2=Header_line2+"\t"+str(self.yUnit)
            Header_line1=Header_line1+"\t"+str(self.yValues[i])
        if data is None:
            data=self.data

        if os.path.isfile(path+filename):
            UserInput=raw_input("Filename already exists! Overwrite File? [y/n]")
            if UserInput=='y':
                np.savetxt(path+filename,np.vstack((Xaxis,np.array(data).T)).T,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
                print "Data saved"
            else:
                print "Specify other name"    
        else:
            np.savetxt(path+filename,np.vstack((Xaxis,np.array(data).T)).T,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
            print "Data saved"
    
    
    
    def saveParametersTXT(self,data=None,error=None,path=None,filename=None,filename_extension="Fit"):
        # data Multiindex panda
        Data_present=True
        if data is None:
            if self.FitParameters is not None and self.FitParametersError is not None:
                data=self.FitParameters
                error=self.FitParametersError
            else:
                print "Specify Data to save!"
                Data_present=False
        
        if Data_present:
            
            
            if filename is None:
                FileNameTrue=True
                
            ### LOOP over all (FIT) Parameters in Panda File (Multiindes line 0)
            parameterNumber=0        
            for parameter in data.columns.levels[0]:
                 ###Define Filename  
                if path is None:
                    path=self.path
                if FileNameTrue:
                    filename=self.filename.replace(self.format,"_"+filename_extension+"_"+parameter+".txt")
                else:
                    filename=self.filename+parameter+".txt"
                
                ###Define Header for TXT files    
                        
                Header_line1 = self.yName+"\t"+parameter
                if error is not None:
                    Header_line1=Header_line1+"\t"+error.columns.levels[0][parameterNumber]
                Header_line2 = self.yUnit
                
                Header_line3="Fit Number"
                
                
                for i in range(len(data.columns.levels[1])):
                    Header_line3=Header_line3+"\t"+data.columns.levels[1][i]
                    if error is not None:
                        Header_line3=Header_line3+"\t"+error.columns.levels[1][i]
                
                
                
                ###Define Array to save
                        
                if error is None:
                    for i in range(len(data.columns.levels[1])):
                        if i ==0:
                            saveData=data[parameter].iloc[:,i].values.tolist()
                            saveData=[saveData]
                        else:
                            saveData.append(data[parameter].iloc[:,i].values.tolist())
                else:
                    for i in range(len(data.columns.levels[1])):
                        if i ==0:
                            saveData=data[parameter].iloc[:,i].values.tolist()
                            saveData=[saveData]
                            saveData.append(error[error.columns[parameterNumber][0]].iloc[:,i].values.tolist())
                        else:
                            saveData.append(data[parameter].iloc[:,i].values.tolist())
                            saveData.append(error[error.columns[parameterNumber][0]].iloc[:,i].values.tolist())
                X=np.array(data.index)
                X=X.astype(np.float)
                print saveData
                if os.path.isfile(path+filename):
                    UserInput=raw_input("Filename already exists! Overwrite File? [y/n]")
                    if UserInput=='y':
                        np.savetxt(path+filename,np.vstack((X,np.array(saveData))).T,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
                        print "Parameter Data saved: "+parameter
                    else:
                        print "Specify other name"    
                else:
                    np.savetxt(path+filename,np.vstack((X,np.array(saveData))).T,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
                    print "Parameter Data saved: "+parameter
        
                parameterNumber+=1
                
                
                
                
    def plotParameters(self,data=None,error=None,style="scatter"):
        """
        style:
        scatter,polar,line,bar
        
        """
        Data_present=True
        if data is None:
            if self.FitParameters is not None and self.FitParametersError is not None:
                data=self.FitParameters
                error=self.FitParametersError
            else:
                print "Specify Data to save!"
                Data_present=False
        
        """
        if Data_present:
            data=data.fillna(0)#replace NaN by 0
            if error is not None:
                error=error.fillna(0)        
        """
        fontsize=20
        if Data_present:
            parameterNumber=0
            fig, Plot=plt.subplots(1,len(data.columns.levels[0]),figsize=(6*len(data.columns.levels[0]),8), sharex=True)
            for parameter in data.columns.levels[0]:

                for i in range(len(data.columns.levels[1])):
                    
                    Y=data[parameter].iloc[:,i].tolist()
                    X=np.array(data.index)
                    X=X.astype(np.float).tolist()
                    if style=="scatter":
                        if len(data.columns.levels[0])==1:                   
                            if error is None:
                                Plot.plot(X,Y,marker="o",linestyle="None",label=data.columns.levels[1][i])
                                Plot.set_xlabel(data.index.name,fontsize=fontsize)
                                Plot.set_ylabel(parameter,fontsize=fontsize)
                                Plot.set_title(parameter,fontsize=fontsize)
                                plt.legend()
                                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
                            else:
                                Yerr=error[error.columns.levels[0][parameterNumber]].iloc[:,i].tolist()
                                Plot.errorbar(X,Y,yerr=Yerr,fmt="o",label=str(i))
                                Plot.set_xlabel(data.index.name,fontsize=fontsize)
                                Plot.set_ylabel(parameter,fontsize=fontsize)
                                Plot.set_title(parameter,fontsize=fontsize)
                                plt.legend()
                                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
                        
                        else:
                            if error is None:
                                Plot[parameterNumber].scatter(X,Y,label=data.columns.levels[1][i])
                                Plot[parameterNumber].set_xlabel(data.index.name,fontsize=fontsize)
                                Plot[parameterNumber].set_ylabel(parameter,fontsize=fontsize)
                                Plot[parameterNumber].set_title(parameter,fontsize=fontsize)
                                plt.legend()
                                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
                            else:
                                Yerr=error[error.columns.levels[0][parameterNumber]].iloc[:,i].tolist()
                                Plot[parameterNumber].errorbar(X,Y,yerr=Yerr,fmt="o",label=str(i))
                                Plot[parameterNumber].set_xlabel(data.index.name,fontsize=fontsize)
                                Plot[parameterNumber].set_ylabel(parameter,fontsize=fontsize)
                                Plot[parameterNumber].set_title(parameter,fontsize=fontsize)
                                plt.legend()
                                plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
                    
        
                parameterNumber+=1
            
class analyseTime(analyseData):
    def __init__(self):
        analyseData.__init__(self)
        self.Decay1=None
    def setZero(self,zero=None):
        if zero==None:
            t0=self.time.iloc[self.data.iloc[:,0].idxmax()]
            self.time=self.time-t0
        elif isinstance(zero,float):
            self.time=self.time-zero
        elif isinstance(zero,int):
            self.time=self.time-zero
        else:
            print "Provide Zero"
            
    
    def fitLifetimeMono(self,tmin=0,tmax=12,tau1=1,Offset=None,data=None):
        if data==None:
            data=self.data
        
        #Create MultiIndex Panda
        if Offset==None:
            iterables=[["Decay Time","Intensity"],["1"]]
            iterablesError=[["Decay Time Error","Intensity Error"],["1"]]
            columns=pd.MultiIndex.from_product(iterables)
            columnsError=pd.MultiIndex.from_product(iterablesError)
            index = [str(i) for i in self.yValuesFloat]
        elif isinstance(Offset, int) or isinstance(Offset, float):
            iterables=[["Decay Time","Intensity","Offset"],["1"]]
            iterablesError=[["Decay Time Error","Intensity Error", "Offset Error"],["1"]]
            columns=pd.MultiIndex.from_product(iterables)
            columnsError=pd.MultiIndex.from_product(iterablesError)
            index = [str(i) for i in self.yValuesFloat]
        else: 
            print " Specify Offset"
        #self.FitParameters=pd.DataFrame(np.zeros((len(peakPositions),len(peakPositions.iloc[0])*4)), index=index, columns=columns)
        self.FitParameters=pd.DataFrame(index=index,columns=columns)
        self.FitParametersError=pd.DataFrame(index=index,columns=columnsError)
        self.FitData=pd.DataFrame(index=self.xValuesFloat,columns=self.yValuesFloat)
        
        for spectrum in range(len(data.columns)):
            
                
            ### Specify Fit Range xmin xmax for fit
            
            xmin=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-tmin))[1]
            xmax=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-tmax))[1]
            
                
            
                
            tminIndex=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-xmin))[0]
            tmaxIndex=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-xmax))[0]
            
            if tminIndex > tmaxIndex:
                a=tminIndex
                tminIndex=tmaxIndex
                tmaxIndex=a
                
            Y=self.data.iloc[tminIndex:tmaxIndex,spectrum].tolist()
            X=self.xValuesFloat.iloc[tminIndex:tmaxIndex].values.tolist()
            
            
            #Fit parameters
            intensity1=1
            if Offset==None:
                parameters=[tau1,intensity1]
                fit, best_parameters,err=MM.MonoExpFit_noOffset(Y,X,parameters)
                
            else:
                parameters=[tau1,intensity1,Offset]
                fit, best_parameters,err=MM.MonoExpFit_wOffset(Y,X,parameters)
                
            
            self.FitParameters.ix[spectrum,("Decay Time","1")]=best_parameters[0]
            self.FitParameters.ix[spectrum,("Intensity","1")]=best_parameters[1]            
            
            self.FitParametersError.ix[spectrum,("Decay Time Error","1")]=err[0]
            self.FitParametersError.ix[spectrum,("Intensity Error","1")]=err[1]
            if Offset!=None:
                self.FitParameters.ix[spectrum,("Offset","1")]=best_parameters[2]
                self.FitParametersError.ix[spectrum,("Offset Error","1")]=err[2]
            
            
            currentFit=pd.DataFrame(fit,index=X,columns=[self.yValuesFloat[spectrum]])
            self.fit=pd.DataFrame(fit,index=X)
            self.FitData.update(currentFit,overwrite=True)
        self.FitParameters.index.name=self.yName+" ("+self.yUnit+")"
                  
        self.plotParameters(data=self.FitParameters,error=self.FitParametersError)
        self.plotData(fits=True)

    def fitLifetimeBi(self,tmin=0,tmax=12,tau1=1,tau2=5,Offset=None,data=None):
        if data==None:        
            data=self.data
        
        #Create MultiIndex Panda
        if Offset==None:
            iterables=[["Decay Time 1","Decay Time 2","Intensity 1","Intensity 2"],["1"]]
            iterablesError=[["Decay Time 1 Error","Decay Time 2 Error","Intensity 1 Error","Intensity 2 Error"],["1"]]
            columns=pd.MultiIndex.from_product(iterables)
            columnsError=pd.MultiIndex.from_product(iterablesError)
            index = [str(i) for i in self.yValuesFloat]
        elif isinstance(Offset, int) or isinstance(Offset, float):
            iterables=[["Decay Time 1","Decay Time 2","Intensity 1","Intensity 2","Offset"],["1"]]
            iterablesError=[["Decay Time 1 Error","Decay Time 2 Error","Intensity 1 Error","Intensity 2 Error", "Offset Error"],["1"]]
            columns=pd.MultiIndex.from_product(iterables)
            columnsError=pd.MultiIndex.from_product(iterablesError)
            index = [str(i) for i in self.yValuesFloat]
        else: 
            print " Specify Offset"
        #self.FitParameters=pd.DataFrame(np.zeros((len(peakPositions),len(peakPositions.iloc[0])*4)), index=index, columns=columns)
        self.FitParameters=pd.DataFrame(index=index,columns=columns)
        self.FitParametersError=pd.DataFrame(index=index,columns=columnsError)
        self.FitData=pd.DataFrame(index=self.xValuesFloat,columns=self.yValuesFloat)
        
        for spectrum in range(len(data.columns)):
            
                
            ### Specify Fit Range xmin xmax for fit
            
            xmin=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-tmin))[1]
            xmax=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-tmax))[1]
            
                
            
                
            tminIndex=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-xmin))[0]
            tmaxIndex=min(enumerate(self.xValuesFloat), key=lambda x: abs(x[1]-xmax))[0]
            
            if tminIndex > tmaxIndex:
                a=tminIndex
                tminIndex=tmaxIndex
                tmaxIndex=a
                
            Y=self.data.iloc[tminIndex:tmaxIndex,spectrum].tolist()
            X=self.xValuesFloat.iloc[tminIndex:tmaxIndex].values.tolist()
            
            
            #Fit parameters
            intensity1=1
            intensity2=1
            if Offset==None:
                parameters=[tau1,tau2,intensity1,intensity2]
                fit, best_parameters,err=MM.BiExpFit_noOffset(Y,X,parameters)
                
            else:
                parameters=[tau1,tau2,intensity1,intensity2,Offset]
                fit, best_parameters,err=MM.BiExpFit_wOffset(Y,X,parameters)
                
            
            self.FitParameters.ix[spectrum,("Decay Time 1","1")]=best_parameters[0]
            self.FitParameters.ix[spectrum,("Decay Time 2","1")]=best_parameters[1]
            self.FitParameters.ix[spectrum,("Intensity 1","1")]=best_parameters[2]
            self.FitParameters.ix[spectrum,("Intensity 2","1")]=best_parameters[3]
            
            
            self.FitParametersError.ix[spectrum,("Decay Time 1 Error","1")]=err[0]
            self.FitParametersError.ix[spectrum,("Decay Time 2 Error","1")]=err[1]
            self.FitParametersError.ix[spectrum,("Intensity 1 Error","1")]=err[2]
            self.FitParametersError.ix[spectrum,("Intensity 1 Error","1")]=err[3]
            if Offset!=None:
                self.FitParameters.ix[spectrum,("Offset","1")]=best_parameters[4]
                self.FitParametersError.ix[spectrum,("Offset Error","1")]=err[4]
            
            
            currentFit=pd.DataFrame(fit,index=X,columns=[self.yValuesFloat[spectrum]])
            self.fit=pd.DataFrame(fit,index=X)
            self.FitData.update(currentFit,overwrite=True)
        self.FitParameters.index.name=self.yName+" ("+self.yUnit+")"
                  
        self.plotParameters(data=self.FitParameters,error=self.FitParametersError)
        self.plotData(fits=True)

      
class analysePeak(analyseData):
    def __init__(self):
        analyseData.__init__(self)
        self.PeakPositionMax=None
        self.PeakPositionMin=None
        self.intCounts=None
        
    def integrateCounts(self,Emin=[],Emax=[],Normalize=True,data=None,yValues=None):
        if data is None:
            data=self.data
        if len(Emin)!=len(Emax):
            print "Size of lists have to match" 
            return
        if not Emin:
            Emin=[0]
            Emax=[10000]
        if yValues is None:
                yValues=self.yValues
                
        index = [str(i) for i in self.yValuesFloat]
        iterables=[["Integrated Counts"],[str(Emin[i])+" - "+str(Emax[i]) for i in range(len(Emin))]]
        columns=pd.MultiIndex.from_product(iterables)
        self.intCounts=pd.DataFrame(index=index,columns=columns)
        
        for i in range(len(Emin)):        
            xlow=min(enumerate(self.energy), key=lambda x: abs(x[1]-Emin[i]))[0]
            xhigh=min(enumerate(self.energy), key=lambda x: abs(x[1]-Emax[i]))[0]
            if xlow>xhigh:
                a=xhigh
                xhigh=xlow
                xlow=a
            
            for j in range(len(yValues)):
                if Normalize:
                    self.intCounts.ix[j,("Integrated Counts",str(Emin[i])+" - "+str(Emax[i]))]=self.data.iloc[xlow:xhigh,j].sum()/float((xhigh-xlow))
                else:
                    self.intCounts.ix[j,("Integrated Counts",str(Emin[i])+" - "+str(Emax[i]))]=self.data.iloc[xlow:xhigh,j].sum()
        self.intCounts.index.name=self.yName+" ("+self.yUnit+")"       

        self.plotParameters(data=self.intCounts)
        
    
    def findPeakPosition(self,sensitivity=50,minPeakDistance=0.2, minPeakHeight=10, E_min=0, E_max=np.inf, data=None):
        maxima=[]
        minima=[]
        if data is None:
            data=self.data
        for i in range(len(self.yValues)):
            maximumEnergy=[]
            minimumEnergy=[]
            maximum, minimum= MM.findpeaks2(self.data.iloc[:,i],lookahead=sensitivity,delta=minPeakDistance,MinPeakheight=minPeakHeight)
                        
            if maximum:            
                for j in range(len(maximum)):
                    if self.energy[maximum[j][0]] < E_max and self.energy[maximum[j][0]] > E_min:
                        maximumEnergy.append(self.energy[maximum[j][0]])
                maxima.append(maximumEnergy)
            if minimum:            
                for j in range(len(minimum)):
                    if self.energy[minimum[j][0]] < E_max and self.energy[minimum[j][0]] > E_min:
                        minimumEnergy.append(self.energy[minimum[j][0]])  
                minima.append(minimumEnergy)
        
        index = [str(i) for i in self.yValuesFloat]
        if maxima:
            iterables=[["PeakPositionMax"],[str(i) for i in range(len(max(maxima,key=len)))]]
            columns=pd.MultiIndex.from_product(iterables)
            self.PeakPositionMax=pd.DataFrame(data=maxima,index=index,columns=columns)
            self.PeakPositionMax.index.name=self.yName+" ("+self.yUnit+")"
            self.plotParameters(self.PeakPositionMax)
        else:
            self.PeakPositionMax=None
        
        

        if minima:
            iterables=[["PeakPositionMin"],[str(i) for i in range(len(max(minima,key=len)))]]
            columns=pd.MultiIndex.from_product(iterables)
            self.PeakPositionMin=pd.DataFrame(data=minima,index=index,columns=columns)
            self.PeakPositionMin.index.name=self.yName+" ("+self.yUnit+")"
        else:
            self.PeakPositionMin=None
    
        
       
        
        
        
        
        
    def fitAllPeaks(self, function="Lorentz",FWHM=0.01,Intensity=100,Offset=0, E_min=0, E_max=np.inf,data=None, peakPositions=None):
        
        if peakPositions is None:
            if self.PeakPositionMax is None:
                self.findPeakPosition(E_min=E_min, E_max=E_max)
            peakPositions=self.PeakPositionMax
        
        if data is None:
            data=self.data
        
        
        
        Emin=min(enumerate(self.energy), key=lambda x: abs(x[1]-E_min))[1]
        
        Emax=min(enumerate(self.energy), key=lambda x: abs(x[1]-E_max))[1]
        
        #Create MultiIndex Panda
        
        iterables=[["Energy","FWHM","Intensity","Offset"],[str(i) for i in range(len(peakPositions.iloc[0]))]]
        iterablesError=[["Energy Error","FWHM Error","Intensity Error","Offset Error"],[str(i) for i in range(len(peakPositions.iloc[0]))]]
        columns=pd.MultiIndex.from_product(iterables)
        columnsError=pd.MultiIndex.from_product(iterablesError)        
        
        index = [str(i) for i in peakPositions.index]
        
        #self.FitParameters=pd.DataFrame(np.zeros((len(peakPositions),len(peakPositions.iloc[0])*4)), index=index, columns=columns)
        self.FitParameters=pd.DataFrame(index=index,columns=columns)
        self.FitParametersError=pd.DataFrame(index=index,columns=columnsError)
        self.FitData=pd.DataFrame(index=self.energy,columns=self.yValues)
        for spectrum in range(len(peakPositions)):
            for peak in range(len(peakPositions.iloc[0])):
                
                
                ### Specify Fit Range xmin xmax for fit
                
                if peak==0:
                    xmin=Emin
                else:
                    if not math.isnan(peakPositions.iloc[spectrum,peak]) and not math.isnan(peakPositions.iloc[spectrum,peak-1]):
                        xmin =(peakPositions.iloc[spectrum,peak]+peakPositions.iloc[spectrum,peak-1])/2.0
                        
                
                
                if peak==len(peakPositions.iloc[0])-1:
                    xmax=Emax
                else:
                    if math.isnan(peakPositions.iloc[spectrum,peak+1]):
                        xmax=Emax
                        
                        
                    else:
                        if not math.isnan(peakPositions.iloc[spectrum,peak]):
                            
                            xmax =(peakPositions.iloc[spectrum,peak+1]+peakPositions.iloc[spectrum,peak])/2.0
                
                if not math.isnan(peakPositions.iloc[spectrum,peak]):
                    
                    EminIndex=min(enumerate(self.energy), key=lambda x: abs(x[1]-xmin))[0]
                    EmaxIndex=min(enumerate(self.energy), key=lambda x: abs(x[1]-xmax))[0]
                    
                    if EminIndex > EmaxIndex:
                        a=EminIndex
                        EminIndex=EmaxIndex
                        EmaxIndex=a
                        
                    Y=self.data.iloc[EminIndex:EmaxIndex,spectrum].tolist()
                    X=self.energy.iloc[EminIndex:EmaxIndex].tolist()
                    
                    #Fit parameters
                    parameters=[FWHM,peakPositions.iloc[spectrum,peak],Intensity,Offset]
                    
                    if function == "Lorentz":
                        fit, best_parameters,err=MM.lorentzfit(Y,X,parameters)
                        
                    elif function == "Gauss":
                        fit, best_parameters,err=MM.gaussfit(Y,X,parameters)
                    
                    self.FitParameters.ix[spectrum,("FWHM",str(peak))]=best_parameters[0]
                    self.FitParameters.ix[spectrum,("Energy",str(peak))]=best_parameters[1]
                    self.FitParameters.ix[spectrum,("Intensity",str(peak))]=best_parameters[2]
                    self.FitParameters.ix[spectrum,("Offset",str(peak))]=best_parameters[3]
                    
                    self.FitParametersError.ix[spectrum,("FWHM Error",str(peak))]=err[0]
                    self.FitParametersError.ix[spectrum,("Energy Error",str(peak))]=err[1]
                    self.FitParametersError.ix[spectrum,("Intensity Error",str(peak))]=err[2]
                    self.FitParametersError.ix[spectrum,("Offset Error",str(peak))]=err[3]
                    currentFit=pd.DataFrame(fit,index=X,columns=[self.yValues[spectrum]])
                    self.fit=pd.DataFrame(fit,index=X)
                    self.FitData.update(currentFit,overwrite=True)
        self.FitParameters.index.name=self.yName+" ("+self.yUnit+")"
                  
        self.plotParameters(data=self.FitParameters,error=self.FitParametersError)   
    
    def fitDefinedPeak2D(self,function="circularGauss",Peak=[0,0],Intensity=1,FWHM=1,Offset=None,data=None,Plot=False):
        if data is None:
            data=np.array(self.data.values.tolist())
        x,y=np.meshgrid(self.xValuesFloat,self.yValuesFloat)
        
        if function == "circularGauss" and Offset is None:
            parameters=[Peak[0],Peak[1],FWHM,Intensity*np.pi*FWHM**2]
            fit, best_parameters,err=MM.circGauss2DFit_noOffset(data,y,x,parameters)
        if function == "circularGauss" and isinstance(Offset, int) or isinstance(Offset, float):
            parameters=[Peak[0],Peak[1],FWHM,Intensity*np.pi*FWHM**2,Offset]
            fit, best_parameters,err=MM.circGauss2DFit_wOffset(data,y,x,parameters)
        
        
        fit=fit.reshape(len(self.xValuesFloat),len(self.yValuesFloat))
       
        currentFit=pd.DataFrame(fit,index=self.xValuesFloat,columns=[self.yValuesFloat])
        self.fit=pd.DataFrame(fit,index=self.xValuesFloat)
        self.FitData.update(currentFit,overwrite=True)
        if Plot:
            self.plotData2D(data=data,fits=True)
        return fit, best_parameters,err
        ############################################################################
        #####    FIT 2D
        ###################################
        
    def fitDefinedPeak2DList(self,function="circularGauss",Peak=[0,0],Intensity=1,FWHM=1,Offset=None,dataList=None,Plot=False, UseFitParameters=False, WeightedsquaredResidual=False):
        if dataList is None:
            dataList=self.dataList
        
        #Create MultiIndex Panda
        if Offset==None:
            iterables=[[self.xName,self.yName,"FWHM","Intensity"],["1"]]
            iterablesError=[[self.xName+" Error",self.yName+" Error","FWHM Error","Intensity Error"],["1"]]
            columns=pd.MultiIndex.from_product(iterables)
            columnsError=pd.MultiIndex.from_product(iterablesError)
            index = [str(i) for i in self.zValuesFloat]
        elif isinstance(Offset, int) or isinstance(Offset, float):
            iterables=[[self.xName,self.yName,"FWHM","Intensity","Offset"],["1"]]
            iterablesError=[[self.xName+" Error",self.yName+" Error","FWHM Error","Intensity Error", "Offset Error"],["1"]]
            columns=pd.MultiIndex.from_product(iterables)
            columnsError=pd.MultiIndex.from_product(iterablesError)
            index = [str(i) for i in self.zValuesFloat]
        else: 
            print " Specify Offset"
        #self.FitParameters=pd.DataFrame(np.zeros((len(peakPositions),len(peakPositions.iloc[0])*4)), index=index, columns=columns)
        self.FitParameters=pd.DataFrame(index=index,columns=columns)
        self.FitParametersError=pd.DataFrame(index=index,columns=columnsError)
        self.FitData=pd.DataFrame(index=self.xValuesFloat,columns=self.yValuesFloat)
        self.FitDataList=[]
        
        Peak0=Peak
        Intensity0=Intensity
        FWHM0=FWHM
        Offset0=Offset
        
        for i in range(len(dataList)):
            data=np.array(dataList[i].values.tolist())
            fit,best_parameters,err = self.fitDefinedPeak2D(function,Peak0,Intensity0,FWHM0,Offset0,data=data,Plot=Plot)
            
            if WeightedsquaredResidual:
                data=np.array(dataList[i].values-self.FitData.values)
                #data[data<0.01]=0.01
                data=np.square(data)
                data[data<1]=1
                data=1./data              
                data=np.multiply(dataList[i].values,data)+dataList[i].values
                MaxPeak=np.unravel_index(data.argmax(),data.shape)
                data=np.array(data.tolist())
                fit,best_parameters,err = self.fitDefinedPeak2D(function,[self.xValuesFloat[MaxPeak[0]],self.yValuesFloat[MaxPeak[1]]],data[MaxPeak[0],MaxPeak[1]],FWHM0,Offset0,data=data,Plot=Plot)
                
            self.FitParameters.ix[i,(self.xName,"1")]=best_parameters[0]
            self.FitParameters.ix[i,(self.yName,"1")]=best_parameters[1]
            self.FitParameters.ix[i,("FWHM","1")]=np.abs(best_parameters[2])
            self.FitParameters.ix[i,("Intensity","1")]=best_parameters[3]
            
            
            self.FitParametersError.ix[i,(self.xName+" Error","1")]=err[0]
            self.FitParametersError.ix[i,(self.yName+" Error","1")]=err[1]
            self.FitParametersError.ix[i,("FWHM Error","1")]=err[2]
            self.FitParametersError.ix[i,("Intensity Error","1")]=err[3]
            
            if UseFitParameters:
                Peak0=[best_parameters[0],best_parameters[1]]
                Intensity0=best_parameters[3]
                FWHM0=best_parameters[2]
                if Offset!=None:
                    Offset0=best_parameters[4]
            if Offset!=None:
                self.FitParameters.ix[i,("Offset","1")]=best_parameters[4]
                self.FitParametersError.ix[i,("Offset Error","1")]=err[4]
            
            self.FitDataList.append(self.FitData)
            #fit=fit.reshape(len(self.xValuesFloat),len(self.yValuesFloat))
            #currentFit=pd.DataFrame(fit,index=self.xValuesFloat,columns=[self.yValuesFloat])
            #self.fit=pd.DataFrame(fit,index=self.xValuesFloat)
            #self.FitData.update(currentFit,overwrite=True)
        self.FitParameters.index.name=self.zName+" ("+self.zUnit+")"
                  
        self.plotParameters(data=self.FitParameters,error=self.FitParametersError)
        self.plotSimple(xValues=self.FitParameters.iloc[:,0].tolist(),yValues=self.FitParameters.iloc[:,1].tolist(),zValues=self.zValuesFloat,xName="x",yName="y",zName="Time (s)",title="x-y (Time)")
        self.plotSimple(xValues=self.FitParameters.iloc[:,0].tolist(),yValues=self.FitParameters.iloc[:,1].tolist(),zValues=self.FitParameters.iloc[:,3],xName="x",yName="y",zName="Intensity",title="x-y (Time)")
        self.plotSimple(xValues=self.FitParameters.iloc[:,0].tolist(),yValues=self.FitParameters.iloc[:,1].tolist(),xName="x",yName="y",title="x-y (Time)")
        
            
    def fitDefinedPeaks(self,function="Lorentz",Energy=[],FWHM=[],Intensity=[],Offset=[], E_min=[], E_max=[],data=None,yValues=None):
        if data is None:
            data=self.data
        if yValues is None:
            yValues=self.yValues   
        if not len(Energy)==len(FWHM) and not len(Energy)==len(Intensity) and not len(Energy)==len(Offset):
            print "Specify all Parameters correctly"
        
            
        
        #Create MultiIndex Panda
        
        iterables=[["Energy","FWHM","Intensity","Offset"],[str(i) for i in range(len(Energy))]]
        iterablesError=[["Energy Error","FWHM Error","Intensity Error","Offset Error"],[str(i) for i in range(len(Energy))]]
        columns=pd.MultiIndex.from_product(iterables)
        columnsError=pd.MultiIndex.from_product(iterablesError)        
        
        index = [str(i) for i in self.yValuesFloat]
        self.FitDataList=[]
        #self.FitParameters=pd.DataFrame(np.zeros((len(peakPositions),len(peakPositions.iloc[0])*4)), index=index, columns=columns)
        self.FitParameters=pd.DataFrame(index=index,columns=columns)
        self.FitParametersError=pd.DataFrame(index=index,columns=columnsError)
        self.FitData=pd.DataFrame(index=self.energy,columns=self.yValuesFloat)
        for spectrum in range(len(data.columns)):
            for peak in range(len(Energy)):
                
                ### Specify Fit Range xmin xmax for fit
                if len(Energy)==len(E_min) and len(Energy)==len(E_max):
                    xmin=min(enumerate(self.energy), key=lambda x: abs(x[1]-E_min[peak]))[1]
                    xmax=min(enumerate(self.energy), key=lambda x: abs(x[1]-E_max[peak]))[1]
                
                
                elif len(E_min)==1 and len(E_max)==1:
                    Emin=min(enumerate(self.energy), key=lambda x: abs(x[1]-E_min[peak]))[1]
                    Emax=min(enumerate(self.energy), key=lambda x: abs(x[1]-E_max[peak]))[1]
                    if peak==0:
                        xmin=Emin
                    else:
                        xmin =(Energy[peak]+Energy[peak-1])/2.0
                            
                    
                    
                    if peak==len(Energy)-1:
                        xmax=Emax
                    else:
                        xmax =(Energy[peak+1]+Energy[peak])/2.0
                        
                                
                                
                else:
                    print "Enter E_min and E_max correctly"
                    
                
                    
                EminIndex=min(enumerate(self.energy), key=lambda x: abs(x[1]-xmin))[0]
                EmaxIndex=min(enumerate(self.energy), key=lambda x: abs(x[1]-xmax))[0]
                
                if EminIndex > EmaxIndex:
                    a=EminIndex
                    EminIndex=EmaxIndex
                    EmaxIndex=a
                    
                Y=self.data.iloc[EminIndex:EmaxIndex,spectrum].tolist()
                X=self.energy.iloc[EminIndex:EmaxIndex].tolist()
                
                #Fit parameters
                parameters=[FWHM[peak],Energy[peak],Intensity[peak],Offset[peak]]
                
                if function == "Lorentz":
                    fit, best_parameters,err=MM.lorentzfit(Y,X,parameters)
                    
                elif function == "Gauss":
                    fit, best_parameters,err=MM.gaussfit(Y,X,parameters)
                
                self.FitParameters.ix[spectrum,("FWHM",str(peak))]=best_parameters[0]
                self.FitParameters.ix[spectrum,("Energy",str(peak))]=best_parameters[1]
                self.FitParameters.ix[spectrum,("Intensity",str(peak))]=best_parameters[2]
                self.FitParameters.ix[spectrum,("Offset",str(peak))]=best_parameters[3]
                
                self.FitParametersError.ix[spectrum,("FWHM Error",str(peak))]=err[0]
                self.FitParametersError.ix[spectrum,("Energy Error",str(peak))]=err[1]
                self.FitParametersError.ix[spectrum,("Intensity Error",str(peak))]=err[2]
                self.FitParametersError.ix[spectrum,("Offset Error",str(peak))]=err[3]
                
                
                currentFit=pd.DataFrame(fit,index=X,columns=[self.yValuesFloat[spectrum]])
                self.fit=pd.DataFrame(fit,index=X)
                self.FitData.update(currentFit,overwrite=True)
                self.FitDataList.append(self.data)
        self.FitParameters.index.name=self.yName+" ("+self.yUnit+")"
                 
        self.plotParameters(data=self.FitParameters,error=self.FitParametersError)
           
class analyse2D(analysePeak,analyseTime):
    
    def __init__(self):
        analysePeak.__init__(self)
        analyseTime.__init__(self)
        self.averageSubtracted=None
        self.dataFFT=None
        #self.header=[]
        return            
    def plotSubtractAverage(self,xRange=[],yRange=[],zRange=[]):
        self.averageSubtracted=self.data.subtract(self.data.T.mean(),axis='index')
        self.plotData2D(xRange,yRange,zRange,Colormap="Seismic",data=self.averageSubtracted)
    
    def saveSubtractAverage(self,path=None,filename=None,filename_extension="_AverageSubtracted"):
        if self.averageSubtracted is None:        
            self.averageSubtracted=self.data.subtract(self.data.T.mean(),axis='index')
        self.saveDataTXT(self.averageSubtracted,path,filename,filename_extension)
    
    def fft2DAbs(self,data=None,convertAxes=True):
        if data is None:
            data=self.data
        self.data=pd.DataFrame(np.abs(np.fft.fft2(self.data)))
        if convertAxes:
            self.xValuesFloat=2*np.pi/self.xValuesFloat
            self.yValuesFloat=2*np.pi/self.yValuesFloat
        return self.data   
        
    def fft2DAbsList(self,dataList=None):
        if dataList is None:
            dataList=self.dataList
        for i in range(len(dataList)):
            if i==0:
                dataList[i]=self.fft2DAbs(data=dataList[i],convertAxes=True)
            else:
                dataList[i]=self.fft2DAbs(data=dataList[i],convertAxes=False)
        self.dataList=dataList
               


class analyseFitParameters(analyse2D):
    def __init__(self):
        analyse2D.__init__(self)
        self.SmoothedFitParameters=None
    def smoothParameters(self,smoothfilter="ma",windowsize=10,FitParameter=1,SubtractSmooth=False):
    
        #Error???
               
            #if self.SmoothedFitParameters is None:            
            #    self.SmoothedFitParameters=self.FitParameters
        SmoothedData,smooth=MM.smooth(smoothfilter,self.FitParameters.iloc[:,FitParameter].tolist(),windowsize)
        #if SubtractSmooth:
        SmoothedData_Subtract=np.subtract(self.FitParameters.iloc[:,FitParameter].tolist(),SmoothedData.tolist())
        #else:
        SmoothedData=SmoothedData.tolist()
         #   self.SmoothedFitParameters.iloc[:,FitParameter]=SmoothedData
        return SmoothedData,SmoothedData_Subtract
    
    
    def smoothParameters2D(self,smoothfilter="ma",windowsize=10,x=0,y=1):
        self.SmoothedFitParameters=self.FitParameters
        
        SmoothX,SmoothX_sub=self.smoothParameters(smoothfilter,windowsize,x,True)
        SmoothY,SmoothY_sub=self.smoothParameters(smoothfilter,windowsize,y,True)
        
        if self.zValuesFloat is not None:        
            self.plotSimple(xValues=SmoothX,yValues=SmoothY,zValues=self.zValuesFloat,xName="x",yName="y",zName="Time",title=smoothfilter)
        elif self.SmoothedFitParameters.iloc[:,3] is not None:
            self.plotSimple(xValues=SmoothX,yValues=SmoothY,zValues=self.SmoothedFitParameters.iloc[:,3],xName="x",yName="y",zName="Intensity",title=smoothfilter)
        else:
            self.plotSimple(xValues=SmoothX,yValues=SmoothY,xName="x",yName="y",title=smoothfilter)
        
        #self.SmoothedFitParameters=None
        
        #self.smoothParameters(smoothfilter,windowsize,x,True)
        #self.smoothParameters(smoothfilter,windowsize,y,True)
        if self.zValuesFloat is not None:        
            self.plotSimple(xValues=SmoothX_sub,yValues=SmoothY_sub,zValues=self.zValuesFloat,xName="x",yName="y",zName="Time",title=smoothfilter)
        if self.SmoothedFitParameters.iloc[:,3] is not None:
            self.plotSimple(xValues=SmoothX_sub,yValues=SmoothY_sub,zValues=self.SmoothedFitParameters.iloc[:,3],xName="x",yName="y",zName="Intensity",title=smoothfilter)
        else:
            self.plotSimple(xValues=SmoothX_sub,yValues=SmoothY_sub,xName="x",yName="y",title=smoothfilter)
##################            
#   SetxValues   #
##################       




#data = analyse2D()
data = analyseFitParameters()
#filename="QD03_Polarization_Det_Exc90Deg_400nm_40MHz_006nW_5K.csv"
#filename="QD01_Spectrum_Exc400nm_40MHz_001nW_G1800_005K.csv"
#filename="TRPL_QD01_400nm_25GHz_5K_5nW.csv"
#filename="Dot04_TimeTrace_Exc400nm_40MHz_090nW_010ms_bin_5K.csv"
#filename='C:/Users/bec/Desktop/temp/test/FA/'
filename="PL_Power_Series_NPL_CsPbBrCL_400nm_18nW_5K_spot3_single.csv"


#data.importDataTIFF(path,filename, sortByDate=False, Offset=32768, ScalingX=1,ScalingY=1,RescaleFactor=5)
"""
path="Z:/groups/photonics/Projects/QuantumPhotonics/Perovskite QDs/MichaelBecker/20170803_singleNPL/QD01/"
data.importDataTIFF(path, sortByDate=True, Offset=32768, ScalingX=16e3/50.0,ScalingY=16e3/50.0,RescaleFactor=5)
data.setBoundariesDataList([8.8e4,9.15e4],[5.75e4,6.1e4],[10,200])
data.fitDefinedPeak2DList(function="circularGauss",Peak=[8.965e4,5.87e4],Intensity=40,FWHM=500,Plot=False,UseFitParameters=True,WeightedsquaredResidual=False)
data.smoothParameters2D("ga",5,0,1)

path="Z:/groups/photonics/Projects/QuantumPhotonics/Perovskite QDs/MichaelBecker/20170728_singleNPL/QD02/"
data.importDataTIFF(path, sortByDate=True, Offset=32768, ScalingX=16e3/50.0,ScalingY=16e3/50.0,RescaleFactor=5)
data.setBoundariesDataList([8.8e4,9.14e4],[5.65e4,6.0e4],[1,200])
data.fitDefinedPeak2DList(function="circularGauss",Peak=[8.965e4,5.80e4],Intensity=150,FWHM=500,Plot=False,UseFitParameters=True,WeightedsquaredResidual=False)


path="Z:/groups/photonics/Projects/QuantumPhotonics/Perovskite QDs/MichaelBecker/20170803_singleNPL/QD03_2/"
data.importDataTIFF(path, sortByDate=True, Offset=32768, ScalingX=16e3/50.0,ScalingY=16e3/50.0,RescaleFactor=5)
data.setBoundariesDataList([9.1e4,9.4e4],[5.8e4,6.1e4],[0,600])
data.fitDefinedPeak2DList(function="circularGauss",Peak=[9.2e4,5.9e4],Intensity=130,FWHM=500,Plot=False,UseFitParameters=True,WeightedsquaredResidual=False)
data.smoothParameters2D("ga",10,0,1)

path="Z:/groups/photonics/Projects/QuantumPhotonics/Perovskite QDs/MichaelBecker/20170803_singleNPL/QD03_2/"
data.importDataTIFF(path, sortByDate=True, Offset=32768, ScalingX=16e3/50.0,ScalingY=16e3/50.0,RescaleFactor=5)
data.setBoundariesDataList([8.95e4,9.23e4],[5.6e4,5.85e4],[0,600])
data.fitDefinedPeak2DList(function="circularGauss",Peak=[9.08e4,5.73e4],Intensity=130,FWHM=500,Plot=False,UseFitParameters=True,WeightedsquaredResidual=False)
data.smoothParameters2D("ga",10,0,1)
"""
#data.setBoundariesDataList([9.1e4,9.28e4],[5.8e4,6.0e4],[0,600])
#data.setBoundariesDataList([275*16e-6/50.0,285*16e-6/50.0],[177*16e-6/50.0,187*16e-6/50.0],[0,140])
#data.setBoundariesDataList([8.86e-5,9.25e-5],[5.65e-5,6.05e-5],[300,400])

#data.setBoundariesData([546,590],[485,530])
#data.setBoundariesDataList([546,590],[485,530],[50,250])
#Z:\groups\photonics\Projects\QuantumPhotonics\Perovskite QDs\MichaelBecker\20161005_Supercrystals_Amplifier\Powerdependence\1\Lifetime_Sample81_1
##################data.importDataCSV('Z:/groups/photonics/Projects/QuantumPhotonics/Perovskite QDs/MichaelBecker/20170927_singleCsPbBr2Cl_NPL/',filename=filename,sortByDate=False)
#data.set_yValues([5,15,30,50,75,105,130,160,200,250,278,300],"Temperature","K")
#data.setBoundariesData(xRange=[0,25],yRange=[-1,500])
#data.normalizeData()
#data.setZero()
#data.fitLifetimeBi(0.05,20,1,5)
#data.saveParametersTXT(filename_extension="BiExpFit_noOffset")
#data.saveDataTXT()

#data.importExportRandomCSV('Z:/groups/photonics/Projects/QuantumPhotonics/Perovskite QDs/MichaelBecker/20170927_singleCsPbBr2Cl_NPL/')
#data.importDataCSV('C:/Users/bec/Documents/_Documents/Analysis/20170314_SetupCalibration/20170314_SetupCalibration/',sortByDate=False)
#data.plotData2D()
#data.saveDataTXT()
#data.integrateCounts([2.5],[2.57])
# data.saveParametersTXT(data.intCounts)

#data.load_yValues('Z:/groups/photonics/Projects/QuantumPhotonics/Perovskite QDs/MichaelBecker/20161007_Supercrystals_Amplifier_lowT/#69/2/Summary/Crystal02_#71_Powersseries_Exc400nm_1kHz_Summary.csv',yName="Power",yUnit="W")
#data.set_yValues(np.linspace(0,3,62).tolist(),yName="Z-axis",yUnit="micro m")
#data.set_yValues(np.linspace(0,2.72,137).tolist(),yName="Z-axis",yUnit="micro m")
#data.set_yValues([5,10,20,30,40,60,80,100,120,140,160,180,240,280],yName="Temperature",yUnit="K")
#data.findPeakPosition(100,0.2,30,0,2.48)
#data.fitAllPeaks("Lorentz",0.01,100,5,0,2.48)
#data.fitDefinedPeaks("Lorentz",[1.76,1.88,1.97],[0.1,0.1,0.1],[100,100,100],[0,0,0],[0,1.81,1.96],[1.81,1.96,2.5])

#data.saveDataTXT()
#data.plotData2D()
#data.fitDefinedPeaks("Lorentz",[2.506],[0.1],[100],[0],[2.5],[2.51])
#data.saveParametersTXT(data.intCounts)
#data.saveParametersTXT()
    
def BoltzmannTemp(Delta,IR):
    return (-1)*1/8.617330E-5*Delta/(np.log(IR)), Delta, IR

def FWMSignal(t12,t13,Delta,T2,omega):
    return np.exp(-4*Delta**2*(1-np.cos(omega*t12))*(1-np.cos(omega*t13)))*np.exp(-4*t12/T2)
    
    
def opticalDensity(a,b,P,R):
    return P/(np.pi*a*b*R*3.06*1.6022e-19)