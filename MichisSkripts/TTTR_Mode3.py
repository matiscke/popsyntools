# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 09:55:28 2016

@author: bec
"""

import os, numpy as np, csv, matplotlib.pyplot as plt, scipy.optimize as opt, math, struct, binascii, gc, time, random
import multiprocessing
from operator import sub
from joblib import Parallel, delayed
import MikeMath as MM
import pandas as pd

#I get errors unless the function that the cores execute in parallel is defined outside the class function 
def hist2(x,y,bins):
    store = np.zeros(len(bins)-1,dtype='float');
    for i in x:
        res = y.searchsorted(bins+i)
        store += res[1:]-res[:-1]
    return store


class datalog:
    def __init__(self):
        self.data=None
        self.length=None
        
    def appendData(self,a):
        self.data.append(a)
    def getlength(self):
        return len(self.data)
    
class LogLifetime:
    
    def __init__(self):
        self.cycle=65536
        self.ts0 = None
        self.nsync = None
        self.fastCorr = None
        self.syncdiv = None
        self.dt = None
        self.resolution = None
        self.bins = None
        self.hist = None
        self.filename = None
        self.crossCorr = None
        self.autoCorr0 = None
        self.autoCorr1 = None
        
        self.binlen = None
        
        self.norm = None
        
        return
    
    def importDataArrays(self,ts0,ts1,dt):
        
        return
    
    def importDataFile(self,filename,nrrec=None):
        self.filename=filename
            
        with open(filename, "rb") as f:

            while True:
                if f.read(1) == b'H':
                    if f.read(13) == b'WSync_Divider': # recognize time unit entry
                        break
            f.read(26) # rest of the tag
            self.syncdiv = struct.unpack('i',f.read(4))[0]
            print('Sync Divider:', self.syncdiv)            
            while True:
                if f.read(1) == b'M':
                    if f.read(18) == b'easDesc_Resolution':
                        break
            f.read(21) # rest of the tag
            self.resolution = struct.unpack('d',f.read(8))[0]
            print('Resolution (s):', self.resolution)
            while True:
                if f.read(1) == b'M':
                    if f.read(24) == b'easDesc_GlobalResolution': # recognize time unit entry
                        break
            f.read(15) # rest of the tag
            self.dt = struct.unpack('d',f.read(8))[0]
            print('Time unit:', self.dt)

            while True:
                if f.read(1) == b'T':
                    if f.read(23) == b'TResult_NumberOfRecords': # recognize number of records entry
                        break
            f.read(16) # rest of the tag
            if nrrec is None:
                nrrec = struct.unpack('q',f.read(8))[0] # extract number of records
                print('Number of records in file:', nrrec)

            while True:
                if f.read(1) == b'H':
                    if f.read(9) == b'eader_End':
                        print('Header_End found')
                        break
            f.read(38) # rest of Header_End
            
            counter0 = 0; sync = 0; counter15 = 0
            start0 = 0; start1 = 0
            macrotime = 0
            
            self.ts0 = np.zeros(nrrec,dtype='int64'); #self.ts0.fill(0)
            self.nsync = np.zeros(nrrec,dtype='int64');
            #self.fastCorr = np.zeros(nrrec,dtype='int32');
            
            for i in range(nrrec):
                entry = f.read(4)
                channel = struct.unpack("I",entry)[0] >> 28 # read channel number, first 4 bits
                
                #totaltime = (struct.unpack("I",entry)[0] & 0xFFFFFFF) + macrotime  
                
                
                
                
                if channel == 0 or channel ==1:
                    # read true arrival time
                    totaltime = (struct.unpack("I",entry)[0] & 0xFFFF) +macrotime
                    self.ts0[counter0] = totaltime
                    
                    
                    
                    sync = (struct.unpack("I",entry)[0] & 0xFFFFFFF)
                     #store absolute time in array ts0
                    self.nsync[counter0]= sync
                    
                    #sync = ((struct.unpack("I",entry)[0]>>16) & 4095)
                    
                     #store absolute time in array ts0
                    self.nsync[counter0]= sync
                    counter0 += 1
                    
                elif channel == 15:
                    #macrotime += 210698240 #freddy
                    macrotime +=self.cycle #40MHz
                    #macrotime = (struct.unpack("I",entry)[0] & 0xFFFFFFF)
                    counter15 += 1
            print counter0,counter15
            print "Time: ",totaltime,sync 
            self.ts0 = self.ts0[0:counter0]
            self.nsync = self.nsync[0:counter0]            
            """
            self.ts0 = self.ts0[0:counter0]
            self.ts1 = self.ts1[0:counter1]
            self.fastCorr = self.fastCorr[0:counter0+counter1]
            
            print('File import time:', (time.time() - calctime0), 's')
            """
            print('Number of events on channel 0:', counter0)
            print('Number of overflows:', counter15)
            print('Total experiment time:', self.ts0[-1]*self.dt, 's')

            #print('ts0:', self.ts0.nbytes/1000000.0, 'MB')
            #print('ts1:', self.ts1.nbytes/1000000.0, 'MB')
    def setTimeBoundaries(self,tmin=0,tmax=1000):
        tminIndex=min(enumerate(self.ts0*self.dt), key=lambda x: abs(x[1]-tmin))[0]
        tmaxIndex=min(enumerate(self.ts0*self.dt), key=lambda x: abs(x[1]-tmax))[0]
        
        self.ts0=self.ts0[tminIndex:tmaxIndex]
        self.nsync=self.nsync[tminIndex:tmaxIndex]
            
    def plotLifetime(self,bins=1000):
        plt.figure(figsize=[16,4])
        legendentries = []
        if self.ts0 is not None:
            trace = np.histogram(self.nsync,bins)

            #plt.plot(trace[1][:-1]*self.dt,trace[0]/trace[0].max(), color='b')
            plt.plot(trace[1][1:-1]/65536*self.resolution*1e9,trace[0][1:])
            
            legendentries.append('channel 0')

        plt.legend(legendentries)
        plt.xlabel('Time (ns)')
        plt.ylabel('Intensity (s)')
        #plt.set_yscale("log")
        #plt.axis([0,max(self.ts0.max(),self.ts1.max())*self.dt,0,2.1])
        ax = plt.gca()
        ax.get_yaxis().set_tick_params(direction='in')
        ax.get_xaxis().set_tick_params(direction='in')
        return trace   
    
    def plotLifetimeRange(self,bins=1000, binsLifetime=1000, Range=[[0,np.inf]]):
        fig, Lifetime=plt.subplots(1,1,figsize=[16,4])
        legendentries = []
        self.LifetimeTraces=[]
        for t in range(len(Range)):
            if self.ts0 is not None:
                timetrace = np.histogram(self.ts0,bins)
                
                self.Rangesync=[]
                runner=0
                for i in range(bins):
                    if  timetrace[0][i] >=Range[t][0] and timetrace[0][i] <=Range[t][1]:
                       for k in range(timetrace[0][i]):
                           self.Rangesync.append(self.nsync[runner+k])
                           
                    runner+=timetrace[0][i]
            if self.Rangesync is not None:
                
                LifetimeTrace = np.histogram(self.Rangesync,binsLifetime)
    
                #plt.plot(trace[1][:-1]*self.dt,trace[0]/trace[0].max(), color='b')
                Lifetime.plot(LifetimeTrace[1][1:-1]/self.cycle*self.resolution*1e9,np.divide(LifetimeTrace[0][1:],float(LifetimeTrace[0][1:].max())))
                legendentries.append('Range '+str(Range[t][0])+" - "+str(Range[t][1]))
            self.LifetimeTraces.append(LifetimeTrace)
    
        Lifetime.legend(legendentries)
        Lifetime.set_yscale("log")
        Lifetime.set_xlabel('Time (s)')
        Lifetime.set_ylabel('Intensity (s)')
        #plt.set_yscale("log")
        #plt.axis([0,max(self.ts0.max(),self.ts1.max())*self.dt,0,2.1])
        ax = plt.gca()
        ax.get_yaxis().set_tick_params(direction='in')
        ax.get_xaxis().set_tick_params(direction='in')        
    
    def saveLifetimeRange(self,bins=1000, binsLifetime=1000, Range=[[0,np.inf]],filename_extension="LifetimeRange"):
        self.LifetimeRange=[]
        for t in range(len(Range)):
            if self.ts0 is not None:
                timetrace = np.histogram(self.ts0,bins)
                
                self.Rangesync=[]
                runner=0
                for i in range(bins):
                    if  timetrace[0][i] >=Range[t][0] and timetrace[0][i] <=Range[t][1]:
                       for k in range(timetrace[0][i]):
                           self.Rangesync.append(self.nsync[runner+k])
                           
                    runner+=timetrace[0][i]
            if self.Rangesync is not None:
                
                LifetimeTrace = np.histogram(self.Rangesync,binsLifetime)
            self.LifetimeRange.append(LifetimeTrace)
    
        if self.LifetimeRange:      
            Header_line1 = "Time \t Intensity "
            Header_line2 = "ns  \t arb. units"
            Header_line3="Bins Trace =  "+str(bins)+", Bins Lifetime: "+str(binsLifetime)
            for y in range(len(Range)):
                Header_line3=Header_line3+"\t Range: "+str(Range[y][0])+" - "+str(Range[y][1])+""
            extension="_"+filename_extension+".txt"
            data=[]
            [data.append(self.LifetimeRange[z][0]) for z in range(len(self.LifetimeRange))]
            saveFile=np.vstack((np.multiply(np.array(self.LifetimeTraces[0][1][1:]),self.resolution*1e9/self.cycle),np.array(data))).T
            np.savetxt(self.filename.replace(".ptu",extension),saveFile,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
                        
                                
            
    def saveLifetime(self,bins=1000,filename_extension="Lifetime"):
        if self.nsync is not None:
            trace = np.histogram(self.nsync,bins)

            traceX=np.multiply(trace[1][1:],self.resolution/65536.0*1e9)
        
        Header_line1 = "Time \t Intensity \t Norm. Intensity"
        Header_line2 = "ns  \t arb. units \t arb. units"
        Header_line3=" "
        extension="_"+filename_extension+".txt"

        np.savetxt(self.filename.replace(".ptu",extension),np.vstack((np.array(traceX),np.array(trace[0]),np.divide(np.array(trace[0]),float(trace[0].max())))).T,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
        
    
        
    def plotTimeTrace(self,bins=1000):
        plt.figure(figsize=[16,4])
        legendentries = []
        if self.ts0 is not None:
            trace = np.histogram(self.ts0,bins)
            plt.plot(trace[1][:-1]*self.dt,trace[0])
            legendentries.append('channel 0')
        """
        if self.ts1 is not None:
            trace = np.histogram(self.ts1,bins)
            plt.plot(trace[1][:-1]*self.dt,trace[0]/trace[0].max()+1, color='r')
            legendentries.append('channel 1')
        """
        plt.legend(legendentries)
        plt.xlabel('Time (s)')
        plt.ylabel('Intensity (s)')
        #plt.axis([0,max(self.ts0.max(),self.ts1.max())*self.dt,0,2.1])
        ax = plt.gca()
        ax.get_yaxis().set_tick_params(direction='out')
        ax.get_xaxis().set_tick_params(direction='out')
        return trace
    
    def saveTimeTrace(self,bins=1000,filename_extension="TimeTrace"):
        if self.ts0 is not None:
            trace = np.histogram(self.ts0,bins)

            traceX=np.multiply(trace[1][1:],self.dt)
        
        Header_line1 = "Time \t Intensity \t Norm. Intensity"
        Header_line2 = "s  \t arb. units \t arb. units"
        Header_line3="Bins = "+str(bins)
        extension="_"+filename_extension+".txt"

        np.savetxt(self.filename.replace(".ptu",extension),np.vstack((np.array(traceX),np.array(trace[0]),np.divide(np.array(trace[0]),float(trace[0].max())))).T,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
        
    
    
    def plotFFT(self, bins=10000, degree=True):
        fig, Fourier=plt.subplots(2,2,figsize=(20,12.5), sharex=False)
        #legendentries = []
        fourier_x= []
        fourier_y= []
        fourier_r= []
        fourier_phi= []
        
        if self.ts0 is not None:
            trace = np.histogram(self.ts0,bins)
            traceX=np.multiply(trace[1],self.dt)
            fourier_y=np.fft.rfft(trace[0])
            fourier_x=[x/traceX.max() for x in range(bins/2+1)]
            fourier_r= np.abs(fourier_y)
            fourier_phi= np.angle(fourier_y, deg=degree)
            Fourier[0][0].plot(fourier_x,np.real(fourier_y))
            Fourier[0][0].legend("Real part")
            Fourier[0][0].set_xlabel('Frequency (Hz)')
            Fourier[0][0].set_ylabel('Real part Intensity')
            Fourier[1][0].plot(fourier_x,np.imag(fourier_y))  
            Fourier[1][0].legend("Imaginary part")
            Fourier[1][0].set_xlabel('Frequency (Hz)')
            Fourier[1][0].set_ylabel('Real part Intensity')
            Fourier[0][1].plot(fourier_x,fourier_r)  
            Fourier[0][1].legend("Amplitude")
            Fourier[0][1].set_xlabel('Frequency (Hz)')
            Fourier[0][1].set_xscale("log", nonposx='clip')
            Fourier[0][1].set_yscale("log", nonposx='clip')
            Fourier[0][1].set_ylabel('Amplitude')
            Fourier[1][1].plot(fourier_x,fourier_phi)  
            Fourier[1][1].legend("Phase")
            Fourier[1][1].set_xlabel('Frequency (Hz)')
            Fourier[1][1].set_ylabel('Phase')
             
        return fourier_x, fourier_y,fourier_r,fourier_phi
    
    def saveFFT(self,bins=1000, filename_extension="FFT",degree=True ):
        if self.ts0 is not None:
            trace = np.histogram(self.ts0,bins)
            traceX=np.multiply(trace[1],self.dt)
            fourier_x= []
            fourier_y= []
            fourier_r= []
            fourier_phi= []
            fourier_y=np.fft.rfft(trace[0])
            fourier_x=[x/traceX.max() for x in range(bins/2+1)]
            fourier_r= np.abs(fourier_y)
            fourier_phi= np.angle(fourier_y, deg=degree)
        
        Header_line1 = "Frequeny \t Real Part FFT \t Imaginary Part FFT \t Magnitude \t Phase "
        Header_line2 = "Hz  \t arb. units \t arb. units \t arb. units"
        if degree:
            Header_line2=Header_line2+"\t Degree"
        else:
            Header_line2=Header_line2+"\t Radians"
        Header_line3="Bins = "+str(bins)
        extension="_"+filename_extension+".txt"

        np.savetxt(self.filename.replace(".ptu",extension),np.vstack((np.array(fourier_x),np.array(np.real(fourier_y)),np.array(np.imag(fourier_y)),np.array(fourier_r),np.array(fourier_phi))).T,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
    


    def plotLifetimeTrace(self,bins=1000, binsLifetime=1000,binsHistogram=100,FitRange_low=0, FitRange_high=np.inf, decaytime=0.2, offset=None):
              
        self.LifetimeTracesDecay=[]
        self.LifetimeTracesOffset=[]
        
        
        if self.ts0 is not None:
            timetrace = np.histogram(self.ts0,bins)
            fig, LifetimeT=plt.subplots(2,1,figsize=[16,8])
            
            runner=0
            for i in range(bins):
                BinSync=[]
                for k in range(timetrace[0][i]):
                       BinSync.append(self.nsync[runner+k])
                
                LifetimeTrace = np.histogram(BinSync,binsLifetime)                
                
                
                t0=min(enumerate(LifetimeTrace[1]/self.cycle*self.resolution*1e9), key=lambda x: abs(x[1]-FitRange_low))
                t1=min(enumerate(LifetimeTrace[1]/self.cycle*self.resolution*1e9), key=lambda x: abs(x[1]-FitRange_high))
                # p[2]-- height
                # p[1]-- centre position
                # p[0]-- decay time
                # p[3]-- Offset  
                
                
                
                X=LifetimeTrace[1][t0[0]:t1[0]]/self.cycle*self.resolution*1e9-FitRange_low
    
                Y=np.divide(LifetimeTrace[0][t0[0]:t1[0]],float(LifetimeTrace[0][1:].max()))
                
                #parameters=[decaytime,FitRange_low,1,offset]
                #fit, best_parameters, err=MM.expfit(Y,X,parameters)
                if offset is None:
                    parameters=[decaytime,1]
                    fit, best_parameters, err=MM.MonoExpFit_noOffset(Y,X,parameters)
                    self.LifetimeTracesDecay.append(best_parameters[0])
                    self.LifetimeTracesOffset.append(0)
                else:
                   parameters=[decaytime,1,offset]
                   fit, best_parameters, err=MM.MonoExpFit_wOffset(Y,X,parameters) 
                   self.LifetimeTracesDecay.append(best_parameters[0])
                   self.LifetimeTracesOffset.append(best_parameters[2])
                #parameters2=[decaytime,FitRange_low,1]
                #fit2, best_parameters2, err=MM.expfitfixed(Y,X,parameters2)
                runner+=timetrace[0][i]
                """
                if i<100:
                    print best_parameters
                    fig, test=plt.subplots(1,1,figsize=[8,8])
                    test.plot(X,Y)
                    test.plot(X,fit,label=str(best_parameters[0])+", "+str(best_parameters[1]))
                    plt.legend()
                
                
                if best_parameters[0]>0.5:
                    fig, test=plt.subplots(1,1,figsize=[8,8])
                    test.plot(X,Y)
                    test.plot(X,fit,label=str(best_parameters[0])+", "+str(best_parameters[1])+", "+str(best_parameters[2]))
                    plt.legend()
                """
                   
            trace = np.histogram(self.ts0,bins)
            LifetimeT[0].plot(trace[1][:-1]*self.dt,self.LifetimeTracesDecay)
            LifetimeT[1].plot(trace[1][:-1]*self.dt,trace[0])
            
            #LifetimeT[0].legend("Lifetime Trace-single exp fit")
            LifetimeT[0].set_xlabel('Time (s)')
            LifetimeT[0].set_ylabel('Decay Time (ns)')
            #LifetimeT[0].legend() 
            #LifetimeT[1].legend("Lifetime Trace-single exp fit")
            LifetimeT[1].set_xlabel('Time (s)')
            LifetimeT[1].set_ylabel('Intensity per bin')
            #LifetimeT[1].legend() 
            ax = plt.gca()
            ax.get_yaxis().set_tick_params(direction='in')
            ax.get_xaxis().set_tick_params(direction='in')
            # 2D Histogram
            self.IntensityLifetimeCorrelation=[]
            self.IntensityLifetimeCorrelation=np.histogram2d(self.LifetimeTracesDecay,trace[0],bins)
            fig, IntensityLifetimeHist=plt.subplots(1,1,figsize=[8,8])            
            IntensityLifetimeHist.hist2d(self.LifetimeTracesDecay,trace[0],bins=binsHistogram)
            IntensityLifetimeHist.legend("Intensity Lifetime Correlation")
            IntensityLifetimeHist.set_xlabel('Decay Time (ns)')
            IntensityLifetimeHist.set_ylabel('Intensity per bin (ns)')
    
    def plotIntHist(self,binsTrace=1000,binsHist=100):
        plt.figure(figsize=[8,8])
        if self.ts0 is not None:
            trace = np.histogram(self.ts0,binsTrace)
            Hist =    np.histogram(trace[0],binsHist)        
            plt.plot(Hist[1][:-1],Hist[0])
        plt.xlabel('Intensity per bin (arb. units)')
        plt.ylabel('Occurance')
        #plt.axis([0,max(self.ts0.max(),self.ts1.max())*self.dt,0,2.1])
        ax = plt.gca()
        ax.get_yaxis().set_tick_params(direction='out')
        ax.get_xaxis().set_tick_params(direction='out')
        return Hist
    def saveIntHist(self,binsTrace=1000,binsHist=100, filename_extension="IntHist"):
        
        if self.ts0 is not None:
            trace = np.histogram(self.ts0,binsTrace)
            Hist =    np.histogram(trace[0],binsHist)        
            
        Header_line1 = "Time \t Intensity \t Norm. Intensity"
        Header_line2 = "s  \t arb. units \t arb. units"
        Header_line3="Bins Trace= "+str(binsTrace)+", Bins Hist= "+str(binsHist)
        extension="_"+filename_extension+".txt"

        np.savetxt(self.filename.replace(".ptu",extension),np.vstack((np.array(Hist[1][:-1]),np.array(Hist[0]))).T,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
    def plotPowerLawBlinking(self,binsTrace=1000,binsOnOff=10,OnOffBoundaries=[]):
        self.plotTimeTrace(binsTrace)
        if self.ts0 is not None:
            trace = np.histogram(self.ts0,binsTrace)
            I0=trace[0][0]
            t0=trace[1][0]
            DurationTime=[[] for i in range(len(OnOffBoundaries))]
            for j in range(len(OnOffBoundaries)):
                if OnOffBoundaries[j][0]<I0<=OnOffBoundaries[j][1]:
                    Range0=j
            for i in range(binsTrace):
                I1=trace[0][i]
                t1=trace[1][i]
                for j in range(len(OnOffBoundaries)):
                    if OnOffBoundaries[j][0]<I1<=OnOffBoundaries[j][1]:
                        Range1=j
                        
                if Range0 != Range1:
                    DurationTime[Range0].append((t1-t0)*self.dt)
                    Range0=Range1
                    I0=I1
                    t0=t1
            
            
            minList=[]
            maxList=[]
            for t in range(len(DurationTime)):
                minList.append(min(DurationTime[t]))
                maxList.append(max(DurationTime[t]))
            tmin=np.log10(min(minList))                
            tmax=np.log10(max(maxList)) 
            
            if tmin<0:
                tmin=-1*math.ceil(abs(tmin))
            else:
                tmin=math.ceil(tmin)
            if tmax<0:
                tmax=math.trunc(tmax)
            else:
                tmax=math.ceil(tmax)
                
            LogHistogram=[]    
            
            
            ##LogSpace##
            LogSpace=np.logspace(tmin,tmax,binsOnOff)
            r=1
            for t in range(len(LogSpace)):
                if t>0:
                    if LogSpace[t]<=max(self.ts0*self.dt)/binsTrace*t:
                        LogSpace[t]=r*max(self.ts0*self.dt)/binsTrace
                        r+=1
                    
            
            
            for k in range(len(DurationTime)):
                LogHistogram.append(np.histogram(DurationTime[k],LogSpace))
            
            
            cmap = plt.cm.jet
            plt.figure(figsize=[8,8])
            ax = plt.gca()
            for m in range(len(LogHistogram)):
                ax.scatter(LogHistogram[m][1][:-1],LogHistogram[m][0],marker="o",s=60,alpha=0.4,color=cmap(m / float(len(LogHistogram))))
            ax.axis([10**tmin,10**tmax,0.9,max([max(LogHistogram[n][0]) for n in range(len(LogHistogram))])])
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlabel("Duration (s)")
            ax.set_ylabel("Occurrence")
            return LogHistogram
    def savePowerLawBlinking(self,binsTrace=1000,binsOnOff=10,OnOffBoundaries=[],filename_extension="PowerLawBlinking"):
        
        if self.ts0 is not None:
            trace = np.histogram(self.ts0,binsTrace)
            I0=trace[0][0]
            t0=trace[1][0]
            DurationTime=[[] for i in range(len(OnOffBoundaries))]
            for j in range(len(OnOffBoundaries)):
                if OnOffBoundaries[j][0]<I0<=OnOffBoundaries[j][1]:
                    Range0=j
            for i in range(binsTrace):
                I1=trace[0][i]
                t1=trace[1][i]
                for j in range(len(OnOffBoundaries)):
                    if OnOffBoundaries[j][0]<I1<=OnOffBoundaries[j][1]:
                        Range1=j
                        
                if Range0 != Range1:
                    DurationTime[Range0].append((t1-t0)*self.dt)
                    Range0=Range1
                    I0=I1
                    t0=t1
            
            
            minList=[]
            maxList=[]
            for t in range(len(DurationTime)):
                minList.append(min(DurationTime[t]))
                maxList.append(max(DurationTime[t]))
            tmin=np.log10(min(minList))                
            tmax=np.log10(max(maxList)) 
            
            if tmin<0:
                tmin=-1*math.ceil(abs(tmin))
            else:
                tmin=math.ceil(tmin)
            if tmax<0:
                tmax=math.trunc(tmax)
            else:
                tmax=math.ceil(tmax)
                
            LogHistogram=[]   
            
            LogSpace=np.logspace(tmin,tmax,binsOnOff)
            r=1
            for t in range(len(LogSpace)):
                if t>0:
                    if LogSpace[t]<=max(self.ts0*self.dt)/binsTrace*t:
                        LogSpace[t]=r*max(self.ts0*self.dt)/binsTrace
                        r+=1
            for k in range(len(DurationTime)):
                LogHistogram.append(np.histogram(DurationTime[k],LogSpace))                 
            Header_line1 = "Time \t Occurrence "
        Header_line2 = "s  \t arb. units "
        Header_line3="Bins Trace= "+str(binsTrace)+", Bins OnOff= "+str(binsOnOff)
        for t in range(len(LogHistogram)):
            Header_line3=Header_line3+str(OnOffBoundaries[t])
        extension="_"+filename_extension+".txt"
        data=[]
        [data.append(LogHistogram[z][0]) for z in range(len(LogHistogram))]
        saveFile=np.vstack((LogHistogram[0][1][:-1],np.array(data))).T
        np.savetxt(self.filename.replace(".ptu",extension),saveFile,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
        
    def saveLifetimeTrace(self,bins=1000, binsLifetime=1000,binsHistogram=100,FitRange_low=0, FitRange_high=np.inf, decaytime=0.2, offset=1e-3,filename_extension="LifetimeTrace"):
              
        self.LifetimeTracesDecay=[]
        self.LifetimeTracesZero=[]
        self.LifetimeTracesHeight=[]
        self.LifetimeTracesOffset=[]
        self.LifetimeTracesDecayErr=[]
        self.LifetimeTracesZeroErr=[]
        self.LifetimeTracesHeightErr=[]
        self.LifetimeTracesOffsetErr=[]
        
        
        if self.ts0 is not None:
            timetrace = np.histogram(self.ts0,bins)
            
            
            runner=0
            for i in range(bins):
                BinSync=[]
                for k in range(timetrace[0][i]):
                       BinSync.append(self.nsync[runner+k])
                
                LifetimeTrace = np.histogram(BinSync,binsLifetime)                
                
                
                t0=min(enumerate(LifetimeTrace[1]/self.cycle*self.resolution*1e9), key=lambda x: abs(x[1]-FitRange_low))
                t1=min(enumerate(LifetimeTrace[1]/self.cycle*self.resolution*1e9), key=lambda x: abs(x[1]-FitRange_high))
                # p[2]-- height
                # p[1]-- centre position
                # p[0]-- decay time
                # p[3]-- Offset  
                
                
                parameters=[decaytime,FitRange_low,1,offset]
                X=LifetimeTrace[1][t0[0]:t1[0]]/self.cycle*self.resolution*1e9
                
                Y=np.divide(LifetimeTrace[0][t0[0]:t1[0]],float(LifetimeTrace[0][1:].max()))
                
                fit, best_parameters, err=MM.expfit(Y,X,parameters)
                
                self.LifetimeTracesDecay.append(best_parameters[0])
                self.LifetimeTracesZero.append(best_parameters[1])
                self.LifetimeTracesHeight.append(best_parameters[2])
                self.LifetimeTracesOffset.append(best_parameters[3])
                self.LifetimeTracesDecayErr.append(err[0])
                self.LifetimeTracesZeroErr.append(err[1])
                self.LifetimeTracesHeightErr.append(err[2])
                self.LifetimeTracesOffsetErr.append(err[3])
                
                runner+=timetrace[0][i]
            trace = np.histogram(self.ts0,bins)
            
            Header_line1 = "Time \t Lifetime \t Lifetime Error \t Zero\t Zero Error\t Height\t Height Err \t Offset \t Offset Error\t "
            Header_line2 = "s  \t ns \t ns \t s \t s \t arb. units \t arb. units \t arb. units \t arb. units\t"
            Header_line3="Bins = "+str(bins)+"Lifetime Bins = "+str(binsLifetime)
            extension="_"+filename_extension+".txt"
            np.savetxt(self.filename.replace(".ptu",extension),np.vstack((np.array(trace[1][:-1]*self.dt),np.array(self.LifetimeTracesDecay),np.array(self.LifetimeTracesDecayErr),np.array(self.LifetimeTracesZero),np.array(self.LifetimeTracesZeroErr),np.array(self.LifetimeTracesHeight),np.array(self.LifetimeTracesHeightErr),np.array(self.LifetimeTracesOffset),np.array(self.LifetimeTracesOffsetErr))).T,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
            # 2D Histogram
            self.IntensityLifetimeCorrelation=[]
            self.IntensityLifetimeCorrelation=np.histogram2d(self.LifetimeTracesDecay,trace[0],binsHistogram)
            
            Header_line1 = "Decay Time \t Correlation "
            Header_line2 = "ns  \t arb. units "
            Header_line3="Bins = "+str(bins)+"Lifetime Bins = "+str(binsLifetime)
            for i in range(len(self.IntensityLifetimeCorrelation[1])-1):
                Header_line3=Header_line3+"\t"+str(self.IntensityLifetimeCorrelation[2][i])
            extension2="_"+filename_extension+"_LifetimeIntensityCorr.txt"
            np.savetxt(self.filename.replace(".ptu",extension2),np.vstack((np.array(self.IntensityLifetimeCorrelation[1][:-1]),np.array(self.IntensityLifetimeCorrelation[0]).T)).T,fmt="%.6e\t",header='\n'.join([Header_line1, Header_line2,Header_line3]),comments='')
            
            self.saveTimeTrace(bins,filename_extension+"TimeTrace")
                       

    """  
    def calcAutoCorr0(self,bins):
        self.autoCorr0 = self.correlate(self.ts0,self.ts0,bins)
        
    def correlate(self,a,b,bins):
        calctime0 = time.time()

        counterA = 0
        counterB = 0

        hist = np.empty(len(bins)-1,dtype='float'); hist.fill(0)
        
        while counterA < len(a):
            res = b.searchsorted(bins+a[counterA])
            hist += res[1:]-res[:-1]
            counterA += 1
            #lastmin = res[0]
        
        print('Total analysis time:', (time.time() - calctime0), 's')
        
        tmax = max(a[-1],b[-1])
        norm = 0.5*(bins[:-1] - bins[1:])*(bins[:-1] + bins[1:] - 2*tmax - 1)
        Iavgsq = len(a)*len(b)/tmax
        
        self.norm = norm/norm[0]*Iavgsq
        
        return np.divide(hist,norm/norm[0]*Iavgsq) # normalize
        #loghistnormunfiltered = loghistnorm
        
    """   
    def blinkingStats(self,bins=1000, threshold_low=1, threshold_high=1):
        state = 0
        t = 0
        OFFtimes = []
        ONtimes = []
        trace=[]
        if self.ts0 is not None:
            trace = np.histogram(self.ts0,bins)
        for i in range(len(trace[0])): # run over all bins
            if trace[0][i] < threshold_low: # if new bin is OFF sate
                if state == 0: # if QD was OFF already
                    t += 1
                elif state == 2: # if QD was ON
                    ONtimes.append(t)
                    t = 1; state = 0
                else: # if QD was between ON and OFF
                    t = 1; state = 0
            elif trace[0][i] < threshold_high: # if new bin is between ON and OFF
                if state == 0: # if QD was OFF
                    OFFtimes.append(t)
                    t = 1; state = 1
                elif state == 2: # if QD was ON
                    ONtimes.append(t)
                    t = 1; state = 1
                else: # if QD was between ON and OFF already
                    t += 1
            else: # if new bin is ON
                if state == 2: # if QD was ON already
                    t += 1
                elif state == 0: # if QD was OFF
                    OFFtimes.append(t)
                    t = 1; state = 2
                else: # if QD was between ON and OFF
                    t = 1; state = 1
               
        OFFhist = np.histogram(OFFtimes,int(1e2/(self.ts0[-1]*self.dt/bins))-1,[1,int(1e2/(self.ts0[-1]*self.dt/bins))])
        ONhist = np.histogram(ONtimes,int(1e2/(self.ts0[-1]*self.dt/bins))-1,[1,int(1e2/(self.ts0[-1]*self.dt/bins))])
        print ONtimes,OFFtimes
        # Below, the correction to estimate the probability of rare events from the finite measurement
        # See Fig 4 from Kuno et al., J. Chem. Phys. 115, 1028-1040 (2001)
        OFFhist2 = np.zeros(len(OFFhist[0])); ONhist2 = np.zeros(len(ONhist[0]))
        for i in range(max(OFFtimes)):
            lower = -1; higher = -1;
            for j in range(i-1,-1,-1):
                if OFFhist[0][j] > 0:
                    lower = j
                    break
            for j in range(i+1,max(OFFtimes)+1,1):
                if OFFhist[0][j] > 0:
                    higher = j
                    break
            OFFhist2[i] = (2*OFFhist[0][i]/(higher-lower))
            
        for i in range(max(ONtimes)):
            lower = -1; higher = -1;
            for j in range(i-1,-1,-1):
                if ONhist[0][j] > 0:
                    lower = j
                    break
            for j in range(i+1,max(ONtimes)+1,1):
                if ONhist[0][j] > 0:
                    higher = j
                    break
            ONhist2[i] = (2*ONhist[0][i]/(higher-lower))
        
        self.OFFstats = [OFFhist[1][:-1],OFFhist2]
        self.ONstats = [ONhist[1][:-1],ONhist2]





"""
filename = 'C:/Users/bec/Desktop/temp/test2.ptu'

corr1 = LogCorrelation()
corr1.importDataFile(filename)
trace=corr1.plotTimeTraces(10000)
Lifetime=corr1.plotLifetime(10000)

filename2 = 'C:/Users/bec/Desktop/temp/test.ptu'

corr2 = LogCorrelation()
corr2.importDataFile(filename2)
trace2=corr2.plotTimeTraces(10000)
Lifetime2=corr2.plotLifetime(10000)

filename3 = 'C:/Users/bec/Desktop/temp/test3.ptu'
"""
corr = LogLifetime()
#filename='C:/Users/bec/Desktop/temp/test3.ptu'
filename='Z:/groups/photonics/Projects/QuantumPhotonics/Perovskite QDs/MichaelBecker/20161221_CsPbI3/QD03_057nW.ptu'

corr.importDataFile(filename)
corr.plotLifetime()
"""
corr.saveLifetime(1000,"Lifetime_1000bins")

corr.plotLifetimeTrace(5000,500,50,2.47,20,0.2,1e-3)
corr.saveLifetimeTrace(5000,500,50,2.47,20,0.2,1e-3,"LifetimeTrace_200msbin")
#corr.saveLifetimeTrace(100000,1000,50,2.47,20,0.2,1e-3,"LifetimeTrace_10msbin")
corr.saveTimeTrace(100000,"TimeTrace_10msbin")
corr.plotFFT(10000)
corr.saveFFT(10000,"FFT_100msbin")


trace3=corr3.plotTimeTraces(10000)

Z:/groups/photonics/Projects/QuantumPhotonics/CdSCdSe_QD/BEC/20160817_Ilaria_CddSeCdS/QD03_TTTR3_97nW.ptu
Lifetime3=corr3.plotLifetime(1000)
fx,fy,fr,fphi=corr3.plotFFT(1000)
"""