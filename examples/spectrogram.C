R__LOAD_LIBRARY(libRootFftwWrapper.so)
#include "FFTtools.h"
#include <iostream>
#include <vector>
#include <string.h>
#include <typeinfo>
/*
spectrograme.C V1.1
This package makes a canavs with 24 spectragrams on them with incressing channel number reading from left to right 
in increassing channel number

how to use this program:
call spectGMain(const char* inputFolder,const char* inputTriggerPick,int station,int run)

the agrument inputFolder is a input string of where the folder that has the data is.  
in the standered data structure rnog data should be stored ex
/data/desy-mirror/inbox/station23/run23/combined.root
The input for the functtion here than should be "/data/desy-mirror/inbox"

the next argument is a string for the trigger filter eneter in the aberivation for the filter expaned on comment for the triggerReader function 

the next agment is the station number and the argument after that is the run number so 
binTime is the last arg in min
running spectGMain("/data/desy-mirror/inbox","",13,142,5) would apply no trigger filter on the 13th station on run number 142 with 5min bins
*/
double timeRange(int,int,const char*);//finds max time for between seperation in a event  
bool triggerReader(const char* , mattak::Dataset*);//triggerReader takes a pointer of the data set and gives back a bool on weather the selected trigger is triggered for no selection pass "", other input options "rf","force","pps","ext","radiant","lt"
void  spectrographGen(int,int,int,TH2*,TH1*,const char*,const char*);//makes individual spectragram 
TH2* spectGSelect(int,int,int ,double,const char*,const char*,double);//passes histograms to spectrographGen
void spectGMain(const char* inputFolder,const char* inputTriggerPick,int station,int run,double binTime){
	TCanvas *ct = new TCanvas("ct","Time on axis",0,0,900,600);
	ct->Divide(4,6);
	TH2* ff[24];
	double timeMax = timeRange(station,run,inputFolder);
	if(timeMax<binTime){
		std::cout<<"ERROR: Binning time is set higher than time in an event  ";
		return;
	}
	timeMax=ceil(timeMax);
	std::cout<<"time"<<timeMax<<typeid(timeMax).name()<<std::endl;
	for(int s=0;s<24;s++){//going though channles
		ct->cd(s+1);
		ff[s]=spectGSelect(station,run,s,binTime,inputFolder,inputTriggerPick,timeMax);
		ff[s]->Draw("colz");
	}
}
TH2* spectGSelect(int stationInput,int runInput,int channelInput,double timeBinSize,const char*  inputFolderFow,const char* triggerPickInput,double histoMax){	
	int binNumber =histoMax/timeBinSize;
	//nameing of histogram relatted to channel
        char * baseName =(char*)"Spectragram Channel:";
        char * numMod   =(char*)Form("%d",channelInput);
        char * fullSName;
        fullSName = (char*)malloc(strlen(baseName)+strlen(numMod));
        sprintf(fullSName,"%s%s",baseName,numMod);

	TH2* qh2SI = new TH2F("h2", fullSName, binNumber, 0.0, histoMax, 205, 0.0, 1600.0);
	TH1* qtimeSum = new  TH1F("ts","tss", binNumber, 0.0, histoMax);
	spectrographGen(stationInput,runInput,channelInput,qh2SI,qtimeSum,inputFolderFow,triggerPickInput);
	return qh2SI;
}
//makes indivdual spectragram
void  spectrographGen(int station,int run,int channel,TH2* h2SI,TH1* timeSum,const char*  folderLocation,const char* triggerPick){
	int numEvents;
	//mattak::Dataset d(13,517,0,"/data/desy-mirror/inbox");
        //mattak::Dataset d(station,run,0,"/data/handcarry22/rootified");
        mattak::Dataset d(station,run,0,folderLocation);
	numEvents =d.N();
	d.setEntry(0);
	const double startTime= d.header()->readout_time;//getting time of first event to find elapes time

	for(int ievent = 0; ievent < numEvents; ievent ++){
     	        d.setEntry(ievent);

		if(!(triggerReader(triggerPick,&d))){//if an event does not satify the trigger it will be skiped
			continue;
		}
		double absTime=d.header()->readout_time;
		double relTime= (absTime-startTime)/60.0;
		bool testB=d.header()->trigger_info.rf_trigger;
		bool sanTest=true;
		TGraph * wf = d.raw()->makeGraph(channel);
		int z = sizeof(wf);
                int waveformN = wf->GetN();
		int newFTSize = waveformN/2 +1;
                double x[waveformN];
                double y[waveformN];
		double freq[newFTSize];
		wf->SetLineColor(2);
		for(int k=0; k<waveformN;k++){
			x[k]=wf->GetPointX(k);
                        y[k]=wf->GetPointY(k);
		}
		FFTWComplex*wfFFT =FFTtools::doFFT(waveformN,y);

		int litN = waveformN/2+1;
		double lit[litN];
		double mag[litN];
		double phase[litN];
                int tbin = h2SI->GetXaxis()->FindBin(relTime); 
		for(int a=0; a<litN;a++){
			mag[a]=wfFFT[a].getAbs();
//			double f = a * 1600.0 / (litN);//converts y axis to frequency in MHz
                        double f = a * 1600.0 / (litN);//converts y axis to frequency in MHz
                        int fbin = h2SI->GetYaxis()->FindBin(f); 
                        h2SI->SetBinContent(tbin,fbin,mag[a] + h2SI->GetBinContent(tbin,fbin));


                       // h2SI->Fill(relTime,f,mag[a]);

			timeSum->Fill(relTime);
		}
	}
	int binNumX=h2SI->GetNbinsX();
	int binNumY=h2SI->GetNbinsY();
	for(int i=1;i<=binNumX;i++){
		for(int j=1;j<=binNumY;j++){
			double binPreNorm1= h2SI->GetBinContent(i,j);
			if(binPreNorm1==0){//bins with no content will be set to -18 db
                                double tempZ=-18.0;
				h2SI->SetBinContent(i,j,tempZ);
				continue;
			}
			double binPreNorm2= timeSum->GetBinContent(i);
			double binNorm = binPreNorm1/binPreNorm2;
			double power = 20.0* log10(binNorm);
			h2SI->SetBinContent(i,j,power);//changes powerinto db
		}
	}
}
bool triggerReader(const char* triggerSel,mattak::Dataset *dd){
	if(strncmp(triggerSel,"",2)==0){
                return true;
        }

	else if(strncmp(triggerSel,"rf",2)==0){
		return dd->header()->trigger_info.rf_trigger;
	}
	else if(strncmp(triggerSel,"force",5)==0){
                return dd->header()->trigger_info.force_trigger;
	}
	else if(strncmp(triggerSel,"pps",3)==0){
                return dd->header()->trigger_info.pps_trigger;

	}
	else if(strncmp(triggerSel,"ext",3)==0){
                return dd->header()->trigger_info.ext_trigger;

	}
	else if(strncmp(triggerSel,"radiant",7)==0){
                return dd->header()->trigger_info.radiant_trigger;

	}
	else if(strncmp(triggerSel,"lt",2)==0){
                return dd->header()->trigger_info.lt_trigger;
	}
        //return dd->header()->trigger_info.rf_trigger;

	std::cout<<"ERROR::TRIGGER NOT FOUND DEFULT VALUE RETURNED-----"<<"ERROR IN SPECTRAGRAM CODE CALL"<<std::endl<<"TRIGGER DOES NOT EXIST OR WAS TYPED IN WRONG";
	return false;
}
double timeRange(int stationInput,int runInput,const char*  inputFolderFow){
	//function finds max time across all channles regradless of filters 
	mattak::Dataset f(stationInput,runInput,0,inputFolderFow);
	double maxTime=-1.0;
	int evtMax = f.N()-1;
	for(int y=0;y<23;y++){
		f.setEntry(0);
                double time_0=f.header()->readout_time;
                f.setEntry(evtMax);
                double time_f=f.header()->readout_time;
		double timeDiff=(time_f-time_0)/60.0;
		if(timeDiff>maxTime){
			maxTime=timeDiff;
		}
	}
	return maxTime;
}

