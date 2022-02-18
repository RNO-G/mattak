R__LOAD_LIBRARY(libRootFftwWrapper.so) 
#include "FFTtools.h" 

TGraph * makeAveragePowerSpectrum(int ch, mattak::Dataset & d, bool arithmetic_mean = true, int max_N = 0, bool zero_mean = true)
{
   TGraph * spec_sum = 0;
   int N = d.N();
   if (max_N > 0 && max_N < N) N = max_N; 
   for (int ievent = 0; ievent < N; ievent ++)
   {
      d.setEntry(ievent);
      TGraph * wf = d.raw()->makeGraph(ch);
      if (zero_mean) 
      {
        double mean = wf->GetMean(2); //this 2 means second axis 
        for (int i = 0; i < wf->GetN(); i++) wf->GetY()[i]-=mean; 
      }

      TGraph * spec = arithmetic_mean ? FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(wf) : FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(wf);
      //no longer need wf, let's delete it
      delete wf;
      // if spec_sum hasn't been set up yet, let's just set it to spec
      if (!spec_sum) spec_sum = spec;
      else
       {
          //otherwise, add spec to it and we no longer need spec
           for (int i = 0; i < spec_sum->GetN(); i++) spec_sum->GetY()[i] += spec->GetY()[i];
           delete spec;
       }
   }
   //now divide by N and convert to dB
   for (int i = 0; i < spec_sum->GetN(); i++)
   {
     spec_sum->GetY()[i] /= N; //take average
     if (arithmetic_mean) 
     {
       if (spec_sum->GetY()[i] > 0)   spec_sum->GetY()[i] = 10 * TMath::Log10(spec_sum->GetY()[i]);
       else spec_sum->GetY()[i] = -100; //just some small value to avoid -infinity, adjust as necessary
     }
   }
   spec_sum->SetTitle(Form("Power Spectrum Ch %d", ch)); 
   spec_sum->GetXaxis()->SetTitle("Freq [MHz]"); 
   spec_sum->GetYaxis()->SetTitle("dB [arb]"); 

   return spec_sum;
}
