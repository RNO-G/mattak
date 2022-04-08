R__LOAD_LIBRARY(libRootFftwWrapper.so) 
#include "FFTtools.h" 

std::vector<TGraph *>  makeAveragePowerSpectra(mattak::Dataset & d, uint32_t channel_mask = 0xfffff, bool arithmetic_mean = true, int max_N = 0, bool zero_mean = true)
{
   std::vector<TGraph*> spec_sum (24); 
   int N = d.N();
   if (max_N > 0 && max_N < N) N = max_N; 
   for (int ievent = 0; ievent < N; ievent ++)
   {
      d.setEntry(ievent);
      for (int ch = 0; ch < 24; ch++) 
      {
        if (!(channel_mask & (1 << ch)) ) continue; 

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
        if (!spec_sum[ch]) spec_sum [ch]= spec;
        else
         {
            //otherwise, add spec to it and we no longer need spec
             for (int i = 0; i < spec_sum[ch]->GetN(); i++) spec_sum[ch]->GetY()[i] += spec->GetY()[i];
             delete spec;
         }
      }
   }

   for (int ch = 0; ch < 24; ch++) 
   {

     if (!(channel_mask & (1 << ch)) ) continue; 

     //now divide by N and convert to dB
     for (int i = 0; i < spec_sum[ch]->GetN(); i++)
     {
       spec_sum[ch]->GetY()[i] /= N; //take average
       if (arithmetic_mean) 
       {
         if (spec_sum[ch]->GetY()[i] > 0)   spec_sum[ch]->GetY()[i] = 10 * TMath::Log10(spec_sum[ch]->GetY()[i]);
         else spec_sum[ch]->GetY()[i] = -100; //just some small value to avoid -infinity, adjust as necessary
       }
     }
     spec_sum[ch]->SetTitle(Form("Power Spectrum Ch %d", ch)); 
     spec_sum[ch]->GetXaxis()->SetTitle("Freq [MHz]"); 
     spec_sum[ch]->GetYaxis()->SetTitle("dB [arb]"); 
   }

   return spec_sum;
}

TGraph * makeAveragePowerSpectrum(int ch, mattak::Dataset & d, bool arithmetic_mean = true, int max_N = 0, bool zero_mean = true)
{
   std::vector<TGraph*> spec = makeAveragePowerSpectra(d, 1 << ch, arithmetic_mean, max_N, zero_mean); 
   return spec[ch];
}
