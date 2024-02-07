R__LOAD_LIBRARY(libRootFftwWrapper.so) 
#include "FFTtools.h" 

typedef bool (*select_fn)(const mattak::Header * h); 

bool select_force_triggers(const mattak::Header *h) 
{
  return h->trigger_info.force_trigger; 
}

std::vector<TGraph *>  makeAveragePowerSpectra(mattak::Dataset & d, uint32_t channel_mask = 0xfffff, int max_N = 0, bool zero_mean = true, select_fn select = nullptr)
{
   std::vector<TGraph*> spec_sum (24); 
   int N = d.N();
   int total = 0;
   if (max_N > 0 && max_N < N) N = max_N; 
   for (int ievent = 0; ievent < N; ievent ++)
   {
      d.setEntry(ievent);
      if (select && !select(d.header()))
        continue; 

      total++;
      for (int ch = 0; ch < 24; ch++) 
      {
        if (!(channel_mask & (1 << ch)) ) continue; 

        TGraph * wf = d.raw()->makeGraph(ch);
        if (zero_mean) 
        {
          double mean = wf->GetMean(2); //this 2 means second axis 
          for (int i = 0; i < wf->GetN(); i++) wf->GetY()[i]-=mean; 
        }
        for (int i = 0; i < wf->GetN(); i++) wf->GetY()[i]*= (2500./4095);  //approximate adc to mV conversion. 

        TGraph * spec = FFTtools::makePowerSpectrum(wf);

        for (int i = 0; i < spec->GetN(); i++) 
        {
          spec->GetX()[i] *= 1e3; // convert to mW/MHz (divide by 50 ohms, convert to mW, note we had mV)
        }

        for (int i = 0; i < spec->GetN(); i++) 
        {
          spec->GetY()[i] /= (50.*1e3 * spec->GetX()[1] * wf->GetN()); // convert to mW/MHz (divide by 50 ohms, convert to mW, note we had mV)
        }
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
       if (total) 
         spec_sum[ch]->GetY()[i] /= total; //take average
       if (spec_sum[ch]->GetY()[i] > 0)   
       {
         spec_sum[ch]->GetY()[i] = 10 * TMath::Log10(spec_sum[ch]->GetY()[i]);
       }
       else spec_sum[ch]->GetY()[i] = -100; //just some small value to avoid -infinity, adjust as necessary
     }
     spec_sum[ch]->SetTitle(Form("Power Spectrum Ch %d", ch)); 
     spec_sum[ch]->GetXaxis()->SetTitle("Freq [MHz]"); 
     spec_sum[ch]->GetYaxis()->SetTitle("Power at digitizer [dBm/MHz]"); 
   }

   return spec_sum;
}

TGraph * makeAveragePowerSpectrum(int ch, mattak::Dataset & d, int max_N = 0, bool zero_mean = true, select_fn select = nullptr)
{
   std::vector<TGraph*> spec = makeAveragePowerSpectra(d, 1 << ch, max_N, zero_mean, select); 
   return spec[ch];
}

TGraph * makeAverageForcedPowerSpectrum(int ch, mattak::Dataset & d, int max_N = 0, bool zero_mean = true)
{
   std::vector<TGraph*> spec = makeAveragePowerSpectra(d, 1 << ch, max_N, zero_mean, select_force_triggers); 
   return spec[ch];
}
