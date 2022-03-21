#include "mattak/Pedestals.h" 
#include <iostream> 
#include "TH1.h" 


ClassImp(mattak::Pedestals); 


mattak::Pedestals::Pedestals(const rno_g_pedestal_t * peds ) 
  : Pedestals() 
{

#ifndef LIBRNO_G_SUPPORT
  std::cerr << "Not compiled with librno-g support. "<< std::endl;
  (void) peds; 
#else

  doInit(peds); 
#endif
}

void mattak::Pedestals::doInit(const rno_g_pedestal_t * peds)
{
#ifndef LIBRNO_G_SUPPORT

#else
  this->when = peds->when; 
  this->nevents = peds->nevents; 
  this->mask = peds->mask; 
  this->station_number = peds->station; 
  this->flags = peds->flags; 
  this->vbias[0] = peds->vbias[0] < 0 ? -1 : 3.3*peds->vbias[0]/4095; 
  this->vbias[1] = peds->vbias[1] < 0 ? -1 : 3.3*peds->vbias[1]/4095; 
  for (unsigned i = 0; i < mattak::k::num_radiant_channels; i++) 
  {
    memcpy(this->pedestals[i], peds->pedestals[i], sizeof(int16_t) * mattak::k::num_lab4_samples); 
  }
#endif
}

mattak::Pedestals::Pedestals(const char * pedfile) 
{
#ifndef LIBRNO_G_SUPPORT
  std::cerr << "Not compiled with librno-g support. "<< std::endl;
  (void) pedfile; 
#else

  rno_g_pedestal_t peds; 
  rno_g_file_handle_t h; 
  rno_g_init_handle(&h, pedfile,"r"); 
  rno_g_pedestal_read(h, &peds); 
  doInit(&peds); 
#endif

}

TH1 * mattak::Pedestals::pedestalHist(int chan, const char * name) const
{
  if (chan < 0 || chan >= mattak::k::num_radiant_channels) return 0; 

  TString sname;
  if (name) sname = name; 
  else sname.Form("peds%d",chan);

  TString stitle; 
  stitle.Form("Pedestals, t=%d, ch %d ; sample ; pedestal", when, chan); 
  TH1 * hist = new TH1I(sname.Data(), stitle.Data(), mattak::k::num_lab4_samples, 0, mattak::k::num_lab4_samples); 

  for (int i = 1; i<= hist->GetNbinsX(); i++) 
  {
    hist->SetBinContent(i, pedestals[chan][i-1]); 
  }
  hist->SetEntries(hist->GetNbinsX()); 
  hist->SetStats(0); 
  return hist; 
}
