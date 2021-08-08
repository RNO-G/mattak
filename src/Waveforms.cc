#include "mattak/Waveforms.h" 
#include <iostream> 


ClassImp(mattak::Waveforms); 


mattak::Waveforms::Waveforms(const rno_g_waveform_t * wf ) 
  : Waveforms() 
{

#ifndef LIBRNO_G_SUPPORT
  std::cerr << "Not compiled with librno-g support. "<< std::endl;
  (void) wf; 
#else

  this->run_number = wf->run_number; 
  this->event_number = wf->event_number; 
  this->buffer_length = wf->radiant_nsamples; 
  this->station_number = wf->station; 
  for (unsigned i = 0; i < mattak::k::num_radiant_channels; i++) 
  {
    memcpy(this->radiant_data[i], wf->radiant_waveforms[i], sizeof(int16_t) * mattak::k::num_radiant_samples); 
  }

#endif

}
