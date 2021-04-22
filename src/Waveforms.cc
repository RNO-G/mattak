#include "mattak/Waveforms.h" 
#include <iostream> 


ClassImp(mattak::Waveforms); 


mattak::Waveforms::Waveforms(const rno_g_waveform_t * wf) 
  : Waveforms() 
{

#ifndef LIBRNO_G_SUPPORT
  std::cerr << "Not compiled with librno-g support. "<< std::endl;
  (void) wf; 
#else




#endif

}
