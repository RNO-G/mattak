#include "mattak/DAQStatus.h" 
#include <iostream> 


ClassImp(mattak::DAQStatus); 


mattak::DAQStatus::DAQStatus(const rno_g_daqstatus_t * status) 
  : DAQStatus() 
{

#ifndef LIBRNO_G_SUPPORT
  std::cerr << "Not compiled with librno-g support. "<< std::endl;
  (void) status; 
#else





#endif

}
