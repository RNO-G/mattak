#include "mattak/Header.h" 
#include <iostream> 


ClassImp(mattak::Header); 


mattak::Header::Header(const rno_g_header_t * head) 
  : mattak::Header() 
{

#ifndef LIBRNO_G_SUPPORT
  std::cerr << "Not compiled with librno-g support. "<< std::endl;
  (void) head; 
#else





#endif

}
