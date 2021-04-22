#include "mattak/Sensors.h" 
#include <iostream> 


ClassImp(mattak::Sensors); 


mattak::Sensors::Sensors(const rno_g_sensors_t * sensors) 
  : Sensors() 
{

#ifndef LIBRNO_G_SUPPORT
  std::cerr << "Not compiled with librno-g support. "<< std::endl;
  (void) sensors; 
#else





#endif

}
