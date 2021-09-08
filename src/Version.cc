#include "mattak/Version.h" 


const char * the_version = GIT_COMMIT_HASH; 

const char * mattak::version() { return the_version; } 

