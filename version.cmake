
execute_process ( COMMAND git describe --tags --always --dirty OUTPUT_VARIABLE GIT_VERSION ERROR_QUIET)

string(STRIP "${GIT_VERSION}" GIT_VERSION) 
SET(GIT_VERSION "\"${GIT_VERSION}\"")

message ( "GIT VERSION IS ${GIT_VERSION}" )
configure_file("Version.cc.in" "Version.cc")


