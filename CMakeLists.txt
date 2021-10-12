cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(mattak)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "RelWithDebInfo")

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED)



set(libname mattak)
set(DICTNAME G__${libname})


if(DEFINED ENV{RNO_G_INSTALL_DIR})
    message("RNO_G_INSTALL_DIR is set to $ENV{RNO_G_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{RNO_G_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{RNO_G_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{RNO_G_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{RNO_G_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{RNO_G_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{RNO_G_INSTALL_DIR}/include)
else()
    message("Will install to /usr/local  RNO_G_INSTALL_DIR not set")
    set(UTIL_LIB_DIR /usr/local/lib)
    set(UTIL_INC_DIR /usr/local/include)
    set(UTIL_BIN_DIR /usr/local/bin)
    set(UTIL_SHARE_DIR /usr/local/share)
    set(LD_UTIL /usr/local/lib)
    set(INC_UTIL /usr/local/include)
endif()


include(${ROOT_USE_FILE} )
include_directories( ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/src ${UTIL_INC_DIR})

option (LIBRNO_G_SUPPORT "Support for reading / converting from raw files." OFF) 

set (LIBRNO_G_PATH  "../librno-g" CACHE STRING "Additional path  to search for librno-g (only necessary if not in RNO_G_INSTALL_DIR/lib)" ) 


set(H_FILES 
  mattak/Header.h 
  mattak/Constants.h 
  mattak/DAQStatus.h
  mattak/Sensors.h
  mattak/TriggerInfo.h
  mattak/GNSS.h
  mattak/Sensors.h
  mattak/Waveforms.h
  mattak/Pedestals.h
  mattak/Converter.h
  mattak/RunInfo.h
 ) 



set(CC_FILES 
  src/Header.cc
  src/DAQStatus.cc 
  src/Sensors.cc
  src/TriggerInfo.cc
  src/GNSS.cc
  src/Sensors.cc
  src/Waveforms.cc 
  src/Pedestals.cc 
  src/Converter.cc 
  src/RunInfo.cc 
  Version.cc 
  Empty.cc
  )

set (INSTALL_SCRIPTS
  scripts/rno-g-convert-run
  scripts/rno-g-make-eventlists
)

execute_process ( COMMAND git describe --always --dirty --match 'notatag' OUTPUT_VARIABLE GIT_VERSION ERROR_QUIET)
string(STRIP "${GIT_VERSION}" GIT_VERSION) 
SET(GIT_VERSION "\"${GIT_VERSION}\"")

message ( "GIT VERSION IS ${GIT_VERSION}" )

configure_file("src/Version.cc.in" "Version.cc")


ROOT_GENERATE_DICTIONARY(${DICTNAME} ${H_FILES} LINKDEF src/LinkDef.h) 

add_library(${libname} SHARED ${CC_FILES} ${DICTNAME}.cxx)

install (DIRECTORY src/mattak DESTINATION ${UTIL_INC_DIR})
install (PROGRAMS ${INSTALL_SCRIPTS} DESTINATION ${UTIL_BIN_DIR})

install (TARGETS ${libname}
	ARCHIVE DESTINATION ${UTIL_LIB_DIR}
	LIBRARY DESTINATION ${UTIL_LIB_DIR}
	RUNTIME DESTINATION ${UTIL_BIN_DIR}
  )


if (LIBRNO_G_SUPPORT) 
  add_definitions( -DLIBRNO_G_SUPPORT )
  find_library(RNO_G_LIB rno-g PATHS ${LIBRNO_G_PATH}/build ${UTIL_LIB_DIR})
  find_package(ZLIB REQUIRED)
  include_directories(${UTIL_INC_DIR}) 
  include_directories(${LIBRNO_G_PATH}/src) 
  include_directories(${LIBNUPHASE_PATH}) 
  target_link_libraries(${libname} ${RNO_G_LIB} ${ZLIB_LIBRARIES} ${ROOT_LIBRARIES}) 
else()
  target_link_libraries(${libname} ${ROOT_LIBRARIES}) 
endif() 

set (BINSRCDIR prog) 
macro(add_program binary_name ) 
  add_executable(${binary_name} ${BINSRCDIR}/${binary_name}.cc) 
  target_link_libraries(${binary_name}  ${ROOT_LIBRARIES} ${libname})
  if (LIBNUPHASE_SUPPORT)
    target_link_libraries(${binary_name}  ${ZLIB_LIBRARIES} rno-g)
  endif()
  INSTALL (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${binary_name} DESTINATION ${UTIL_BIN_DIR})
endmacro() 




if (LIBRNO_G_SUPPORT) 
add_program(rno-g-convert) 
endif( )
add_program(rno-g-combine) 
add_program(rno-g-make-eventlist) 
