
.PHONY: all configure clean cleaner install cmake-build cmake-clean cmake-install

CMAKE?=cmake
CCMAKE?=ccmake 


all: cmake-build 
clean: cmake-clean
install: cmake-install 


cmake-build: build/Makefile 
	@touch CMakeLists.txt
	@+make -C  ./build


configure: build/Makefile 
	@$(CCMAKE) . build 

cmake-install: 
	@make -C ./build install 

build/Makefile: 
	@echo "Setting up cmake build."
	@mkdir -p build 
#	rm -rf bin
	@cd build && $(CMAKE) ../ 
#	ln -sf build/bin bin 

distclean: 
	@echo "Removing cmake directory" 
	@rm -rf build 

cmake-clean: build/Makefile 
	@make -C ./build clean 

