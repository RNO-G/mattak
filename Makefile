# This Makefile is just a cmake wrapper

.PHONY: all configure clean cleaner install cmake-build cmake-clean cmake-install cmake-rootless-build cmake-rootless-clean cmake-rootless-configure

CMAKE?=cmake
CCMAKE?=ccmake
builddir=build
rootlessdir=build-noroot
SYSTEMD_DIR=/etc/systemd/system
SERVICE_FILE=./scripts/services/rno-g-autoconverter@.service
TARGET_FILE=./scripts/services/rno-g-autoconverter.target

all: cmake-build #cmake-rootless-build
rootless: cmake-rootless-build
clean: cmake-clean cmake-rootless-clean
install: cmake-install

clean: cmake-clean
install: cmake-install



cmake-build: $(builddir)/Makefile
	@+make -C  ./$(builddir)

cmake-rootless-build: $(rootlessdir)/Makefile
	@+make -C  ./$(rootlessdir)


configure: $(builddir)/Makefile
	@$(CCMAKE) . $(builddir)

rootless-configure: $(rootlessdir)/Makefile
	@$(CCMAKE) . $(rootlessdir)

cmake-install:
	@make -C ./$(builddir) install

$(rootlessdir)/Makefile:
	@echo "Setting up rootless cmake build."
	@mkdir -p $(rootlessdir)
	@cd $(rootlessdir) && $(CMAKE) -DROOTLESS=yes ../


$(builddir)/Makefile:
	@echo "Setting up cmake build."
	@mkdir -p $(builddir)
	@cd $(builddir) && $(CMAKE) ../

distclean:
	@echo "Removing cmake directory"
	@rm -rf $(builddir) $(rootlessdir)

cmake-clean: build/Makefile
	@make -C ./$(builddir) clean

cmake-rootless-clean: build-noroot/Makefile
	@make -C ./$(rootlessdir) clean

# Install systemd service and target files
install-systemd:
	@echo "Installing systemd service and target files..."
	@cp $(SERVICE_FILE) $(SYSTEMD_DIR)/rno-g-autoconverter@.service
	@cp $(TARGET_FILE) $(SYSTEMD_DIR)/rno-g-autoconverter.target
	@systemctl daemon-reload
	@echo "Done. Start with: systemctl start rno-g-autoconverter.target"
