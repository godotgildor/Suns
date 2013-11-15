all: debian zipfile
.PHONY: all

VERSION := 1.0.0

SRC_DIR := src
INSTALL_DIR := installer
BUILD_DIR := build

ZIP_DIR := suns
ZIP_FILE := suns.zip

PLUGIN_INIT   := $(SRC_DIR)/__init__.py
PLUGIN_MOTIFS := $(SRC_DIR)/suns_aa_motifs.py
PLUGIN_PATH   := usr/lib/python2.7/dist-packages/pmg_tk/startup

DEB_DIR := debian
DEB_BUILD_DIR := $(DEB_DIR)/build
DEB_NAME := pymol-suns-search
DEB_VERSION := 1
DEB_FULL := $(DEB_NAME)_$(VERSION)-$(DEB_VERSION)_all.deb
DOC_PATH := usr/share/doc/$(DEB_NAME)

debian: $(BUILD_DIR)/$(DEB_FULL)
.PHONY: debian

$(BUILD_DIR)/$(DEB_FULL): \
    $(DEB_DIR)/control \
    $(PLUGIN_INIT) \
    $(PLUGIN_MOTIFS)
	mkdir -p $(BUILD_DIR) \
	 $(DEB_BUILD_DIR)/DEBIAN \
	 $(DEB_BUILD_DIR)/$(PLUGIN_PATH) \
	 $(DEB_BUILD_DIR)/$(DOC_PATH)
	cp $(PLUGIN_INIT) $(DEB_BUILD_DIR)/$(PLUGIN_PATH)/suns-search.py
	cp $(PLUGIN_MOTIFS) $(DEB_BUILD_DIR)/$(PLUGIN_PATH)
	cp $(DEB_DIR)/copyright $(DEB_BUILD_DIR)/$(DOC_PATH)/copyright
	cp $(DEB_DIR)/control $(DEB_BUILD_DIR)/DEBIAN/control
	gzip -9 -c $(DEB_DIR)/changelog > $(BUILD_DIR)/changelog.gz
	chmod 644 $(BUILD_DIR)/changelog.gz
	cp $(BUILD_DIR)/changelog.gz $(DEB_BUILD_DIR)/$(DOC_PATH)
	rm $(BUILD_DIR)/changelog.gz
	ln -fs $(DOC_PATH)/changelog.gz \
	      $(DEB_BUILD_DIR)/$(DOC_PATH)/changelog.Debian.gz
	find $(DEB_BUILD_DIR) -type d | xargs chmod 755
	fakeroot dpkg-deb --build $(DEB_BUILD_DIR) $(BUILD_DIR)/$(DEB_FULL)

zipfile: $(INSTALL_DIR)/$(ZIP_FILE)
.PHONY: zipfile

$(INSTALL_DIR)/$(ZIP_FILE): $(SRC_DIR)/__init__.py $(SRC_DIR)/pika
	ln -s src suns
	zip -r $(INSTALL_DIR)/$(ZIP_FILE) suns
	rm suns

.PHONY: suns_aa_motifs
suns_aa_motifs: 
	python motifs/make_word_dicts.py motifs/ src/suns_aa_motifs.py

.PHONY: clean
clean:
	rm -rf $(INSTALL_DIR)/$(ZIP_FILE)
