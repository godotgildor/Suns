all: zipfile
.PHONY: all

VERSION := 1.0.0

SRC_DIR := src
ZIP_DIR := suns
ZIP_FILE := suns.zip

INSTALL_DIR := installer

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
