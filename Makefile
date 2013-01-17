install_suns.py: \
    install/install_suns_header.template \
    src/suns_search.py                   \
    install/install_suns_middle.template \
    src/suns_plugin.py                   \
    install/install_suns_footer.template
	cat install/install_suns_header.template \
	    src/suns_search.py                   \
            install/install_suns_middle.template \
            src/suns_plugin.py                   \
	    install/install_suns_footer.template \
	    >install/install_suns.py
