install_suns.py: \
    install_suns_header.template \
    src/suns_search.py           \
    install_suns_middle.template \
    src/suns_plugin.py           \
    install_suns_footer.template
	cat install_suns_header.template \
	    src/suns_search.py           \
            install_suns_middle.template \
            src/suns_plugin.py           \
	    install_suns_footer.template \
	    >install_suns.py
