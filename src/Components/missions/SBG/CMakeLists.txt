esma_set_this()


# these are executable scripts
set (PYSCRIPTS
    sbg_vlidort_pyexample.py
    )

set (SCRIPTS
    setup_env
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/missions/${this}
    )

# these are non-executaables
set (RCFILES
    sbg.pcf
    ssd650.pcf
    imager_files.pcf
    imager_sampler.pcf
    )

install (
    FILES ${RCFILES}
    DESTINATION bin/missions/${this}
    )

