esma_set_this()


# these are executable scripts
set (PYSCRIPTS
    sbg_vlidort_pyexample.py
    sbg_vlidort_setup.py
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/missions/${this}
    )

# these are non-executaables
set (FILES
    inputs_vlidort.py
    sbg_vlidort.yaml
    sbg.yaml
    ssd650.yaml
    sbg_files.yaml
    m2_aop.yaml
    sbg_vlidort_run.j
    )

install (
    FILES ${FILES}
    DESTINATION bin/missions/${this}
    )

