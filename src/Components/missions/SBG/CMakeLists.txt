esma_set_this()


# these are executable scripts
set (PYSCRIPTS
    sbg_vlidort_setup.py
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/missions/${this}
    )

# these are non-executaables
set (FILES
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

# these are executable scripts
set (PYSCRIPTS
    benchmarking/sbg_vlidort_pyexample.py
    benchmarking/lamb_vlidort_pyexample.py
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/missions/${this}/benchmarking
    )

# these are non-executaables
set (FILES
    benchmarking/sbg_vlidort_run.j
    benchmarking/sbg_files.yaml
    benchmarking/lamb_vlidort_run.j
    benchmarking/lamb_vlidort.yaml
    benchmarking/lamb_files.yaml
    benchmarking/m2_aop.yaml
    benchmarking/writers.py
    )

install (
    FILES ${FILES}
    DESTINATION bin/missions/${this}/benchmarking
    )

