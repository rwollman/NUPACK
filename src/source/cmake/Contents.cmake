set(LIBNUPACK_FILES
    source/Partition.cc
    source/Model.cc
    source/ComplexSet.cc
    source/State.cc
    source/Loop.cc
    # source/thermo/Thermo.cc source/thermo/Thermo2.cc
    source/Utility.cc
    source/Concentration.cc
    source/Types.cc
    source/Structure.cc
    source/Domain.cc
    source/Config.cc
    source/Constants.cc
    source/Runtime.cc
    source/Costs.cc
)

set(NUPACK_MODULE_FILES
    bind/Constants.cc
    bind/Design.cc
    bind/Document.cc
    bind/Core.cc
    bind/Math.cc
    bind/Model.cc
    bind/Thermo.cc
    bind/Thermo2.cc
)

set(NUPACK_EXECUTABLES complexdefect complexes concentrations count energy mfe pairs pfunc prob sample subopt)

set(NUPACK_PYTHON_FILES
    __init__.py
    core.py concentration.py std.py constants.py model.py rotation.py analysis.py utility.py thermo.py named.py
    design/core.py design/components.py design/objectives.py design/results.py design/weights.py design/__init__.py
    design/trials.py design/parsers.py
)


set(NUPACK_PYTHON_TEST_FILES
    __init__.py
    test_concentration.py
    test_model.py
    test_analysis.py
    test_design.py
    test_notebooks.py
    test_utilities.py
)

