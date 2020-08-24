

np_add_executable(count count.c)
np_add_executable(complexdefect defect.c)
np_add_executable(energy energy.c)
np_add_executable(mfe mfe.c)
np_add_executable(pairs pairs.c)
np_add_executable(pfunc pfunc.c)
np_add_executable(prob prob.c)
np_add_executable(subopt subopt.c)

if(SAMPLE)
    np_add_executable(sample sample.c)
endif(SAMPLE)
