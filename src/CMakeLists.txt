set(SAS_SRC
        AliasTable.h
        AliasTable.cpp
        SystematicAliasSampler.h)

set(LIB_TYPE STATIC)

add_library(sas ${LIB_TYPE} ${SAS_SRC})
