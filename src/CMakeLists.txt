
SET(OTOC_LIB_SRCS
	BitUtility.cpp
	SingleTrace.cpp
	TraceState.cpp
	StateId.cpp
	StateCollection.cpp
	StateGenerator.cpp
	MixState.cpp
	HamOperator.cpp
	PowerSeries.cpp
	ScriptGenerator.cpp
	Spectrum.cpp
    NormCalculator.cpp
    MixedTrace.cpp
    otoc.cpp
    BlockMatrix.cpp
)

IF(WIN32)
	SET(OTOC_LIB_SRCS 
		${OTOC_LIB_SRCS}
        BitUtility.h
        SingleTrace.h
        TraceState.h
        StateId.h
        StateCollection.h
        StateGenerator.h
        MixState.h
        HamOperator.h
        PowerSeries.h
        ScriptGenerator.h
        Spectrum.h
        NormCalculator.h
        MixedTrace.h
        otoc.h
        BlockMatrix.h
	)	
	source_group("Header Files" FILES *.h)
ENDIF(WIN32)

SET(OTOC_MAIN_SRCS
	main.cpp)

#I like to keep test files in a separate source directory called test
file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)

INCLUDE_DIRECTORIES (${CMAKE_CURRENT_SOURCE_DIR}/../eigen/)

#find_package(Boost 1.70 REQUIRED COMPONENTS
#             unit_test_framework)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_library (otoclib STATIC ${OTOC_LIB_SRCS})
    add_executable (otoc ${OTOC_MAIN_SRCS})
    target_link_libraries(otoc ${Boost_LIBRARIES} otoclib)
    
    if (NOT WIN32)
        add_executable(otoctest ${TEST_SRCS})
        target_link_libraries(otoctest ${Boost_LIBRARIES} otoclib)
        add_test(NAME otoctest COMMAND otoctest)
     endif()
endif()
