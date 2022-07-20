
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/src/Utilities)

add_library(modelLibrary
	model.cpp
	../Utilities/iio.c
	../Utilities/LibImages.cpp
	../Utilities/comparators.cpp
	../Utilities/mt19937ar.c
	../Utilities/PatchManager/imagePatchManager.cpp
	../Utilities/PatchManager/databasePatchManager.cpp
	../Utilities/PartitionTree/vptree.cpp
	../Utilities/PartitionTree/forestManager.cpp
	../Utilities/PatchSearch/localRefinement.cpp
	../Utilities/nfa.cpp
	)


set_target_properties(modelLibrary PROPERTIES
	SOVERSION ${PROJECT_VERSION})

install(TARGETS modelLibrary
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    COMPONENT libs)
