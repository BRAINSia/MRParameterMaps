cmake_minimum_required(VERSION 2.8)

project(MRParameterSuite)

#The following lines are required to use Dart
enable_testing()

find_package(ITK 4.7 REQUIRED)
include(${ITK_USE_FILE})

add_library(Bruker2DSEQImageIO
  itkBruker2DSEQImageIO.cxx
  itkBruker2DSEQImageIOFactory.cxx
)
target_link_libraries(Bruker2DSEQImageIO ${ITK_LIBRARIES})
add_library(PhilipsRECImageIO
  itkPhilipsRECImageIO.cxx itkPhilipsRECImageIOFactory.cxx itkPhilipsPAR.cxx
)
target_link_libraries(PhilipsRECImageIO ${ITK_LIBRARIES})

set(CUSTOM_IO_LIBRARIES PhilipsRECImageIO Bruker2DSEQImageIO )

set(CurrentExe ImageCompare)
add_executable(${CurrentExe} ImageCompare.cxx)
set_target_properties(${CurrentExe} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${CurrentExe} ${CUSTOM_IO_LIBRARIES})

set(CurrentExe BrukerT2Map)
set(CurrentSrc
  BrukerT2Map.cxx
  )
add_executable(${CurrentExe} ${CurrentSrc})

target_link_libraries(${CurrentExe} ${CUSTOM_IO_LIBRARIES})


set(CurrentExe BrukerT1Map)
set(CurrentSrc
  BrukerT1Map.cxx
  )
add_executable(${CurrentExe} ${CurrentSrc})
target_link_libraries(${CurrentExe} ${CUSTOM_IO_LIBRARIES})

set(CurrentExe PhilipsT2Map)
set(CurrentSrc
  PhilipsT2Map.cxx
  )

add_executable(${CurrentExe} ${CurrentSrc})
target_link_libraries(${CurrentExe} ${CUSTOM_IO_LIBRARIES})

set(CurrentExe PhilipsT1Map)
set(CurrentSrc
  PhilipsT1Map.cxx
  itkPhilipsPAR.cxx
  )
add_executable(${CurrentExe} ${CurrentSrc})
target_link_libraries(${CurrentExe} ${CUSTOM_IO_LIBRARIES})


find_package(SlicerExecutionModel REQUIRED GenerateCLP)
include(${GenerateCLP_USE_FILE})

set(CurrentExe T1rhoMap)
set(CurrentSrc T1rhoMap.cxx
  itkFDFCommonImageIO.cxx
  itkFDFImageIO.cxx
  itkFDFImageIOFactory.cxx
)
GenerateCLP(CurrentSrc T1rhoMap.xml)
add_executable(${CurrentExe} ${CurrentSrc})
target_link_libraries(${CurrentExe} ${CUSTOM_IO_LIBRARIES})

set(CurrentParam
  ${CMAKE_SOURCE_DIR}/DCB021304.oj1/3/pdata/1/2dseq
  ${CMAKE_CURRENT_BINARY_DIR}/brkt2test.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/brkt2exptest.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/brkt2consttest.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/brkt2rsqrtest.mhd
  0
  0
  "2.0"
  1000
  )

add_test(NAME BrukerT2MapTest
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND $<TARGET_FILE:BrukerT2Map> ${CurrentParam})
add_test(NAME CompareImage1
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND $<TARGET_FILE:ImageCompare> ${CMAKE_SOURCE_DIR}/brukert2test.mhd ${CMAKE_CURRENT_BINARY_DIR}/brkt2test.mhd)
add_test(NAME CompareImage2
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND $<TARGET_FILE:ImageCompare> ${CMAKE_SOURCE_DIR}/brukert2exptest.mhd ${CMAKE_CURRENT_BINARY_DIR}/brkt2exptest.mhd)
add_test(NAME CompareImage3
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND $<TARGET_FILE:ImageCompare> ${CMAKE_SOURCE_DIR}/brukert2consttest.mhd ${CMAKE_CURRENT_BINARY_DIR}/brkt2consttest.mhd)
add_test(NAME CompareImage4
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND $<TARGET_FILE:ImageCompare> ${CMAKE_SOURCE_DIR}/brukert2rsqrtest.mhd ${CMAKE_CURRENT_BINARY_DIR}/brkt2rsqrtest.mhd)

set(CurrentParam
  ${CMAKE_SOURCE_DIR}/DCB021304.oj1/5/pdata/1/2dseq
  ${CMAKE_CURRENT_BINARY_DIR}/brkt1test.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/brkt1exptest.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/brkt1consttest.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/brkt1rsqrtest.mhd
  0
  5
  "5.0"
  20546109
  )
add_test(NAME BrukerT1MapTest
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND $<TARGET_FILE:BrukerT1Map> ${CurrentParam})

add_test(NAME CompareImage5
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND $<TARGET_FILE:ImageCompare> ${CMAKE_SOURCE_DIR}/brukert1test.mhd ${CMAKE_CURRENT_BINARY_DIR}/brkt1test.mhd)

set(CurrentParam
  ${CMAKE_SOURCE_DIR}/T2_Map_4_3.PAR
  ${CMAKE_CURRENT_BINARY_DIR}/philt2test.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/philt2exptest.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/philt2consttest.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/philt2rsqrtest.mhd
  0
  1
  "2.0"
  14
  )
add_test(NAME PhilipsT2MapTest
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND $<TARGET_FILE:PhilipsT2Map> ${CurrentParam})

add_test(NAME CompareImage6
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND $<TARGET_FILE:ImageCompare> ${CMAKE_SOURCE_DIR}/philipst2test.mhd ${CMAKE_CURRENT_BINARY_DIR}/philt2test.mhd)

set(CurrentParam
  ${CMAKE_SOURCE_DIR}/T1_LL_10_1.PAR
  ${CMAKE_CURRENT_BINARY_DIR}/philt1test.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/philt1exptest.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/philt1consttest.mhd
  ${CMAKE_CURRENT_BINARY_DIR}/philt1rsqrtest.mhd
  0
  1
  "5.0"
  50
  )
add_test(NAME PhilipsT1MapTest
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND $<TARGET_FILE:PhilipsT1Map> ${CurrentParam})
add_test(NAME CompareImage7
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND $<TARGET_FILE:ImageCompare> ${CMAKE_SOURCE_DIR}/philipst1test.mhd ${CMAKE_CURRENT_BINARY_DIR}/philt1test.mhd)
