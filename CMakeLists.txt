INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

#Change PROJECT_NAME to the name of your project
PROJECT(MRParameterSuite)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
SET(Required_Packages
ITK3.6
)

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
ITKIO
ITKCommon
)

#the following block of code is an example of how to build an executable in
#cmake.  Unmodified, it will add an executable called "MyExe" to the project.
#MyExe will be built using the files MyClass.h and MyClass.cxx, and it will
#be linked to all the libraries you specified above. 
#You can build more than one executable per project

#Add a definition for making programs single threaded for ITK < 3.
IF ( ${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR} LESS 3.0 )
  ADD_DEFINITIONS(-DNO_MULTI_THREADING)
ENDIF(${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR} LESS 3.0 )

SET(CurrentExe "ImageCompare")
ADD_EXECUTABLE(${CurrentExe} ImageCompare.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
SET(CurrentExe "BrukerT2Map")
SET(CurrentSrc
BrukerT2Map.cxx
itkBruker2DSEQImageIO.cxx
itkBruker2DSEQImageIOFactory.cxx
)
ADD_EXECUTABLE(${CurrentExe} ${CurrentSrc}) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

SET(CurrentExe "BrukerT1Map")
SET(CurrentSrc
BrukerT1Map.cxx
itkBruker2DSEQImageIO.cxx
itkBruker2DSEQImageIOFactory.cxx
)
ADD_EXECUTABLE(${CurrentExe} ${CurrentSrc}) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

SET(CurrentExe "PhilipsT2Map")
SET(CurrentSrc
PhilipsT2Map.cxx
itkPhilipsPAR.cxx
itkPhilipsRECImageIO.cxx
itkPhilipsRECImageIOFactory.cxx
)
ADD_EXECUTABLE(${CurrentExe} ${CurrentSrc}) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

SET(CurrentExe "PhilipsT1Map")
SET(CurrentSrc
PhilipsT1Map.cxx
itkPhilipsPAR.cxx
itkPhilipsRECImageIO.cxx
itkPhilipsRECImageIOFactory.cxx
)
ADD_EXECUTABLE(${CurrentExe} ${CurrentSrc}) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.
SET(CurrentParam
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
ADD_TEST(BrukerT2MapTest BrukerT2Map ${CurrentParam})
ADD_TEST(CompareImage1 ImageCompare ${CMAKE_SOURCE_DIR}/brukert2test.mhd ${CMAKE_CURRENT_BINARY_DIR}/brkt2test.mhd)
ADD_TEST(CompareImage2 ImageCompare ${CMAKE_SOURCE_DIR}/brukert2exptest.mhd ${CMAKE_CURRENT_BINARY_DIR}/brkt2exptest.mhd)
ADD_TEST(CompareImage3 ImageCompare ${CMAKE_SOURCE_DIR}/brukert2consttest.mhd ${CMAKE_CURRENT_BINARY_DIR}/brkt2consttest.mhd)
ADD_TEST(CompareImage4 ImageCompare ${CMAKE_SOURCE_DIR}/brukert2rsqrtest.mhd ${CMAKE_CURRENT_BINARY_DIR}/brkt2rsqrtest.mhd)

SET(CurrentParam
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
ADD_TEST(BrukerT1MapTest BrukerT1Map ${CurrentParam})
ADD_TEST(CompareImage5 ImageCompare ${CMAKE_SOURCE_DIR}/brukert1test.mhd ${CMAKE_CURRENT_BINARY_DIR}/brkt1test.mhd)

SET(CurrentParam
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
ADD_TEST(PhilipsT2MapTest PhilipsT2Map ${CurrentParam})
ADD_TEST(CompareImage6 ImageCompare ${CMAKE_SOURCE_DIR}/philipst2test.mhd ${CMAKE_CURRENT_BINARY_DIR}/philt2test.mhd)

SET(CurrentParam
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
ADD_TEST(PhilipsT1MapTest PhilipsT1Map ${CurrentParam})
ADD_TEST(CompareImage7 ImageCompare ${CMAKE_SOURCE_DIR}/philipst1test.mhd ${CMAKE_CURRENT_BINARY_DIR}/philt1test.mhd)
