/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: PhilipsT2Map.cxx,v $
  Language:  C++
  Date:      $Date: 2008/06/05 13:44:35 $
  Version:   $Revision: 1.00 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkWin32Header.h"
#include <iostream>
#include <fstream>
#include "itkNumericTraits.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPhilipsRECImageIO.h"
#include "itkMRT2ParameterMap3DImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkShiftScaleInPlaceImageFilter.h"
#include "itkMetaDataObject.h"
#include "itkVectorContainer.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "vnl/vnl_vector_fixed.h"
#include <algorithm>

int main(int argc, char **argv)
{
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImage outputT2Image outputExpConstImage";
    std::cerr << " outputConstImage outputRSquaredImage [T2map=0,R2map=1]";
    std::cerr << " [algorithm] [maxT2Time] [threshold]" << std::endl;
    return 1;
    }

  const char * inputFilename   = argv[1];
  const char * outputT2Filename = argv[2];
  const char * outputExpConstFilename = argv[3];
  const char * outputConstFilename = argv[4];
  const char * outputRSquaredFilename = argv[5];

  // Philips PAR/REC images are signed 16 bit and the
  // reader supports 4D images.
  typedef itk::Image< short, 4 > PhilipsImageType;
  typedef itk::Image< float, 4 > ImageType4D;
  typedef itk::Image< float, 3 > ImageType;

  typedef itk::PhilipsRECImageIO
    PhilipsRECImageIOType;
  typedef itk::MRT2ParameterMap3DImageFilter< ImageType::PixelType >
    MRT2ParameterMap3DImageFilterType;
  typedef itk::VectorIndexSelectionCastImageFilter< 
    MRT2ParameterMap3DImageFilterType::OutputImageType, ImageType >
    VectorIndexSelectionCastImageFilterType;
  typedef itk::ImageFileReader< PhilipsImageType >
    ReaderType;
  typedef itk::ImageFileWriter< ImageType >
    WriterType;
  typedef itk::ThresholdImageFilter< PhilipsImageType >
    ThresholdImageFilterType;
  typedef itk::ShiftScaleImageFilter< PhilipsImageType, ImageType4D >
    ShiftScaleImageFilterType;
  typedef itk::ShiftScaleInPlaceImageFilter< ImageType4D >
    ShiftScaleInPlaceImageFilterType;
  typedef itk::ExtractImageFilter< ImageType4D, ImageType >
    ExtractImageFilterType;
  typedef itk::VectorContainer< unsigned int, ExtractImageFilterType::Pointer >
    ExtractImageFilterContainerType;

  bool r2Mapping = false;
  int algorithm = MRT2ParameterMap3DImageFilterType::LINEAR;
  double maxT2Time = 10.0f;
  short threshold = 0;

  if( argc > 6 )
    {
    if( atoi( argv[6] ) )
      {
      r2Mapping = true;
      }
    }

  if( argc > 7 )
    {
    algorithm = atoi( argv[7] );
    }

  if( argc > 8 )
    {
    maxT2Time = atof( argv[8] );
    }

  if( argc > 9 )
    {
    threshold = (short)atoi( argv[9] );
    }

  int dims[4] = {0};
  // Create Philips PAR/REC reader and check the file if it can be read.
  PhilipsRECImageIOType::Pointer imageIO = PhilipsRECImageIOType::New();
  if( !imageIO->CanReadFile(inputFilename) )
    {
    std::cerr << "Could not read Philips PAR/REC file" << std::endl;
    return 1;
    }

  // Read the image information.
  imageIO->SetFileName(inputFilename);
  try
    {
    imageIO->ReadImageInformation();
    }
  catch( itk::ExceptionObject &err )
    {
    std::cerr << "ExceptionObject caught";
    std::cerr << " : "  << err.GetDescription();
    return 1;
    }

  unsigned int numberOfEchoTimes = 0;
  int tempNum = 0;
  // Make sure that the number of images
  // stored in the REC file matches the number of echoes.
  if( !itk::ExposeMetaData<int>(imageIO->GetMetaDataDictionary(),
    itk::PAR_MaxNumberOfEchoes,tempNum) )
    {
    std::cerr << "Could not determine the number of echoes" << std::endl;
    return 1;
    }
  numberOfEchoTimes = static_cast<unsigned int>(tempNum);

  // Get the image dimensions and make sure that
  // there exists at least numberOfEchoTimes image
  // volumes.  It's possible to have more if the 
  // REC file contains more than one image type
  // (i.e. magnitude, phase, real, imaginary, etc.)
  dims[0] = imageIO->GetDimensions(0);
  dims[1] = imageIO->GetDimensions(1);
  dims[2] = imageIO->GetDimensions(2);
  dims[3] = imageIO->GetDimensions(3);
  if( numberOfEchoTimes > (unsigned int)dims[3] )
    {
    std::cerr << "The number of echoes is larger than the number of image ";
    std::cerr << "blocks" << std::endl;
    return 1;
    }

  // Check to see if the REC file has at least the magnitude image.
  int haveMagnitude = 0;
  int numberOfImageTypes = 0;
  PhilipsRECImageIOType::ImageTypesType imageTypes;
  if( !itk::ExposeMetaData<int>(imageIO->GetMetaDataDictionary(),
    itk::PAR_NumberOfImageTypes,numberOfImageTypes) )
    {
    std::cerr << "Could not determine the number of image types" << std::endl;
    return 1;
    }
  if( !itk::ExposeMetaData<PhilipsRECImageIOType::ImageTypesType>
    (imageIO->GetMetaDataDictionary(), itk::PAR_ImageTypes,imageTypes) )
    {
    std::cerr << "Could not get the image types vector" << std::endl;
    return 1;
    }
  for(int j=0; j<numberOfImageTypes; j++)
    {
    if( imageTypes[j] == 0 )
      {
      haveMagnitude = 1;
      }
    }
  if( !haveMagnitude )
    {
    std::cerr << "Magnitude image type not found in REC file" << std::endl;
    return 1;
    }

  // Check to make sure that there is only one scanning sequence, 
  // otherwise the T2 map cannot be processed.
  int numberOfScanningSequences = 0;
  if( !itk::ExposeMetaData<int>(imageIO->GetMetaDataDictionary(),
    itk::PAR_NumberOfScanningSequences,numberOfScanningSequences) )
    {
    std::cerr << "Could not determine the number of scanning sequences"; 
    std::cerr << std::endl;
    return 1;
    }
  if( numberOfScanningSequences > 1 )
    {
    std::cerr << "Cannot process a T2 map when the number of scanning";
    std::cerr << " sequences is greater than 1" << std::endl;
    return 1;
    }

  // Get rescale values for converting the 16 bit image to floating point.
  PhilipsRECImageIOType::ScanningSequenceImageTypeRescaleValuesContainerType
    ::Pointer scanSequenceImageTypeRescaleValues = NULL;
  if( !itk::ExposeMetaData<PhilipsRECImageIOType
    ::ScanningSequenceImageTypeRescaleValuesContainerType::Pointer>
    (imageIO->GetMetaDataDictionary(),
    itk::PAR_ScanningSequenceImageTypeRescaleValues,
    scanSequenceImageTypeRescaleValues) )
    {
    std::cerr << "Could not get the rescale values for each scanning";
    std::cerr << " sequence and image type" << std::endl;
    return 1;
    }
  if( !scanSequenceImageTypeRescaleValues )
    {
    std::cerr << "Received NULL scanning sequence/image types vector"; 
    std::cerr << " pointer from meta dictionary" << std::endl;
    return 1;
    }
  PhilipsRECImageIOType::ImageTypeRescaleValuesContainerType::Pointer 
    rescaleValueVector = 
    scanSequenceImageTypeRescaleValues->ElementAt(0); //Only 1 scanning sequence.
  if( !rescaleValueVector )
    {
    std::cerr << "Received NULL rescale values vector pointer from"; 
    std::cerr << " meta dictionary" << std::endl;
    return 1;
    }

  // Read the image.
  ReaderType::Pointer baselineReader = ReaderType::New();
  baselineReader->SetFileName(inputFilename);
  baselineReader->SetImageIO(imageIO);
  try
    {
    baselineReader->UpdateLargestPossibleRegion();
    }
  catch( itk::ExceptionObject& err )
    {
    std::cerr << "ExceptionObject caught";
    std::cerr << " : "  << err.GetDescription();
    return 1;
    }

  // Threshold the image.
  ThresholdImageFilterType::Pointer t2Mask = ThresholdImageFilterType::New();
  t2Mask->SetOutsideValue(0);
  t2Mask->SetInput(baselineReader->GetOutput());
  t2Mask->ThresholdBelow(threshold);

  // Change image to floating point value.
  ShiftScaleImageFilterType::Pointer scaleOnly = NULL;
  ShiftScaleInPlaceImageFilterType::Pointer shiftAndScale = NULL;
  PhilipsRECImageIOType::ImageTypeRescaleValuesType rescaleValues = 
    rescaleValueVector->ElementAt(0); // Magnitude image is the first element.
  if( (rescaleValues[2] != 0) && // scale slope (SS)
    (rescaleValues[1] != 0) ) // rescale slope (RS)
    {
    scaleOnly = ShiftScaleImageFilterType::New();
    scaleOnly->SetInput(t2Mask->GetOutput());
    scaleOnly->SetScale(rescaleValues[1]); //RS
    shiftAndScale = ShiftScaleInPlaceImageFilterType::New();
    shiftAndScale->SetInput(scaleOnly->GetOutput());
    shiftAndScale->SetShift(rescaleValues[0]); //rescale intercept (RI)
    shiftAndScale->SetScale(1.0/(rescaleValues[2]*rescaleValues[1])); //1/(SS*RS)
    }
  else
    {
    std::cerr << "Invalid rescale values" << std::endl;
    return 1;
    }

  // Create T2 mapping class.
  MRT2ParameterMap3DImageFilterType::Pointer t2Map = 
    MRT2ParameterMap3DImageFilterType::New();
  // Select the fit type.
  switch(algorithm)
    {
    case MRT2ParameterMap3DImageFilterType::LINEAR:
    t2Map->SetAlgorithm(MRT2ParameterMap3DImageFilterType::LINEAR);
      break;
    case MRT2ParameterMap3DImageFilterType::NON_LINEAR:
    t2Map->SetAlgorithm(MRT2ParameterMap3DImageFilterType::NON_LINEAR);
      break;
    case MRT2ParameterMap3DImageFilterType::NON_LINEAR_WITH_CONSTANT:
    t2Map->SetAlgorithm(
      MRT2ParameterMap3DImageFilterType::NON_LINEAR_WITH_CONSTANT);
      break;
    default:
    std::cerr << "In valid algorithm = " << algorithm << std::endl;
      return 1;
    }
  t2Map->SetMaxT2Time(maxT2Time);
  if( r2Mapping )
    {
    t2Map->PerformR2MappingOn();
    }
#ifdef NO_MULTI_THREADING
  t2Map->SetNumberOfThreads(1);
#endif

  // Extract the volumes and echoes and set to T2 map filter.
  // Typically a multi-echo spin-echo image is used to meaure T2.
  // A spin-echo image is acquired using a 90 degree rf-pulse
  // followed by a 180 degree rf-pulse.  To acquire a multi-echo 
  // sequence a series of 180 degree rf-pulses follow the first 90 
  // and 180 degree pulses.  The application of a perfectly homogenous 
  // 180 rf-pulse is difficult and imperfections in the flip-angle 
  // lead to what are called stimulated echoes for the echo images
  // acquired after the first 180 degree rf-pulse.  Therefore in order
  // to measure T2 the first echo image without the stimulated echo
  // signal is usually thrown out.  The following code below will
  // extract the echo image volumes from the 4D image.  It will also
  // get the echo times from the imageIO class and convert the times
  // to seconds.
  ExtractImageFilterContainerType::Pointer extractVOI = 
    ExtractImageFilterContainerType::New();
  extractVOI->resize(numberOfEchoTimes-1);
  ExtractImageFilterType::InputImageRegionType extractionRegion;
  ExtractImageFilterType::InputImageSizeType extractionSize;
  extractionSize[0] = dims[0];
  extractionSize[1] = dims[1];
  extractionSize[2] = dims[2];
  extractionSize[3] = 0;
  ExtractImageFilterType::InputImageIndexType extractionIndex;
  extractionIndex[0] = 0;
  extractionIndex[1] = 0;
  extractionIndex[2] = 0;
  extractionRegion.SetSize(extractionSize);
  PhilipsRECImageIOType::EchoTimesContainerType::Pointer ptrToEchoes = NULL;
  if(!itk::ExposeMetaData<PhilipsRECImageIOType::EchoTimesContainerType::Pointer>
    (imageIO->GetMetaDataDictionary(), itk::PAR_EchoTimes,ptrToEchoes) )
    {
    std::cerr << "Could not get the echo times" << std::endl;
    return 1;
    }
  if( !ptrToEchoes )
    {
    std::cerr << "Received NULL echo times pointer from";
    std::cerr << " meta dictionary" << std::endl;
    return 1;
    }
  if( ptrToEchoes->size() != numberOfEchoTimes )
    {
    std::cerr << "The size of the echo times vector does";
    std::cerr << " not match the number of echoes listed in";
    std::cerr << " the PAR file" << std::endl;
    return 1;
    }
  for(unsigned int i=1; i<numberOfEchoTimes; i++)
    {
    extractVOI->SetElement(i-1, ExtractImageFilterType::New());
    extractVOI->ElementAt(i-1)->SetInput(shiftAndScale->GetOutput());
    extractionIndex[3] = i;
    extractionRegion.SetIndex(extractionIndex);
    extractVOI->ElementAt(i-1)->SetExtractionRegion(extractionRegion);
    t2Map->AddMREchoImage(ptrToEchoes->ElementAt(i)/1000.0f, //convert to seconds
      extractVOI->ElementAt(i-1)->GetOutput());
    }

  // Extract each output component and write to disk.
  VectorIndexSelectionCastImageFilterType::Pointer extractComp = 
    VectorIndexSelectionCastImageFilterType::New();
  extractComp->SetInput(t2Map->GetOutput());
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(extractComp->GetOutput());

  // T2/R2 map.
  writer->SetFileName( outputT2Filename );
  extractComp->SetIndex(0);
  try
    {
    writer->Update();
    }
  catch(...)
    {
    std::cerr << "Error during write of " << outputT2Filename << std::endl;
    return 1;
    }
  // Exponent Constant map.
  extractComp->SetIndex(1);
  writer->SetFileName( outputExpConstFilename );
  try
    {
    writer->Update();
    }
  catch(...)
    {
    std::cerr << "Error during write of " << outputExpConstFilename << std::endl;
    return 1;
    }
  // Constant map.
  extractComp->SetIndex(2);
  writer->SetFileName( outputConstFilename );
  try
    {
    writer->Update();
    }
  catch(...)
    {
    std::cerr << "Error during write of " << outputConstFilename << std::endl;
    return 1;
    }
  // Rsquared map.
  extractComp->SetIndex(3);
  writer->SetFileName( outputRSquaredFilename );
  try
    {
    writer->Update();
    }
  catch(...)
    {
    std::cerr << "Error during write of " << outputRSquaredFilename << std::endl;
    return 1;
    }

  return 0;
}
