/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: PhilipsT1Map.cxx,v $
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
#include "itkMRT1ParameterMap3DImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkShiftScaleImageFilter.h"
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
    std::cerr << argv[0] << " inputImage outputT1Image outputExpConstImage";
    std::cerr << " outputConstImage outputRSquaredImage [T1map=0,R1map=1]";
    std::cerr << " [algorithm] [maxT1Time] [threshold]" << std::endl;
    return 1;
    }

  const char * inputFilename   = argv[1];
  const char * outputT1Filename = argv[2];
  const char * outputExpConstFilename = argv[3];
  const char * outputConstFilename = argv[4];
  const char * outputRSquaredFilename = argv[5];

  // Philips PAR/REC images are signed 16 bit and the
  // reader supports 4D images.
  typedef itk::Image< short, 4 > PhilipsImageType4D;
  typedef itk::Image< short, 3 > PhilipsImageType;
  typedef itk::Image< float, 4 > ImageType4D;
  typedef itk::Image< float, 3 > ImageType;

  typedef itk::PhilipsRECImageIO
    PhilipsRECImageIOType;
  typedef itk::MRT1ParameterMap3DImageFilter< ImageType::PixelType >
    MRT1ParameterMap3DImageFilterType;
  typedef itk::VectorIndexSelectionCastImageFilter<
    MRT1ParameterMap3DImageFilterType::OutputImageType, ImageType >
    VectorIndexSelectionCastImageFilterType;
  typedef itk::ImageFileReader< PhilipsImageType4D >
    ReaderType;
  typedef itk::ImageFileWriter< ImageType >
    WriterType;
  typedef itk::ThresholdImageFilter< PhilipsImageType >
    ThresholdImageFilterType;
  typedef itk::ShiftScaleImageFilter< PhilipsImageType4D, ImageType4D >
    ShiftScaleImageFilterType;
  typedef itk::ShiftScaleImageFilter< ImageType4D, ImageType4D >
    ShiftScaleInPlaceImageFilterType;
  typedef itk::ExtractImageFilter< ImageType4D, ImageType >
    ExtractImageFilterType;
  typedef itk::VectorContainer< unsigned int, ExtractImageFilterType::Pointer >
    ExtractImageFilterContainerType;
  typedef itk::MaskImageFilter< ImageType, PhilipsImageType ,ImageType >
    MaskImageFilterType;
  typedef itk::VectorContainer< unsigned int, MaskImageFilterType::Pointer >
    MaskImageFilterContainerType;
  typedef itk::ExtractImageFilter< PhilipsImageType4D, PhilipsImageType >
    ExtractPhilipsImageFilterType;

  bool r1Mapping = false;
  int algorithm = MRT1ParameterMap3DImageFilterType::INVERSION_RECOVERY;
  double maxT1Time = 10.0f;
  short threshold = 0;

  if( argc > 6 )
    {
    if( atoi( argv[6] ) )
      {
      r1Mapping = true;
      }
    }

  if( argc > 7 )
    {
    algorithm = atoi( argv[7] );
    }

  if((algorithm == MRT1ParameterMap3DImageFilterType::IDEAL_STEADY_STATE) ||
    (algorithm == MRT1ParameterMap3DImageFilterType::HYBRID_STEADY_STATE_3PARAM))
    {
    std::cerr << "The IDEAL_STEADY_STATE and HYBRID_STEADY_STATE_3PARAM";
    std::cerr << " algorithms are not supported" << std::endl;
    return 1;
    }

  if( argc > 8 )
    {
    maxT1Time = atof( argv[8] );
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

  unsigned int numberOfTimePoints = 0;
  int tempNum = 0;
  // Find the maximum number of cardiac phases.
  if( !itk::ExposeMetaData<int>(imageIO->GetMetaDataDictionary(),
    itk::PAR_MaxNumberOfCardiacPhases,tempNum) )
    {
    std::cerr << "Could not determine the number of cardiac";
    std::cerr << " phases" << std::endl;
    return 1;
    }
  numberOfTimePoints = static_cast<unsigned int>(tempNum);
  if( numberOfTimePoints < 2 )
    {
    std::cerr << "Must have at least 2 cardiac phase images to";
    std::cerr << " calculate T1" << std::endl;
    return 1;
    }

  // Get the image dimensions and make sure that
  // there exists at least numberOfTimePoints image
  // volumes.  It's possible to have more if the
  // REC file contains more than one image type
  // (i.e. magnitude, phase, real, imaginary, etc.)
  dims[0] = imageIO->GetDimensions(0);
  dims[1] = imageIO->GetDimensions(1);
  dims[2] = imageIO->GetDimensions(2);
  dims[3] = imageIO->GetDimensions(3);
  if( numberOfTimePoints > (unsigned int)dims[3] )
    {
    std::cerr << "The number of time points is larger than the number";
    std::cerr << " of image blocks" << std::endl;
    return 1;
    }

  // Check to see if we have the correct image types.
  int haveCorrectImage = 0;
  int numberOfImageTypes = 0;
  PhilipsRECImageIOType::ImageTypesType imageTypes;
  if( !itk::ExposeMetaData<int>(imageIO->GetMetaDataDictionary(),
    itk::PAR_NumberOfImageTypes,numberOfImageTypes) )
    {
    std::cerr << "Could not determine the number of image types" << std::endl;
    return 1;
    }
  if( !itk::ExposeMetaData<PhilipsRECImageIOType::ImageTypesType>
    (imageIO->GetMetaDataDictionary(),itk::PAR_ImageTypes,imageTypes) )
    {
    std::cerr << "Could not get the image types vector" << std::endl;
    return 1;
    }
  // Need image 4 (corrected real image) for inversion recovery/Look-Locker
  // fit type.
  if( (algorithm == MRT1ParameterMap3DImageFilterType::INVERSION_RECOVERY) ||
    (algorithm ==
      MRT1ParameterMap3DImageFilterType::INVERSION_RECOVERY_3PARAM) ||
    (algorithm == MRT1ParameterMap3DImageFilterType::LOOK_LOCKER) )
    {
    for(int j=0; j<numberOfImageTypes; j++)
      {
      if( imageTypes[j] == 4 )
        {
        haveCorrectImage = 1;
        }
      }
    if( !haveCorrectImage )
      {
      std::cerr << "No corrected real image type detected in PAR";
      std::cerr << " file" << std::endl;
      return 1;
      }
    }
  // Now check for existence of the magnitude image.
  haveCorrectImage = 0;
  for(int j=0; j<numberOfImageTypes; j++)
    {
    if( imageTypes[j] == 0 )
      {
      haveCorrectImage = 1;
      }
    }
  if( !haveCorrectImage )
    {
    std::cerr << "No magnitude image type detected in PAR file" << std::endl;
    return 1;
    }

  // Get rescale values.
  PhilipsRECImageIOType::ScanningSequenceImageTypeRescaleValuesContainerType
    ::Pointer scanSequenceImageTypeRescaleValues = NULL;
  if( !itk::ExposeMetaData<PhilipsRECImageIOType::
    ScanningSequenceImageTypeRescaleValuesContainerType::Pointer>
    (imageIO->GetMetaDataDictionary(),
    itk::PAR_ScanningSequenceImageTypeRescaleValues,
    scanSequenceImageTypeRescaleValues))
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
  // All of the images we need should be located in the first scanning
  // sequence or scanning sequence 0.
  PhilipsRECImageIOType::ImageTypeRescaleValuesContainerType::Pointer
    rescaleValueVector = scanSequenceImageTypeRescaleValues->ElementAt(0);
  if( !rescaleValueVector )
    {
    std::cerr << "Received NULL rescale values vector pointer from";
    std::cerr << " meta dictionary" << std::endl;
    return 1;
    }

  // Get trigger times.
  PhilipsRECImageIOType::TriggerTimesContainerType::Pointer ptrToTimePoints =
    NULL;
  if(!itk::ExposeMetaData<
    PhilipsRECImageIOType::TriggerTimesContainerType::Pointer>
    (imageIO->GetMetaDataDictionary(),itk::PAR_TriggerTimes,ptrToTimePoints) )
    {
    std::cerr << "Could not get the trigger times" << std::endl;
    return 1;
    }
  if( !ptrToTimePoints )
    {
    std::cerr << "Received NULL trigger times pointer from meta";
    std::cerr << " dictionary" << std::endl;
    return 1;
    }
  if( ptrToTimePoints->size() != numberOfTimePoints )
    {
    std::cerr << "The size of the time points vector does not match the";
    std::cerr << " number of cardiac phases listed in the PAR file" << std::endl;
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

  // Change image to floating point value.
  ShiftScaleImageFilterType::Pointer scaleOnly = NULL;
  ShiftScaleInPlaceImageFilterType::Pointer shiftAndScale = NULL;
  PhilipsRECImageIOType::ImageTypeRescaleValuesType rescaleValues;
  if((algorithm == MRT1ParameterMap3DImageFilterType::INVERSION_RECOVERY) ||
    (algorithm == MRT1ParameterMap3DImageFilterType::INVERSION_RECOVERY_3PARAM)
    || (algorithm == MRT1ParameterMap3DImageFilterType::LOOK_LOCKER) )
    {
    // Corrected image should be at location 4, based
    // on the image type, which is 4.
    rescaleValues = rescaleValueVector->ElementAt(4);
    }
  else
    {
    // Magnitude image will be the first image volume.
    rescaleValues = rescaleValueVector->ElementAt(0);
    }
  if( (rescaleValues[2] != 0) && //scale slope (SS)
    (rescaleValues[1] != 0) ) //rescale slope (RS)
    {
    scaleOnly = ShiftScaleImageFilterType::New();
    scaleOnly->SetInput(baselineReader->GetOutput());
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

  // Create T1 mapping class.
  MRT1ParameterMap3DImageFilterType::Pointer t1Map =
    MRT1ParameterMap3DImageFilterType::New();
  // Select the fit type.
  switch(algorithm)
    {
    case MRT1ParameterMap3DImageFilterType::IDEAL_STEADY_STATE:
    t1Map->SetAlgorithm(MRT1ParameterMap3DImageFilterType::IDEAL_STEADY_STATE);
      break;
    case MRT1ParameterMap3DImageFilterType::INVERSION_RECOVERY:
    t1Map->SetAlgorithm(MRT1ParameterMap3DImageFilterType::INVERSION_RECOVERY);
      break;
    case MRT1ParameterMap3DImageFilterType::ABSOLUTE_INVERSION_RECOVERY:
    t1Map->SetAlgorithm(
      MRT1ParameterMap3DImageFilterType::ABSOLUTE_INVERSION_RECOVERY);
      break;
    case MRT1ParameterMap3DImageFilterType::LOOK_LOCKER:
    t1Map->SetAlgorithm(MRT1ParameterMap3DImageFilterType::LOOK_LOCKER);
      break;
    case MRT1ParameterMap3DImageFilterType::ABSOLUTE_LOOK_LOCKER:
    t1Map->SetAlgorithm(MRT1ParameterMap3DImageFilterType::ABSOLUTE_LOOK_LOCKER);
      break;
    case MRT1ParameterMap3DImageFilterType::HYBRID_STEADY_STATE_3PARAM:
    t1Map->SetAlgorithm(
      MRT1ParameterMap3DImageFilterType::HYBRID_STEADY_STATE_3PARAM);
      break;
    case MRT1ParameterMap3DImageFilterType::INVERSION_RECOVERY_3PARAM:
    t1Map->SetAlgorithm(
      MRT1ParameterMap3DImageFilterType::INVERSION_RECOVERY_3PARAM);
      break;
    case MRT1ParameterMap3DImageFilterType::ABSOLUTE_INVERSION_RECOVERY_3PARAM:
    t1Map->SetAlgorithm(
      MRT1ParameterMap3DImageFilterType::ABSOLUTE_INVERSION_RECOVERY_3PARAM);
      break;
    default:
    std::cerr << "In valid algorithm = " << algorithm << std::endl;
      return 1;
    }
  t1Map->SetMaxT1Time(maxT1Time);
  if( r1Mapping )
    {
    t1Map->PerformR1MappingOn();
    }
#ifdef NO_MULTI_THREADING
  t1Map->SetNumberOfThreads(1);
#endif

  // Extract the volumes and set to T1 map filter.
  // Convert time values to seconds.
  ExtractImageFilterContainerType::Pointer extractVOI =
    ExtractImageFilterContainerType::New();
  extractVOI->resize(numberOfTimePoints);
  MaskImageFilterContainerType::Pointer maskFilterContainer
    = MaskImageFilterContainerType::New();
  maskFilterContainer->resize(numberOfTimePoints);
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
  extractionIndex[3] = 0;
  extractionRegion.SetSize(extractionSize);
  // Generate a mask using the magnitude image.
  ExtractPhilipsImageFilterType::Pointer magnitudeImage =
    ExtractPhilipsImageFilterType::New();
  ThresholdImageFilterType::Pointer magnitudeMask
    = ThresholdImageFilterType::New();
  magnitudeImage->SetInput(baselineReader->GetOutput());
  extractionRegion.SetIndex(extractionIndex);
  magnitudeImage->SetExtractionRegion(extractionRegion);
  magnitudeMask->SetInput(magnitudeImage->GetOutput());
  magnitudeMask->SetOutsideValue(0);
  magnitudeMask->ThresholdBelow(threshold);
  try
    {
    magnitudeMask->UpdateLargestPossibleRegion();
    }
  catch( itk::ExceptionObject &err )
    {
    std::cerr << "ExceptionObject caught";
    std::cerr << " : "  << err.GetDescription();
    return 1;
    }

  // Extract volumes according to algorithm type.
  if( (algorithm == MRT1ParameterMap3DImageFilterType::INVERSION_RECOVERY) ||
    (algorithm == MRT1ParameterMap3DImageFilterType::INVERSION_RECOVERY_3PARAM)
    || (algorithm == MRT1ParameterMap3DImageFilterType::LOOK_LOCKER) )
    {
    // Corrected real image is after the magnitude and real images in first
    // scanning sequence.
    for(unsigned int i=0; i<numberOfTimePoints; i++)
      {
      extractVOI->SetElement(i, ExtractImageFilterType::New());
      extractVOI->ElementAt(i)->SetInput(shiftAndScale->GetOutput());
      extractionIndex[3] = 2*numberOfTimePoints + i;
      extractionRegion.SetIndex(extractionIndex);
      extractVOI->ElementAt(i)->SetExtractionRegion(extractionRegion);
      maskFilterContainer->SetElement(i,MaskImageFilterType::New());
      maskFilterContainer->ElementAt(i)->SetInput1(
        extractVOI->ElementAt(i)->GetOutput());
      maskFilterContainer->ElementAt(i)->SetInput2(magnitudeMask->GetOutput());
      // convert to seconds
      t1Map->AddMRImage(ptrToTimePoints->ElementAt(i)/1000.0f,
        maskFilterContainer->ElementAt(i)->GetOutput());
      }
    }
  else
    {
    // Magnitude image is at front.
    for(unsigned int i=0; i<numberOfTimePoints; i++)
      {
      extractVOI->SetElement(i, ExtractImageFilterType::New());
      extractVOI->ElementAt(i)->SetInput(shiftAndScale->GetOutput());
      maskFilterContainer->SetElement(i,MaskImageFilterType::New());
      maskFilterContainer->ElementAt(i)->SetInput1(
        extractVOI->ElementAt(i)->GetOutput());
      maskFilterContainer->ElementAt(i)->SetInput2(magnitudeMask->GetOutput());
      extractionIndex[3] = i;
      extractionRegion.SetIndex(extractionIndex);
      extractVOI->ElementAt(i)->SetExtractionRegion(extractionRegion);
      // convert to seconds
      t1Map->AddMRImage(ptrToTimePoints->ElementAt(i)/1000.0f,
        maskFilterContainer->ElementAt(i)->GetOutput());
      }
    }

  // Extract each output component and write to disk.
  VectorIndexSelectionCastImageFilterType::Pointer extractComp =
    VectorIndexSelectionCastImageFilterType::New();
  extractComp->SetInput(t1Map->GetOutput());
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(extractComp->GetOutput());

  // T1/R1 map.
  writer->SetFileName( outputT1Filename );
  extractComp->SetIndex(0);
  try
    {
    writer->Update();
    }
  catch(...)
    {
    std::cerr << "Error during write of " << outputT1Filename << std::endl;
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
