/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: BrukerT1Map.cxx,v $
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
#include "itkBruker2DSEQImageIO.h"
#include "itkMRT1ParameterMap3DImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkVectorImage.h"
#include "itkVariableLengthVector.h"
#include "itkMetaDataObject.h"
#include "vnl/vnl_vector_fixed.h"

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

  const unsigned int Dimension = 3;

  // Bruker 2dseq images are usually signed 16 bit or
  // 32 bit integers.  Just use 32 bit floating point.
  typedef itk::Image< float, Dimension >  ImageType;

  typedef itk::Bruker2DSEQImageIO
    Bruker2DSEQImageIOType;
  typedef itk::MRT1ParameterMap3DImageFilter<ImageType::PixelType>
    MRT1ParameterMap3DImageFilterType;
  typedef itk::VectorImage< float, Dimension >
    VectorImageType;
  typedef itk::VectorIndexSelectionCastImageFilter<
    MRT1ParameterMap3DImageFilterType::OutputImageType, ImageType>
    VectorIndexSelectionCastImageFilterType;
  typedef itk::ThresholdImageFilter< ImageType >
    ThresholdImageFilterType;
  typedef itk::ImageFileReader< ImageType >
    ReaderType;
  typedef itk::ImageFileWriter< ImageType >
    WriterType;

  bool r1Mapping = false;
  int algorithm = MRT1ParameterMap3DImageFilterType::IDEAL_STEADY_STATE;
  double maxT1Time = 10.0f;
  float threshold = 0.0f;

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

  if( argc > 8 )
    {
    maxT1Time = atof( argv[8] );
    }

  if( argc > 9 )
    {
    threshold = atof( argv[9] );
    }

  int dims[3] = {0};
  // Create 2DSEQ reader and check the file if it can be read.
  Bruker2DSEQImageIOType::Pointer imageIO = Bruker2DSEQImageIOType::New();
  if( !imageIO->CanReadFile(inputFilename) )
    {
    std::cerr << "Could not read 2dseq file" << std::endl;
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

  // Get the image dimensions
  dims[0] = imageIO->GetDimensions(0);
  dims[1] = imageIO->GetDimensions(1);
  dims[2] = imageIO->GetDimensions(2);

  // Extract repetition/inversion times depending on the algorithm.
  Bruker2DSEQImageIOType::ACQRepetitionTimeContainerType::Pointer ptrToTimePoints
    = NULL;
  if((algorithm == MRT1ParameterMap3DImageFilterType::IDEAL_STEADY_STATE) ||
    (algorithm == MRT1ParameterMap3DImageFilterType::HYBRID_STEADY_STATE_3PARAM))
    {
    // Repetition times for the saturation recovery fit type.
    if(!itk::ExposeMetaData<
      Bruker2DSEQImageIOType::ACQRepetitionTimeContainerType::Pointer>
      (imageIO->GetMetaDataDictionary(),itk::ACQ_REPETITION_TIME,
      ptrToTimePoints))
      {
      std::cerr << "Could not get the repetition times" << std::endl;
      return 1;
      }
    }
  else
    {
    // Inversion times for all others.
    if(!itk::ExposeMetaData<
      Bruker2DSEQImageIOType::ACQRepetitionTimeContainerType::Pointer>
      (imageIO->GetMetaDataDictionary(),itk::ACQ_INVERSION_TIME,ptrToTimePoints))
      {
      std::cerr << "Could not get the inversion times" << std::endl;
      return 1;
      }
    }
  if( !ptrToTimePoints )
    {
    std::cerr << "Received NULL repetition/inversion times pointer";
    std::cerr << " from meta dictionary" << std::endl;
    return 1;
    }
  unsigned int numberOfTimePoints = ptrToTimePoints->Size();
  if( numberOfTimePoints < 2 )
    {
    std::cerr << "Must have at least 2 images to calculate T1" << std::endl;
    return 1;
    }

  // Convert the times to seconds.
  for(unsigned int i=0; i<(unsigned int)numberOfTimePoints; i++)
    {
    // convert to seconds
    ptrToTimePoints->SetElement(i,ptrToTimePoints->ElementAt(i)/1000.0f); 
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

  // Mask the images if the mask is set.
  ThresholdImageFilterType::Pointer t1Mask = ThresholdImageFilterType::New();
  t1Mask->SetOutsideValue(0);
  t1Mask->SetInput(baselineReader->GetOutput());
  t1Mask->ThresholdBelow(threshold);
  try
    {
    t1Mask->UpdateLargestPossibleRegion();
    }
  catch( itk::ExceptionObject &err )
    {
    std::cerr << "ExceptionObject caught";
    std::cerr << " : "  << err.GetDescription();
    return 1;
    }

  // Get real number of slices and create vector image.
  // Also threshold the image at the same time.
  int realSlices = dims[2]/numberOfTimePoints;
  VectorImageType::RegionType region;
  VectorImageType::SizeType size;
  size[0] = dims[0];
  size[1] = dims[1];
  size[2] = realSlices;
  VectorImageType::IndexType index;
  index[0] = 0;
  index[1] = 0;
  index[2] = 0;
  region.SetSize(size);
  region.SetIndex(index);
  VectorImageType::Pointer vectorImage = VectorImageType::New();
  vectorImage->SetRegions(region);
  vectorImage->SetVectorLength(numberOfTimePoints);
  VectorImageType::PointType origin = baselineReader->GetOutput()->GetOrigin();
  origin[2] = -baselineReader->GetOutput()->GetSpacing()[2]*size[2]/2.0f;
  vectorImage->SetOrigin(origin);
  vectorImage->SetSpacing(baselineReader->GetOutput()->GetSpacing());
  vectorImage->SetDirection(baselineReader->GetOutput()->GetDirection());
  vectorImage->Allocate();
  ImageType::IndexType imageTimePointIndex;
  imageTimePointIndex.Fill(0);
  for(index[0]=0,imageTimePointIndex[0]=0;index[0]<(int)size[0]; 
    index[0]++,imageTimePointIndex[0]++)
    {
    for(index[1]=0,imageTimePointIndex[1]=0;index[1]<(int)size[1]; 
      index[1]++,imageTimePointIndex[1]++)
      {
      for(index[2]=0;index[2]<(int)size[2];index[2]++ )
      {
      // Multi-echo inversion recovery Bruker images are stored as they are 
      // acquired.  This means that the images are not stored as volumes.  We 
      // need to put each echo from each slice into the vector as follows:
      if((algorithm != MRT1ParameterMap3DImageFilterType::IDEAL_STEADY_STATE) &&
        (algorithm != 
          MRT1ParameterMap3DImageFilterType::HYBRID_STEADY_STATE_3PARAM))
        {
        // Skip to next slice
        imageTimePointIndex[2] = index[2]*numberOfTimePoints; 
        }
      VectorImageType::PixelType timePointVector(numberOfTimePoints);
      for(unsigned int timePoint=0; timePoint<numberOfTimePoints; timePoint++)
        {
        // Images are stored as volumes in the saturation recovery case, so we 
        // need to skip to each new repetition while filling the vector.
        if((algorithm == MRT1ParameterMap3DImageFilterType::IDEAL_STEADY_STATE) 
          || (algorithm == 
            MRT1ParameterMap3DImageFilterType::HYBRID_STEADY_STATE_3PARAM) )
          {
          imageTimePointIndex[2] = index[2] + (realSlices*timePoint);
          }
          ImageType::PixelType pixelVal = 
            baselineReader->GetOutput()->GetPixel(imageTimePointIndex);
          ImageType::PixelType maskVal = 0;
          if((algorithm == MRT1ParameterMap3DImageFilterType::IDEAL_STEADY_STATE)
            || (algorithm == 
              MRT1ParameterMap3DImageFilterType::HYBRID_STEADY_STATE_3PARAM))
            {
            maskVal = t1Mask->GetOutput()->GetPixel(index);
            }
          else
            {
            VectorImageType::IndexType tempIndex = index;
            tempIndex[2] = tempIndex[2]*numberOfTimePoints; // Skip to next slice
            maskVal = t1Mask->GetOutput()->GetPixel(tempIndex);
            }
          timePointVector[timePoint] = (maskVal==0)?0:pixelVal;
          if((algorithm != MRT1ParameterMap3DImageFilterType::IDEAL_STEADY_STATE)
            && (algorithm != 
              MRT1ParameterMap3DImageFilterType::HYBRID_STEADY_STATE_3PARAM) )
            {
            ++imageTimePointIndex[2];
            }
          }
        vectorImage->SetPixel(index, timePointVector);
        }
      }
    }
  // Not needed anymore.
  t1Mask = NULL;
  baselineReader = NULL;

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
  t1Map->SetMRImage(ptrToTimePoints, vectorImage);
  if( r1Mapping )
    {
    t1Map->PerformR1MappingOn();
    }
#ifdef NO_MULTI_THREADING
  t1Map->SetNumberOfThreads(1);
#endif

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
