/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: BrukerT2Map.cxx,v $
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
#include "itkMRT2ParameterMap3DImageFilter.h"
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

  const unsigned int Dimension = 3;

  // Bruker 2dseq images are usually signed 16 bit or
  // 32 bit integers.  Just use 32 bit floating point.
  typedef itk::Image< float, Dimension > ImageType;

  typedef itk::Bruker2DSEQImageIO
    Bruker2DSEQImageIOType;
  typedef itk::MRT2ParameterMap3DImageFilter<ImageType::PixelType>
    MRT2ParameterMap3DImageFilterType;
  typedef itk::VectorImage< float, Dimension >
    VectorImageType;
  typedef itk::VectorIndexSelectionCastImageFilter<
    MRT2ParameterMap3DImageFilterType::OutputImageType, ImageType>
    VectorIndexSelectionCastImageFilterType;
  typedef itk::ImageFileReader< ImageType >
    ReaderType;
  typedef itk::ImageFileWriter< ImageType >
    WriterType;

  bool r2Mapping = false;
  int algorithm = MRT2ParameterMap3DImageFilterType::LINEAR;
  double maxT2Time = 10.0f;
  float threshold = 0.0f;

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

  // Get the echo times in ms.
  Bruker2DSEQImageIOType::ACQEchoTimeContainerType::Pointer ptrToEchoes = NULL;
  if(!itk::ExposeMetaData<
    Bruker2DSEQImageIOType::ACQEchoTimeContainerType::Pointer>
   (imageIO->GetMetaDataDictionary(),itk::ACQ_ECHO_TIME,ptrToEchoes))
    {
    std::cerr << "Could not get the echo times" << std::endl;
    return 1;
    }
  if( !ptrToEchoes )
    {
    std::cerr << "Received NULL echo times pointer from meta dictionary"
      << std::endl;
    return 1;
    }
  unsigned int numberOfEchoTimes = ptrToEchoes->Size();

  // Typically a multi-echo spin-echo image is used to measure T2.
  // A spin-echo image is acquired using a 90 degree rf-pulse
  // followed by a 180 degree rf-pulse.  To acquire a multi-echo
  // sequence a series of 180 degree rf-pulses follow the first 90
  // and 180 degree pulses.  The application of a perfectly homogenous
  // 180 rf-pulse is difficult and imperfections in the flip-angle
  // lead to what are called stimulated echoes for the echo images
  // acquired after the first 180 degree rf-pulse.  Therefore in order
  // to measure T2 the first echo image without the stimulated echo
  // signal is usually thrown out.  The following for loop below will
  // do that as well as convert the echo times to seconds.
  for(unsigned int i=1; i<(unsigned int)numberOfEchoTimes; i++)
    {
    // convert to seconds
    ptrToEchoes->SetElement(i-1,ptrToEchoes->ElementAt(i)/1000.0f);
    }
  ptrToEchoes->CastToSTLContainer().pop_back();
  --numberOfEchoTimes;
  if( numberOfEchoTimes < 2 )
    {
    std::cerr << "Must have at least 2 echo images to calculate T2" << std::endl;
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

  // Get real number of slices and create vector image.
  // Also threshold the image at the same time.
  int realSlices = dims[2]/(numberOfEchoTimes+1);
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
  vectorImage->SetVectorLength(numberOfEchoTimes);
  VectorImageType::PointType origin = baselineReader->GetOutput()->GetOrigin();
  origin[2] = -baselineReader->GetOutput()->GetSpacing()[2]*size[2]/2.0f;
  vectorImage->SetOrigin(origin);
  vectorImage->SetSpacing(baselineReader->GetOutput()->GetSpacing());
  vectorImage->SetDirection(baselineReader->GetOutput()->GetDirection());
  vectorImage->Allocate();
  ImageType::IndexType echoIndex;
  for(index[0]=0,echoIndex[0]=0;index[0]<(int)size[0];index[0]++,echoIndex[0]++)
    {
    for(index[1]=0,echoIndex[1]=0;index[1]<(int)size[1];
      index[1]++,echoIndex[1]++)
      {
      for(index[2]=0;index[2]<(int)size[2];index[2]++)
        {
        // Multi-echo Bruker images are stored as they are acquired.  This
        // means that the images are not stored as volumes.  We need to put
        // each echo from each slice into the vector as follows:
        // Skip to next slice ignoring the first echo.
        echoIndex[2] = index[2]*(numberOfEchoTimes+1) + 1;
        VectorImageType::PixelType echoVector(numberOfEchoTimes);
        for(unsigned int echo=0; echo<numberOfEchoTimes; echo++)
          {
          ImageType::PixelType pixelVal
            = baselineReader->GetOutput()->GetPixel(echoIndex);
          echoVector[echo] = (pixelVal < threshold)?0:pixelVal;
          ++echoIndex[2];
          }
        vectorImage->SetPixel(index, echoVector);
        }
      }
    }
  baselineReader = NULL; // Not needed anymore.

  // Create T2 mapping class.
  MRT2ParameterMap3DImageFilterType::Pointer t2Map
    = MRT2ParameterMap3DImageFilterType::New();
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
  t2Map->SetMREchoImage(ptrToEchoes, vectorImage);
  if( r2Mapping )
    {
    t2Map->PerformR2MappingOn();
    }
#ifdef NO_MULTI_THREADING
  t2Map->SetNumberOfThreads(1);
#endif

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
