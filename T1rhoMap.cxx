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
#include "itkMetaDataObject.h"
#include "itkVectorContainer.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "vnl/vnl_vector_fixed.h"
#include <algorithm>

#include "itkFDFImageIOFactory.h"
#include "itkFDFImageIO.h"


#include "T1rhoMapCLP.h"
 

int main(int argc, char **argv)
{
  PARSE_ARGS;
  
  
  typedef itk::Image< float, 3 > ImageType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::ThresholdImageFilter< ImageType > ThresholdImageFilterType;
  typedef itk::MRT2ParameterMap3DImageFilter< ImageType::PixelType >
    MRT2ParameterMap3DImageFilterType;
  typedef itk::VectorIndexSelectionCastImageFilter< 
    MRT2ParameterMap3DImageFilterType::OutputImageType, ImageType >
    VectorIndexSelectionCastImageFilterType;
      
  // Create T2 mapping class.
  
  MRT2ParameterMap3DImageFilterType::Pointer t2Map = 
    MRT2ParameterMap3DImageFilterType::New();
    
  int numberOfImages = inputVolumes.size();
  
  for (int i=0;i<numberOfImages;i++)
  {
    // Read the image.
    ReaderType::Pointer imageReader = ReaderType::New();
    imageReader->SetFileName(inputVolumes[i]);
    try
      {
      imageReader->UpdateLargestPossibleRegion();
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
    t2Mask->SetInput( imageReader->GetOutput() );
    t2Mask->ThresholdBelow( threshold );
    t2Mask->Update( );
    
    t2Map->AddMREchoImage(T1rhoTimes[i], t2Mask->GetOutput());
  }
  
  
  if (mappingAlgorithm == "Linear")
    {
    t2Map->SetAlgorithm(MRT2ParameterMap3DImageFilterType::LINEAR);
    }
  else if (mappingAlgorithm == "Nonlinear") 
    {
    t2Map->SetAlgorithm(MRT2ParameterMap3DImageFilterType::NON_LINEAR);
    }
  else if (mappingAlgorithm == "Nonlinear With Constant") 
    {
    t2Map->SetAlgorithm(MRT2ParameterMap3DImageFilterType::NON_LINEAR_WITH_CONSTANT);
    }
  else
    {
    std::cerr << "In valid algorithm = " << mappingAlgorithm << std::endl;
    return 1;
    }
  t2Map->SetMaxT2Time( maxTime );
  
#ifdef NO_MULTI_THREADING
  t2Map->SetNumberOfThreads(1);
#endif

  
  // Extract each output component and write to disk.
  VectorIndexSelectionCastImageFilterType::Pointer extractComp = 
    VectorIndexSelectionCastImageFilterType::New();
  extractComp->SetInput(t2Map->GetOutput());
  WriterType::Pointer writer = WriterType::New();
  extractComp->SetIndex(0);
  writer->SetInput(extractComp->GetOutput());

  // T2/R2 map.
  if ( outputFilename.length() > 0)
  {
    writer->SetFileName( outputFilename );
    extractComp->SetIndex(0);
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & err )
    {
      std::cerr << "Error during write of " << outputFilename << std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return 1; 
    }
  }
  
  // Exponent Constant map.
  if ( outputExpConstFilename.length() > 0)
  {
    extractComp->SetIndex(1);
    writer->SetFileName( outputExpConstFilename );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & err )
    {
      std::cerr << "Error during write of " << outputExpConstFilename << std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return 1; 
    }
  }
  
  // Constant map.
  if ( outputConstFilename.length() > 0)
  {
    extractComp->SetIndex(2);
    writer->SetFileName( outputConstFilename );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & err )
    {
      std::cerr << "Error during write of " << outputConstFilename << std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return 1; 
    }
  }
  
  // Rsquared map.
  if ( outputRSquaredFilename.length() > 0)
  {
    extractComp->SetIndex(3);
    writer->SetFileName( outputRSquaredFilename );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & err )
    {
      std::cerr << "Error during write of " << outputRSquaredFilename << std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return 1; 
    }
  }

  return 0;
}
