/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMRT1ParameterMap3DImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/01/13 02:14:11 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
 
  Contributor(s):

=========================================================================*/
#ifndef __itkMRT1ParameterMap3DImageFilter_h
#define __itkMRT1ParameterMap3DImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"
#include "vnl/vnl_vector.h"

namespace itk
{

/** \class MRT1ParameterMap3DImageFilter
 * \brief Generate a 3D MR T1 parameter map using multiple images.
 *
 * This filter is templated over the input pixel type and the output pixel
 * type. The input T1-weighted MR images must all have the same size and 
 * dimensions.  The 3D \c VectorImage output will have four components.  
 * \c VectorImageToImageAdaptor should be used to extract the various 
 * components for viewing or saving to file.
 *
 * Given multiple T1-weighted MR images and a vector of TR (saturation
 * recovery) times or TI times (inversion recovery), generate the T1 or 
 * R1 parameter map.  The exponential fitting function may be changed 
 * using the \c SetAlgorithm function.
 *
 * \par Inputs and Usage
 * There are two ways to use this class. When you have multiple images 
 * you would use the class as
 * \code
 *       filter->AddMRImage( time1, image1 );
 *       filter->AddMRImage( time2, image2 );
 *   ...
 * \endcode
 *
 * \par
 * When you have the 'n' MR images in a single multi-component image 
 * (VectorImage), you can specify the images simply as
 * \code
 *       filter->SetMRImage( timeContainer, vectorImage );
 * \endcode
 *
 * \par Outputs
 * The output image is a vector image containing 4 values:
 * The first component of the output will be the T1 time in seconds (R1 in Hz), 
 * the second component will be the constant A as shown in the functions below,
 * and the fourth component will be the r-squared value from the curve fitting.
 * The third component will vary depending on the type of T1 fitting selected.  
 * For all of the 2 parameter models below the third component will be zero.  
 * For the remaining models the third component will be the value B as shown 
 * below. 
 * 
 * IDEAL_STEADY_STATE (Non-linear least squares using Levenberg-Marquardt):
 * S = A(1-e^-(TR/T1))
 *
 * HYBRID_STEADY_STATE_3PARAM (Non-linear least squares using 
 * Levenberg-Marquardt):
 * S = A(B-e^-(TR/T1))
 *
 * INVERSION_RECOVERY (Non-linear least squares using Levenberg-Marquardt):
 * S = A(1-2e^-(TI/T1))
 *
 * INVERSION_RECOVERY_3PARAM (Non-linear least squares using 
 * Levenberg-Marquardt):
 * S = A(1-Be^-(TI/T1))
 *
 * ABSOLUTE_INVERSION_RECOVERY (Non-linear least squares using 
 * Levenberg-Marquardt):
 * S = abs(A(1-2e^-(TI/T1)))
 *
 * ABSOLUTE_INVERSION_RECOVERY_3PARAM (Non-linear least squares using 
 * Levenberg-Marquardt):
 * S = abs(A(1-Be^-(TI/T1)))
 *
 * LOOK_LOCKER (Non-linear least squares using Levenberg-Marquardt):
 * S = A(1-Be^-(TI/T1*)), T1=T1*(B-1)
 * Deichmann R, Haase A. Quantification of T1 valeus by snapshot
 * FLASH NMR imaging. J Magn Reson 1992;96:608–612.
 *
 * ABSOLUTE_LOOK_LOCKER (Non-linear least squares using Levenberg-Marquardt):
 * S = abs(A(1-Be^-(TI/T1*))), T1=T1*(B-1)
 *
 * \code
 *       VectorImage< TMRParameterMapImagePixelType, 3 >
 * \endcode
 *
 * \par Parameters
 * \li Algorithm -  Set/Get the T1 fitting algorithm to use.
 * \li MaxT1Time - Set the maximum T1 time (T1 times greater than or equal to
 * to this value will be set to this value).
 * \li PerformR1Mapping - If On R1 (in Hz) will be calculated instead of T1.
 *
 * \ingroup Multithreaded
 *
 * \author Don Bigler, Center for Nuclear Magnetic Resonance Research, 
 * Penn State Milton S. Hershey Medical Center
 *
 */
template< class TMRImagePixelType, class TMRParameterMapImagePixelType=double >
class ITK_EXPORT MRT1ParameterMap3DImageFilter:
    public ImageToImageFilter<Image< TMRImagePixelType, 3 >, 
                              VectorImage< TMRParameterMapImagePixelType, 3 > >
{
public:
  /** Standard class typedefs. */
  typedef MRT1ParameterMap3DImageFilter              Self;
  typedef ImageToImageFilter<Image< TMRImagePixelType, 3 >, 
          VectorImage< TMRParameterMapImagePixelType, 3 > >  
                                                     Superclass;
  typedef SmartPointer<Self>                         Pointer;
  typedef SmartPointer<const Self>                   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MRT1ParameterMap3DImageFilter, ImageToImageFilter);
  
  
  typedef TMRImagePixelType                          MRPixelType;
  typedef typename VectorImage< TMRParameterMapImagePixelType, 3 >::PixelType
                                                     MRParameterMapPixelType;
  typedef TMRParameterMapImagePixelType              MRParameterPixelType;

  /** T1-weighted MR image data. */
  typedef typename Superclass::InputImageType        MRImageType;
  
  /** An alternative typedef defining one (of the many) images. 
   * It will be assumed that the vectorImage has a vector length 
   * parameter of \c n (number of time points) */
  typedef VectorImage< MRPixelType, 3 >              MRImagesType;
  
  typedef typename Superclass::OutputImageType       MRParameterMapImageType;
  typedef MRParameterMapImageType                    OutputImageType;
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

  /** Data type of echo times */
  typedef double                                     TimeType;
  typedef vnl_vector<TimeType>                       ExponentialFitType;

  /** Container to hold echo times of the 'n' MR image measurements */
  typedef VectorContainer< unsigned int,TimeType >   TimeContainerType;

  /** Set method to add the time (in seconds) and its corresponding image. */
  void AddMRImage( TimeType, const MRImageType *image);

  /** Another set method to add a time point (in seconds) and its corresponding
   * image. The image here is a VectorImage. The user is expected to pass the 
   * times in a container. The ith element of the container corresponds 
   * to the time of the ith component image of the VectorImage. */
  void SetMRImage( TimeContainerType *, const MRImagesType *image);

  /** define values used to determine which algorithm to use */
  static const int IDEAL_STEADY_STATE = 0;
  static const int INVERSION_RECOVERY = 1;
  static const int ABSOLUTE_INVERSION_RECOVERY = 2;
  static const int LOOK_LOCKER = 3;
  static const int ABSOLUTE_LOOK_LOCKER = 4;
  static const int HYBRID_STEADY_STATE_3PARAM = 5;
  static const int INVERSION_RECOVERY_3PARAM = 6;
  static const int ABSOLUTE_INVERSION_RECOVERY_3PARAM = 7;
  
  /** Set/Get the T1 fitting algorithm (default is IDEAL_STEADY_STATE). */
  itkSetMacro(Algorithm, int);
  itkGetMacro(Algorithm, int);
  
  /** Set/Get the maximum T1 time (default is 10.0). */
  itkSetMacro(MaxT1Time, TimeType);
  itkGetMacro(MaxT1Time, TimeType);
    
  /** Set/Get whether to perform R1 mapping instead of T21 mapping (default is 
  /* false). */
  itkSetMacro( PerformR1Mapping, bool );
  itkGetConstReferenceMacro( PerformR1Mapping, bool );
  itkBooleanMacro( PerformR1Mapping );
  
#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(MREchoPixelConvertibleToDoubleCheck,
    (Concept::Convertible<MRPixelType, double>));
  itkConceptMacro(MRParameterMapPixelConvertibleToDoubleCheck,
    (Concept::Convertible<MRParameterPixelType, double>));
 /** End concept checking */
#endif
  
protected:
  MRT1ParameterMap3DImageFilter();
  ~MRT1ParameterMap3DImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  void FitIdealSteadyState(ExponentialFitType X, ExponentialFitType Y, 
      unsigned int num, MRParameterMapPixelType &output);
  void FitHybridSteadyState3Param(ExponentialFitType X, ExponentialFitType Y, 
      unsigned int num, MRParameterMapPixelType &output);
  void FitInversionRecovery(ExponentialFitType X, ExponentialFitType Y, 
      unsigned int num, MRParameterMapPixelType &output);
  void FitInversionRecovery3Param(ExponentialFitType X, ExponentialFitType Y, 
      unsigned int num, MRParameterMapPixelType &output);
  void FitAbsoluteInversionRecovery(ExponentialFitType X, ExponentialFitType Y, 
      unsigned int num, MRParameterMapPixelType &output);
  void FitAbsoluteInversionRecovery3Param(ExponentialFitType X, 
      ExponentialFitType Y, unsigned int num, MRParameterMapPixelType &output);
  void FitLookLocker(ExponentialFitType X, ExponentialFitType Y, 
      unsigned int num, MRParameterMapPixelType &output);
  void FitAbsoluteLookLocker(ExponentialFitType X, ExponentialFitType Y, 
      unsigned int num, MRParameterMapPixelType &output);
  
  /** Setup vector image vector length. */
  virtual void GenerateOutputInformation();
  
  void BeforeThreadedGenerateData();
  void ThreadedGenerateData( const 
      OutputImageRegionType &outputRegionForThread, int);
  
  /** enum to indicate if the MR echo image is specified as a single multi-
   * component image or as several separate images */
  typedef enum
    {
    IsInASingleImage = 1,
    IsInManyImages,
    Else
    } MRImageTypeEnumeration;

private:
  MRT1ParameterMap3DImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  int      m_Algorithm;
  
  TimeType m_MaxT1Time;
  
  bool     m_PerformR1Mapping;
  
  /** container to hold the time points */
  TimeContainerType::Pointer m_TimeContainer;

  /** Number of images */
  unsigned int               m_NumberOfImages;

  /** MR image was specified in a single image or in multiple images */
  MRImageTypeEnumeration     m_MRImageTypeEnumeration;

};

  
} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMRT1ParameterMap3DImageFilter.txx"
#endif
  
#endif
