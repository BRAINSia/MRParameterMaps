/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMRT1ParameterMap3DImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008/01/13 02:14:11 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMRT1ParameterMap3DImageFilter_txx
#define __itkMRT1ParameterMap3DImageFilter_txx

#include "itkMRT1ParameterMap3DImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "itkArray.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_least_squares_function.h"
#include "vnl/algo/vnl_levenberg_marquardt.h"
#include <math.h>

namespace itk
{

class vnl_ideal_steady_state_function : public vnl_least_squares_function
{
public:
    vnl_ideal_steady_state_function(bool with_grad, unsigned int n, 
      vnl_vector<double> t, vnl_vector<double> s)
    : vnl_least_squares_function(2,n,with_grad ? use_gradient : no_gradient) 
    {this->m_NumSignals=n; this->m_Time=t; this->m_Signal=s;}

  static double compute(double t, double a, double b) {
    return (a * ( 1.0 - exp( -t * b )));
    }
  static double compute_a(double t, double /* a */, double b) {
    return (1.0 - exp( -t * b ));
    }
  static double compute_b(double t, double a, double b) {
    return (t * a * exp( -t * b ));
    }

  void f(vnl_vector<double> const& x, vnl_vector<double>& y) {
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      y[i] = compute(this->m_Time[i], x(0), x(1) ) - this->m_Signal[i];
      }
    }

  void gradf(vnl_vector<double> const& x, vnl_matrix<double> &J) {
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      J(i,0) = compute_a(this->m_Time[i], x(0), x(1) );
      }
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      J(i,1) = compute_b(this->m_Time[i], x(0), x(1) );
      }
    }
  
private:
    unsigned int        m_NumSignals;
    vnl_vector<double>  m_Time;
    vnl_vector<double>  m_Signal;
};

class vnl_hybrid_steady_state_3param_function : public vnl_least_squares_function
{
public:
    vnl_hybrid_steady_state_3param_function(bool with_grad, unsigned int n, 
      vnl_vector<double> t, vnl_vector<double> s)
    : vnl_least_squares_function(3,n,with_grad ? use_gradient : no_gradient) 
    {this->m_NumSignals=n; this->m_Time=t; this->m_Signal=s;}

  static double compute(double t, double a, double b, double c) {
    return (a * ( b - exp( -t * c )));
    }
  static double compute_a(double t, double /* a */, double b, double c) {
    return (b - exp( -t * c ));
    }
  static double compute_b(double /* t */, double a, double /* b */, double /* c */) {
    return ( a );
    }
  static double compute_c(double t, double a, double /* b */, double c) {
    return (t * a * exp( -t * c ));
    }

  void f(vnl_vector<double> const& x, vnl_vector<double>& y) {
    for (unsigned int i=0; i<this->m_NumSignals; ++i)       {
      y[i] = compute(this->m_Time[i], x(0), x(1), x(2) ) - this->m_Signal[i];
      }
    }

  void gradf(vnl_vector<double> const& x, vnl_matrix<double> &J) {
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      J(i,0) = compute_a(this->m_Time[i], x(0), x(1), x(2) );
      }
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      J(i,1) = compute_b(this->m_Time[i], x(0), x(1), x(2) );
      }
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      J(i,2) = compute_c(this->m_Time[i], x(0), x(1), x(2) );
      }
    }
  
private:
    unsigned int        m_NumSignals;
    vnl_vector<double>  m_Time;
    vnl_vector<double>  m_Signal;
};

class vnl_inversion_recovery_function : public vnl_least_squares_function
{
public:
  vnl_inversion_recovery_function(bool with_grad, unsigned int n,
    vnl_vector<double> t, vnl_vector<double> s)
  : vnl_least_squares_function(2, n,with_grad ? use_gradient : no_gradient) 
  {this->m_NumSignals=n; this->m_Time=t; this->m_Signal=s;}

  static double compute(double t, double a, double b) {
    return (a * ( 1.0 - 2.0 * exp( -t * b )));
    }
  static double compute_a(double t, double /* a */, double b) {
    return (1.0 - 2.0 * exp( -t * b ));
    }
  static double compute_b(double t, double a, double b) {
    return (t * a * 2.0 * exp( -t * b ));
    }

  void f(vnl_vector<double> const& x, vnl_vector<double>& y) {
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      y[i] = compute(this->m_Time[i], x(0), x(1) ) - this->m_Signal[i];
      }
    }

  void gradf(vnl_vector<double> const& x, vnl_matrix<double> &J) {
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      J(i,0) = compute_a(this->m_Time[i], x(0), x(1) );
      }
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      J(i,1) = compute_b(this->m_Time[i], x(0), x(1) );
      }
    }
  
private:
    unsigned int        m_NumSignals;
    vnl_vector<double>  m_Time;
    vnl_vector<double>  m_Signal;
};

class vnl_inversion_recovery_3param_function : public vnl_least_squares_function
{
public:
  vnl_inversion_recovery_3param_function(bool with_grad, unsigned int n,
    vnl_vector<double> t, vnl_vector<double> s)
  : vnl_least_squares_function(3, n,with_grad ? use_gradient : no_gradient) 
  {this->m_NumSignals=n; this->m_Time=t; this->m_Signal=s;}

  static double compute(double t, double a, double b, double c) {
    return (a * ( 1.0 - b * exp( -t * c )));
    }
  static double compute_a(double t, double /* a */, double b, double c) {
    return (1.0 - b * exp( -t * c ));
    }
  static double compute_b(double t, double a, double /* b */, double c) {
    return ( - a * exp( -t * c ));
    }
  static double compute_c(double t, double a, double b, double c) {
    return (t * a * b * exp( -t * c ));
    }

  void f(vnl_vector<double> const& x, vnl_vector<double>& y) {
    for (unsigned int i=0; i<this->m_NumSignals; ++i)       {
      y[i] = compute(this->m_Time[i], x(0), x(1), x(2) ) - this->m_Signal[i];
      }
    }

  void gradf(vnl_vector<double> const& x, vnl_matrix<double> &J) {
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      J(i,0) = compute_a(this->m_Time[i], x(0), x(1), x(2) );
      }
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      J(i,1) = compute_b(this->m_Time[i], x(0), x(1), x(2) );
      }
    for (unsigned int i=0; i<this->m_NumSignals; ++i) {
      J(i,2) = compute_c(this->m_Time[i], x(0), x(1), x(2) );
      }
    }
  
private:
    unsigned int        m_NumSignals;
    vnl_vector<double>  m_Time;
    vnl_vector<double>  m_Signal;
};

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::MRT1ParameterMap3DImageFilter() 
{
  // At least 1 input is necessary for a vector image.
  // For images added one at a time we need at least 2
  this->SetNumberOfRequiredInputs( 1 ); 
  this->m_NumberOfImages = 0;
  this->m_MaxT1Time = 10.0f;
  this->m_PerformR1Mapping = false;
  this->m_MRImageTypeEnumeration = Else;
  this->m_TimeContainer = NULL;
  this->m_Algorithm = IDEAL_STEADY_STATE;
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void 
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "PerformR1Mapping: " << (this->m_PerformR1Mapping ? "On" 
    : "Off") << std::endl;
  if ( this->m_TimeContainer )
    {
    os << indent << "TimeContainer: "
      << this->m_TimeContainer << std::endl;
    }
  else
    {
    os << indent << "TimeContainer: (Times not set)" << std::endl;
    }
  os << indent << "NumberOfImages: " << this->m_NumberOfImages << std::endl;
  os << indent << "Maximum T1 time: " << this->m_MaxT1Time << std::endl;
  if ( this->m_MRImageTypeEnumeration == IsInASingleImage )
    {
    os << indent << "Multiple images haven been supplied " << std::endl;
    }
  else if ( this->m_MRImageTypeEnumeration == IsInManyImages )
    {
    os << indent << "A multicomponent image has been supplied" << std::endl;
    }
  if ( this->m_Algorithm == IDEAL_STEADY_STATE )
    {
    os << indent << "The IDEAL_STEADY_STATE algorithm is being used for the T1 "
      "fitting" << std::endl;
    }
  else if ( this->m_Algorithm == HYBRID_STEADY_STATE_3PARAM )
    {
    os << indent << "The HYBRID_STEADY_STATE_3PARAM algorithm is being used for "
      "the T1 fitting" << std::endl;
    }
   else if ( this->m_Algorithm == INVERSION_RECOVERY )
    {
    os << indent << "The INVERSION_RECOVERY algorithm is being used for the T1 "
      "fitting" << std::endl;
    }
  else if ( this->m_Algorithm == INVERSION_RECOVERY_3PARAM )
    {
    os << indent << "The INVERSION_RECOVERY_3PARAM algorithm is being used for "
      "the T1 fitting" << std::endl;
    }
  else if ( this->m_Algorithm == ABSOLUTE_INVERSION_RECOVERY )
    {
    os << indent << "The ABSOLUTE_INVERSION_RECOVERY algorithm is being used for"
      " the T1 fitting" << std::endl;
    }
  else if ( this->m_Algorithm == ABSOLUTE_INVERSION_RECOVERY_3PARAM )
    {
    os << indent << "The ABSOLUTE_INVERSION_RECOVERY_3PARAM algorithm is being "
      "used for the T1 fitting" << std::endl;
    }
  else if ( this->m_Algorithm == LOOK_LOCKER )
    {
    os << indent << "The LOOK_LOCKER algorithm is being used for the T1 fitting" 
      << std::endl;
    }
  else
    {
    os << indent << "The ABSOLUTE_LOOK_LOCKER algorithm is being used for the T1"
      " fitting" << std::endl;
    }
}

//----------------------------------------------------------------------------
template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void 
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::GenerateOutputInformation()
{
  // Override the method to set vector length
  Superclass::GenerateOutputInformation();

  typename Superclass::OutputImagePointer output = this->GetOutput();
  if( !output )
    {
    return;
    }
  // Vector length is always 4.
  output->SetVectorLength(4);
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void 
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::BeforeThreadedGenerateData()
{
  const unsigned int numberOfInputs = this->GetNumberOfInputs();

  // There need to be at least 2 images to be able to compute the T1 map.
  if( this->m_NumberOfImages < 2 )
    {
    itkExceptionMacro( << "At least 2 images are required" );
    }
    
  // If there is only 1 image, it must be an itk::VectorImage. Otherwise we must 
  // have a container of (numberOfInputs-1) itk::Image. Check to make sure
  if ( numberOfInputs == 1
      && this->m_MRImageTypeEnumeration != IsInASingleImage )
    {
    std::string imageClassName(
        this->ProcessObject::GetInput(0)->GetNameOfClass());
    if ( strcmp(imageClassName.c_str(),"VectorImage") != 0 )
      {
      itkExceptionMacro( << 
          "There is only one image. It should be a VectorImage. "
          << "But its of type: " << imageClassName );
      }
    }
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::FitIdealSteadyState(ExponentialFitType X, ExponentialFitType Y, 
  unsigned int num, MRParameterMapPixelType &output)
{
  ExponentialFitType temp;
  TimeType averageY = 0;
  TimeType SSE = 0;
  TimeType SST = 0;

  // Do non-linear fitting now.
  vnl_ideal_steady_state_function f(true,num,X,Y);

  vnl_levenberg_marquardt lm(f);

  // Find initial estimate first before doing nonlinear fitting.
  vnl_vector<TimeType> x1(2);
  x1[0] = Y[num-1]/0.95f; // Assume last input is 95% recovered signal.
  x1[1] = 1.0f; // Assume T1 is 1.0 seconds
  if (f.has_gradient())
    {
    lm.minimize_using_gradient(x1);
    }
  else
    {
    lm.minimize_without_gradient(x1);
    }
  //lm.diagnose_outcome(std::cout);
  output[0] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[1]);
  output[1] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[0]);
  output[2] = 0; // Not used for this fitting type.

  // Calculate R-squared.
  for(unsigned int i=0; i<num; i++)
    {
    double err = vnl_ideal_steady_state_function::compute(X[i],x1[0],x1[1])-Y[i];
    SSE += (err*err);
    averageY += Y[i];
    }
  averageY /= (double)num;
  temp = Y - averageY;
  SST = temp.squared_magnitude();
  output[3] = static_cast<typename MRParameterMapPixelType::ValueType>((SST != 0)
    ?fabs(1.0-(SSE/SST)):0);
  if( output[3] > 1 )
    {
    output[3] = 0;
    }
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::FitHybridSteadyState3Param(ExponentialFitType X, ExponentialFitType Y, 
  unsigned int num, MRParameterMapPixelType &output)
{
  ExponentialFitType temp;
  TimeType averageY = 0;
  TimeType SSE = 0;
  TimeType SST = 0;

  // Do non-linear fitting now.
  vnl_hybrid_steady_state_3param_function f(true,num,X,Y);

  vnl_levenberg_marquardt lm(f);

  // Find initial estimate first before doing nonlinear fitting.
  vnl_vector<TimeType> x1(3);
  x1[0] = Y[num-1]/0.95f; // Assume last input is 95% recovered signal.
  x1[1] = 1.0f; // Should always be close to 1 for steady state.
  x1[2] = 1.0f; // Assume T1 is 1.0 seconds
  if (f.has_gradient())
    {
    lm.minimize_using_gradient(x1);
    }
  else
    {
    lm.minimize_without_gradient(x1);
    }
  //lm.diagnose_outcome(std::cout);
  output[0] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[2]);
  output[1] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[0]);
  output[2] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[1]);

  // Calculate R-squared.
  for(unsigned int i=0; i<num; i++)
    {
    double err = vnl_hybrid_steady_state_3param_function::compute(X[i],x1[0],
      x1[1],x1[2])-Y[i];
    SSE += (err*err);
    averageY += Y[i];
    }
  averageY /= (double)num;
  temp = Y - averageY;
  SST = temp.squared_magnitude();
  output[3] = static_cast<typename MRParameterMapPixelType::ValueType>((SST != 0)
    ?fabs(1.0-(SSE/SST)):0);
  if( output[3] > 1 )
    {
    output[3] = 0;
    }
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::FitInversionRecovery(ExponentialFitType X, ExponentialFitType Y, 
  unsigned int num, MRParameterMapPixelType &output)
{
  ExponentialFitType temp;
  TimeType averageY = 0;
  TimeType SSE = 0;
  TimeType SST = 0;

  // Do non-linear fitting now.
  vnl_inversion_recovery_function f(true,num,X,Y);

  vnl_levenberg_marquardt lm(f);

  // Find initial estimate first before doing nonlinear fitting.
  vnl_vector<TimeType> x1(2);
  x1[0] = Y[num-1]/0.95f; // Assume last input is 95% recovered signal.
  x1[1] = 1.0f; // Assume T1 is 1.0 seconds
  if (f.has_gradient())
    {
    lm.minimize_using_gradient(x1);
    }
  else
    {
    lm.minimize_without_gradient(x1);
    }
  //lm.diagnose_outcome(std::cout);
  output[1] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[0]);
  output[0] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[1]);
  output[2] = 0; // Not used for this fitting type.

  // Calculate R-squared.
  for(unsigned int i=0; i<num; i++)
    {
    double err = vnl_inversion_recovery_function::compute(X[i],x1[0],x1[1])-Y[i];
    SSE += (err*err);
    averageY += Y[i];
    }
  averageY /= (double)num;
  temp = Y - averageY;
  SST = temp.squared_magnitude();
  output[3] = static_cast<typename MRParameterMapPixelType::ValueType>((SST != 0)
    ?fabs(1.0-(SSE/SST)):0);
  if( output[3] > 1 )
    {
    output[3] = 0;
    }
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::FitInversionRecovery3Param(ExponentialFitType X, ExponentialFitType Y, 
  unsigned int num, MRParameterMapPixelType &output)
{
  ExponentialFitType temp;
  TimeType averageY = 0;
  TimeType SSE = 0;
  TimeType SST = 0;

  // Do non-linear fitting now.
  vnl_inversion_recovery_3param_function f(true,num,X,Y);

  vnl_levenberg_marquardt lm(f);

  // Find initial estimate first before doing nonlinear fitting.
  vnl_vector<TimeType> x1(3);
  x1[0] = Y[num-1]/0.95f; // Assume last input is 95% recovered signal.
  x1[1] = 2.0f; // Should always be close to 2 for inversion recovery. 
  x1[2] = 1.0f; // Assume T1 is 1.0 seconds
  if (f.has_gradient())
    {
    lm.minimize_using_gradient(x1);
    }
  else
    {
    lm.minimize_without_gradient(x1);
    }
  //lm.diagnose_outcome(std::cout);
  output[0] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[2]);
  output[1] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[0]);
  output[2] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[1]);

  // Calculate R-squared.
  for(unsigned int i=0; i<num; i++)
    {
    double err = vnl_inversion_recovery_3param_function::compute(X[i],x1[0],x1[1]
      ,x1[2])-Y[i];
    SSE += (err*err);
    averageY += Y[i];
    }
  averageY /= (double)num;
  temp = Y - averageY;
  SST = temp.squared_magnitude();
  output[3] = static_cast<typename MRParameterMapPixelType::ValueType>((SST != 0)
    ?fabs(1.0-(SSE/SST)):0);
  if( output[3] > 1 )
    {
    output[3] = 0;
    }
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::FitAbsoluteInversionRecovery(ExponentialFitType X, ExponentialFitType Y, 
  unsigned int num, MRParameterMapPixelType &output)
{
  ExponentialFitType temp, Yfixed;
  TimeType averageY = 0;
  TimeType SSE = 0;
  TimeType SST = 0;

  // Make function non absolute again.  This is done by making sure the
  // Y values are always increasing.
  Yfixed = Y;
  for(unsigned int i=0; i<num-1; i++)
    {
    // Negative slope
    if( Y[i+1]-Y[i] < 0 )
      {
      Yfixed[i] = -Yfixed[i]; 
      }
    // Positive slope
    else
      {
      break;
      }
    }

  // Do non-linear fitting now.
  vnl_inversion_recovery_function f(true,num,X,Yfixed);

  vnl_levenberg_marquardt lm(f);

  // Find initial estimate first before doing nonlinear fitting.
  vnl_vector<TimeType> x1(2);
  x1[0] = Yfixed[num-1]/0.95f; // Assume last input is 75% recovered signal.
  x1[1] = 1.0f; // Assume T1 is 1.0 seconds
  if (f.has_gradient())
    {
    lm.minimize_using_gradient(x1);
    }
  else
    {
    lm.minimize_without_gradient(x1);
    }
  //lm.diagnose_outcome(std::cout);
  output[1] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[0]);
  output[0] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[1]);
  output[2] = 0; // Not used for this fitting type.

  // Calculate R-squared.
  for(unsigned int i=0; i<num; i++)
    {
    double err = vnl_inversion_recovery_function::compute(X[i],x1[0],x1[1])
      -Yfixed[i];
    SSE += (err*err);
    averageY += Yfixed[i];
    }
  averageY /= (double)num;
  temp = Yfixed - averageY;
  SST = temp.squared_magnitude();
  output[3] = static_cast<typename MRParameterMapPixelType::ValueType>((SST != 0)
    ?fabs(1.0-(SSE/SST)):0);
  if( output[3] > 1 )
    {
    output[3] = 0;
    }
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::FitAbsoluteInversionRecovery3Param(ExponentialFitType X, ExponentialFitType Y, 
  unsigned int num, MRParameterMapPixelType &output)
{
  ExponentialFitType temp, Yfixed;
  TimeType averageY = 0;
  TimeType SSE = 0;
  TimeType SST = 0;

  // Make function non absolute again.  This is done by making sure the
  // Y values are always increasing.
  Yfixed = Y;
  for(unsigned int i=0; i<num-1; i++)
    {
    // Negative slope
    if( Y[i+1]-Y[i] < 0 )
      {
      Yfixed[i] = -Yfixed[i]; 
      }
    // Positive slope
    else
      {
      break;
      }
    }

  // Do non-linear fitting now.
  vnl_inversion_recovery_3param_function f(true,num,X,Yfixed);

  vnl_levenberg_marquardt lm(f);

  // Find initial estimate first before doing nonlinear fitting.
  vnl_vector<TimeType> x1(3);
  x1[0] = Yfixed[num-1]/0.95f; // Assume last input is 95% recovered signal.
  x1[1] = 2.0f; // Should always be close to 2 for inversion recovery.
  x1[2] = 1.0f; // Assume T1 is 1.0 seconds
  if (f.has_gradient())
    {
    lm.minimize_using_gradient(x1);
    }
  else
    {
    lm.minimize_without_gradient(x1);
    }
  //lm.diagnose_outcome(std::cout);
  output[0] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[2]);
  output[1] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[0]);
  output[2] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[1]);

  // Calculate R-squared.
  for(unsigned int i=0; i<num; i++)
    {
    double err = vnl_inversion_recovery_3param_function::compute(X[i],x1[0],x1[1]
      ,x1[2])-Yfixed[i];
    SSE += (err*err);
    averageY += Yfixed[i];
    }
  averageY /= (double)num;
  temp = Yfixed - averageY;
  SST = temp.squared_magnitude();
  output[3] = static_cast<typename MRParameterMapPixelType::ValueType>((SST != 0)
    ?fabs(1.0-(SSE/SST)):0);
  if( output[3] > 1 )
    {
    output[3] = 0;
    }
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::FitLookLocker(ExponentialFitType X, ExponentialFitType Y, unsigned int num, 
  MRParameterMapPixelType &output)
{
  ExponentialFitType temp;
  TimeType averageY = 0;
  TimeType SSE = 0;
  TimeType SST = 0;

  // Do non-linear fitting now.
  vnl_inversion_recovery_3param_function f(true,num,X,Y);

  vnl_levenberg_marquardt lm(f);

  // Find initial estimate first before doing nonlinear fitting.
  vnl_vector<TimeType> x1(3);
  x1[0] = Y[num-1]/0.95f; // Assume last input is 95% recovered signal.
  x1[1] = 2.0f; // Should always be close to 2 for Look-Locker. 
  x1[2] = 1.0f; // Assume T1 is 1.0 seconds
  if (f.has_gradient())
    {
    lm.minimize_using_gradient(x1);
    }
  else
    {
    lm.minimize_without_gradient(x1);
    }
  //lm.diagnose_outcome(std::cout);
  output[0] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[2]);
  output[1] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[0]);
  output[2] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[1]);

  // Calculate R-squared.
  for(unsigned int i=0; i<num; i++)
    {
    double err = vnl_inversion_recovery_3param_function::compute(X[i],x1[0],x1[1]
      ,x1[2])-Y[i];
    SSE += (err*err);
    averageY += Y[i];
    }
  averageY /= (double)num;
  temp = Y - averageY;
  SST = temp.squared_magnitude();
  output[3] = static_cast<typename MRParameterMapPixelType::ValueType>((SST != 0)
    ?fabs(1.0-(SSE/SST)):0);
  if( output[3] > 1 )
    {
    output[3] = 0;
    }
  // Change T1* to T1
  if( output[1] != 0.0f )
    {
    output[0] = output[0]*(output[2]-1.0f);
    } 
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::FitAbsoluteLookLocker(ExponentialFitType X, ExponentialFitType Y, 
  unsigned int num, MRParameterMapPixelType &output)
{
  ExponentialFitType temp, Yfixed;
  TimeType averageY = 0;
  TimeType SSE = 0;
  TimeType SST = 0;

  // Make function non absolute again.  This is done by making sure the
  // Y values are always increasing.
  Yfixed = Y;
  for(unsigned int i=0; i<num-1; i++)
    {
    // Negative slope
    if( Y[i+1]-Y[i] < 0 )
      {
      Yfixed[i] = -Yfixed[i]; 
      }
    // Positive slope
    else
      {
      break;
      }
    }

  // Do non-linear fitting now.
  vnl_inversion_recovery_3param_function f(true,num,X,Yfixed);

  vnl_levenberg_marquardt lm(f);

  // Find initial estimate first before doing nonlinear fitting.
  vnl_vector<TimeType> x1(3);
  x1[0] = Yfixed[num-1]/0.95f; // Assume last input is 95% recovered signal.
  x1[1] = 2.0f; // Should always be close to 2 for Look-Locker.
  x1[2] = 1.0f; // Assume T1 is 1.0 seconds
  if (f.has_gradient())
    {
    lm.minimize_using_gradient(x1);
    }
  else
    {
    lm.minimize_without_gradient(x1);
    }
  //lm.diagnose_outcome(std::cout);
  output[0] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[2]);
  output[1] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[0]);
  output[2] = static_cast<typename MRParameterMapPixelType::ValueType>(x1[1]);

  // Calculate R-squared.
  for(unsigned int i=0; i<num; i++)
    {
    double err = vnl_inversion_recovery_3param_function::compute(X[i],x1[0],x1[1]
      ,x1[2])-Yfixed[i];
    SSE += (err*err);
    averageY += Yfixed[i];
    }
  averageY /= (double)num;
  temp = Yfixed - averageY;
  SST = temp.squared_magnitude();
  output[3] = static_cast<typename MRParameterMapPixelType::ValueType>((SST != 0)
    ?fabs(1.0-(SSE/SST)):0);
  if( output[3] > 1 )
    {
    output[3] = 0;
    }
  // Change T1* to T1
  if( output[1] != 0.0f )
    {
    output[0] = output[0]*(output[2]-1.0f);
    } 
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
  ThreadIdType threadId)
{
  typename OutputImageType::Pointer outputImage = 
    static_cast< OutputImageType * >(this->ProcessObject::GetOutput(0));
  
  ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
  oit.GoToBegin();
  
  ProgressReporter progress(this, threadId,
    outputRegionForThread.GetNumberOfPixels(), 100);
    
  // Two cases here .
  // 1. If the images have been specified in multiple images, we will create
  // 'n' iterators for each of the images and fit the T1 curve for each voxel. 
  // 2. If the images have been specified in a single multi-component image,
  // one iterator will suffice to do the same.

  if( this->m_MRImageTypeEnumeration == IsInManyImages )
    {

    typedef ImageRegionConstIterator< MRImageType > MRImageIteratorType;
    std::vector< MRImageIteratorType * > mrImageItContainer;
    int nonzeroCount = 0;

    for( unsigned int i = 0; i< this->m_NumberOfImages; i++ )
      {
      typename MRImageType::Pointer mrImagePointer = NULL;
      
      mrImagePointer = static_cast< MRImageType * >( 
        this->ProcessObject::GetInput(i) );

      MRImageIteratorType *eit = new MRImageIteratorType( 
                          mrImagePointer, outputRegionForThread );
      eit->GoToBegin();
      mrImageItContainer.push_back(eit);
      }
    
    // Iterate over the images and fit the T1 curve.
    while( !oit.IsAtEnd() )
      {
      
      MRParameterMapPixelType map;
      map.SetSize(4);
      map.Fill(0);
      // Create array for T1 calculation.
      nonzeroCount = -1;
      vnl_vector<TimeType> imageValues(this->m_NumberOfImages);
      vnl_vector<TimeType> timePoints(this->m_NumberOfImages);
      for( unsigned int i = 0; i< this->m_NumberOfImages; i++ )
        {
        imageValues[i] = static_cast<TimeType>(mrImageItContainer[i]->Get());
        timePoints[i] = this->m_TimeContainer->ElementAt(i);
        // The calculations above don't like zeros.  This is a problem for
        // inversion recovery where a zero is possible.  However, the 
        // assumption will be made that an exact zero is not likely and is
        // due to a masked or thresholded voxel.
        if( (imageValues[i] == 0) && (nonzeroCount < 0) )
          {
          nonzeroCount = i;
          }
        ++(*mrImageItContainer[i]);
        }
      if( nonzeroCount < 0 )
        {
        nonzeroCount = this->m_NumberOfImages;
        }

      // Only do the calculation if we have at least 2 contiguous non-zero 
      // values.
      if( nonzeroCount >= 2 )
        {
        switch( this->m_Algorithm )
          {
          case IDEAL_STEADY_STATE:
            FitIdealSteadyState(timePoints,imageValues,nonzeroCount,map);
            break;
          case INVERSION_RECOVERY:
            FitInversionRecovery(timePoints,imageValues,nonzeroCount,map);
            break;
          case ABSOLUTE_INVERSION_RECOVERY:
            FitAbsoluteInversionRecovery(timePoints,imageValues,
              nonzeroCount,map);
            break;
          case LOOK_LOCKER:
            FitLookLocker(timePoints,imageValues,nonzeroCount,map);
            break;
          case ABSOLUTE_LOOK_LOCKER:
            FitAbsoluteLookLocker(timePoints,imageValues,nonzeroCount,map);
            break;
          case HYBRID_STEADY_STATE_3PARAM:
            FitHybridSteadyState3Param(timePoints,imageValues,nonzeroCount,map);
            break;
          case INVERSION_RECOVERY_3PARAM:
            FitInversionRecovery3Param(timePoints,imageValues,nonzeroCount,map);
            break;
          case ABSOLUTE_INVERSION_RECOVERY_3PARAM:
            FitAbsoluteInversionRecovery3Param(timePoints,imageValues,
              nonzeroCount,map);
            break;
          default:
            itkExceptionMacro( << "Unknown fit type: " << this->m_Algorithm );
          }
        // Do inverse of R1 if not
        // performing R1 mapping.
        if( !this->m_PerformR1Mapping )
          {
          if( map[0] != 0 )
            {
            map[0] = 1.0f/map[0];
            }
          if( static_cast<TimeType>(map[0]) > this->m_MaxT1Time )
            {
            map[0] = static_cast<MRParameterPixelType>(this->m_MaxT1Time);
            }
          }
        else
          {
          if( this->m_MaxT1Time == 0 )
            {
            this->m_MaxT1Time = 1.0f/NumericTraits< TimeType >::max();
            }
          if( static_cast<TimeType>(map[0]) < (1.0f/this->m_MaxT1Time) )
            {
              map[0] = static_cast<MRParameterPixelType>(1.0f/this->m_MaxT1Time);
            }
          }

        // Should never have a value less than zero,
        // so automatically filter that out.
        if( map[0] < 0 )
          {
          map.Fill(0);
          }
      }

      oit.Set( map );
      ++oit;
      progress.CompletedPixel();
      }

      for( unsigned int i = 0; i< mrImageItContainer.size(); i++ )
        {
        delete mrImageItContainer[i];
        }
    }
  // The images are specified in a single multi-component image
  else if( this->m_MRImageTypeEnumeration == IsInASingleImage )
    {
    typedef ImageRegionConstIterator< MRImagesType > MRImageIteratorType;
    typedef typename MRImagesType::PixelType         MREchoVectorType;
    typename MRImagesType::Pointer mrImagePointer = NULL;
    int nonzeroCount = 0;

    mrImagePointer = static_cast< MRImagesType * >(
      this->ProcessObject::GetInput(0) );

    MRImageIteratorType eit(mrImagePointer, outputRegionForThread );
    eit.GoToBegin();

    while( !eit.IsAtEnd() )
      {
      MREchoVectorType echoValues = eit.Get();

      MRParameterMapPixelType map;
      map.SetSize(4);
      map.Fill(0);
      // Create array for T1 calculation.
      nonzeroCount = -1;
      vnl_vector<TimeType> imageValues(this->m_NumberOfImages);
      vnl_vector<TimeType> timePoints(this->m_NumberOfImages);
      for( unsigned int i = 0; i< this->m_NumberOfImages; i++ )
        {
        imageValues[i] = static_cast<TimeType>(echoValues[i]);
        timePoints[i] = this->m_TimeContainer->ElementAt(i);
        // The calculations above don't like zeros.  This is a problem for
        // inversion recovery where a zero is possible.  However, the 
        // assumption will be made that an exact zero is not likely and is
        // due to a masked or thresholded voxel.
        if( (imageValues[i] == 0) && (nonzeroCount < 0) )
          {
          nonzeroCount = i;
          }
        }
      if( nonzeroCount < 0 )
        {
        nonzeroCount = this->m_NumberOfImages;
        }

      // Only do the calculation if we have at least 2 contiguous non-zero 
      // values.
      if( nonzeroCount >= 2 )
        {
        switch( this->m_Algorithm )
          {
          case IDEAL_STEADY_STATE:
          FitIdealSteadyState(timePoints,imageValues,nonzeroCount,map);
            break;
          case INVERSION_RECOVERY:
          FitInversionRecovery(timePoints,imageValues,nonzeroCount,map);
            break;
          case ABSOLUTE_INVERSION_RECOVERY:
          FitAbsoluteInversionRecovery(timePoints,imageValues,nonzeroCount,map);
            break;
          case LOOK_LOCKER:
          FitLookLocker(timePoints,imageValues,nonzeroCount,map);
            break;
          case ABSOLUTE_LOOK_LOCKER:
          FitAbsoluteLookLocker(timePoints,imageValues,nonzeroCount,map);
            break;
          case HYBRID_STEADY_STATE_3PARAM:
          FitHybridSteadyState3Param(timePoints,imageValues,nonzeroCount,map);
            break;
          case INVERSION_RECOVERY_3PARAM:
          FitInversionRecovery3Param(timePoints,imageValues,nonzeroCount,map);
            break;
          case ABSOLUTE_INVERSION_RECOVERY_3PARAM:
          FitAbsoluteInversionRecovery3Param(timePoints,imageValues,nonzeroCount,
            map);
            break;
          default:
            itkExceptionMacro( << "Unknown fit type: " << this->m_Algorithm );
          }

        // Do inverse of R1 if not
        // performing R1 mapping.
        if( !this->m_PerformR1Mapping )
          {
          if( map[0] != 0 )
            {
            map[0] = 1.0f/map[0];
            }
          if( static_cast<TimeType>(map[0]) > this->m_MaxT1Time )
            {
            map[0] = static_cast<MRParameterPixelType>(this->m_MaxT1Time);
            }
          }
        else
          {
          if( this->m_MaxT1Time == 0 )
            {
            this->m_MaxT1Time = 1.0f/NumericTraits< TimeType >::max();
            }
          if( static_cast<TimeType>(map[0]) < (1.0f/this->m_MaxT1Time) )
            {
            map[0] = static_cast<MRParameterPixelType>(1.0f/this->m_MaxT1Time);
            }
          }
      
        // Should never have a value less than zero,
        // so automatically filter that out.
        if( map[0] < 0 )
          {
          map.Fill(0);
          }
      }

    oit.Set( map );
    ++oit; // Output (fitted image parameters) iterator
    ++eit; // Echo image iterator
    progress.CompletedPixel();
    }
  }
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::AddMRImage( TimeType timePoint, const MRImageType *image)
{
  // Make sure crazy users did not call both AddMRImage and 
  // SetMRImage
  if( this->m_MRImageTypeEnumeration == IsInASingleImage)
    {
    itkExceptionMacro( << "Cannot call both methods:" 
    << "AddMRImage and SetMRImage. Please call only one of them.");
    }

  // If the container to hold the echo times hasn't been allocated
  // yet, allocate it.
  if( !this->m_TimeContainer )
    {
    this->m_TimeContainer = TimeContainerType::New();
    }
    
  m_TimeContainer->InsertElement( this->m_NumberOfImages, timePoint );
  this->ProcessObject::SetNthInput( this->m_NumberOfImages, 
    const_cast< MRImageType* >(image) );
  ++this->m_NumberOfImages;
  this->m_MRImageTypeEnumeration = IsInManyImages;
}

template< class TMRImagePixelType, class TMRParameterMapImagePixelType >
void
MRT1ParameterMap3DImageFilter<TMRImagePixelType,TMRParameterMapImagePixelType>
::SetMRImage( TimeContainerType *timePointContainer, const MRImagesType *image )
{
  // Make sure crazy users did not call both AddMRImage and 
  // SetMRImage
  if( this->m_MRImageTypeEnumeration == IsInManyImages )
    {
    itkExceptionMacro( << "Cannot call both methods:" 
    << "AddMRImage and SetMRImage. Please call only one of them.");
    }

  this->m_TimeContainer = timePointContainer;

  this->m_NumberOfImages = timePointContainer->Size();

  // ensure that the image we received has as many components as 
  // the number of time points
  if( image->GetVectorLength() != this->m_NumberOfImages )
    {
    itkExceptionMacro( << this->m_NumberOfImages << 
      " time points specified but image has " << image->GetVectorLength()
      << " components.");
    }
  
  this->ProcessObject::SetNthInput( 0, const_cast< MRImagesType* >(image) );
  this->m_MRImageTypeEnumeration = IsInASingleImage;
}


} // end namespace itk

#endif
