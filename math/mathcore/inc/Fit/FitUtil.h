// @(#)root/mathcore:$Id$
// Author: L. Moneta Tue Nov 28 10:52:47 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class FitUtil

#ifndef ROOT_Fit_FitUtil
#define ROOT_Fit_FitUtil

#include "Math/IParamFunctionfwd.h"
#include "Math/IParamFunction.h"
#include "ROOT/TThreadExecutor.hxx"
#include "ROOT/TProcessExecutor.hxx"

#include "Fit/BinData.h"
#include "Fit/UnBinData.h"

#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"

#include "TError.h"
#include "TSystem.h"

namespace ROOT {

   namespace Fit {

/**
   namespace defining utility free functions using in Fit for evaluating the various fit method
   functions (chi2, likelihood, etc..)  given the data and the model function

   @ingroup FitMain
*/
namespace FitUtil {

  typedef  ROOT::Math::IParamMultiFunction IModelFunction;
  typedef  ROOT::Math::IParamMultiGradFunction IGradModelFunction;

  template<class T>
  using IModelFunctionTempl = ROOT::Math::IParamMultiFunctionTempl<T>;

   //internal class defining
         template<class T>
         class LikelihoodAux{
          public:

          LikelihoodAux(T logv={}, T w={}, T w2={}): logvalue(logv), weight(w), weight2(w2){
          }

           LikelihoodAux operator +( const LikelihoodAux & l) const{
              return LikelihoodAux<T>(logvalue + l.logvalue, weight  + l.weight, weight2 + l.weight2);
           }

           LikelihoodAux &operator +=(const LikelihoodAux & l){
              logvalue += l.logvalue;
              weight  += l.weight;
              weight2 += l.weight2;
              return *this;
           }

            T logvalue;
            T weight;
            T weight2;

         };

         template<>
         class LikelihoodAux<double>{
          public:

          LikelihoodAux(double logv =0.0, double w = 0.0, double w2 = 0.0):logvalue(logv), weight(w), weight2(w2){};

           LikelihoodAux operator +( const LikelihoodAux & l) const{
              return LikelihoodAux<double>(logvalue + l.logvalue, weight  + l.weight, weight2 + l.weight2);
           }

           LikelihoodAux &operator +=(const LikelihoodAux & l){
              logvalue += l.logvalue;
              weight  += l.weight;
              weight2 += l.weight2;
              return *this;
           }

            double logvalue;
            double weight;
            double weight2;
         };

         // internal class to evaluate the function or the integral
         // and cached internal integration details
         // if useIntegral is false no allocation is done
         // and this is a dummy class
         // class is templated on any parametric functor implementing operator()(x,p) and NDim()
         // contains a constant pointer to the function

         template <class ParamFunc = ROOT::Math::IParamMultiFunctionTempl<double>>
         class IntegralEvaluator {

         public:

            IntegralEvaluator(const ParamFunc & func, const double * p, bool useIntegral = true) :
               fDim(0),
               fParams(0),
               fFunc(0),
               fIg1Dim(0),
               fIgNDim(0),
               fFunc1Dim(0),
               fFuncNDim(0)
            {
               if (useIntegral) {
                  SetFunction(func, p);
               }
            }

            void SetFunction(const ParamFunc & func, const double * p = 0) {
               // set the integrand function and create required wrapper
               // to perform integral in (x) of a generic  f(x,p)
               fParams = p;
               fDim = func.NDim();
               // copy the function object to be able to modify the parameters
               //fFunc = dynamic_cast<ROOT::Math::IParamMultiFunction *>( func.Clone() );
               fFunc = &func;
               assert(fFunc != 0);
               // set parameters in function
               //fFunc->SetParameters(p);
               if (fDim == 1) {
                  fFunc1Dim = new ROOT::Math::WrappedMemFunction< IntegralEvaluator, double (IntegralEvaluator::*)(double ) const > (*this, &IntegralEvaluator::F1);
                  fIg1Dim = new ROOT::Math::IntegratorOneDim();
                  //fIg1Dim->SetFunction( static_cast<const ROOT::Math::IMultiGenFunction & >(*fFunc),false);
                  fIg1Dim->SetFunction( static_cast<const ROOT::Math::IGenFunction &>(*fFunc1Dim) );
               }
               else if (fDim > 1) {
                  fFuncNDim = new ROOT::Math::WrappedMemMultiFunction< IntegralEvaluator, double (IntegralEvaluator::*)(const double *) const >  (*this, &IntegralEvaluator::FN, fDim);
                  fIgNDim = new ROOT::Math::IntegratorMultiDim();
                  fIgNDim->SetFunction(*fFuncNDim);
               }
               else
                  assert(fDim > 0);
            }


            void SetParameters(const double *p) {
               // copy just the pointer
               fParams = p;
            }

            ~IntegralEvaluator() {
               if (fIg1Dim) delete fIg1Dim;
               if (fIgNDim) delete fIgNDim;
               if (fFunc1Dim) delete fFunc1Dim;
               if (fFuncNDim) delete fFuncNDim;
               //if (fFunc) delete fFunc;
            }

            // evaluation of integrand function (one-dim)
            double F1 (double x) const {
               double xx= x;
               return ExecFunc(fFunc, &xx, fParams);
            }
            // evaluation of integrand function (multi-dim)
            double FN(const double * x) const {
               return ExecFunc(fFunc, x, fParams);
            }

            double Integral(const double *x1, const double * x2) {
               // return unormalized integral
               return (fIg1Dim) ? fIg1Dim->Integral( *x1, *x2) : fIgNDim->Integral( x1, x2);
            }

            double operator()(const double *x1, const double * x2) {
               // return normalized integral, divided by bin volume (dx1*dx...*dxn)
               if (fIg1Dim) {
                  double dV = *x2 - *x1;
                  return fIg1Dim->Integral( *x1, *x2)/dV;
               }
               else if (fIgNDim) {
                  double dV = 1;
                  for (unsigned int i = 0; i < fDim; ++i)
                     dV *= ( x2[i] - x1[i] );
                  return fIgNDim->Integral( x1, x2)/dV;
//                   std::cout << " do integral btw x " << x1[0] << "  " << x2[0] << " y " << x1[1] << "  " << x2[1] << " dV = " << dV << " result = " << result << std::endl;
//                   return result;
               }
               else
                  assert(1.); // should never be here
               return 0;
            }

         private:

            template<class T>
            inline double ExecFunc(T *f, const double *x, const double *p) const{
                return (*f)(x, p);
            }

            inline double ExecFunc(const IModelFunctionTempl<Double_v> *f, const double *x, const double *p) const{
                const Double_v xx = vecCore::FromPtr<Double_v>(x);
                const double *p0 = p;
                auto res =  (*f)( &xx, (const double *)p0);
                return res[0];
            }

            // objects of this class are not meant to be copied / assigned
            IntegralEvaluator(const IntegralEvaluator& rhs);
            IntegralEvaluator& operator=(const IntegralEvaluator& rhs);

            unsigned int fDim;
            const double * fParams;
            //ROOT::Math::IParamMultiFunction * fFunc;  // copy of function in order to be able to change parameters
            // const ParamFunc * fFunc;       //  reference to a generic parametric function
            const ParamFunc * fFunc;
            ROOT::Math::IntegratorOneDim * fIg1Dim;
            ROOT::Math::IntegratorMultiDim * fIgNDim;
            ROOT::Math::IGenFunction * fFunc1Dim;
            ROOT::Math::IMultiGenFunction * fFuncNDim;
         };

   /** Chi2 Functions */

   /**
       evaluate the Chi2 given a model function and the data at the point x.
       return also nPoints as the effective number of used points in the Chi2 evaluation
   */
   double EvaluateChi2(const IModelFunction & func, const BinData & data, const double * x, unsigned int & nPoints, const unsigned int &executionPolicy, unsigned nChunks=0);

   /**
       evaluate the effective Chi2 given a model function and the data at the point x.
       The effective chi2 uses the errors on the coordinates : W = 1/(sigma_y**2 + ( sigma_x_i * df/dx_i )**2 )
       return also nPoints as the effective number of used points in the Chi2 evaluation
   */
   double EvaluateChi2Effective(const IModelFunction & func, const BinData & data, const double * x, unsigned int & nPoints);

   /**
       evaluate the Chi2 gradient given a model function and the data at the point x.
       return also nPoints as the effective number of used points in the Chi2 evaluation
   */
   void EvaluateChi2Gradient(const IModelFunction & func, const BinData & data, const double * x, double * grad, unsigned int & nPoints);

   /**
       evaluate the LogL given a model function and the data at the point x.
       return also nPoints as the effective number of used points in the LogL evaluation
   */
   template<class T>
   double EvaluateLogL(const IModelFunctionTempl<T> & func, const UnBinData & data, const double * p, int iWeight, bool extended, unsigned int & nPoints, const unsigned int &executionPolicy, unsigned nChunks=0);

   double EvaluateLogL(const IModelFunction & func, const UnBinData & data, const double * p, int iWeight, bool extended, unsigned int & nPoints, const unsigned int &executionPolicy, unsigned nChunks=0);

   /**
       evaluate the LogL gradient given a model function and the data at the point x.
       return also nPoints as the effective number of used points in the LogL evaluation
   */
   void EvaluateLogLGradient(const IModelFunction & func, const UnBinData & data, const double * x, double * grad, unsigned int & nPoints);
   void EvaluateLogLGradient(const IModelFunctionTempl<Double_v> &, const UnBinData &, const double *, double *, unsigned int & ) ;

   /**
       evaluate the Poisson LogL given a model function and the data at the point x.
       return also nPoints as the effective number of used points in the LogL evaluation
       By default is extended, pass extedend to false if want to be not extended (MultiNomial)
   */
   double EvaluatePoissonLogL(const IModelFunction & func, const BinData & data, const double * x, int iWeight, bool extended, unsigned int & nPoints);

   /**
       evaluate the Poisson LogL given a model function and the data at the point x.
       return also nPoints as the effective number of used points in the LogL evaluation
   */
   void EvaluatePoissonLogLGradient(const IModelFunction & func, const BinData & data, const double * x, double * grad);

//    /**
//        Parallel evaluate the Chi2 given a model function and the data at the point x.
//        return also nPoints as the effective number of used points in the Chi2 evaluation
//    */
//    double ParallelEvalChi2(const IModelFunction & func, const BinData & data, const double * x, unsigned int & nPoints);

   // methods required by dedicate minimizer like Fumili

   /**
       evaluate the residual contribution to the Chi2 given a model function and the BinPoint data
       and if the pointer g is not null evaluate also the gradient of the residual.
       If the function provides parameter derivatives they are used otherwise a simple derivative calculation
       is used
   */
   double EvaluateChi2Residual(const IModelFunction & func, const BinData & data, const double * x, unsigned int ipoint, double *g = 0);

   /**
       evaluate the pdf contribution to the LogL given a model function and the BinPoint data.
       If the pointer g is not null evaluate also the gradient of the pdf.
       If the function provides parameter derivatives they are used otherwise a simple derivative calculation
       is used
   */
   double EvaluatePdf(const IModelFunction & func, const UnBinData & data, const double * x, unsigned int ipoint, double * g = 0);
   double EvaluatePdf(const IModelFunctionTempl<Double_v> & func, const UnBinData & data, const double * p, unsigned int i, double *);
   /**
       evaluate the pdf contribution to the Poisson LogL given a model function and the BinPoint data.
       If the pointer g is not null evaluate also the gradient of the Poisson pdf.
       If the function provides parameter derivatives they are used otherwise a simple derivative calculation
       is used
   */
   double EvaluatePoissonBinPdf(const IModelFunction & func, const BinData & data, const double * x, unsigned int ipoint, double * g = 0);

   unsigned setAutomaticChunking(unsigned nEvents);

   template<class T>
   struct Evaluate{
      static double EvalChi2(const IModelFunctionTempl<T> &func, const BinData & data, const double * p, unsigned int &nPoints, const unsigned int &executionPolicy, unsigned nChunks = 0)
      {
         // evaluate the chi2 given a  vectorized function reference  , the data and returns the value and also in nPoints
         // the actual number of used points
         // normal chi2 using only error on values (from fitting histogram)
         // optionally the integral of function in the bin is used

         unsigned int n = data.Size();

         nPoints = 0; // count the effective non-zero points
         // set parameters of the function to cache integral value
#ifdef USE_PARAMCACHE
         (const_cast<IModelFunctionTempl<T> &>(func)).SetParameters(p);
#endif
         // do not cache parameter values (it is not thread safe)
         //func.SetParameters(p);


         // get fit option and check case if using integral of bins
         const DataOptions &fitOpt = data.Opt();
         if (fitOpt.fExpErrors || fitOpt.fIntegral || fitOpt.fExpErrors)
            Error("FitUtil::EvaluateChi2", "The vectorized implementation doesn't support Integrals, BinVolume or ExpErrors\n. Aborting operation.");

         (const_cast<IModelFunctionTempl<T> &>(func)).SetParameters(p);

         double maxResValue = std::numeric_limits<double>::max() / n;
         std::vector<double> ones{1, 1, 1, 1};
         auto vecSize = vecCore::VectorSize<T>();

         auto mapFunction = [&](unsigned int i) {
            // in case of no error in y invError=1 is returned
            T x, y, invErrorVec;
            vecCore::Load<T>(x, data.GetCoordComponent(i * vecSize, 0));
            vecCore::Load<T>(y, data.ValuePtr(i * vecSize));
            const auto invError = data.ErrorPtr(i * vecSize);
            auto invErrorptr = (invError != nullptr) ? invError : &ones.front();
            vecCore::Load<T>(invErrorVec, invErrorptr);

            T fval{};

#ifdef USE_PARAMCACHE
            fval = func(&x);
#else
            fval = func(&x, p);
#endif
            nPoints++;

            T tmp = (y - fval) * invErrorVec;
            T chi2 = tmp * tmp;


            // avoid inifinity or nan in chi2 values due to wrong function values
            auto m = vecCore::Mask_v<T>(chi2 > maxResValue);

            vecCore::MaskedAssign<T>(chi2, m, maxResValue);

            return chi2;
         };

         auto redFunction = [](const std::vector<T> &objs) {
            return std::accumulate(objs.begin(), objs.end(), T{});
         };

         T res{};
         if (executionPolicy == 0) {
            for (unsigned int i = 0; i < (data.Size() / vecSize); i++) {
               res += mapFunction(i);
            }
         } else if (executionPolicy == 1) {
            auto chunks = nChunks != 0 ? nChunks : setAutomaticChunking(data.Size() / vecSize);
            ROOT::TThreadExecutor pool;
            res = pool.MapReduce(mapFunction, ROOT::TSeq<unsigned>(0, data.Size() / vecSize), redFunction, chunks);
            // } else if(executionPolicy == 2){
            //   ROOT::TProcessExecutor pool;
            //   res = pool.MapReduce(mapFunction, ROOT::TSeq<unsigned>(0, data.Size()/vecSize), redFunction);
         } else {
            Error("FitUtil::EvaluateChi2", "Execution policy unknown. Avalaible choices:\n 0: Serial (default)\n 1: MultiThread\n 2: MultiProcess");
         }
         nPoints = n;

#ifdef DEBUG
         std::cout << "chi2 = " << chi2 << " n = " << nPoints  /*<< " rejected = " << nRejected */ << std::endl;
#endif

         return res.sum();
      }

      static double EvalLogL(const IModelFunctionTempl<T> &func, const UnBinData & data, const double * const p, int iWeight,
                             bool extended, unsigned int &nPoints, const unsigned int &executionPolicy, unsigned nChunks = 0)
      {
         // evaluate the LogLikelihood
         unsigned int n = data.Size();

         //unsigned int nRejected = 0;

         // set parameters of the function to cache integral value
#ifdef USE_PARAMCACHE
         (const_cast<IModelFunctionTempl<T> &>(func)).SetParameters(p);
#endif

         // this is needed if function must be normalized
         bool normalizeFunc = false;
         double norm = 1.0;
         if (normalizeFunc) {
            // compute integral of the function
            std::vector<double> xmin(data.NDim());
            std::vector<double> xmax(data.NDim());
            IntegralEvaluator<IModelFunctionTempl<T>> igEval(func, p, true);
            // compute integral in the ranges where is defined
            if (data.Range().Size() > 0) {
               norm = 0;
               for (unsigned int ir = 0; ir < data.Range().Size(); ++ir) {
                  data.Range().GetRange(&xmin[0], &xmax[0], ir);
                  norm += igEval.Integral(xmin.data(), xmax.data());
               }
            } else {
               // use (-inf +inf)
               data.Range().GetRange(&xmin[0], &xmax[0]);
               // check if funcition is zero at +- inf
               const auto xmin_v = vecCore::FromPtr<T>(xmin.data());
               const auto xmax_v = vecCore::FromPtr<T>(xmax.data());
               if (func(&xmin_v, p).sum() != 0 || func(&xmax_v, p).sum() != 0) {
                  MATH_ERROR_MSG("FitUtil::EvaluateLogLikelihood", "A range has not been set and the function is not zero at +/- inf");
                  return 0;
               }
               norm = igEval.Integral(&xmin[0], &xmax[0]);
            }
         }

         // needed to compute effective global weight in case of extended likelihood

         auto vecSize = vecCore::VectorSize<T>();

         auto mapFunction = [ &, p](const unsigned i) {
            T W{};
            T W2{};
            const auto x = vecCore::FromPtr<T>(data.GetCoordComponent(i * vecSize, 0));

#ifdef USE_PARAMCACHE
            auto fval = func(&x);
#else
            auto fval = func(&x, p);
#endif
            if (normalizeFunc) fval = fval * (1 / norm);

            // function EvalLog protects against negative or too small values of fval
            auto logval =  ROOT::Math::Util::EvalLog(fval);
            if (iWeight > 0) {
               auto weight = data.WeightsPtr(i) == nullptr ? 1 : vecCore::FromPtr<T>(data.WeightsPtr(i));
               logval *= weight;
               if (iWeight == 2) {
                  logval *= weight; // use square of weights in likelihood
                  if (!extended) {
                     // needed sum of weights and sum of weight square if likelkihood is extended
                     W  = weight;
                     W2 = weight * weight;
                  }
               }
            }
            nPoints++;
            return LikelihoodAux<T>(logval, W, W2);
         };

         auto redFunction = [](const std::vector<LikelihoodAux<T>> &objs) {
            return std::accumulate(objs.begin(), objs.end(), LikelihoodAux<T>(),
            [](const LikelihoodAux<T> &l1, const LikelihoodAux<T> &l2) {
               return l1 + l2;
            });
         };

         T logl_v{};
         T sumW_v{};
         T sumW2_v{};

         if (executionPolicy == 0) {
            for (unsigned int i = 0; i < n / vecSize; ++i) {
               auto resArray = mapFunction(i);
               logl_v += resArray.logvalue;
               sumW_v += resArray.weight;
               sumW2_v += resArray.weight2;
            }
         } else if (executionPolicy == 1) {
            auto chunks = nChunks != 0 ? nChunks : setAutomaticChunking(data.Size() / vecSize);
            ROOT::TThreadExecutor pool(4);
            auto resArray = pool.MapReduce(mapFunction, ROOT::TSeq<unsigned>(0, data.Size() / vecSize), redFunction, chunks);
            logl_v = resArray.logvalue;
            sumW_v = resArray.weight;
            sumW2_v = resArray.weight2;
         } else if (executionPolicy == 2) {
            // ROOT::TProcessExecutor pool;
            // res = pool.MapReduce(mapFunction, ROOT::TSeq<unsigned>(0, n), redFunction);
         } else {
            Error("FitUtil::EvaluateLogL", "Execution policy unknown. Avalaible choices:\n 0: Serial (default)\n 1: MultiThread\n 2: MultiProcess");
         }

         //reduce vector type to double.
         double logl  = 0.;
         double sumW  = 0.;
         double sumW2 = 0;;

         for (unsigned vIt = 0; vIt < vecSize; vIt++) {
            logl += logl_v[vIt];
            sumW += sumW_v[vIt];
            sumW2 += sumW2_v[vIt];
         }

         if (extended) {
            // add Poisson extended term
            double extendedTerm = 0; // extended term in likelihood
            double nuTot = 0;
            // nuTot is integral of function in the range
            // if function has been normalized integral has been already computed
            if (!normalizeFunc) {
               IntegralEvaluator<IModelFunctionTempl<T>> igEval(func, p, true);
               std::vector<double> xmin(data.NDim());
               std::vector<double> xmax(data.NDim());

               // compute integral in the ranges where is defined
               if (data.Range().Size() > 0) {
                  nuTot = 0;
                  for (unsigned int ir = 0; ir < data.Range().Size(); ++ir) {
                     data.Range().GetRange(&xmin[0], &xmax[0], ir);
                     nuTot += igEval.Integral(xmin.data(), xmax.data());
                  }
               } else {
                  // use (-inf +inf)
                  data.Range().GetRange(&xmin[0], &xmax[0]);
                  // check if funcition is zero at +- inf
                  const auto xmin_v = vecCore::FromPtr<T>(xmin.data());
                  const auto xmax_v = vecCore::FromPtr<T>(xmax.data());
                  if (func(&xmin_v, p).sum() != 0 || func(&xmax_v, p).sum() != 0) {
                     MATH_ERROR_MSG("FitUtil::EvaluateLogLikelihood", "A range has not been set and the function is not zero at +/- inf");
                     return 0;
                  }
                  nuTot = igEval.Integral(&xmin[0], &xmax[0]);
               }

               // force to be last parameter value
               //nutot = p[func.NDim()-1];
               if (iWeight != 2)
                  extendedTerm = - nuTot;  // no need to add in this case n log(nu) since is already computed before
               else {
                  // case use weight square in likelihood : compute total effective weight = sw2/sw
                  // ignore for the moment case when sumW is zero
                  extendedTerm = - (sumW2 / sumW) * nuTot;
               }

            } else {
               nuTot = norm;
               extendedTerm = - nuTot + double(n) *  ROOT::Math::Util::EvalLog(nuTot);
               // in case of weights need to use here sum of weights (to be done)
            }

            logl += extendedTerm;
         }

         // reset the number of fitting data points
         //  nPoints = n;
// std::cout<<", n: "<<nPoints<<std::endl;
         nPoints = 0;
         return -logl;

      }

      static double EvalChi2Effective(const IModelFunctionTempl<T> &func, const BinData & data, const double * p, unsigned int &nPoints)
      {
         Error("FitUtil::Evaluate<T>::EvalChi2Effective", "The vectorized evaluation of the Chi2 with coordinate errors is still not supported");
         return -1.;
      }
   };

   template<>
   struct Evaluate<double>{
      static double EvalChi2(const IModelFunction & func, const BinData & data, const double * p, unsigned int &nPoints, const unsigned int &executionPolicy, unsigned nChunks = 0)
      {
         // evaluate the chi2 given a  function reference, the data and returns the value and also in nPoints
         // the actual number of used points
         // normal chi2 using only error on values (from fitting histogram)
         // optionally the integral of function in the bin is used
         return FitUtil::EvaluateChi2(func, data, p, nPoints, executionPolicy, nChunks);
      }
      static double EvalLogL(const IModelFunctionTempl<double> &func, const UnBinData & data, const double * p, int iWeight,
      bool extended, unsigned int &nPoints, const unsigned int &executionPolicy, unsigned nChunks = 0)
      {
         return FitUtil::EvaluateLogL(func, data, p, iWeight, extended, nPoints, executionPolicy, nChunks);
      }
      static double EvalChi2Effective(const IModelFunctionTempl<double> &func, const BinData & data, const double * p, unsigned int &nPoints)
      {
         return FitUtil::EvaluateChi2Effective(func, data, p, nPoints);
      }
   };


} // end namespace FitUtil

   } // end namespace Fit

} // end namespace ROOT

Vc_DECLARE_ALLOCATOR(ROOT::Fit::FitUtil::LikelihoodAux<Double_v>);

#endif /* ROOT_Fit_FitUtil */
