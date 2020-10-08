// @(#)root/thread:$Id$
// Author: Xavier Valls September 2020

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TExecutor
#define ROOT_TExecutor

#include <ROOT/RMakeUnique.hxx>
#include "ROOT/TExecutorCRTP.hxx"
#include "ROOT/TSequentialExecutor.hxx"
#ifdef R__USE_IMT
#include "ROOT/TThreadExecutor.hxx"
#endif
#include "ROOT/TProcessExecutor.hxx"
#include "TROOT.h"
#include "ExecutionPolicy.hxx"
#include <memory>
#include <stdexcept>
#include <thread>

//////////////////////////////////////////////////////////////////////////
///
/// \class ROOT::Internal::TExecutor
/// \brief This class defines an interface to execute the same task
/// multiple times, sequentially or in parallel depending on the execution policy passed
/// as a first parameter on construction, and possibly with different arguments every time.
/// The classes implementing it mimic the behaviour of python's pool.Map method.
///
/// ###ROOT::Internal::TExecutor::Map
/// The two possible usages of the Map method are:\n
/// * `Map(F func, unsigned nTimes)`: func is executed nTimes with no arguments
/// * `Map(F func, T& args)`: func is executed on each element of the collection of arguments args
///
/// For either signature, func is executed as many times as needed by a pool of
/// nThreads threads; It defaults to the number of cores.\n
/// A collection containing the result of each execution is returned.\n
/// **Note:** the user is responsible for the deletion of any object that might
/// be created upon execution of func, returned objects included. ROOT::::Internal::TExecutor never
/// deletes what it returns, it simply forgets it.\n
///
/// \param func
/// \parblock
/// a lambda expression, an std::function, a loaded macro, a
/// functor class or a function that takes zero arguments (for the first signature)
/// or one (for the second signature).
/// \endparblock
/// \param args
/// \parblock
/// a standard vector, a ROOT::TSeq of integer type or an initializer list for the second signature.
/// An integer only for the first.\n
/// \endparblock
///
/// **Note:** in cases where the function to be executed takes more than
/// zero/one argument but all are fixed except zero/one, the function can be wrapped
/// in a lambda or via std::bind to give it the right signature.\n
///
/// #### Return value:
/// An std::vector. The elements in the container
/// will be the objects returned by func.
///
/// ### ROOT::Internal::TExecutor::MapReduce
/// This set of methods behaves exactly like Map, but takes an additional
/// function as a third argument. This function is applied to the set of
/// objects returned by the corresponding Map execution to "squash" them
/// to a single object. This function should be independent of the size of
/// the vector returned by Map due to optimization of the number of chunks.
///
/// #### Examples:
/// ~~~{.cpp}
/// root[] ROOT::Internal::TExecutor pool; auto ten = pool.MapReduce([]() { return 1; }, 10, [](std::vector<int> v) { return std::accumulate(v.begin(), v.end(), 0); })
/// root[] ROOT::Internal::TExecutor pool(ROOT::Internal::ExecutionPolicy::kMultiprocess); auto hist = pool.MapReduce(CreateAndFillHists, 10, PoolUtils::ReduceObjects);
/// ~~~
///
//////////////////////////////////////////////////////////////////////////


namespace ROOT{

namespace Internal{
class TExecutor: public TExecutorCRTP<TExecutor> {
public:

   /// \brief Class constructor. Sets the default execution policy and initializes the corresponding executor.
   /// Defaults to multithreaded execution policy if ROOT is compiled with IMT=ON and IsImplicitMTEnabled. Otherwise it defaults to a serial execution policy
   /// \param nProcessingUnits [optional] Number of parallel processing units, only taken into account if the execution policy is kMultithread
   explicit TExecutor(unsigned nProcessingUnits = 0) :
      TExecutor(ROOT::IsImplicitMTEnabled() ? ROOT::Internal::ExecutionPolicy::kMultithread : ROOT::Internal::ExecutionPolicy::kSerial, nProcessingUnits) {}

   /// \brief Class constructor. Sets the execution policy and initializes the corresponding executor.
   /// \param execPolicy Execution policy(kMultithread, kMultiprocess, kSerial) to process the data
   /// \param nProcessingUnits [optional] Number of parallel processing units, only taken into account if the execution policy is kMultithread
   explicit TExecutor(ROOT::Internal::ExecutionPolicy execPolicy, unsigned nProcessingUnits = 0) : fExecPolicy(execPolicy) {
      fExecPolicy = execPolicy;
      switch(fExecPolicy) {
         case ROOT::Internal::ExecutionPolicy::kSerial:
            fSeqPool = std::make_unique<ROOT::TSequentialExecutor>();
            break;
#ifdef R__USE_IMT
         case ROOT::Internal::ExecutionPolicy::kMultithread:
            fThreadPool = std::make_unique<ROOT::TThreadExecutor>(nProcessingUnits);
            break;
#endif
         case ROOT::Internal::ExecutionPolicy::kMultiprocess:
            fProcPool = std::make_unique<ROOT::TProcessExecutor>(nProcessingUnits);
            break;
         default:
            throw std::invalid_argument("kMultithread policy not available when ROOT is compiled with IMT=OFF.");
      }
   }

TExecutor(const TExecutor &) = delete;
TExecutor &operator=(const TExecutor &) = delete;

/// Return the execution policy the executor is set to
ROOT::Internal::ExecutionPolicy Policy(){ return fExecPolicy; }

using TExecutorCRTP<TExecutor>::Map;
template<class F, class Cond = noReferenceCond<F>>
auto Map(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type>;
template<class F, class INTEGER, class Cond = noReferenceCond<F, INTEGER>>
auto Map(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type>;
template<class F, class T, class Cond = noReferenceCond<F, T>>
auto Map(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>;
template<class F, class T, class Cond = noReferenceCond<F, T>>
auto Map(F func, const std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type>;

// // MapReduce
// // the late return types also check at compile-time whether redfunc is compatible with func,
// // other than checking that func is compatible with the type of arguments.
// // a static_assert check in TExecutor::Reduce is used to check that redfunc is compatible with the type returned by func
using TExecutorCRTP<TExecutor>::MapReduce;
template<class F, class R, class Cond = noReferenceCond<F>>
auto MapReduce(F func, unsigned nTimes, R redfunc, unsigned nChunks) -> typename std::result_of<F()>::type;
template<class F, class INTEGER, class R, class Cond = noReferenceCond<F, INTEGER>>
auto MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(INTEGER)>::type;
template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
auto MapReduce(F func, std::initializer_list<T> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type;
template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
auto MapReduce(F func, std::vector<T> &args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type;
template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
auto MapReduce(F func, const std::vector<T> &args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type;

using TExecutorCRTP<TExecutor>::Reduce;

unsigned GetPoolSize();

protected:
template<class F, class R, class Cond = noReferenceCond<F>>
auto Map(F func, unsigned nTimes, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F()>::type>;
template<class F, class INTEGER, class R, class Cond = noReferenceCond<F, INTEGER>>
auto Map(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(INTEGER)>::type>;
template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
auto Map(F func, std::vector<T> &args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(T)>::type>;
template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
auto Map(F func, const std::vector<T> &args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(T)>::type>;
template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
auto Map(F func, std::initializer_list<T> args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(T)>::type>;

private:
   ROOT::Internal::ExecutionPolicy fExecPolicy;
#ifdef R__USE_IMT
   std::unique_ptr<ROOT::TThreadExecutor> fThreadPool;
#endif
   std::unique_ptr<ROOT::TProcessExecutor> fProcPool;
   std::unique_ptr<ROOT::TSequentialExecutor> fSeqPool;
};


//////////////////////////////////////////////////////////////////////////
/// \copydoc TExecutorCRTP::Map(F func,unsigned nTimes)
template<class F, class Cond>
auto TExecutor::Map(F func, unsigned nTimes) -> std::vector<typename std::result_of<F()>::type> {
   using retType = decltype(func());
   std::vector<retType> res;
   switch(fExecPolicy){
      case ROOT::Internal::ExecutionPolicy::kSerial:
         res = fSeqPool->Map(func, nTimes);
         break;
#ifdef R__USE_IMT
      case ROOT::Internal::ExecutionPolicy::kMultithread:
         res = fThreadPool->Map(func, nTimes);
         break;
#endif
      case ROOT::Internal::ExecutionPolicy::kMultiprocess:
         res = fProcPool->Map(func, nTimes);
         break;
      default:
         break;
   }
   return res;
}

//////////////////////////////////////////////////////////////////////////
/// \copydoc TExecutorCRTP::Map(F func,ROOT::TSeq<INTEGER> args)
template<class F, class INTEGER, class Cond>
auto TExecutor::Map(F func, ROOT::TSeq<INTEGER> args) -> std::vector<typename std::result_of<F(INTEGER)>::type> {
   using retType = decltype(func(args.front()));
   std::vector<retType> res;

   switch(fExecPolicy){
      case ROOT::Internal::ExecutionPolicy::kSerial:
         res = fSeqPool->Map(func, args);
         break;
#ifdef R__USE_IMT
      case ROOT::Internal::ExecutionPolicy::kMultithread:
         res = fThreadPool->Map(func, args);
         break;
#endif
      case ROOT::Internal::ExecutionPolicy::kMultiprocess:
         res = fProcPool->Map(func, args);
         break;
      default:
         break;
   }
   return res;
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute func (with no arguments) nTimes, dividing the execution in nChunks and providing a result per chunk if
/// the execution policy is multithreaded. Otherwise, it ignores the two last arguments and performs a normal Map operation.
///
/// \param func Function to be executed.
/// \param nTimes Number of times function should be called.
/// \param redfunc Reduction function, used both to generate the partial results and the end result. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A vector with the results of the function calls.
template<class F, class R, class Cond>
auto TExecutor::Map(F func, unsigned nTimes, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F()>::type> {
#ifdef R__USE_IMT
   if (fExecPolicy == ROOT::Internal::ExecutionPolicy::kMultithread) {
      return fThreadPool->Map(func, nTimes, redfunc, nChunks);
   }
#endif
   return Map(func, nTimes);
}

//////////////////////////////////////////////////////////////////////////
/// \copydoc TExecutorCRTP::Map(F func,std::vector<T> &args)
template<class F, class T, class Cond>
auto TExecutor::Map(F func, std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type> {
   // //check whether func is callable
   using retType = decltype(func(args.front()));
   std::vector<retType> res;
   switch(fExecPolicy){
      case ROOT::Internal::ExecutionPolicy::kSerial:
         res = fSeqPool->Map(func, args);
         break;
#ifdef R__USE_IMT
      case ROOT::Internal::ExecutionPolicy::kMultithread:
         res = fThreadPool->Map(func, args);
         break;
#endif
      case ROOT::Internal::ExecutionPolicy::kMultiprocess:
         res = fProcPool->Map(func, args);
         break;
      default:
         break;
   }
   return res;
}

//////////////////////////////////////////////////////////////////////////
/// \copydoc TExecutorCRTP::Map(F func,const std::vector<T> &args)
template<class F, class T, class Cond>
auto TExecutor::Map(F func, const std::vector<T> &args) -> std::vector<typename std::result_of<F(T)>::type> {
   // //check whether func is callable
   using retType = decltype(func(args.front()));
   std::vector<retType> res;
   switch(fExecPolicy){
      case ROOT::Internal::ExecutionPolicy::kSerial:
         res = fSeqPool->Map(func, args);
         break;
#ifdef R__USE_IMT
      case ROOT::Internal::ExecutionPolicy::kMultithread:
         res = fThreadPool->Map(func, args);
         break;
#endif
      case ROOT::Internal::ExecutionPolicy::kMultiprocess:
         res = fProcPool->Map(func, args);
         break;
      default:
         break;
   }
   return res;
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over a sequence of indexes, dividing the execution in nChunks and providing a result per chunk if
/// the execution policy is multithreaded. Otherwise, it ignores the two last arguments and performs a normal Map operation.
///
/// \param func Function to be executed. Must take an element of the sequence passed assecond argument as a parameter.
/// \param args Sequence of indexes to execute `func` on.
/// \param redfunc Reduction function, used to combine the results of the calls to `func` into partial results. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A vector with the results of the function calls.
template<class F, class INTEGER, class R, class Cond>
auto TExecutor::Map(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(INTEGER)>::type> {
#ifdef R__USE_IMT
   if (fExecPolicy == ROOT::Internal::ExecutionPolicy::kMultithread) {
      return fThreadPool->Map(func, args, redfunc, nChunks);
   }
#endif
   return Map(func, args);
}



//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of a vector, dividing the execution in nChunks and providing a result per chunk if
/// the execution policy is multithreaded. Otherwise, it ignores the two last arguments and performs a normal Map operation.
///
/// \param func Function to be executed on the elements of the vector passed as second parameter.
/// \param args Vector of elements passed as an argument to `func`.
/// \param redfunc Reduction function, used to combine the results of the calls to `func` into partial results. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A vector with the results of the function calls.
template<class F, class T, class R, class Cond>
auto TExecutor::Map(F func, std::vector<T> &args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(T)>::type> {
#ifdef R__USE_IMT
   if (fExecPolicy == ROOT::Internal::ExecutionPolicy::kMultithread) {
      return fThreadPool->Map(func, args, redfunc, nChunks);
   }
#endif
   return Map(func, args);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of an immutable vector, dividing the execution in nChunks and providing a result per chunk if
/// the execution policy is multithreaded. Otherwise, it ignores the two last arguments and performs a normal Map operation.
///
/// \param func Function to be executed on the elements of the vector passed as second parameter.
/// \param args Immutable vector of elements passed as an argument to `func`.
/// \param redfunc Reduction function, used to combine the results of the calls to `func` into partial results. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A vector with the results of the function calls.
template<class F, class T, class R, class Cond>
auto TExecutor::Map(F func, const std::vector<T> &args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(T)>::type> {
#ifdef R__USE_IMT
   if (fExecPolicy == ROOT::Internal::ExecutionPolicy::kMultithread) {
      return fThreadPool->Map(func, args, redfunc, nChunks);
   }
#endif
   return Map(func, args);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of an initializer_list, dividing the execution in nChunks and providing a result per chunk if
/// the execution policy is multithreaded. Otherwise, it ignores the two last arguments and performs a normal Map operation.
///
/// \param func Function to be executed on the elements of the initializer_list passed as second parameter.
/// \param args initializer_list for a vector to apply `func` on.
/// \param redfunc Reduction function, used to combine the results of the calls to `func` into partial results. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A vector with the results of the function calls.
template<class F, class T, class R, class Cond>
auto TExecutor::Map(F func, std::initializer_list<T> args, R redfunc, unsigned nChunks) -> std::vector<typename std::result_of<F(T)>::type> {
   std::vector<T> vargs(std::move(args));
   const auto &reslist = Map(func, vargs, redfunc, nChunks);
   return reslist;
}


//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function `nTimes` (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results if the execution policy is multithreaded.
/// Otherwise, it ignores the two last arguments and performs a normal Map operation.
///
/// \param func Function to be executed. Must take an element of the sequence passed as second argument as a parameter.
/// \param nTimes Number of times function should be called.
/// \param redfunc Reduction function to combine the results of the calls to `func` into partial results, and these
/// into a final result. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A value result of "reducing" the vector returned by the Map operation into a single object.
template<class F, class R, class Cond>
auto TExecutor::MapReduce(F func, unsigned nTimes, R redfunc, unsigned nChunks) -> typename std::result_of<F()>::type {
   return Reduce(Map(func, nTimes, redfunc, nChunks), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over a sequence of indexes (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results if the execution policy is multithreaded.
/// Otherwise, it ignores the two last arguments and performs a normal Map operation.
///
/// \param func Function to be executed. Must take an element of the sequence passed assecond argument as a parameter.
/// \param args Sequence of indexes to execute `func` on.
/// \param redfunc Reduction function to combine the results of the calls to `func` into partial results, and these
/// into a final result. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A value result of "reducing" the vector returned by the Map operation into a single object.
template<class F, class INTEGER, class R, class Cond>
auto TExecutor::MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(INTEGER)>::type {
   return Reduce(Map(func, args, redfunc, nChunks), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of an initializer_list (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results if the execution policy is multithreaded.
/// Otherwise, it ignores the two last arguments and performs a normal Map operation.
///
/// \param func Function to be executed. Must take an element of the sequence passed assecond argument as a parameter.
/// \param args initializer_list for a vector to apply `func` on.
/// \param redfunc Reduction function to combine the results of the calls to `func` into partial results, and these
/// into a final result. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A value result of "reducing" the vector returned by the Map operation into a single object.
template<class F, class T, class R, class Cond>
auto TExecutor::MapReduce(F func, std::initializer_list<T> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type {
   return Reduce(Map(func, args, redfunc, nChunks), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of a vector (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results if the execution policy is multithreaded.
/// Otherwise, it ignores the two last arguments and performs a normal Map operation.
///
/// \param func Function to be executed. Must take an element of the sequence passed assecond argument as a parameter.
/// \param args Vector of elements passed as an argument to `func`.
/// \param redfunc Reduction function to combine the results of the calls to `func` into partial results, and these
/// into a final result. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A value result of "reducing" the vector returned by the Map operation into a single object.
template<class F, class T, class R, class Cond>
auto TExecutor::MapReduce(F func, std::vector<T> &args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type {
   return Reduce(Map(func, args, redfunc, nChunks), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Execute a function over the elements of an immutable vector (Map) and accumulate the results into a single value (Reduce).
/// Benefits from partial reduction into `nChunks` intermediate results if the execution policy is multithreaded.
/// Otherwise, it ignores the two last arguments and performs a normal Map operation.
///
/// \param func Function to be executed. Must take an element of the sequence passed assecond argument as a parameter.
/// \param args Immutable vector, whose elements are passed as an argument to `func`.
/// \param redfunc Reduction function to combine the results of the calls to `func` into partial results, and these
/// into a final result. Must return the same type as `func`.
/// \param nChunks Number of chunks to split the input data for processing.
/// \return A value result of "reducing" the vector returned by the Map operation into a single object.
template<class F, class T, class R, class Cond>
auto TExecutor::MapReduce(F func, const std::vector<T> &args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type {
   return Reduce(Map(func, args, redfunc, nChunks), redfunc);
}

//////////////////////////////////////////////////////////////////////////
/// \brief Return the number of pooled workers.
///
/// \return The number of workers in the pool in the executor used as a backend.

unsigned TExecutor::GetPoolSize()
{
   unsigned poolSize{0u};
   switch(fExecPolicy){
      case ROOT::Internal::ExecutionPolicy::kSerial:
         poolSize = fSeqPool->GetPoolSize();
         break;
#ifdef R__USE_IMT
      case ROOT::Internal::ExecutionPolicy::kMultithread:
         poolSize = fThreadPool->GetPoolSize();
         break;
#endif
      case ROOT::Internal::ExecutionPolicy::kMultiprocess:
         poolSize = fProcPool->GetPoolSize();
         break;
      default:
         break;
   }
   return poolSize;
}

} // namespace Internal
} // namespace ROOT

#endif
