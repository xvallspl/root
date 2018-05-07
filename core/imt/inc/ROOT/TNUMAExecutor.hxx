#ifndef ROOT_TNUMAExecutor
#define ROOT_TNUMAExecutor

#include "ROOT/TProcessExecutor.hxx"
#include "ROOT/TThreadExecutorImpl.hxx"
#include "ROOT/RArrayView.hxx" // std::array_view

#include <numa.h>
#include <algorithm> // std::min, std::max
#include <thread> //std::ma
namespace ROOT {
namespace Experimental {

class TNUMAExecutor {
public:

   template< class F, class... T>
   using noReferenceCond = typename std::enable_if<"Function can't return a reference" && !(std::is_reference<typename std::result_of<F(T...)>::type>::value)>::type;

   explicit TNUMAExecutor():TNUMAExecutor(ROOT::Internal::TPoolManager::GetPoolSize()==0? std::thread::hardware_concurrency() : ROOT::Internal::TPoolManager::GetPoolSize()){}
   explicit TNUMAExecutor(unsigned nThreads):fNDomains(numa_max_node()+1)
   {
      fDomainNThreads = nThreads/fNDomains;
   }

   unsigned GetNUMADomains(){
      return fNDomains;
   }

   template<class F>
   void Foreach(F func, unsigned nTimes);
   template<class F, class INTEGER>
   void Foreach(F func, ROOT::TSeq<INTEGER> args);
   /// \cond
   template<class F, class T>
   void Foreach(F func, std::initializer_list<T> args);
   /// \endcond
   template<class F, class T>
   void Foreach(F func, std::vector<T> &args);

   template<class F, class R, class Cond = noReferenceCond<F>>
   auto MapReduce(F func, unsigned nTimes, R redfunc, unsigned nChunks = 0) -> typename std::result_of<F()>::type;
   template<class F, class INTEGER, class R, class Cond = noReferenceCond<F, INTEGER>>
   auto MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks = 0) -> typename std::result_of<F(INTEGER)>::type;
   template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
   auto MapReduce(F func, std::initializer_list<T> args, R redfunc, unsigned nChunks = 0) -> typename std::result_of<F(T)>::type
   {
      std::vector<T> vargs(std::move(args));
      return MapReduce(func, vargs, redfunc, nChunks);
   }
   template<class F, class T, class R, class Cond = noReferenceCond<F, T>>
   auto MapReduce(F func, std::vector<T> &args, R redfunc, unsigned nChunks = 0) -> typename std::result_of<F(T)>::type;

private:

   template<class T>
   std::vector<std::array_view<T>> splitData(std::vector<T> &vec);

   unsigned fNDomains{};
   unsigned fDomainNThreads{};
};


template<class T>
std::vector<std::array_view<T>> TNUMAExecutor::splitData(std::vector<T> &vec)
{
   unsigned int nToProcess = vec.size();
   unsigned stride = (nToProcess + fNDomains - 1) / fNDomains; //ceiling the division
   auto av = std::make_view(vec);
   std::vector<std::array_view<T>> v;
   unsigned i;
   for(i=0; i*stride<av.size()-stride; i++) {
      v.emplace_back(av.slice(av.begin()+ i*stride, av.begin()+(i+1)*stride));
   }
   v.emplace_back(av.slice(av.begin() + i*stride, av.end()));

   return v;
}


template<class F>
void TNUMAExecutor::Foreach(F func, unsigned nTimes) {
      TThreadExecutorImpl pool(fDomainNThreads);
      pool.Foreach(func, nTimes);
}

//////////////////////////////////////////////////////////////////////////
/// Execute func in parallel, taking an element of a
/// sequence as argument.
template<class F, class INTEGER>
void TNUMAExecutor::Foreach(F func, ROOT::TSeq<INTEGER> args) {
      TThreadExecutorImpl pool{fDomainNThreads};
      pool.Foreach(func, args);
}

/// \cond
//////////////////////////////////////////////////////////////////////////
/// Execute func in parallel, taking an element of a
/// initializer_list as argument.
template<class F, class T>
void TNUMAExecutor::Foreach(F func, std::initializer_list<T> args) {
    TThreadExecutorImpl pool{fDomainNThreads};
    std::vector<T> vargs(std::move(args));
    pool.Foreach(func, vargs);
}
/// \endcond

//////////////////////////////////////////////////////////////////////////
/// Execute func in parallel, taking an element of an
/// std::vector as argument.
template<class F, class T>
void TNUMAExecutor::Foreach(F func, std::vector<T> &args) {
      TThreadExecutorImpl pool{fDomainNThreads};
      pool.Foreach(func, args);
}

template<class F, class R, class Cond>
auto TNUMAExecutor::MapReduce(F func, unsigned nTimes, R redfunc, unsigned nChunks) -> typename std::result_of<F()>::type
{
   auto runOnNode = [&](unsigned int i) {
      numa_run_on_node(i);
      ROOT::TThreadExecutorImpl pool{fDomainNThreads};
      auto res = nChunks? pool.MapReduce(func, nTimes, redfunc, nChunks/fNDomains) : pool.MapReduce(func, nTimes, redfunc); 
      numa_run_on_node_mask(numa_all_nodes_ptr);
      return res;
   };

   ROOT::TProcessExecutor proc(fNDomains);
   return proc.MapReduce(runOnNode, ROOT::TSeq<unsigned>(fNDomains), redfunc);
}

template<class F, class T, class R, class Cond>
auto TNUMAExecutor::MapReduce(F func, std::vector<T> &args, R redfunc, unsigned nChunks) -> typename std::result_of<F(T)>::type
{
   auto dataRanges = splitData(args);
   auto runOnNode = [&](unsigned int i) {
      numa_run_on_node(i);
      ROOT::TThreadExecutorImpl pool{fDomainNThreads};
      auto res = nChunks? pool.MapReduce(func, dataRanges[i], redfunc, nChunks/fNDomains) : pool.MapReduce(func, dataRanges[i], redfunc); 
      numa_run_on_node_mask(numa_all_nodes_ptr);
      return res;
   };

   ROOT::TProcessExecutor proc(fNDomains);
   return proc.MapReduce(runOnNode, ROOT::TSeq<unsigned>(fNDomains), redfunc);
}



template<class F, class INTEGER, class R, class Cond>
auto TNUMAExecutor::MapReduce(F func, ROOT::TSeq<INTEGER> args, R redfunc, unsigned nChunks) -> typename std::result_of<F(INTEGER)>::type
{
   unsigned stride = (*args.end() - *args.begin() + fNDomains - 1) / fNDomains; //ceiling the division
   auto runOnNode = [&](unsigned int i) {
      numa_run_on_node(i);
      ROOT::TThreadExecutorImpl pool{fDomainNThreads};
      ROOT::TSeq<unsigned> sequence(std::max(*args.begin(), i*stride), std::min((i+1)*stride, *args.end()));
      auto res =  pool.MapReduce(func, sequence, redfunc, nChunks/fNDomains);
      numa_run_on_node_mask(numa_all_nodes_ptr);
      return res;
   };

   ROOT::TThreadExecutorImpl proc(fNDomains);
   return proc.MapReduce(runOnNode, ROOT::TSeq<unsigned>(fNDomains), redfunc);
}

} // namespace Experimental
} // namespace ROOT

#endif
