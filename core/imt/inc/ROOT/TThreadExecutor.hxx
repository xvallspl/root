// @(#)root/thread:$Id$
// Author: Xavier Valls March 2016

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TThreadExecutor
#define ROOT_TThreadExecutor

#include "RConfigure.h"

// exclude in case ROOT does not have IMT support
#ifndef R__USE_IMT
// No need to error out for dictionaries.
# if !defined(__ROOTCLING__) && !defined(G__DICTIONARY)
#  error "Cannot use ROOT::TThreadExecutorImpl without defining R__USE_IMT."
# endif
#else

#include "ROOT/TNUMAExecutor.hxx"
#include "ROOT/TThreadExecutorImpl.hxx"
#define NUMAExecutor
namespace ROOT {

   #ifndef NUMAExecutor
      using TThreadExecutor = TThreadExecutorImpl;
   #else
      using TThreadExecutor = ::ROOT::Experimental::TNUMAExecutor;
   #endif
} // namespace ROOT

#endif   // R__USE_IMT
#endif
