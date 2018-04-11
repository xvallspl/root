############################################################################
# CMakeLists.txt file for building ROOT core/imt package
############################################################################

set(headers ROOT/TTaskGroup.hxx)
set(sources base.cxx TTaskGroup.cxx)

if (imt)
  set(headers ${headers} ROOT/TPoolManager.hxx ROOT/TThreadExecutor.hxx ROOT/TThreadExecutorImpl.hxx ROOT/TFuture.hxx ROOT/TNUMAExecutor.hxx)
  set(sources ${sources} TImplicitMT.cxx TThreadExecutorImpl.cxx TPoolManager.cxx G__Imt.cxx)
  ROOT_GENERATE_DICTIONARY(G__Imt ${headers} STAGE1
                           MODULE Imt LINKDEF LinkDef.h
                           DEPENDENCIES Core Thread BUILTINS TBB)
endif()

include_directories(SYSTEM ${TBB_INCLUDE_DIRS})

ROOT_LINKER_LIBRARY(Imt ${sources}
                    LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ${TBB_LIBRARIES}
                    DEPENDENCIES Core Thread numa MultiProc BUILTINS TBB)

ROOT_INSTALL_HEADERS(${installoptions})

if(imt)
    ROOT_ADD_TEST_SUBDIRECTORY(test)
endif()
