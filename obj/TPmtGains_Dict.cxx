// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TPmtGains_Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "TPmtGains.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TPmtGains(void *p = 0);
   static void *newArray_TPmtGains(Long_t size, void *p);
   static void delete_TPmtGains(void *p);
   static void deleteArray_TPmtGains(void *p);
   static void destruct_TPmtGains(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TPmtGains*)
   {
      ::TPmtGains *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TPmtGains >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TPmtGains", ::TPmtGains::Class_Version(), "TPmtGains.hxx", 14,
                  typeid(::TPmtGains), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TPmtGains::Dictionary, isa_proxy, 4,
                  sizeof(::TPmtGains) );
      instance.SetNew(&new_TPmtGains);
      instance.SetNewArray(&newArray_TPmtGains);
      instance.SetDelete(&delete_TPmtGains);
      instance.SetDeleteArray(&deleteArray_TPmtGains);
      instance.SetDestructor(&destruct_TPmtGains);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TPmtGains*)
   {
      return GenerateInitInstanceLocal((::TPmtGains*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TPmtGains*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TPmtGains::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TPmtGains::Class_Name()
{
   return "TPmtGains";
}

//______________________________________________________________________________
const char *TPmtGains::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPmtGains*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TPmtGains::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPmtGains*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TPmtGains::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPmtGains*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TPmtGains::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPmtGains*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TPmtGains::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPmtGains.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TPmtGains::Class(),this);
   } else {
      R__b.WriteClassBuffer(TPmtGains::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TPmtGains(void *p) {
      return  p ? new(p) ::TPmtGains : new ::TPmtGains;
   }
   static void *newArray_TPmtGains(Long_t nElements, void *p) {
      return p ? new(p) ::TPmtGains[nElements] : new ::TPmtGains[nElements];
   }
   // Wrapper around operator delete
   static void delete_TPmtGains(void *p) {
      delete ((::TPmtGains*)p);
   }
   static void deleteArray_TPmtGains(void *p) {
      delete [] ((::TPmtGains*)p);
   }
   static void destruct_TPmtGains(void *p) {
      typedef ::TPmtGains current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TPmtGains

namespace {
  void TriggerDictionaryInitialization_TPmtGains_Dict_Impl() {
    static const char* headers[] = {
"TPmtGains.hxx",
0
    };
    static const char* includePaths[] = {
"/usr/local/root/include",
"/.",
"/home/admin/root-6.08.00/include",
"/home/gold/captain/pds/pdsAnalysis/obj/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TPmtGains_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TPmtGains.hxx")))  TPmtGains;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TPmtGains_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TPmtGains.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TPmtGains", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TPmtGains_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TPmtGains_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TPmtGains_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TPmtGains_Dict() {
  TriggerDictionaryInitialization_TPmtGains_Dict_Impl();
}
