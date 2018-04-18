// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TPdsSummary_Dict

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
#include "TPdsSummary.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TPdsSummary(void *p = 0);
   static void *newArray_TPdsSummary(Long_t size, void *p);
   static void delete_TPdsSummary(void *p);
   static void deleteArray_TPdsSummary(void *p);
   static void destruct_TPdsSummary(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TPdsSummary*)
   {
      ::TPdsSummary *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TPdsSummary >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TPdsSummary", ::TPdsSummary::Class_Version(), "TPdsSummary.hxx", 28,
                  typeid(::TPdsSummary), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TPdsSummary::Dictionary, isa_proxy, 4,
                  sizeof(::TPdsSummary) );
      instance.SetNew(&new_TPdsSummary);
      instance.SetNewArray(&newArray_TPdsSummary);
      instance.SetDelete(&delete_TPdsSummary);
      instance.SetDeleteArray(&deleteArray_TPdsSummary);
      instance.SetDestructor(&destruct_TPdsSummary);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TPdsSummary*)
   {
      return GenerateInitInstanceLocal((::TPdsSummary*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TPdsSummary*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TPdsSummary::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TPdsSummary::Class_Name()
{
   return "TPdsSummary";
}

//______________________________________________________________________________
const char *TPdsSummary::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPdsSummary*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TPdsSummary::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPdsSummary*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TPdsSummary::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPdsSummary*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TPdsSummary::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPdsSummary*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TPdsSummary::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPdsSummary.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TPdsSummary::Class(),this);
   } else {
      R__b.WriteClassBuffer(TPdsSummary::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TPdsSummary(void *p) {
      return  p ? new(p) ::TPdsSummary : new ::TPdsSummary;
   }
   static void *newArray_TPdsSummary(Long_t nElements, void *p) {
      return p ? new(p) ::TPdsSummary[nElements] : new ::TPdsSummary[nElements];
   }
   // Wrapper around operator delete
   static void delete_TPdsSummary(void *p) {
      delete ((::TPdsSummary*)p);
   }
   static void deleteArray_TPdsSummary(void *p) {
      delete [] ((::TPdsSummary*)p);
   }
   static void destruct_TPdsSummary(void *p) {
      typedef ::TPdsSummary current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TPdsSummary

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 210,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace {
  void TriggerDictionaryInitialization_TPdsSummary_Dict_Impl() {
    static const char* headers[] = {
"TPdsSummary.hxx",
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
#line 1 "TPdsSummary_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <class _CharT> struct __attribute__((annotate("$clingAutoload$string")))  char_traits;
}
namespace std{template <typename > class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$TPdsSummary.hxx")))  TPdsSummary;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TPdsSummary_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TPdsSummary.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TPdsSummary", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TPdsSummary_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TPdsSummary_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TPdsSummary_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TPdsSummary_Dict() {
  TriggerDictionaryInitialization_TPdsSummary_Dict_Impl();
}
