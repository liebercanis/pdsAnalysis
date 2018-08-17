// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TPmtEvent_Dict

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
#include "TPmtEvent.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TPmtEvent(void *p = 0);
   static void *newArray_TPmtEvent(Long_t size, void *p);
   static void delete_TPmtEvent(void *p);
   static void deleteArray_TPmtEvent(void *p);
   static void destruct_TPmtEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TPmtEvent*)
   {
      ::TPmtEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TPmtEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TPmtEvent", ::TPmtEvent::Class_Version(), "TPmtEvent.hxx", 16,
                  typeid(::TPmtEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TPmtEvent::Dictionary, isa_proxy, 4,
                  sizeof(::TPmtEvent) );
      instance.SetNew(&new_TPmtEvent);
      instance.SetNewArray(&newArray_TPmtEvent);
      instance.SetDelete(&delete_TPmtEvent);
      instance.SetDeleteArray(&deleteArray_TPmtEvent);
      instance.SetDestructor(&destruct_TPmtEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TPmtEvent*)
   {
      return GenerateInitInstanceLocal((::TPmtEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TPmtEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TPmtEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TPmtEvent::Class_Name()
{
   return "TPmtEvent";
}

//______________________________________________________________________________
const char *TPmtEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPmtEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TPmtEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPmtEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TPmtEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPmtEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TPmtEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPmtEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TPmtEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPmtEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TPmtEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(TPmtEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TPmtEvent(void *p) {
      return  p ? new(p) ::TPmtEvent : new ::TPmtEvent;
   }
   static void *newArray_TPmtEvent(Long_t nElements, void *p) {
      return p ? new(p) ::TPmtEvent[nElements] : new ::TPmtEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_TPmtEvent(void *p) {
      delete ((::TPmtEvent*)p);
   }
   static void deleteArray_TPmtEvent(void *p) {
      delete [] ((::TPmtEvent*)p);
   }
   static void destruct_TPmtEvent(void *p) {
      typedef ::TPmtEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TPmtEvent

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 210,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 210,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlETPmtHitgR_Dictionary();
   static void vectorlETPmtHitgR_TClassManip(TClass*);
   static void *new_vectorlETPmtHitgR(void *p = 0);
   static void *newArray_vectorlETPmtHitgR(Long_t size, void *p);
   static void delete_vectorlETPmtHitgR(void *p);
   static void deleteArray_vectorlETPmtHitgR(void *p);
   static void destruct_vectorlETPmtHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TPmtHit>*)
   {
      vector<TPmtHit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TPmtHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TPmtHit>", -2, "vector", 210,
                  typeid(vector<TPmtHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETPmtHitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TPmtHit>) );
      instance.SetNew(&new_vectorlETPmtHitgR);
      instance.SetNewArray(&newArray_vectorlETPmtHitgR);
      instance.SetDelete(&delete_vectorlETPmtHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlETPmtHitgR);
      instance.SetDestructor(&destruct_vectorlETPmtHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TPmtHit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<TPmtHit>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETPmtHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TPmtHit>*)0x0)->GetClass();
      vectorlETPmtHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETPmtHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETPmtHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TPmtHit> : new vector<TPmtHit>;
   }
   static void *newArray_vectorlETPmtHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TPmtHit>[nElements] : new vector<TPmtHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETPmtHitgR(void *p) {
      delete ((vector<TPmtHit>*)p);
   }
   static void deleteArray_vectorlETPmtHitgR(void *p) {
      delete [] ((vector<TPmtHit>*)p);
   }
   static void destruct_vectorlETPmtHitgR(void *p) {
      typedef vector<TPmtHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TPmtHit>

namespace {
  void TriggerDictionaryInitialization_TPmtEvent_Dict_Impl() {
    static const char* headers[] = {
"TPmtEvent.hxx",
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
#line 1 "TPmtEvent_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$TPmtEvent.hxx")))  TPmtHit;
class __attribute__((annotate("$clingAutoload$TPmtEvent.hxx")))  TPmtEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TPmtEvent_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TPmtEvent.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TPmtEvent", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TPmtEvent_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TPmtEvent_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TPmtEvent_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TPmtEvent_Dict() {
  TriggerDictionaryInitialization_TPmtEvent_Dict_Impl();
}
