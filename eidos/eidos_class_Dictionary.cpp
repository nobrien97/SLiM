//
//  eidos_class_Dictionary.cpp
//  Eidos
//
//  Created by Ben Haller on 9/30/20.
//  Copyright (c) 2020 Philipp Messer.  All rights reserved.
//	A product of the Messer Lab, http://messerlab.org/slim/
//

//	This file is part of Eidos.
//
//	Eidos is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
//	Eidos is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License along with Eidos.  If not, see <http://www.gnu.org/licenses/>.

#include "eidos_class_Dictionary.h"
#include "eidos_property_signature.h"
#include "eidos_call_signature.h"

#include <utility>
#include <algorithm>
#include <vector>


#pragma mark -
#pragma mark EidosDictionaryUnretained
#pragma mark -

EidosDictionaryUnretained::EidosDictionaryUnretained(__attribute__((unused)) const EidosDictionaryUnretained &p_original) : EidosObjectElement()
{
	// copy hash_symbols_ from p_original; I think this is called only when a Substitution is created from a Mutation
	if (p_original.hash_symbols_)
	{
		hash_symbols_ = new std::unordered_map<std::string, EidosValue_SP>;
		
		*hash_symbols_ = *(p_original.hash_symbols_);
	}
}


//
// Eidos support
//
#pragma mark -
#pragma mark Eidos support
#pragma mark -

const EidosObjectClass *EidosDictionaryUnretained::Class(void) const
{
	return gEidos_EidosDictionaryUnretained_Class;
}

EidosValue_SP EidosDictionaryUnretained::ExecuteInstanceMethod(EidosGlobalStringID p_method_id, const std::vector<EidosValue_SP> &p_arguments, EidosInterpreter &p_interpreter)
{
	switch (p_method_id)
	{
		case gEidosID_getValue:		return ExecuteMethod_getValue(p_method_id, p_arguments, p_interpreter);
		//case gEidosID_setValue:	return ExecuteMethod_Accelerated_setValue(p_method_id, p_arguments, p_interpreter);
		default:					return EidosObjectElement::ExecuteInstanceMethod(p_method_id, p_arguments, p_interpreter);
	}
}

//	*********************	- (*)getValue(string $key)
//
EidosValue_SP EidosDictionaryUnretained::ExecuteMethod_getValue(EidosGlobalStringID p_method_id, const std::vector<EidosValue_SP> &p_arguments, EidosInterpreter &p_interpreter)
{
#pragma unused (p_method_id, p_arguments, p_interpreter)
	EidosValue *key_value = p_arguments[0].get();
	
	std::string key = key_value->StringAtIndex(0, nullptr);
	
	if (!hash_symbols_)
		return gStaticEidosValueNULL;
	
	auto found_iter = hash_symbols_->find(key);
	
	if (found_iter == hash_symbols_->end())
	{
		return gStaticEidosValueNULL;
	}
	else
	{
		return found_iter->second;
	}
}

//	*********************	- (void)setValue(string $key, * value)
//
EidosValue_SP EidosDictionaryUnretained::ExecuteMethod_Accelerated_setValue(EidosObjectElement **p_elements, size_t p_elements_size, EidosGlobalStringID p_method_id, const std::vector<EidosValue_SP> &p_arguments, EidosInterpreter &p_interpreter)
{
#pragma unused (p_method_id, p_arguments, p_interpreter)
	EidosValue *key_value = p_arguments[0].get();
	std::string key = key_value->StringAtIndex(0, nullptr);
	EidosValue_SP value = p_arguments[1];
	EidosValueType value_type = value->Type();
	
	// Object values can only be remembered if their class is under retain/release, so that we have control over the object lifetime
	// See also Eidos_ExecuteFunction_defineConstant(), which enforces the same rule
	if (value_type == EidosValueType::kValueObject)
	{
		const EidosObjectClass *value_class = ((EidosValue_Object *)value.get())->Class();
		
		if (!value_class->UsesRetainRelease())
			EIDOS_TERMINATION << "ERROR (EidosDictionaryUnretained::ExecuteMethod_Accelerated_setValue): setValue() can only accept object classes that are under retain/release memory management internally; class " << value_class->ElementType() << "is not.  This restriction is necessary in order to guarantee that the kept object elements remain valid." << EidosTerminate(nullptr);
	}
	
	if (value_type == EidosValueType::kValueNULL)
	{
		// Setting a key to NULL removes it from the map
		for (size_t element_index = 0; element_index < p_elements_size; ++element_index)
		{
			EidosDictionaryUnretained *element = (EidosDictionaryUnretained *)(p_elements[element_index]);
			
			if (element->hash_symbols_)
				element->hash_symbols_->erase(key);
		}
	}
	else
	{
		// BCH: Copy values just as EidosSymbolTable does, to prevent them from being modified underneath us etc.
		// Note that when setting a value across multiple object targets, however, they all receive the same copy.
		// That should be safe; there should be no way for that value to get modified after we have copied it.  1/28/2019
		// If we have the only reference to the value, we don't need to copy it; otherwise we copy, since we don't want to hold
		// onto a reference that somebody else might modify under us (or that we might modify under them, with syntaxes like
		// x[2]=...; and x=x+1;). If the value is invisible then we copy it, since the symbol table never stores invisible values.
		if ((value->UseCount() != 1) || value->Invisible())
			value = value->CopyValues();
		
		for (size_t element_index = 0; element_index < p_elements_size; ++element_index)
		{
			EidosDictionaryUnretained *element = (EidosDictionaryUnretained *)(p_elements[element_index]);
			
			if (!element->hash_symbols_)
				element->hash_symbols_ = new std::unordered_map<std::string, EidosValue_SP>;
			
			(*element->hash_symbols_)[key] = value;
		}
	}
	
	return gStaticEidosValueVOID;
}


//
//	EidosDictionaryUnretained_Class
//
#pragma mark -
#pragma mark EidosDictionaryUnretained_Class
#pragma mark -

EidosObjectClass *gEidos_EidosDictionaryUnretained_Class = new EidosDictionaryUnretained_Class();


const std::string &EidosDictionaryUnretained_Class::ElementType(void) const
{
	return gEidosStr_EidosDictionaryUnretained;		// Note that this class name is not visible to users at this point; "EidosDictionary" is EidosDictionaryRetained
}

const std::vector<EidosMethodSignature_CSP> *EidosDictionaryUnretained_Class::Methods(void) const
{
	static std::vector<EidosMethodSignature_CSP> *methods = nullptr;
	
	if (!methods)
	{
		methods = new std::vector<EidosMethodSignature_CSP>(*EidosObjectClass::Methods());
		
		methods->emplace_back((EidosInstanceMethodSignature *)(new EidosInstanceMethodSignature(gEidosStr_getValue, kEidosValueMaskAny))->AddString_S("key"));
		methods->emplace_back(((EidosInstanceMethodSignature *)(new EidosInstanceMethodSignature(gEidosStr_setValue, kEidosValueMaskVOID))->AddString_S("key")->AddAny("value"))->DeclareAcceleratedImp(EidosDictionaryUnretained::ExecuteMethod_Accelerated_setValue));
		
		std::sort(methods->begin(), methods->end(), CompareEidosCallSignatures);
	}
	
	return methods;
}


//
//	EidosDictionaryRetained_Class
//
#pragma mark -
#pragma mark EidosDictionaryRetained_Class
#pragma mark -

EidosObjectClass *gEidos_EidosDictionaryRetained_Class = new EidosDictionaryRetained_Class();


const std::string &EidosDictionaryRetained_Class::ElementType(void) const
{
	return gEidosStr_Dictionary;
}

bool EidosDictionaryRetained_Class::UsesRetainRelease(void) const
{
	return true;
}









































