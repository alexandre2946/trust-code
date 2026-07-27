/****************************************************************************
* Copyright (c) 2026, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/

#ifndef Entree_included
#define Entree_included

#include <Input.h>

class Objet_U;

class Entree: virtual public Input {
	public:
		using Input::Input;

		Entree(const Input& other): Input(other) {}
		Entree(Input& other): Input(other) {}


	    #pragma GCC diagnostic push
	    #pragma GCC diagnostic ignored "-Wextra"

	    Entree(const Entree& other) : Input(other.rdbuf()) {}
	    Entree(Entree& other) : Input(other.rdbuf()) {}

	    #pragma GCC diagnostic pop
	
		// forward operator>> to Input implementation (need to do one explicitly for lvalues)
	    Entree& operator>>(Objet_U& ob) override {
	        Input::operator>>(ob);
	        return *this;
	    }
	
	    template <class Type>
	    Entree& operator>>(Type& value) & {
	        Input::operator>>(value);
	        return *this;
	    }
	
	    template <class Type>
	    Entree& operator>>(Type&) && = delete;

		void set_check_types(int) {}

		virtual void set_bin(int value) { set_bin(!!value); }
		virtual void set_bin(bool value) { set_binary(value); }
		virtual void set_64b(bool value) { set_64_bits(value); }
		virtual void set_diffuse(bool value) { diffuse_ = true; } // true by default and cannot be change
		virtual bool get_diffuse() const { return diffuse_; }

		void set_istream(std::istream* stream_pointer);
		void set_istream(std::istream& stream);
		virtual const std::istream& get_istream() const;
		virtual std::istream& get_istream();
  		virtual int jumpOfLines();

	public:

		enum Error_Action { ERROR_EXIT, ERROR_CONTINUE, ERROR_EXCEPTION };
		virtual void set_error_action(Error_Action);
		Error_Action get_error_action();
	
	protected:
		
		bool diffuse_ = true; // this flag as no real meaning in TRUST, but some classe use it
};

// for storing opened binary files
class Nom;
int is_a_binary_file(Nom&);

void convert_to(const char *s, True_int& ob);
void convert_to(const char *s, long& ob);
void convert_to(const char *s, long long& ob);
void convert_to(const char *s, float& ob);
void convert_to(const char *s, double& ob);

#endif
