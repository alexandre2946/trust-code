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
#ifndef InputConcat_included
#define InputConcat_included

#include <sstream>

#include <Input.h>
#include <Nom.h>

class InputConcat : virtual public Input {
		
		/////////////
		// classes //
		/////////////
	
	private:
		
		class ConcatBuffer : public std::streambuf {
			
				//////////////////
				// constructors //
				//////////////////

			public:
	
				ConcatBuffer(std::streambuf* a, std::streambuf* b): first(a), second(b), current(a) {}

				/////////////
				// methods //
				/////////////
		
			protected:

				int_type underflow() override {
					if (!current)
						return traits_type::eof();
		
					int_type c = current->sgetc();
					if (traits_type::eq_int_type(c, traits_type::eof())) {
						if (current == first) {
							current = second;
							return underflow();
						}
						return traits_type::eof();
					}
					return c;
				}
		
				int_type uflow() override {
					int_type c = underflow();
					if (!traits_type::eq_int_type(c, traits_type::eof()))
						current->sbumpc();
					return c;
				}

				////////////////
				// attributes //
				////////////////
		
			private:

				std::streambuf* first;
				std::streambuf* second;
				std::streambuf* current;
		};

		//////////////////
		// constructors //
		//////////////////

	public:

		InputConcat(std::streambuf* a, std::streambuf* b);
		InputConcat(const Nom& string_, Input& stream);

		////////////////
		// attributes //
		////////////////
	
	protected:
	
		std::stringbuf string; // when used with an input string and a stream
		ConcatBuffer buffer;
};

#endif
