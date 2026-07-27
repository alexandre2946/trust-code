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
#ifndef Input_included
#define Input_included

#include <iostream>
#include <cstdint>
#include <limits>
#include <cassert>

#include <Process.h>
#include <InputOutputException.h>

class Objet_U;
template <typename T> class TRUST_Ref;
class TRUST_Ref_Objet_U;
class Entree;

#include <arch.h>

class Input: public std::istream {
		
		///////////
		// enums //
		///////////
	
	public:
		
		enum class ReadMode {
			ASCII,
			Bin32, // FIXME should be remove at some point
			Bin64
		};

		enum class CommunicationMode {
			Isolated,
			Send,
			Receive 
		};

		enum class ErrorMode {
			Continue,
			Exit,
			Throw
		};
		
		//////////////////
		// constructors //
		//////////////////

		using std::istream::istream;

		Input(): std::istream(nullptr) {}

    	#pragma GCC diagnostic push
    	#pragma GCC diagnostic ignored "-Wextra"

		Input(const Input& other): std::istream(other.rdbuf()) // attach same buffer
		{
		    this->clear(other.rdstate()); // copy stream state flags
			copy_modes(other);
		}

		Input(const std::istream& stream): std::istream(stream.rdbuf()) {}
		Input(const std::streambuf& buffer): Input(const_cast<std::streambuf&>(buffer)) {}

		Input(Input& other): std::istream(other.rdbuf()) // attach same buffer
		{
		    this->clear(other.rdstate()); // copy stream state flags
			copy_modes(other);
		}


		Input(std::istream& stream): std::istream(stream.rdbuf()) {}
		Input(std::streambuf& buffer): std::istream(&buffer) {}

    	#pragma GCC diagnostic pop

		///////////////
		// operators //
		///////////////

		virtual Input& operator>>(Objet_U& ob); // could be set to final once Entree is remove

		template <class T>
		Input& operator>>(TRUST_Ref<T>&) { std::cerr << __func__ << " :: SHOULD NOT BE CALLED ! Use -> !! " << std::endl ; throw; }

		Input& operator>>(TRUST_Ref_Objet_U&) { std::cerr << __func__ << " :: SHOULD NOT BE CALLED ! Use -> !! " << std::endl ; throw; }

		template <class Type, class = std::enable_if_t<!std::is_base_of<Objet_U, Type>::value>>
		Input& operator>>(Type& value) {
			get<Type>(value);
			return *this;
		}

		///////////////
		// operators //
		///////////////

		// idealy we would prevent dynamic allocation (it could mess the internal buffers)
		// however it is use for TRUST_REF so I comment the deletes
		//void* operator new(std::size_t) = delete;
		//void* operator new[](std::size_t) = delete;

		/////////////
		// methods //
		/////////////

		void attach(std::istream& stream) {
			this->rdbuf(stream.rdbuf());
			this->clear(stream.rdstate()); // copy state
		}

		void detach() {
			this->rdbuf(nullptr);
			this->setstate(std::ios::eofbit | std::ios::badbit);
		}

		//using std::istream::get;

		// get a single value
		template <class Type>
		int get(Type& value) {
			int ok = 0;

			//static_assert(std::is_trivially_copyable<Type>::value, "Can only read trivially copyable objects (for MPI support).");
			if (communication_mode != CommunicationMode::Isolated && !std::is_trivially_copyable<Type>::value) { // FIXME add a proper serializer, not just one working with Objets_U
				std::cerr << "Can only read trivially copyable objects on MPI communications (for MPI support)." << std::endl;
				Process::exit(1);
			}

			if (communication_mode == CommunicationMode::Isolated || communication_mode == CommunicationMode::Send) {
				switch (read_mode) {
					case ReadMode::ASCII:
						ok = read_ascii(value);
						break;
					case ReadMode::Bin32:
						if (!std::is_integral<Type>::value && (sizeof(Type) == 8)) // allows to read 8 bits for double
							ok = read_bin<uint64_t>(value);
						else
							ok = read_bin<uint32_t>(value);
						break;
					case ReadMode::Bin64:
						ok = read_bin<uint64_t>(value);
						break;
					default:
						std::cerr << "Invalid reading mode encounter while reading stream." << std::endl;
						Process::exit(1);
						ok = 1;
				}
			}

#ifndef LATATOOLS
			if (communication_mode == CommunicationMode::Send || communication_mode == CommunicationMode::Receive) {
				share(ok, reinterpret_cast<char*>(&value), sizeof(Type));
			}
#endif

			return ok;
		}

		template <class Type>
		int read_ascii(Type& value) {
			//std::istream::operator>>(value);
			static_cast<std::istream&>(*this) >> value; // use cast to avoid ambiguity
			return error_handler();
		}

		template <class Buffer, class Type>
		int read_bin(Type& value) {
			// the Type is non arithmetic (integer or float), we use the operator>> method to redirect to correct method
#if __cplusplus >= 201703L
			if constexpr (!std::is_arithmetic<Type>::value) { // constexpr required here if non-integer types are used, otherwise it will try to compile read_bin<Buffer>(Type* array, std::streamsize size) with the type and will fail when trying to compile the underlying get method
#else
			if (!std::is_arithmetic<Type>::value) {
#endif
				operator>>(value);
				return error_handler();
			} else {
				// if it is an arithmetic type, we can use the method for arrays
				return read_bin<Buffer>(&value, 1);
			}
		}

		// get an array
		template <class Type>
		int get(Type* array, std::streamsize size) {
			static_assert(std::is_trivially_copyable<Type>::value, "Can only read trivially copyable objects (for MPI support).");

			int ok = 1;

/*
#ifndef LATATOOLS
			if (communication_mode == CommunicationMode::Send || communication_mode == CommunicationMode::Receive) {
				if constexpr (std::is_same<Type, char>::value) {
					// use strlen to get the array size (zero-terminated) // FIXME should probably remove this, it seems strange that only arrays of char have this behavior (occure in Nom::readOn for example), furthermore, it cannot work with lata tools
					if (Process::je_suis_maitre())
						size = std::strlen(array);

					share(ok, reinterpret_cast<char*>(&size), sizeof(std::streamsize));
				}
			}

			if (!ok)
				return ok;
#endif
*/

			if (communication_mode == CommunicationMode::Isolated || communication_mode == CommunicationMode::Send) {
				// for text we use a special method
				if constexpr (std::is_same<Type, char>::value) {
					ok = read_text(reinterpret_cast<char*>(array), size);
				} else {
					switch (read_mode) {
						case ReadMode::ASCII:
							ok = read_ascii(array, size);
							break;
						case ReadMode::Bin32:
							if (!std::is_integral<Type>::value && (sizeof(Type) == 8)) // allows to read 8 bits for double
								ok = read_bin<uint64_t>(array, size);
							else
								ok = read_bin<uint32_t>(array, size);
							break;
						case ReadMode::Bin64:
							ok = read_bin<uint64_t>(array, size);
							break;
						default:
							std::cerr << "Invalid reading mode encounter while reading stream." << std::endl;
							Process::exit(1);
							ok = 0;
					}
				}
			}

/*
			if (!ok)
				return ok;
*/

#ifndef LATATOOLS
			if (communication_mode == CommunicationMode::Send || communication_mode == CommunicationMode::Receive) {
				share(ok, reinterpret_cast<char*>(array), size * sizeof(Type));
			}
#endif

			return ok;
		}

		int read_text(char* array, std::streamsize size);

		template <class Type>
		int read_ascii(Type* array, std::streamsize size) {
			for (std::streamsize index = 0; index < size; index++)
				read_ascii(array[index]);
			return error_handler();
		}

		template <class Buffer, class Type>
		int read_bin(Type* array, std::streamsize size) {
			// if the Buffer size is the same as type size and type is an integer, or if the Type is floating point, we simply read the array
			if (avoid_conversion || (std::is_integral<Type>::value && (sizeof(Buffer) == sizeof(Type))) || !std::is_integral<Type>::value) { // could use constexpr but for some reason lata_tools compile in c++14 and constexpr is thus not available
				this->read(reinterpret_cast<char*>(array), size * sizeof(Type));
				return error_handler();
			}

			// if Type is not an integer, we call the get on every element
			if (!std::is_integral<Type>::value) { 
				for (std::streamsize index = 0; index < size; index++)
					read_bin<Buffer, Type>(array[index]);

				return error_handler();
			}

			// now we only deal with integer types
			// there is two cases: Type can be signed or unsigned
			if (std::is_signed<Type>::value) {
				// read a buffer or signed version of Buffer
				using SBuffer = typename std::make_signed<Buffer>::type;
				return read_bin_integer<SBuffer, Type>(array, size);
			} else {
				using UBuffer = typename std::make_unsigned<Buffer>::type;
				return read_bin_integer<UBuffer, Type>(array, size);
			}

			return error_handler();
		}

		template <class Buffer, class Type>
		int read_bin_integer(Type* array, std::streamsize size) {
			// we read a buffer of type Buffer then convert while checking boundaries
			assert(std::is_integral<Type>::value && std::is_integral<Buffer>::value); // Type and Buffer should both be integers.;
			assert((std::is_signed<Type>::value ^ std::is_signed<Buffer>::value) == 0); // Type and Buffer types should both be signed or both be unsigned.;
			assert(sizeof(Buffer) != sizeof(Type)); // this method can be call in that case but it slower

			int ok = 1;

			// get
			Buffer* buffer = new Buffer[size];
			this->read(reinterpret_cast<char*>(buffer), size * sizeof(Buffer));
			ok &= error_handler();

			// convert
			if (sizeof(Buffer) <= sizeof(Type)) {
				// if the size of Buffer smaller than Type, than a Buffer value always fit in a Type
				for (std::streamsize index = 0; index < size; index++) {
					array[index] = static_cast<Type>(buffer[index]);
				}
			} else {
				// if the size of Type is smaller than Buffer, there might be some overflow
				for (std::streamsize index = 0; index < size; index++) {
					// we check if the value fits (we can cast the min and max because we know that sizeof(Buffer) is greater than sizeof(Type) here)
					if (buffer[index] >= static_cast<Buffer>(std::numeric_limits<Type>::min()) && buffer[index] <= static_cast<Buffer>(std::numeric_limits<Type>::max())) {
						array[index] = static_cast<Type>(buffer[index]);
					} else {
						std::cerr << "Error while reading stream: the target variable will overflow." << std::endl;
						Process::exit(1);
						ok = 0;
					}
				}
			}

			delete[] buffer;

			return ok;
		}
		
	private:
		
		void update_mode();
		int error_handler();

#ifndef LATATOOLS

	protected:

		void share(int& ok, char* array, std::streamsize size);
		void share_state();
#endif

	protected:

		void copy_modes(const Input& other);

	public:

		void set_binary(bool value);
		void set_64_bits(bool value);
		void set_diffusion(bool value);
		bool get_diffusion() const;
		void set_avoid_conversion(bool value);
		void set_error_mode(ErrorMode error_mode);
		const ErrorMode& get_error_mode() const;

		////////////////
		// attributes //
		////////////////
	
	public:

		ErrorMode error_mode = ErrorMode::Continue;
	
	protected:
		
		ReadMode read_mode = ReadMode::ASCII;
		CommunicationMode communication_mode = CommunicationMode::Isolated;

		bool is_bin = false;

#ifdef INT_is_64_
  		bool is_64_bits = true;
#else
  		bool is_64_bits = false;
#endif
	
		bool avoid_conversion = false;
};


#endif
