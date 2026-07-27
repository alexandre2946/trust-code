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
#ifndef InputFile_included
#define InputFile_included

#include <fstream>
#include <string>

#if __cplusplus >= 201703L
#include <filesystem>
#endif

#include <Input.h>
#include <Objet_U.h>

class InputFile: virtual public Input, public std::ifstream, public Objet_U {

		//////////////////
		// constructors //
		//////////////////
	
	public:
		
		// we hide the std::ifstream constructors and use our own open method
		InputFile(): Input(), std::ifstream() {}
		explicit InputFile(const char* file_name, std::ios_base::openmode mode = std::ios_base::in): Input(), std::ifstream() {
			open(file_name, mode);
		}

		explicit InputFile(const std::string& file_name, std::ios_base::openmode mode = std::ios_base::in):
			InputFile(file_name.c_str(), mode) {}

    	#pragma GCC diagnostic push
    	#pragma GCC diagnostic ignored "-Wextra"
		InputFile(const InputFile& other): Input(dynamic_cast<const Input&>(other)) {}
		InputFile(InputFile& other): Input(dynamic_cast<Input&>(other)) {}
    	#pragma GCC diagnostic pop

		/////////////
		// methods //
		/////////////
	
	public:
		
		void open(const char* file_name, std::ios_base::openmode mode = std::ios_base::in);
		void open(const std::string& file_name, std::ios_base::openmode mode = std::ios_base::in)                       { open(file_name.c_str(), mode); }
#if __cplusplus >= 201703L
		void open(const std::filesystem::path& file_name, std::ios_base::openmode mode = std::ios_base::in)             { open(file_name.c_str(), mode); }
#endif

		void close();

	private:
		
		void share_buffer();
};

#endif
