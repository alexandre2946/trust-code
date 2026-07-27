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
#include <InputFile.h>

#ifndef LATATOOLS
#include <communications.h>
#endif


void InputFile::open(const char* file_name, std::ios_base::openmode mode) {
	// open the file only on process 0 if shared
	if (communication_mode == CommunicationMode::Isolated || communication_mode == CommunicationMode::Send) {
		std::ifstream::open(file_name, mode);
	}

	// shared the buffer accross process (not usefull but allow to check the buffer state on all process
	if (communication_mode == CommunicationMode::Send || communication_mode == CommunicationMode::Receive) {
		share_buffer();
	}
}

void InputFile::close() {
	// close on process 0 if shared
	if (communication_mode == CommunicationMode::Isolated || communication_mode == CommunicationMode::Send) {
		std::ifstream::close();
	}

	// share buffer (none now) (it also share state)
	if (communication_mode == CommunicationMode::Send || communication_mode == CommunicationMode::Receive) {
		share_buffer();
	}
}

void InputFile::share_buffer() {
#ifndef LATATOOLS
	/*
	// get the buffer on process 0
	std::streambuf* buffer = nullptr;
	if (Process::me() == 0)
		buffer = this->rdbuf();

	// share the pointer to buffer on all process
	char address[sizeof(std::streambuf*)];
	std::memcpy(address, &buffer, sizeof(std::streambuf*));
	envoyer_broadcast_array(address, sizeof(std::streambuf*), 0);
	std::memcpy(&buffer, address, sizeof(std::streambuf*));
	*/
	
	// share the state of the istream
	share_state();
#endif
}
