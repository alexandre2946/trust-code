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

#include <cstring>

#include <Input.h>
#include <Objet_U.h>
#include <Process.h>
#include <InputOutputException.h>

#ifndef LATATOOLS
#include <communications.h>
#endif

Input& Input::operator>>(Objet_U& ob) {
	return ob.readOn(*this);
}

#ifndef LATATOOLS
void Input::share(int& ok, char* array, std::streamsize size) {
	assert(get_diffusion()); // this method call only make sens in diffusion mode

	assert(size < std::numeric_limits<std::streamsize>::max());

	envoyer_broadcast(ok, 0);

	if (ok)
        envoyer_broadcast_array(array, static_cast<int>(size), 0); // send from processor 0 to other processors
	
    // broadcast stream state
	share_state();
}

void Input::share_state() {
	// get the state for the process
    std::ios::iostate state = std::ios::goodbit;
    if (Process::me() == 0) {
        state = this->rdstate();
    }

    // broadcast stream state of process 0
    int state_int = static_cast<int>(state);
    envoyer_broadcast(state_int, 0);

    // apply state locally
	this->clear();
	state = static_cast<std::ios::iostate>(state_int);
    this->setstate(state);
}
#endif

void Input::update_mode() {
	if (is_bin) {
		if (is_64_bits) 
			read_mode = Input::ReadMode::Bin64;
		else
			read_mode = Input::ReadMode::Bin32;
	} else {
		read_mode = Input::ReadMode::ASCII;
	}
}
	
void Input::set_binary(bool value) {
	is_bin = value;
	update_mode();
}

void Input::set_64_bits(bool value) {
	is_64_bits = value;
	update_mode();
}

void Input::set_diffusion(bool value) {
	if (Process::nproc() > 1) {
		if (value) {
			if (Process::me() == 0)
				communication_mode = CommunicationMode::Send;
			else
				communication_mode = CommunicationMode::Receive;
		} else {
			communication_mode = CommunicationMode::Isolated;
		}
	} else { // no need for diffusion operations
		communication_mode = CommunicationMode::Isolated;
	}
}

bool Input::get_diffusion() const {
	return communication_mode != CommunicationMode::Isolated;
}

void Input::set_avoid_conversion(bool value) {
	avoid_conversion = value;
}

int Input::read_text(char* ob, std::streamsize bufsize) {

	assert(rdbuf() != 0);
	assert(bufsize > 0);
	ob[bufsize-1] = 1;


	//std::cout << "read text start (bin: " << is_bin << ", size: " << bufsize << ")" << std::endl;
	//ssize_t position = tellg();

	if (is_bin) {
	    // En binaire, on lit jusqu'au prochain caractere 0
	    // (on lit les espaces, retours a la ligne etc)
	    std::streamsize i;
	    for (i = 0; i < bufsize-1; i++)
	      {
	        (*this).read(ob+i, sizeof(char));
			//std::cout << "read " << i << ": " << ob << std::endl;
	        if (fail()) {
	          ob[i] = 0;
			  //std::cout << "fail" << std::endl;
			}
	        if (ob[i] == 0) {
			  //std::cout << "eos, peak: " << std::string(1, char(this->peek())) << std::endl;
			  /*
			  if (i == 0)
			  	Process::exit("Empty");
				*/
	          break;
			 }
	      }
	    ob[i] = 0;
	} else {
	    // L'appel suivant permet a l'operateur>> de limiter le nombre
	    // de caracteres lus. On lira au maximum bufsize-1 caracteres.
	    std::istream& is = dynamic_cast<std::istream&>(*this);

		is >> std::ws; // skip leading spaces

		std::streamsize i = 0;
		for (; i < bufsize - 1; ++i) {
			int c = is.peek();
			if (c == EOF || std::isspace(static_cast<unsigned char>(c)))
				break;

			ob[i] = static_cast<char>(is.get()); // advance in the stream
		}
		ob[i] = '\0';

		//is >> std::ws; // skip remaining spaces

	    if (fail())
	      ob[0] = 0;
	}

	if (ob[bufsize-1] == 0) {
		// Note Benoit Mathieu:
		// Si on a rempli le buffer jusqu'au bout, c'est qu'il est probablement
		// trop petit. Lire la suite est dangereux car en ascii on ne sait pas
		// si on a pu lire pile poil la chaine (donc c'est ok), ou si la chaine
		// n'a pas ete lue en int. Il faut tester tres serieusement la stl et
		// ca depend sans doute de l'implementation. Donc, si le buffer est plein,
		// on plante le code et tant pis.
		Cerr << "Error in Entree::lire(char* ob, int bufsize) : buffer too small" << finl;
		Process::exit();
	}

	//std::cout << "read '" << ob << "' at " << position << std::endl;

	return error_handler();
}

void Input::set_error_mode(ErrorMode error_mode_) {
	error_mode = error_mode_;
}

int Input::error_handler() {
	// if fail bit not set, its okay
	if (!fail())
		return 1;
	
	switch (error_mode) {
		case ErrorMode::Continue: break; // do nothing
		case ErrorMode::Exit:
			Cerr << "Error while reading in Entree object. Exiting.\n";
			if (eof())
				Cerr << " End of file reached." << finl;
			else
				Cerr << " IO error (not an EOF error)." << finl;

			Process::exit();
			break;
		case ErrorMode::Throw:
			if (eof()) {
				InputOutputException e("End of file reached.");
				throw(e);
			} else {
				InputOutputException e("Unknown error.");
				throw(e);
			}
			break;
	}
	return 0;
}


void Input::copy_modes(const Input& other) {
	read_mode = other.read_mode;
	error_mode = other.error_mode;
	communication_mode = other.communication_mode;

	is_bin = other.is_bin;
  	is_64_bits = other.is_64_bits;
	avoid_conversion = other.avoid_conversion;
}
