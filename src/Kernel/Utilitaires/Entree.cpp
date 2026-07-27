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

#include <fstream>

#include <Entree.h>
#include <Objet_U.h>
#include <Nom.h>

int Entree::jumpOfLines()
{
  if(rdbuf() != nullptr)
    {
      int jol = 0;
      char tmp=(char)this->peek();
      while(isspace(tmp))  //tmp=='\n')
        {
          if(tmp=='\n')
            jol++;
          std::istream::get(tmp);
          tmp=(char)this->peek();
        }
      return jol;
    }
  return -1;
}

int is_a_binary_file(Nom& filename)
{
  // On parcourt les 1000 premiers octets
  // Detection tres imparfaite donc limite
  // aux fichiers geom de TRUST
  int n=0;
  int c;
  std::ifstream fic(filename.getChar());
  // Si on rencontre un caractere ASCII>127
  // alors le fichier est de type binaire
  while((c = fic.get()) != EOF && n++<1000)
    if ((c>127)||(c<9))
      return 1;
  // GF sous windows les caracteres binaire sont surtout <9
  //  else printf("ici %d %c \n",c,c);
  return 0;
}

void error_convert(const char * s, const char * type)
{
  Cerr << "Error converting a string to type " <<  type << " : string = " << s << finl;
  Process::exit();
}

/*! @brief methode de conversion
 *
 */
void convert_to(const char *s, int& ob)
{
  errno = 0;
  char * errorptr = 0;
  ob = (int)strtol(s, &errorptr, 0 /* base 10 par defaut */);
  if (errno || *errorptr != 0) error_convert(s,"int");
}

void convert_to(const char *s, long& ob)
{
  errno = 0;
  char * errorptr = 0;
  ob = strtol(s, &errorptr, 0 /* base 10 par defaut */);
  if (errno || *errorptr != 0)  error_convert(s,"long");
}

void convert_to(const char *s, long long& ob)
{
  errno = 0;
  char * errorptr = 0;
#ifdef HPPA_11 /* NO_PROCESS */
  ob = strtol(s, &errorptr, 0 /* base 10 par defaut */);
#else /* NO_PROCESS */
#ifdef MICROSOFT
  ob = _strtoi64(s, &errorptr, 0 /* base 10 par defaut */);
#else
  ob = strtoll(s, &errorptr, 0 /* base 10 par defaut */);
#endif
#endif /* NO_PROCESS */
  if (errno || *errorptr != 0)  error_convert(s,"long long");
}

void convert_to(const char *s, float& ob)
{
  errno = 0;
  char * errorptr = 0;
  ob = strtof(s, &errorptr);
  if (errno || *errorptr != 0)  error_convert(s,"float");
}

void convert_to(const char *s, double& ob)
{
  errno = 0;
  char * errorptr = 0;
  ob = strtod(s, &errorptr);
  if (errno || *errorptr != 0)  error_convert(s,"double");
}

/*! @brief renvoie error_action_ pour cette entree (permet de la modifier et de restaurer ensuite la valeur anterieure)
 *
 */
Entree::Error_Action Entree::get_error_action()
{
	switch (error_mode) {
		case ErrorMode::Continue: return ERROR_CONTINUE;
		case ErrorMode::Exit:     return ERROR_EXIT;
		case ErrorMode::Throw:    return ERROR_EXCEPTION;
		default:
			Cerr << "Invalid mode." << finl;
			Process::exit();
			return ERROR_CONTINUE;
	}
}

/*! @brief Change le comportement en cas d'erreur de l'entree, voir error_handle_() et get_error_action()
 *
 */
void Entree::set_error_action(Entree::Error_Action action)
{
	switch (action) {
		case ERROR_CONTINUE:  set_error_mode(ErrorMode::Continue); break;
		case ERROR_EXIT:      set_error_mode(ErrorMode::Exit);     break;
		case ERROR_EXCEPTION: set_error_mode(ErrorMode::Throw);    break;
	}
}

void Entree::set_istream(std::istream* stream_pointer) {
	if (stream_pointer == nullptr)
		detach();
	else
		attach(*stream_pointer);
}

void Entree::set_istream(std::istream& stream) {
	attach(stream);
}

const std::istream& Entree::get_istream() const {
	return *this;
}

std::istream& Entree::get_istream() {
	return *this;
}
