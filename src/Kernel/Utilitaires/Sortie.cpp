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

#include <Objet_U.h>
#include <Sortie.h>
#include <Separateur.h>
#include <iostream>
#include <iomanip>
#include <cstring>

const Separateur finl(Separateur::ENDL);
const Separateur tspace(Separateur::SPACE);

Sortie::Sortie() :
  AbstractIO(),
  col_width_(-1)
{
  // Constructor does **not** instanciate ostream_ - typically done in derived classes
}


void Sortie::setf(IOS_FORMAT code)
{
  if(ostream_)
    ostream_->setf(code);
}

void Sortie::precision(int pre)
{
  if(ostream_)
    ostream_->precision(pre);
}

Sortie::Sortie(ostream& os) : Sortie()
{
  if(os.rdbuf())
    ostream_ = std::make_unique<ostream>(os.rdbuf());
  else
    Process::exit();
}

Sortie::Sortie(const Sortie& os) : Sortie()
{
  if (os.has_ostream())
    {
      Cerr<<"we try to copy a Sortie with a non-empty std::ostream!!!"<<finl;
      abort();
    }
}

// Assignment operators
Sortie& Sortie::operator=(ostream& os)
{
  // Make a new ostream_:
  ostream_ = std::make_unique<ostream>(os.rdbuf());
  return *this;
}

Sortie& Sortie::operator=(Sortie& os)
{
  if (os.has_ostream())
    {
      Cerr<<"we try to copy a Sortie with a non-empty std::stream!!!"<<finl;
      abort();
    }
  else  // We are copying from a void Sortie - make sure we clean our side:
    ostream_ = nullptr;
  return *this;
}

int Sortie::add_col(const double ob)
{
  if (bin_ or col_width_ == -1)
    abort();
  else
    (*ostream_) << std::right << std::setw(col_width_) << ob;
  return ostream_->good();
}

int Sortie::add_col(const char * ob)
{
  if (bin_ or col_width_ == -1)
    abort();
  else
    (*ostream_) << std::right << std::setw(col_width_) << ob;
  return ostream_->good();
}

int Sortie::put(const unsigned * ob, std::streamsize n, std::streamsize nb_col) { return put_template<unsigned>(ob,n,nb_col); }
int Sortie::put(const int * ob, std::streamsize n, std::streamsize nb_col) { return put_template<int>(ob,n,nb_col); }
int Sortie::put(const long * ob, std::streamsize n, std::streamsize nb_col) { return put_template<long>(ob,n,nb_col); }
int Sortie::put(const long long * ob, std::streamsize n, std::streamsize nb_col) { return put_template<long long>(ob,n,nb_col); }
int Sortie::put(const float * ob, std::streamsize n, std::streamsize nb_col) { return put_template<float>(ob,n,nb_col); }
int Sortie::put(const double * ob, std::streamsize n, std::streamsize nb_col) { return put_template<double>(ob,n,nb_col); }

Sortie& Sortie::operator<<(const unsigned ob) { return operator_template<unsigned>(ob); }
Sortie& Sortie::operator<<(const int ob) { return operator_template<int>(ob); }
Sortie& Sortie::operator<<(const float ob) { return operator_template<float>(ob); }
Sortie& Sortie::operator<<(const double ob) { return operator_template<double>(ob); }
Sortie& Sortie::operator<<(const long ob) { return operator_template<long>(ob); }
Sortie& Sortie::operator<<(const long long ob) { return operator_template<long long>(ob); }
Sortie& Sortie::operator<<(const unsigned long ob) { return operator_template<unsigned long>(ob); }

Sortie& Sortie::operator <<(ostream& (*f)(ostream&))
{
  // This is not ideal, but it is needed to allow "<< endl"
  // Problem: it does not work if ostream_ == 0.
  if(ostream_ && !bin_)
    (*f)(*ostream_);
  return *this;
}

Sortie& Sortie::operator <<(Sortie& (*f)(Sortie&))
{
  (*f)(*this);
  return *this;
}

Sortie& Sortie::operator <<(ios& (*f)(ios&))
{
  if (ostream_)
    (*f)(*ostream_);
  return *this;
}

Sortie& Sortie::flush()
{
  ostream_->flush();
  return *this;
}

Sortie& Sortie::operator<<(const Separateur& ob)
{
  if (bin_)
    {
      // In binary mode, separators are not written
    }
  else
    {
      switch (ob.get_type())
        {
        case Separateur::ENDL:
          // endl = '\n' + flush
#if defined(__CYGWIN__) || defined(MICROSOFT)
          // GF issue under Windows with old line
          (*ostream_)<<endl;
#else
          (*ostream_) << '\n';
#endif
          // Optional flush (can be overridden by Sortie_Fichier_base::flush())
          flush();
          break;
        case Separateur::SPACE:
          (*ostream_) << ' ';
          break;
        }
    }
  return *this;
}

Sortie& Sortie::lockfile()
{
  Cerr<<"Sortie::lockfile() : we should not go through it !!"<<finl;
  Process::exit();
  return *this;
}

Sortie& Sortie::unlockfile()
{
  Cerr<<"Sortie::unlockfile() : we should not go through it !!"<<finl;
  Process::exit();
  return *this;
}

Sortie& Sortie::syncfile()
{
  Cerr<<"Sortie::syncfile() : we should not go through it !!"<<finl;
  Process::exit();
  return *this;
}

/*! @brief Writes an object or variable.
 *
 * In this implementation (and in most derived classes), ob.printOn is simply called
 *   (except in Sortie_Nulle). Note: if the stream must be usable in either ASCII or BINARY mode,
 *   "<< space <<" or "<< finl <<" must be inserted between objects.
 *
 */
Sortie& Sortie::operator<<(const Objet_U& ob)
{
  ob.printOn(*this);
  return *this;
}

/*! @brief Writes a character string.
 *
 * Note: to correctly re-read the string in ASCII mode, it must not contain
 *  any separator (neither space nor newline, ...).
 *
 */
Sortie& Sortie::operator <<(const char* ob)
{
  if(bin_)
    {
      // This is unfortunate: in LIST there is "<< blanc <<"
      // which forces this test:
      // One should use "<< space <<" but that would change the binary .Zones files...
      if (strcmp(ob, " "))
        {
          const int n = (int)strlen(ob) + 1;
          ostream_->write((char *) ob, n * sizeof(char));
        }
    }
  else
    {
      (*ostream_) << ob;
    }
  // B.Mathieu, 7/10/2004: removing the space after ecrire(ob).
  // Warning, now one must do
  //  fichier << "chaine" << space << nombre;
  // to be able to re-read the file transparently with
  //  fichier >> motcle >> nombre;
  // Ancien code:
  //   if(!bin_)
  //     ecrire(" ");
  return *this;
}

Sortie& Sortie::operator <<(const std::string& str) { return (*this) << str.c_str(); }


/*! @brief Changes the write mode of the file.
 *
 * This method can be called at any time. However, for Ecrire_Fichier_Partage files,
 *   the change must be made only at the beginning of a block write,
 *   immediately after syncfile() (otherwise, newlines are incorrectly translated
 *   during the next syncfile() call).
 *
 */
void Sortie::set_bin(bool bin)
{
  bin_ = bin;
  if (ostream_)
    {
      Cerr<<"Error you cant change binary format after open "<<finl;
      assert(0);
      Process::exit();
    }
}

/*! @brief Low-level method to write an array of ints or floats to the stream.
 *
 * In the base class implementation, writes to ostream_.
 *   In binary mode uses ostream::write(); in ASCII mode uses ostream::operator<<().
 *   In ASCII mode, a newline is inserted after every "nb_col" values and at the end of the array.
 *   Return value: ostream_->good().
 *
 */
template<typename _TYPE_>
int Sortie::put_template(const _TYPE_ *ob, std::streamsize n, std::streamsize nb_col)
{
  assert(n >= 0);
  if (bin_)
    {
      if (must_convert<_TYPE_>())
        {
          // Need to cast, use '>>' operator  - see doc in operator_template<>
          for (int i = 0; i < n; i++) (*this) << ob[i];
        }
      else
        {
          // In binary, optimized block writing:
          std::streamsize sz = sizeof(_TYPE_);
          sz *= n;
          ostream_->write((const char*) ob, sz);
        }
    }
  else
    {
      std::streamsize j = nb_col;
      for (std::streamsize i = 0; i < n; i++)
        {
          (*ostream_) << (ob[i]) << (' ');
          j--;
          if (j <= 0)
            {
              (*ostream_) << (endl);
              j = nb_col;
            }
        }
      // If we did not end with a newline, add one
      if (j != nb_col && n > 0) (*ostream_) << (endl);

      ostream_->flush();
    }
  return ostream_->good();
}

// Explicit instanciations
template int Sortie::put_template(const unsigned *ob, std::streamsize n, std::streamsize nb_col);
template int Sortie::put_template(const int *ob, std::streamsize n, std::streamsize nb_col);
template int Sortie::put_template(const long *ob, std::streamsize n, std::streamsize nb_col);
template int Sortie::put_template(const long long *ob, std::streamsize n, std::streamsize nb_col);
template int Sortie::put_template(const float *ob, std::streamsize n, std::streamsize nb_col);
template int Sortie::put_template(const double *ob, std::streamsize n, std::streamsize nb_col);

/*! @brief Low-level method to write an int or float to the stream.
 *
 * In the base class implementation, writes to ostream_.
 *   In binary mode uses ostream::write(); in ASCII mode uses ostream::operator<<().
 *
 */
template<typename _TYPE_>
Sortie& Sortie::operator_template(const _TYPE_ &ob)
{
  if (bin_)
    {
      if (this->must_convert<_TYPE_>())
        {
          trustIdType val = static_cast<trustIdType>(ob);
          ostream_->write((char*) &val, sizeof(trustIdType));
        }
      else
        ostream_->write((char*) &ob, sizeof(_TYPE_));
    }
  else
    (*ostream_) << ob;
  return *this;
}

// Explicit instanciations
template Sortie& Sortie::operator_template(const unsigned& ob);
template Sortie& Sortie::operator_template(const int& ob);
template Sortie& Sortie::operator_template(const long& ob);
template Sortie& Sortie::operator_template(const long long& ob);
template Sortie& Sortie::operator_template(const float& ob);
template Sortie& Sortie::operator_template(const double& ob);
template Sortie& Sortie::operator_template(const unsigned long& ob);
