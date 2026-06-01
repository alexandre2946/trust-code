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

#ifndef EcrFicPartage_included
#define EcrFicPartage_included

#include <SFichier.h>
#include <OBuffer.h>
#include <Nom.h>

/*! @brief Writing to a shared file. This class implements the operators and virtual methods of the SFichier class as follows:
 *
 *     A file physically located on the disk of the machine hosting the master task of
 *     the Trio-U application (the process with rank 0 in the "all" group).
 *     Each process writes its data one after the other, in ascending order of rank in the "all" group.
 *     Synchronization is achieved through a semaphore system.
 *     This type of file is notably used for the creation of a single post-processing file.
 *
 */

class OBuffer;

// Each PE writes to the same() file
class EcrFicPartage : public SFichier
{
  Declare_instanciable_sans_constructeur_ni_destructeur(EcrFicPartage);
public:
  EcrFicPartage();
  EcrFicPartage(const char* name,IOS_OPEN_MODE mode=ios::out);
  int ouvrir(const char* name,IOS_OPEN_MODE mode=ios::out) override;
  ~EcrFicPartage() override;
  void close();
  Sortie& lockfile() override;
  Sortie& unlockfile() override;
  Sortie& syncfile() override;
  void precision(int) override;
  int get_precision() override;

  Sortie& flush() override;
  void set_bin(bool bin) override;
  void set_64b(bool is64) override;

  Sortie& operator <<(const char* ob) override;
  Sortie& operator <<(const Objet_U& ob) override;
  Sortie& operator <<(const std::string& ob) override;
  Sortie& operator <<(const Separateur& ) override;

  Sortie& operator <<(const int ob) override;
  Sortie& operator <<(const unsigned ob) override;
  Sortie& operator <<(const long ob) override;
  Sortie& operator <<(const long long ob) override;
  Sortie& operator <<(const unsigned long ob) override;
  Sortie& operator <<(const float ob) override;
  Sortie& operator <<(const double ob) override;

  int put(const unsigned* ob, std::streamsize n, std::streamsize pas) override;
  int put(const int* ob, std::streamsize n, std::streamsize pas) override;
  int put(const long* ob, std::streamsize n, std::streamsize pas) override;
  int put(const long long* ob, std::streamsize n, std::streamsize pas) override;
  int put(const float* ob, std::streamsize n, std::streamsize pas) override;
  int put(const double* ob, std::streamsize n, std::streamsize pas) override;

protected:
  Nom nom_fic_;

private:
  OBuffer& get_obuffer();
  OBuffer * obuffer_ptr_ = nullptr; // Pointer: avoids including OBuffer.h

  template <typename _TYPE_>
  int put_template(const _TYPE_* ob, std::streamsize n, std::streamsize pas)
  {
    get_obuffer().put(ob,n,pas);
    return 1;
  }

  template <typename _TYPE_>
  Sortie& operator_template(const _TYPE_& ob)
  {
    get_obuffer() << ob;
    return *this;
  }
};

#endif /* EcrFicPartage_included */
