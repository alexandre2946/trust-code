/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef LecFicDistribue_sansnum_included
#define LecFicDistribue_sansnum_included

#include <Separateur.h>
#include <EFichier.h>

class Objet_U;

/*! @brief This class implements the operators and virtual methods of the EFichier class as follows: there are as many files as there are processes, physically located on the disk of the machine hosting the master task of the Trio-U application (the process of rank 0 in the "all" group).
 *
 *     The master process reads one item at a time from each file and sends it to the corresponding process.
 *     The same applies to the methods for inspecting the state of a file.
 *
 */
class LecFicDistribue_sansnum : public EFichier
{
  // the master reads the file and propagates the information
private :
  LecFicDistribue_sansnum(int);
public:
  LecFicDistribue_sansnum();
  LecFicDistribue_sansnum(const char* name,IOS_OPEN_MODE mode=ios::in);

  int ouvrir(const char* name,IOS_OPEN_MODE mode=ios::in) override;

  ~LecFicDistribue_sansnum() override;

protected:

private:

};

#endif
