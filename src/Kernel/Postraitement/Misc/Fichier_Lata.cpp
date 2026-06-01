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

#include <EcrFicPartageMPIIO.h>
#include <EcrFicPartageBin.h>
#include <Fichier_Lata.h>

/*! @brief Builds a file of type EcrFicPartage(Bin) or EcrFicPrive(Bin), binary or not depending on the "format" parameter.
 *
 *   If parallel==MULTIPLE_FILES, the file is of type EcrFicPrive(Bin).
 *     In this case, each processor opens a different file named
 *     "basename_XXXXXextension", where XXXXX equals Process::me().
 *     All processors will return is_master() == 1.
 *   If parallel==SINGLE_FILE is non-zero, the file is of type EcrFicPartage(Bin).
 *     Only the master processor opens the file; the file name is
 *     "basenameextension".
 *     is_master() returns 1 on the master, 0 on other processors.
 *
 * @param (basename) beginning of the file name
 * @param (extension) end of the file name
 * @param (mode_append) If mode_append==ERASE, opens in write mode; if mode_append==APPEND, opens in append mode.
 * @param (format) Determines whether to open in binary mode or not. (possible values: Format_Post_Lata::ASCII or Format_Post_Lata::BINAIRE)
 * @param (parallel) single shared file or multiple private files...
 */
Fichier_Lata::Fichier_Lata(const char * basename, const char * extension,
                           Mode mode_append,
                           Format_Post_Lata::Format format,
                           Format_Post_Lata::Options_Para parallel) :
  filename_(""),
  fichier_(0),
  is_parallel_(0)
{
  char s[20] = "";

  switch(parallel)
    {
    case Format_Post_Lata::SINGLE_FILE_MPIIO:
    case Format_Post_Lata::SINGLE_FILE:
      {
        is_parallel_ = 1;
        filename_ = basename;
        filename_ += extension;
        // For sequential execution, open a SFichier
        // to avoid buffering in memory
        if  (Process::is_sequential())
          fichier_ = new SFichier;
        else
          {
            if (format == Format_Post_Lata::BINAIRE && parallel == Format_Post_Lata::SINGLE_FILE_MPIIO)
              fichier_ = new EcrFicPartageMPIIO;
            else
              fichier_ = new EcrFicPartage;
          }
        break;
      }
    case Format_Post_Lata::MULTIPLE_FILES:
      {
        is_parallel_ = 0;
        fichier_ = new SFichier;
        const int moi = Process::me();
        snprintf(s, 20, "_%05d", (int)moi);
        break;
      }
    default:
      Cerr << "Fichier_Lata::Fichier_Lata: parallel option not supported " << (int)parallel << finl;
      Process::exit();
    }
  filename_ = basename;
  filename_ += s;
  filename_ += extension;

  switch(format)
    {
    case Format_Post_Lata::ASCII:
      fichier_->set_bin(0);
      break;
    case Format_Post_Lata::BINAIRE:
      fichier_->set_bin(1);
      break;
    default:
      Cerr << "Fichier_Lata::Fichier_Lata: format not supported " << (int)format << finl;
      Process::exit();
    }


  IOS_OPEN_MODE mode = ios::out;
  switch(mode_append)
    {
    case ERASE:
      mode = ios::out;
      break;
    case APPEND:
      mode = ios::out | ios::app;
      break;
    default:
      Cerr << "Fichier_Lata::Fichier_Lata: open mode not supported " << (int)mode_append << finl;
      Process::exit();
    }
  //if (Process::je_suis_maitre() || parallel==Format_Post_Lata::MULTIPLE_FILES)
  {
    int ok = fichier_->ouvrir(filename_, mode);
    if (!ok)
      {
        Cerr << "Error in Fichier_Lata::Fichier_Lata\n"
             << " Error while opening file : " << filename_
             << finl;
        Process::exit();
      }
    fichier_->setf(ios::scientific);
    fichier_->precision(8);
  }
}

Fichier_Lata::~Fichier_Lata()
{
  if (fichier_)
    {
      delete fichier_;
      fichier_ = 0;
    }
}

SFichier& Fichier_Lata::get_SFichier()
{
  assert(fichier_);
  return *fichier_;
}

/*! @brief Returns the file name with its path.
 *
 */
const Nom& Fichier_Lata::get_filename() const
{
  return filename_;
}

/*! @brief If the file is of shared type, returns 1 if me() equals the group master, 0 otherwise.
 *
 *   If the file is private, returns 1 on all processors.
 *
 */
int Fichier_Lata::is_master() const
{
  int resu = 0;
  if (is_parallel_ == 0)
    {
      // Sequential execution, private files: each processor is master
      resu = 1;
    }
  else
    {
      // Parallel execution: a single master
      if (Process::je_suis_maitre())
        resu = 1;
    }
  return resu;
}

/*! @brief If the file is of shared type, calls the syncfile() method; otherwise does nothing.
 *
 */
void Fichier_Lata::syncfile()
{
  if (is_parallel_ && Process::is_parallel())
    fichier_->syncfile();
}

Fichier_Lata_maitre::Fichier_Lata_maitre(const char * basename,
                                         const char * extension,
                                         Mode mode_append,
                                         Format_Post_Lata::Options_Para parallel) :
  Fichier_Lata(basename, extension,
               mode_append, Format_Post_Lata::ASCII, parallel)
{
  fichier_->setf(ios::scientific);
  // The precision of the master file can be changed here:
  fichier_->precision(8);
}
