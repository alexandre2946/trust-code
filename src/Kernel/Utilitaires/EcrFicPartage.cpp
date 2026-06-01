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

#include <EcrFicPartage.h>
#include <PE_Groups.h>
#include <Comm_Group.h>
#include <communications.h>
#include <Perf_counters.h>

Implemente_instanciable_sans_constructeur_ni_destructeur(EcrFicPartage,"EcrFicPartage",SFichier);
Entree& EcrFicPartage::readOn(Entree& s) { throw; }
Sortie& EcrFicPartage::printOn(Sortie& s) const { throw; }

EcrFicPartage::EcrFicPartage() : SFichier()
{
  obuffer_ptr_ = new OBuffer;
  set_bin(false);
}

/*! @brief Opens the file with the given mode and prot parameters. These parameters are the parameters of the standard open method
 *
 */
EcrFicPartage::EcrFicPartage(const char* name,IOS_OPEN_MODE mode)
{
  obuffer_ptr_ = new OBuffer;
  set_bin(false);
  obuffer_ptr_->set_64b(this->is_64b());

  ouvrir(name, mode);
}

inline OBuffer& EcrFicPartage::get_obuffer()
{
  assert(obuffer_ptr_);
  return *obuffer_ptr_;
}

/*! @brief Opens the file with the given mode and prot parameters. These parameters are the parameters of the standard open method
 *
 */
int EcrFicPartage::ouvrir(const char* name,IOS_OPEN_MODE mode)
{
  // Sanity check: are all processors present?
  barrier();

  int ok = 1;
  if(je_suis_maitre())
    {
#ifdef FILESYSTEM_NON_GLOBAL
      nom_fic_ = pwd();
      nom_fic_ += "/";
      nom_fic_ += name;
#else
      nom_fic_ = name;
#endif
      // Only the master opens the file
      ok = SFichier::ouvrir((const char *)nom_fic_, mode);
    }
  syncfile();

  // Modif B.Math. 22/09/2004: all processors go through the buffer,
  //  including the master.
  get_obuffer().new_buffer();
  return ok;
}


/*! @brief Closes the file
 *
 */
EcrFicPartage::~EcrFicPartage()
{
  close();
  delete obuffer_ptr_;
  obuffer_ptr_ = 0;
}

void EcrFicPartage::close()
{
  // Sanity check: is everyone present?
  barrier();
  const int buflen = get_obuffer().len();
  if(buflen > 0)
    {
      Cerr << "***** WARNING : EcrFicPartage::close() ******* "<<nom_fic_
           << "\non PE " << me() << " the buffer is not empty\n"
           << "  (Missing syncfile) : one makes a last syncfile" << finl;
      Cerr<<get_obuffer().str()<<finl;
    }
#ifndef NDEBUG
  // Is there a processor on which data still remains
  const trustIdType maxbuflen = mp_sum(buflen);
  if (maxbuflen > 0)
    syncfile();
#endif
  if(je_suis_maitre())
    SFichier::close();
#ifndef NDEBUG
  if (maxbuflen > 0)
    {
      Cerr<<"It missed a syncfile somewhere"<<finl;
      Cerr<<"Indeed, maxbuflen =" << maxbuflen << finl;
      Cerr<<"GF prefers to stop the calculation to correct "<<finl;
      exit();
    }
#endif
}
/*! @brief Allows the calling process to block waiting for the common resource shared by all processes, which is the shared file.
 *
 * If the calling process is not the first, it waits for the previous process to indicate the position in the file where it should write next.
 *     This method is systematically called before any new write to the file.
 *
 * @return (Sortie&) *this
 */
Sortie& EcrFicPartage::lockfile()
{
  return *this;
}


/*! @brief Releases the critical resource for the next process.
 *
 * The calling process, unless it is the first process of the group, sends
 *     its current position to the next process in the group. This method should
 *     be called after each write to the file.
 *
 * @return (Sortie&) *this
 */
Sortie& EcrFicPartage::unlockfile()
{
  return *this;
}

/*! @brief Triggers writing to disk of the data accumulated on the different processors since the last call to syncfile().
 *
 *   Data is written in ascending processor order.
 *   This function must be called the same number of times on all processors!
 *
 *   Example:
 *    processor 0:                           processor 1:
 *     file << "pe0 : 1" << finl;            file << "pe1 : 1" << finl;
 *     file << "pe0 : 2" << finl;            file << "pe1 : 2" << finl;
 *     file.syncfile();                      file.syncfile();
 *     file << "pe0 : 3" << finl;            file << "pe1 : 3" << finl;
 *     file << "pe0 : 4" << finl;
 *     file.syncfile();                      file.syncfile();
 *     file << "pe0 : end" << finl;          // processor 1 writes no data
 *     file.syncfile();                      file.syncfile();
 *
 *   File contents:
 *     pe0 : 1
 *     pe0 : 2
 *     pe1 : 1
 *     pe1 : 2
 *     pe0 : 3
 *     pe0 : 4
 *     pe1 : 3
 *     pe0 : end
 *
 */
Sortie& EcrFicPartage::syncfile()
{
  // In ASCII mode, data is converted to ASCII when written to the buffer.
  //  A second conversion occurs when writing to disk
  //  (depending on the operating system, '\n' may be encoded differently for example).
  // So, in ASCII mode, we use:
  //   file << buffer;
  // and in binary mode:
  //   file.write(buffer, size);
  // Because file << buffer determines the buffer length by looking for the '\0' character,
  // in ASCII mode the buffer must end with a '\0' character,
  // which is added here:

  if (! bin_)
    {
      get_obuffer().put_null_char();
    }

  const Comm_Group& group = PE_Groups::current_group();
  if(je_suis_maitre())
    {
      int p;
      const int nb_proc = nproc();
      for(p=0; p<nb_proc; p++)
        {
          const char * buffer_data = 0;
          char * allocated_buffer = 0;
          int buf_size;

          // We retrieve data from processor p, either directly (p==me()),
          //  or via communication:
          if (p == me())
            {
              // Writing my own data: I take it from the buffer.
              // This pointer may be null:
              buffer_data = get_obuffer().str();
              buf_size = get_obuffer().len();
            }
          else
            {
              // Writing data from another processor; retrieve it.
              int dummy = 0;
              envoyer(dummy, p, 100); // Signal processor p so it sends its data
              recevoir(buf_size, p, 100);
              if (buf_size > 0)
                {
                  buffer_data = allocated_buffer = new char[buf_size];
                  group.recv(p, allocated_buffer, buf_size, 100);
                }
            }
          // Write to disk file
          if (buf_size > 0)
            {
              assert(buffer_data);
              ostream& os = get_ostream();
              if (bin_)
                {
                  // Binary write without conversion:
                  statistics().begin_count(STD_COUNTERS::IO_EcrireFicPartageBin,statistics().get_last_opened_counter_level()+1);
                  os.write(buffer_data, buf_size);
                  statistics().end_count(STD_COUNTERS::IO_EcrireFicPartageBin,1,buf_size);
                }
              else
                {
                  // Verify that the buffer indeed ends with a 0 character:
                  assert(buffer_data[buf_size-1] == 0);
                  // Write buffer_data as a string
                  // (conversion of \n on certain systems, etc...)
                  os << buffer_data;
                }
            }
          if (allocated_buffer)
            delete[] allocated_buffer;
        }
      // Force everything to be written to disk immediately:
      // (call to the low-level flush function, not trio's).
      get_ostream().flush();
    }
  else
    {
      // Send the buffer to the master processor:
      // We wait for it to request the data to avoid congesting the network:
      // (otherwise all processors send their data simultaneously)
      int dummy;
      recevoir(dummy, 0, 100);
      int buf_size = get_obuffer().len();
      envoyer(buf_size, 0, 100);
      // If the size is non-zero, send the buffer:
      if (buf_size > 0)
        {
          const char * buffer_data = get_obuffer().str();
          group.send(0, buffer_data, buf_size, 100);
        }
    }

  // Empty the buffer:
  get_obuffer().new_buffer();
  return *this;
}

void EcrFicPartage::precision(int i) { get_obuffer().precision(i); }

int EcrFicPartage::get_precision() { return get_obuffer().get_precision(); }

Sortie& EcrFicPartage::operator <<(const char* ob)
{
  get_obuffer() << ob;
  return *this;
}

Sortie& EcrFicPartage::operator <<(const std::string& str) { return operator_template<std::string>(str);}
Sortie& EcrFicPartage::operator <<(const Separateur& s) { return operator_template<Separateur>(s);}
Sortie& EcrFicPartage::operator <<(const Objet_U& ob) { return operator_template<Objet_U>(ob);}
Sortie& EcrFicPartage::operator <<(const int ob) { return operator_template<int>(ob);}
Sortie& EcrFicPartage::operator <<(const unsigned ob) { return operator_template<unsigned>(ob);}
Sortie& EcrFicPartage::operator <<(const float ob) { return operator_template<float>(ob);}
Sortie& EcrFicPartage::operator <<(const double ob) { return operator_template<double>(ob);}
Sortie& EcrFicPartage::operator <<(const long ob) { return operator_template<long>(ob);}
Sortie& EcrFicPartage::operator <<(const long long ob) { return operator_template<long long>(ob);}
Sortie& EcrFicPartage::operator <<(const unsigned long ob) { return operator_template<unsigned long>(ob);}


int EcrFicPartage::put(const unsigned* ob, std::streamsize n, std::streamsize pas) { return put_template<unsigned>(ob,n,pas); }
int EcrFicPartage::put(const int* ob, std::streamsize n, std::streamsize pas) { return put_template<int>(ob,n,pas); }
int EcrFicPartage::put(const long* ob, std::streamsize n, std::streamsize pas) { return put_template<long>(ob,n,pas); }
int EcrFicPartage::put(const long long* ob, std::streamsize n, std::streamsize pas) { return put_template<long long>(ob,n,pas); }
int EcrFicPartage::put(const float* ob, std::streamsize n, std::streamsize pas) { return put_template<float>(ob,n,pas); }
int EcrFicPartage::put(const double* ob, std::streamsize n, std::streamsize pas) { return put_template<double>(ob,n,pas); }


void EcrFicPartage::set_bin(bool bin)
{
  SFichier::set_bin(bin);
  get_obuffer().set_bin(bin_);
}

void EcrFicPartage::set_64b(bool is64)
{
  SFichier::set_64b(is64);
  get_obuffer().set_64b(is64);
}

Sortie& EcrFicPartage::flush() { return (*this); }
