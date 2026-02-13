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

#include <Process.h>
#include <Comm_Group.h>
#include <PE_Groups.h>
#include <communications.h>
#include <Sortie_Nulle.h>
#include <Journal.h>
#include <SFichier.h>
#include <Synonyme_info.h>
#include <petsc_for_kernel.h>
#include <comm_incl.h>
#include <TRUST_Error.h>
#include <Comm_Group_MPI.h>
#include <unistd.h> // sleep() pour certaines machines
#include <SChaine.h>
#include <FichierHDFPar.h>
#include <EChaineJDD.h>
#include <DeviceMemory.h>
#include <kokkos++.h>
#include <fstream>

// Chacun des fichiers Cerr, Cout et Journal(i)
// peut etre redirige vers l'un des quatre fichiers suivants:
// Instance de la Sortie nulle (equivalent de /dev/null)
static Sortie_Nulle  journal_zero_;
// Instance de la Sortie pointant vers cerr
static Sortie        std_err_(cerr);
// Instance de la Sortie pointant vers cout
static Sortie        std_out_(cout);
// Instances du fichier Journal
static SFichier     journal_file_;

static int          journal_file_open_;
static Nom          journal_file_name_;
// Niveau maximal des messages ecrits. La valeur initiale determine
// si les messages ecrits avant l'initialisation du journal sont ecrits
// ou pas.
static int        verbose_level_ = 0;
static int        disable_stop_ = 0;

// Drapeau indiquant si les sorties cerr et cout doivent
// etre redirigees vers le fichier journal
static int        cerr_to_journal_ = 0;
int Process::exception_sur_exit=0;
int Process::multiple_files=5120; // Valeur modifiable avec la variable d'environnement TRUST_MultipleFiles
bool Process::force_single_file(const int ranks, const Nom& filename)
{
  char* theValue = getenv("TRUST_MultipleFiles");
  if (theValue != nullptr) multiple_files=atoi(theValue);
  if (ranks>multiple_files)
    {
      if (Process::je_suis_maitre())   // Attention, necessaire, car appel possible tres tot dans main.cpp alors que Cerr par defini completement sur les processes
        {
          Cerr << "======================================================================================================" << finl;
          Cerr << "Warning! Single file option is forced for " << filename << " above " << multiple_files << " MPI ranks." << finl;
          Cerr << "for I/O performance reasons on cluster and inodes number limitation." << finl;
          Cerr << "If you want to keep multiple files, add at the beginning of your data file to outpass the limitation:" << finl;
          Cerr << "MultipleFiles " << ranks << finl;
          Cerr << "=====================================================================================================" << finl;
        }
      return true;
    }
  else
    return false;
}

/*! @brief renvoie 1 si on est sur le processeur maitre du groupe courant (c'est a dire me() == 0), 0 sinon.
 *
 * Voir Comm_Group::rank()
 *
 */
int Process::je_suis_maitre()
{
  const int r = PE_Groups::current_group().rank();
  return r == 0;
}

/*! @brief renvoie 1 si on est sur le processeur maitre du noeud numa, 0 sinon.
 *
 */
int Process::node_master()
{
  const int r = PE_Groups::get_node_group().rank();
  return r == 0;
}

/*! @brief renvoie le nombre de processeurs dans le groupe courant Voir Comm_Group::nproc() et PE_Groups::current_group()
 *
 */
int Process::nproc()
{
  const int n = PE_Groups::current_group().nproc();
  return n;
}

bool Process::is_parallel()
{
  return Process::nproc() > 1;
}

bool Process::is_sequential()
{
  return Process::nproc() == 1;
}

/*! @brief renvoie mon rang dans le groupe de communication courant.
 *
 * Voir Comm_Group::rank() et PE_Groups::current_group()
 *
 */
int Process::me()
{
  const int r = PE_Groups::current_group().rank();
  return r;
}

/*! @brief Synchronise tous les processeurs du groupe courant (attend que tous les processeurs soient arrives a la barriere)
 *
 *    Instruction a executer sur tous les processeurs du groupe.
 *
 */
void Process::barrier()
{
  PE_Groups::current_group().barrier(0);
}


/*! @brief Calcule la somme de x sur tous les processeurs du groupe courant.
 *
 * @sa mp_max()
 */
double Process::mp_sum(double x)
{
  const Comm_Group& grp = PE_Groups::current_group();
  double y;
  grp.mp_collective_op(&x, &y, 1, Comm_Group::COLL_SUM);
  return y;
}

float Process::mp_sum(float x)
{
  const Comm_Group& grp = PE_Groups::current_group();
  float y;
  grp.mp_collective_op(&x, &y, 1, Comm_Group::COLL_SUM);
  return y;
}

/*! @brief Calcule la somme de x sur tous les processeurs du groupe courant.
 *
 * !!! Note that the sum of many int might result in a long !!!
 *
 * @sa mp_max()
 */
trustIdType Process::mp_sum(trustIdType x)
{
  const Comm_Group& grp = PE_Groups::current_group();
  trustIdType y;
  grp.mp_collective_op(&x, &y, 1, Comm_Group::COLL_SUM);
  return y;
}

template<typename _TYPE_>
void mp_collective_op_arr(TRUSTArray<_TYPE_>& x, Comm_Group::Collective_Op op, int n)
{
  int sz = n==-1 ? x.size_array() : n;
  assert_parallel<_TYPE_>(sz);
  if (sz > 0)
    {
      _TYPE_ *data = x.addr();
      _TYPE_ *tmp = new _TYPE_[sz];
      const Comm_Group& grp = PE_Groups::current_group();
      grp.mp_collective_op(data, tmp, sz, op);
      memcpy(data, tmp, sz * sizeof(_TYPE_));
      delete[] tmp;
    }
}

template<typename _TYPE_>
void Process::mp_sum_for_each_item(TRUSTArray<_TYPE_>& x, int n) { mp_collective_op_arr(x, Comm_Group::COLL_SUM, n); }

template<typename _TYPE_>
void Process::mp_max_for_each_item(TRUSTArray<_TYPE_>& x, int n) { mp_collective_op_arr(x, Comm_Group::COLL_MAX, n); }

template<typename _TYPE_>
void Process::mp_min_for_each_item(TRUSTArray<_TYPE_>& x, int n) { mp_collective_op_arr(x, Comm_Group::COLL_MIN, n); }

/*! @brief C++14 compatible mp_sum_for_each: combine multiple mp_sum calls into one collective operation
 *  Usage: mp_sum_for_each(a, b); mp_sum_for_each(a, b, c); mp_sum_for_each(a, b, c, d); mp_sum_for_each(a, b, c, d, e);
 *  All arguments must be of the same type (double or int) and are modified in place.
 *  Supports 2-5 parameters.
 */
template<typename T>
void Process::mp_sum_for_each(T& arg1, T& arg2)
{
  T data[2] = {arg1, arg2};
  T tmp[2];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 2, Comm_Group::COLL_SUM);
  arg1 = tmp[0];
  arg2 = tmp[1];
}

template<typename T>
void Process::mp_sum_for_each(T& arg1, T& arg2, T& arg3)
{
  T data[3] = {arg1, arg2, arg3};
  T tmp[3];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 3, Comm_Group::COLL_SUM);
  arg1 = tmp[0];
  arg2 = tmp[1];
  arg3 = tmp[2];
}

template<typename T>
void Process::mp_sum_for_each(T& arg1, T& arg2, T& arg3, T& arg4)
{
  T data[4] = {arg1, arg2, arg3, arg4};
  T tmp[4];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 4, Comm_Group::COLL_SUM);
  arg1 = tmp[0];
  arg2 = tmp[1];
  arg3 = tmp[2];
  arg4 = tmp[3];
}

/*! @brief C++14 compatible mp_max_for_each: combine multiple mp_max calls into one collective operation */
template<typename T>
void Process::mp_max_for_each(T& arg1, T& arg2)
{
  T data[2] = {arg1, arg2};
  T tmp[2];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 2, Comm_Group::COLL_MAX);
  arg1 = tmp[0];
  arg2 = tmp[1];
}

template<typename T>
void Process::mp_max_for_each(T& arg1, T& arg2, T& arg3)
{
  T data[3] = {arg1, arg2, arg3};
  T tmp[3];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 3, Comm_Group::COLL_MAX);
  arg1 = tmp[0];
  arg2 = tmp[1];
  arg3 = tmp[2];
}

template<typename T>
void Process::mp_max_for_each(T& arg1, T& arg2, T& arg3, T& arg4)
{
  T data[4] = {arg1, arg2, arg3, arg4};
  T tmp[4];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 4, Comm_Group::COLL_MAX);
  arg1 = tmp[0];
  arg2 = tmp[1];
  arg3 = tmp[2];
  arg4 = tmp[3];
}

/*! @brief C++14 compatible mp_min_for_each: combine multiple mp_min calls into one collective operation */
template<typename T>
void Process::mp_min_for_each(T& arg1, T& arg2)
{
  T data[2] = {arg1, arg2};
  T tmp[2];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 2, Comm_Group::COLL_MIN);
  arg1 = tmp[0];
  arg2 = tmp[1];
}

template<typename T>
void Process::mp_min_for_each(T& arg1, T& arg2, T& arg3)
{
  T data[3] = {arg1, arg2, arg3};
  T tmp[3];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 3, Comm_Group::COLL_MIN);
  arg1 = tmp[0];
  arg2 = tmp[1];
  arg3 = tmp[2];
}

template<typename T>
void Process::mp_min_for_each(T& arg1, T& arg2, T& arg3, T& arg4)
{
  T data[4] = {arg1, arg2, arg3, arg4};
  T tmp[4];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 4, Comm_Group::COLL_MIN);
  arg1 = tmp[0];
  arg2 = tmp[1];
  arg3 = tmp[2];
  arg4 = tmp[3];
}

// 5-parameter versions
template<typename T>
void Process::mp_sum_for_each(T& arg1, T& arg2, T& arg3, T& arg4, T& arg5)
{
  T data[5] = {arg1, arg2, arg3, arg4, arg5};
  T tmp[5];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 5, Comm_Group::COLL_SUM);
  arg1 = tmp[0];
  arg2 = tmp[1];
  arg3 = tmp[2];
  arg4 = tmp[3];
  arg5 = tmp[4];
}

template<typename T>
void Process::mp_max_for_each(T& arg1, T& arg2, T& arg3, T& arg4, T& arg5)
{
  T data[5] = {arg1, arg2, arg3, arg4, arg5};
  T tmp[5];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 5, Comm_Group::COLL_MAX);
  arg1 = tmp[0];
  arg2 = tmp[1];
  arg3 = tmp[2];
  arg4 = tmp[3];
  arg5 = tmp[4];
}

template<typename T>
void Process::mp_min_for_each(T& arg1, T& arg2, T& arg3, T& arg4, T& arg5)
{
  T data[5] = {arg1, arg2, arg3, arg4, arg5};
  T tmp[5];
  const Comm_Group& grp = PE_Groups::current_group();
  grp.mp_collective_op(data, tmp, 5, Comm_Group::COLL_MIN);
  arg1 = tmp[0];
  arg2 = tmp[1];
  arg3 = tmp[2];
  arg4 = tmp[3];
  arg5 = tmp[4];
}

namespace
{

template <typename T>
T mp_operations_commun_(T x, Comm_Group::Collective_Op op)
{
  const Comm_Group& grp = PE_Groups::current_group();
  T y;
  grp.mp_collective_op(&x, &y, 1, op);
  return y;
}
}

/*! @brief renvoie le plus grand int i sur l'ensemble des processeurs du groupe courant.
 *
 */
int Process::mp_max(int x) { return mp_operations_commun_(x,Comm_Group::COLL_MAX); }
double Process::mp_max(double x) { return mp_operations_commun_(x,Comm_Group::COLL_MAX); }
#if INT_is_64_ == 2
trustIdType Process::mp_max(trustIdType x) { return mp_operations_commun_(x,Comm_Group::COLL_MAX); }
#endif


/*! @brief renvoie le plus petit int i sur l'ensemble des processeurs du groupe courant.
 *
 */
int Process::mp_min(int x) { return mp_operations_commun_(x,Comm_Group::COLL_MIN); }
double Process::mp_min(double x) { return mp_operations_commun_(x,Comm_Group::COLL_MIN); }
#if INT_is_64_ == 2
trustIdType Process::mp_min(trustIdType x) { return mp_operations_commun_(x,Comm_Group::COLL_MIN); }
#endif

/*! @brief Calul de la somme partielle de i sur les processeurs 0 a me()-1 (renvoie 0 sur le processeur 0).
 *
 * Voir Comm_Group::mppartial_sum()
 *
 */
trustIdType Process::mppartial_sum(trustIdType x)
{
  const Comm_Group& grp = PE_Groups::current_group();
  trustIdType xx = x;
  trustIdType y;

  grp.mp_collective_op(&xx, &y, 1, Comm_Group::COLL_PARTIAL_SUM);
  return y;
}

/*! @brief Calcule le 'et' logique de b sur tous les processeurs du groupe courant.
 *
 */
bool Process::mp_and(bool b)
{
  const Comm_Group& grp = PE_Groups::current_group();
  int x = b ? 1 : 0;
  int y;
  grp.mp_collective_op(&x, &y, 1, Comm_Group::COLL_MIN);
  return y == 1;
}

bool Process::mp_or(bool b)
{
  const Comm_Group& grp = PE_Groups::current_group();
  int x = b ? 1 : 0;
  int y;
  grp.mp_collective_op(&x, &y, 1, Comm_Group::COLL_MAX);
  return y == 1;
}


int Process::check_int_overflow(trustIdType v)
{
  if (v >= std::numeric_limits<int>::max())
    Process::exit("Value too big - above 32b and can not be converted to int!!");
  return static_cast<int>(v);
}

/*! @brief Routine de sortie de TRUST dans une region Kokkos
 *
 */
/*
KOKKOS_FUNCTION
void Process::Kokkos_exit(const char* str)
{
#ifdef TRUST_USE_GPU
 // ToDo Kokkos: try to exit more properly on device...
 Kokkos::abort(str);
 //Kokkos::finalize();
#else
    Process::exit(str);
#endif
}*/

/*! @brief Routine de sortie de TRUST sur une erreur.
 *
 * Sauvegarde la memoire et de la hierarchie dans les fichiers "memoire.dump" et "hierarchie.dump"
 */
void Process::exit(int i)
{
  Nom message="=========================================\nTRUST has caused an error and will stop.\nUnexpected error during TRUST calculation.";
  std::string jddLine = "\nError triggered at line " + std::to_string(EChaineJDD::file_cur_line_) + " in " + Objet_U::nom_du_cas().getString() + ".data";
  message+=jddLine;
  exit(message,i);
}
void Process::exit(const Nom& message ,int i)
{
  if (exception_sur_exit == 2)
    {
      ::exit(-1); // ND 11/01/23 utilisation d'un second ::exit(-1) dans TRUST car si pas droits d'ecriture appel recursif a Process::exit()
    }

  if(je_suis_maitre())
    {
      Cerr << message << finl;
      Cerr.flush();
      // Utile pour XData et la creation de syno.py
      if (getenv("TRUST_USE_XDATA")!=nullptr)
        {
          SFichier hier("hierarchie.dump");
          hier << "\n             KEYWORDS\n";
          Type_info::hierarchie(hier);
          hier << "\n             SYNONYMS\n";
          Synonyme_info::hierarchie(hier);
        }
      if (!get_disable_stop() && Process::je_suis_maitre())
        {
          Nom nomfic( Objet_U::nom_du_cas() );
          nomfic += ".stop";
          {
            SFichier ficstop( nomfic );
            ficstop <<message<<finl;
          }
        }
    }
  Journal() << message << finl;

  if (exception_sur_exit)
    {
      // Lancement d'une exception (utilise par Execute_parallel)
      throw TRUST_Error("Error in trust ",Process::me());
    }
  else
    {
      int abort=0;
#ifdef MPI_
      if (Process::is_parallel())
        {
          // user defined groups (if any !)
          if (PE_Groups::has_user_defined_group())
            {
              auto& grp = PE_Groups::get_user_defined_group();
              if (sub_type(Comm_Group_MPI,grp))
                ref_cast_non_const(Comm_Group_MPI,grp).free_all(); // free comm + group
            }

          const MPI_Comm& mpi_comm=ref_cast(Comm_Group_MPI,PE_Groups::groupe_TRUST()).get_mpi_comm();
          int tag = 666;
          int buffer[1]= {1};
          MPI_Request request;

          // Envoi non bloquant vers me()+1
          int to_pe = (me()==nproc()-1?0:me()+1);
          MPI_Isend(buffer, 1, MPI_ENTIER, to_pe, tag, mpi_comm, &request);

          // Reception non bloquante depuis me()-1
          int from_pe = (me()==0?nproc()-1:me()-1);
          MPI_Irecv(buffer, 1, MPI_ENTIER, from_pe, tag, mpi_comm, &request);

          // Attente
          sleep(1);

          // Test si me() a recu de me()-1
          True_int ok;
          MPI_Status status;
          MPI_Test(&request,&ok,&status);
          if (!ok)
            abort=1;
        }
#endif
      if (abort)
        {
          if (!je_suis_maitre())
            {
              std_err_ << "!!! TRUST process number " << Process::me() << " exited unexpectedly ! See error message at the end of the file " << journal_file_name_ << " -> Aborting calculation..." << finl;
            }
          PE_Groups::groupe_TRUST().abort();
        }
      else
        {
#ifdef MPI_
          // On MPI_Finalize si MPI_Initialized and not MPI_Finalized
          True_int flag;
          MPI_Initialized(&flag);
          if (flag)
            {
              MPI_Finalized(&flag);
              if (!flag)
                MPI_Finalize();
            }
#endif
          PE_Groups::finalize();
        }
    }
  // Kokkos::finalize();
  // On force exit();
  if (i==0) i=-1;
  ::exit(i); //Seul ::exit utilise dans le code jusqu'a 01/23. second ajoute car appel recursif a Process::exit si droits ecriture dossier etude manquants
}

/*! @brief Routine de sortie de Trio-U sur une erreur abort()
 *
 */
void Process::abort()
{
#ifdef NDEBUG
  // En optimise on sort proprement.
  exit();
#else
  // En debug, on sort brutal pour avoir des infos avec le debugger ?
  ::abort(); //Seul ::abort() utilise dans le code
#endif
}

/*! @brief Renvoie un objet statique de type Sortie qui sert de journal d'evenements.
 *
 * Si message_level <= verbose_level_, on ecrit le message, sinon
 *   on l'envoie sur une Sortie_Nulle.
 *   Si le fichier journal est ouvert, on ecrit dans le fichier, sinon dans stderr.
 *
 */
Sortie& Process::Journal(int message_level)
{
  if (message_level <= verbose_level_ && verbose_level_ > 0)
    {
      if (journal_file_open_)
        return journal_file_;
      else
        return std_err_;
    }
  return journal_zero_;
}

// Renvoie la ram occupee par un processeur
double Process::ram_processeur()
{
#ifdef PETSCKSP_H
  PetscLogDouble memoire;
  PetscMemoryGetCurrentUsage(&memoire);
  return memoire;
#else
  return 0;
#endif
}
#include <sys/resource.h>
double ru_maxrss()
{
  // Best to track OOM
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  // ru_maxrss is in kilobytes
  long rss_kb = usage.ru_maxrss;
  return static_cast<double>(rss_kb*1024);
}

#ifndef __APPLE__
#include <malloc.h>

/*
struct mallinfo2 {
    size_t hblkhd;    // Space in mmapped regions (bytes)
    size_t hblks;     // Number of mmapped regions
    size_t usmblks;   // Always 0 (obsolete field)
    size_t fsmblks;   // Always 0 (obsolete field)
    size_t uordblks;  // Total allocated space (bytes)
    size_t fordblks;  // Total free space (bytes)
    size_t keepcost;  // Top-most, releasable space (bytes)
};
*/
double heap_allocation()
{
  // Best to track memory leak
#if defined(__GLIBC__) && __GLIBC_PREREQ(2, 33)
  struct mallinfo2 info = mallinfo2();
#else
  struct mallinfo info = mallinfo();
#endif
  return (double)info.uordblks;
}

static double heap_allocated_old=0;

#endif /* ndef __APPLE__ */

void Process::imprimer_ram_totale(int all_process)
{
  double memoire;
  //memoire = ram_processeur();
  memoire = ru_maxrss();

#ifndef __APPLE__
  double heap_allocated = heap_allocation();
#endif

  if (memoire)
    {
      //Cout << "RAM provisoire: PETSc " << ram_processeur() << "  ru_maxrss " << memoire << " mallinfo " << heap_allocated << finl;
      int Mo=1024*1024;
      if (all_process) Journal() << (int)(memoire/Mo) << " MBytes of RAM taken by the processor " << Process::me() << finl;
      {
        double max_memoire=Process::mp_max(memoire);
        double total_memoire=Process::mp_sum(memoire);
        Cout << (int)(total_memoire/Mo) << " MBytes of RAM taken by the calculation (max on a rank: "<<(int)(max_memoire/Mo)<<" MB)." << finl;
#ifndef __APPLE__
        Cout << "[RAM] Allocated heap on master rank: " << (int)(heap_allocated/Mo) << " Mbytes";
        double delta = heap_allocated - heap_allocated_old;
        if (delta!=0 && heap_allocated_old>0) Cout << " (" << (delta>0 ? "+" : "") << (long)delta << " bytes)";
        Cout << finl;
        heap_allocated_old = heap_allocated;
#ifdef TRUST_USE_GPU
        int Go = 1024 * Mo;
        double allocated = mp_max((double)DeviceMemory::allocatedBytesOnDevice());
        double total = static_cast<double>(DeviceMemory::deviceMemGetInfo(1));
        Cout << 0.1*(int)(10*allocated/Go) << " GBytes of maximal RAM allocated on a GPU (" <<  int(100 * allocated / total) << "%)" << finl;
#endif
#endif /* ndef __APPLE__ */
      }
#ifdef TRUST_USE_ROCM /* Seulement sur adastra */
      // sUnreclaim sur chaque process:
      std::ifstream meminfo("/proc/meminfo");
      std::string line;
      size_t sunreclaim_kb = 0;
      size_t mem_available_kb = 0;
      size_t mem_total_kb = 0;
      while (std::getline(meminfo, line))
        {
          if (line.substr(0, 9) == "MemTotal:")
            {
              size_t pos = line.find_first_of("0123456789");
              mem_total_kb = std::stoull(line.substr(pos));
            }
          if (line.substr(0, 13) == "MemAvailable:")
            {
              size_t pos = line.find_first_of("0123456789");
              mem_available_kb = std::stoull(line.substr(pos));
            }
          if (line.substr(0, 11) == "SUnreclaim:")
            {
              size_t pos = line.find_first_of("0123456789");
              sunreclaim_kb = std::stoull(line.substr(pos));
              break;
            }
        }
      Process::Journal() << "[RAM] SUnreclaim: " << sunreclaim_kb/1024 << " MB MemAvailable: " << mem_available_kb/1024 << " MB MemTotal: " << mem_total_kb/1024 << " MB " << finl;
#endif
    }
}

/*! @brief Initialisation du journal
 *
 * @param (verbose_level) les messages de niveau <= verbose_level seront affiches, les autres seront mis a la poubelle.
 * @param (file_name) si pointeur nul, tout le monde ecrit dans cerr, sinon c'est le nom du fichier (doit etre different sur chaque processeur)
 * @param (append) indique si on ouvre le fichier en mode append ou pas.
 */
void init_journal_file(int verbose_level, const char * file_name, int append)
{
  end_journal(verbose_level);

  if (verbose_level > 0)
    {
      if (file_name)
        {
          IOS_OPEN_MODE mode = ios::out;
          if (append)
            mode = ios::app;
          if (!journal_file_.ouvrir(file_name, mode))
            {
              Cerr << "Fatal error in init_journal_file: cannot open journal file" << finl;
              Process::exit();
            }

          journal_file_open_ = 1;
          journal_file_name_ = file_name;
        }
    }
  verbose_level_ = verbose_level;
}

void end_journal(int verbose_level)
{
  // Attention: acrobatie pour que ca "plante proprement" si le destructeur
  // ecrit dans le journal !
  journal_file_.close();
  journal_file_open_ = 0;
}

/*! @brief Renvoie l'objet Sortie sur lequel seront redirigees les objets ecrits dans Cerr.
 *
 * Cela peut etre std_err_ ou journal_file_
 *
 */
Sortie& get_Cerr()
{
  if (journal_file_open_ && cerr_to_journal_)
    return journal_file_;
  else
    {
      // dans le cas ou on a pas initialise les groupes
      // on ne peut pas tester si on est maitre
      if (PE_Groups::get_nb_groups()==0)
        return std_err_;
      // Seul le processeur maitre ecrit sur std_err_, les autres ecrivent sur journal_file_
      if (Process::je_suis_maitre())
        return std_err_;
      else if (verbose_level_)
        return journal_file_;
      else
        return journal_zero_;
    }
}

/*! @brief Si on est sur le maitre, on renvoie cout ou le fichier journal sinon journal_zero_.
 *
 *   @sa cerr_to_journal_
 *
 */
Sortie& get_Cout()
{
  if (Process::je_suis_maitre())
    {
      if (journal_file_open_ && cerr_to_journal_)
        return journal_file_;
      else
        return std_out_;
    }
  else
    {
      return journal_zero_;
    }
}

/*! @brief change la destination de Cerr et Cout Si flag=0, c'est stderr et stdout, sinon, si le fichier
 *
 *   journal est ouvert, c'est le journal, sinon c'est
 *   Sortie_Nulle
 *
 */
void set_Cerr_to_journal(int flag)
{
  cerr_to_journal_ = flag;
}

int get_journal_level()
{
  return verbose_level_;
}

void change_journal_level(int level)
{
  verbose_level_ = level;
}

/*! @brief Returns the disable_stop_ flag (Disable or not the writing of the .
 *
 * stop file)
 *
 */
int get_disable_stop()
{
  return disable_stop_;
}

/*! @brief Affects a new value to disable_stop_ flag (Disable or not the writing of the .
 *
 * stop file)
 *
 */
void change_disable_stop(int new_stop)
{
  disable_stop_ = new_stop;
}

// Explicit template instantiations for mp_*_for_each_item
template void Process::mp_sum_for_each_item<double>(TRUSTArray<double>&, int);
template void Process::mp_sum_for_each_item<long>(TRUSTArray<long>&, int);
template void Process::mp_sum_for_each_item<int>(TRUSTArray<int>&, int);
template void Process::mp_max_for_each_item<double>(TRUSTArray<double>&, int);
template void Process::mp_max_for_each_item<int>(TRUSTArray<int>&, int);
template void Process::mp_min_for_each_item<double>(TRUSTArray<double>&, int);
template void Process::mp_min_for_each_item<int>(TRUSTArray<int>&, int);

// Explicit template instantiations for mp_*_for_each (C++14 overloads)
// mp_sum_for_each
template void Process::mp_sum_for_each<double>(double&, double&);
template void Process::mp_sum_for_each<double>(double&, double&, double&);
template void Process::mp_sum_for_each<double>(double&, double&, double&, double&);
template void Process::mp_sum_for_each<double>(double&, double&, double&, double&, double&);
template void Process::mp_sum_for_each<int>(int&, int&);
template void Process::mp_sum_for_each<int>(int&, int&, int&);
template void Process::mp_sum_for_each<int>(int&, int&, int&, int&);
template void Process::mp_sum_for_each<int>(int&, int&, int&, int&, int&);
// mp_max_for_each
template void Process::mp_max_for_each<double>(double&, double&);
template void Process::mp_max_for_each<double>(double&, double&, double&);
template void Process::mp_max_for_each<double>(double&, double&, double&, double&);
template void Process::mp_max_for_each<double>(double&, double&, double&, double&, double&);
template void Process::mp_max_for_each<int>(int&, int&);
template void Process::mp_max_for_each<int>(int&, int&, int&);
template void Process::mp_max_for_each<int>(int&, int&, int&, int&);
template void Process::mp_max_for_each<int>(int&, int&, int&, int&, int&);
// mp_min_for_each
template void Process::mp_min_for_each<double>(double&, double&);
template void Process::mp_min_for_each<double>(double&, double&, double&);
template void Process::mp_min_for_each<double>(double&, double&, double&, double&);
template void Process::mp_min_for_each<double>(double&, double&, double&, double&, double&);
template void Process::mp_min_for_each<int>(int&, int&);
template void Process::mp_min_for_each<int>(int&, int&, int&);
template void Process::mp_min_for_each<int>(int&, int&, int&, int&);
template void Process::mp_min_for_each<int>(int&, int&, int&, int&, int&);
