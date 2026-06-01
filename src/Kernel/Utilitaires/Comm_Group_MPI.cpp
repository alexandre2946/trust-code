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
#include <Comm_Group_MPI.h>
#include <petsc_for_kernel.h>
#include <communications.h>
#include <PE_Groups.h>
#include <vector>
#include <Perf_counters.h>
#ifdef INT_is_64_
#include <algorithm>
#endif

Implemente_instanciable_sans_constructeur_ni_destructeur(Comm_Group_MPI,"Comm_Group_MPI",Comm_Group);


#ifdef MPI_

MPI_Status  * Comm_Group_MPI::mpi_status_ = 0;
MPI_Request * Comm_Group_MPI::mpi_requests_ = 0;
int Comm_Group_MPI::mpi_nrequests_ = -1;
int Comm_Group_MPI::mpi_maxrequests_ = -1;
int Comm_Group_MPI::current_msg_size_;
MPI_Comm Comm_Group_MPI::trio_u_world_ = MPI_COMM_WORLD;
// By default, we initialize mpi at statup (see set_must_mpi_initialize())
bool Comm_Group_MPI::must_mpi_initialize_ = true;

namespace
{
/*! @brief Non-inline part of the MPI error handler.
 *
 * Displays an MPI error code using MPI_Error_string.
 *
 */
void mpi_print_error(int error_code)
{
  Cerr << "mpi_error in Comm_Group_MPI : error_code = " << error_code << finl;
  Process::Journal() << "mpi_error in Comm_Group_MPI : error_code = " << error_code << finl;
  int length = 0;
  char message[MPI_MAX_ERROR_STRING];
  MPI_Error_string(error_code, message, & length);
  if (length > 0)
    {
      Cerr << message << finl;
      Process::Journal() << message << finl;
    }
  // Normally we would have used trio_u_world_, but it is not accessible here.
  // For calling abort, in this case it does not matter, but please do not
  // use this as an example...
  assert(0);
  MPI_Abort(MPI_COMM_WORLD,-1);
  Process::exit();
}

/*! @brief Inline part of the MPI error handler (the test is inlined: except in error cases, there is no additional
 *
 *  function call.
 *
 */
inline void mpi_error(int error_code)
{
  if (error_code != MPI_SUCCESS)
    mpi_print_error(error_code);
}

} // end anonymous NS
#endif


Sortie& Comm_Group_MPI::printOn(Sortie& os) const
{
  exit();
  return os;
}

Entree& Comm_Group_MPI::readOn(Entree& is)
{
  exit();
  return is;
}

/*! @brief Default constructor.
 *
 * You must then call init_group() or init_group_trio() to finish constructing the group.
 *
 */
Comm_Group_MPI::Comm_Group_MPI()
#ifdef MPI_
  : mpi_group_(MPI_GROUP_NULL),
    mpi_comm_(MPI_COMM_NULL),
    must_finalize_(-1) // -1 indicates that the group has not been initialized
#endif
{
}

Comm_Group_MPI::~Comm_Group_MPI()
{
#ifdef MPI_
  // If group not initialized, do nothing:
  // Modified by BM (20/08/2012): only destroy these static members if this is the main group (end of execution)
  if ((mpi_comm_!=MPI_COMM_NULL) && (mpi_comm_ == trio_u_world_))
    {
      delete [] mpi_status_;
      mpi_status_=0;

      for (int r=0; r<mpi_maxrequests_; r++)
        {
          if(mpi_requests_[r]!=MPI_REQUEST_NULL)
            {
              MPI_Request_free(&(mpi_requests_[r]));
            }
        }

      delete [] mpi_requests_;
      mpi_requests_=0;
    }

  else // destroy non-principal groups
    {
      if (mpi_comm_!=MPI_COMM_NULL)
        {
          // destroy the group then the mpi_comm
          mpi_error(MPI_Comm_free(&mpi_comm_));
          assert(mpi_comm_==MPI_COMM_NULL);
        }
      if (mpi_group_!=MPI_GROUP_NULL)
        {
          mpi_error(MPI_Group_free( &mpi_group_));
          assert(mpi_group_==MPI_GROUP_NULL);
        }
    }
#endif
}

/*! @brief Calls MPI_Abort and returns.
 *
 */
void Comm_Group_MPI::abort() const
{
#ifdef MPI_
  MPI_Abort(trio_u_world_,-1);
#endif
}

#ifdef MPI_
template <typename _TYPE_, int TYP_IDX>
void Comm_Group_MPI::mp_collective_op_template(const _TYPE_ *x, _TYPE_ *resu, int n, Comm_Group::Collective_Op op) const
{
  static_assert(TYP_IDX >= 1 && TYP_IDX <= 4, "Invalid type index!");
  MPI_Datatype mpi_typ = TYP_IDX==1 ? MPI_INT : (TYP_IDX==2 ? MPI_LONG : (TYP_IDX==3 ? MPI_DOUBLE : MPI_FLOAT));
  if (n <= 0) return;
  double s = -1;
  bool clock_on = statistics().is_gpu_verbose_on() && Process::je_suis_maitre();
  switch(op)
    {
    case Comm_Group::COLL_SUM:
      statistics().begin_count(STD_COUNTERS::mpi_sumdouble);
      mpi_error(MPI_Allreduce(x, resu, n, mpi_typ, MPI_SUM, mpi_comm_));
      if (clock_on && statistics().is_running(STD_COUNTERS::mpi_sumdouble))
        s = statistics().get_time_since_last_open(STD_COUNTERS::mpi_sumdouble);
      statistics().end_count(STD_COUNTERS::mpi_sumdouble);
      break;
    case Comm_Group::COLL_MIN:
      statistics().begin_count(STD_COUNTERS::mpi_mindouble);
      mpi_error(MPI_Allreduce(x, resu, n, mpi_typ, MPI_MIN, mpi_comm_));
      if (clock_on && statistics().is_running(STD_COUNTERS::mpi_mindouble))
        s = statistics().get_time_since_last_open(STD_COUNTERS::mpi_mindouble);
      statistics().end_count(STD_COUNTERS::mpi_mindouble);
      break;
    case Comm_Group::COLL_MAX:
      statistics().begin_count(STD_COUNTERS::mpi_maxdouble);
      mpi_error(MPI_Allreduce(x, resu, n, mpi_typ, MPI_MAX, mpi_comm_));
      if (clock_on && statistics().is_running(STD_COUNTERS::mpi_maxdouble))
        s = statistics().get_time_since_last_open(STD_COUNTERS::mpi_maxdouble);
      statistics().end_count(STD_COUNTERS::mpi_maxdouble);
      break;
    case Comm_Group::COLL_PARTIAL_SUM:
      internal_collective(x, resu, n, &op, -1 /* only one operation */, 0 /* recursion level */);
      break;
    }
  if (s>0) // Display
    {
      std::string clock(Process::is_parallel() ? "[clock]#" + std::to_string(Process::me()) : "[clock]  ");
      std::string mpi_reduce = "mp_sum";
      if (op==Comm_Group::COLL_MIN)  mpi_reduce = "mp_min";
      else if (op==Comm_Group::COLL_MAX)  mpi_reduce = "mp_max";
      printf("%s %7.3f ms [MPI]    %s\n", clock.c_str(), 0.001 * s, mpi_reduce.c_str());
      fflush(stdout);
    }
}
#endif

void Comm_Group_MPI::mp_collective_op(const double *x, double *resu, int n, Collective_Op op) const
{
#ifdef MPI_
  mp_collective_op_template<double, 3 /*double*/>(x, resu, n, op);
#endif
}

void Comm_Group_MPI::mp_collective_op(const float *x, float *resu, int n, Collective_Op op) const
{
#ifdef MPI_
  mp_collective_op_template<float, 4 /*float*/>(x, resu, n, op);
#endif
}

void Comm_Group_MPI::mp_collective_op(const int *x, int *resu, int n, Collective_Op op) const
{
#ifdef MPI_
  mp_collective_op_template<int, 1 /*int*/>(x, resu, n, op);
#endif
}

#if INT_is_64_ == 2
void Comm_Group_MPI::mp_collective_op(const trustIdType *x, trustIdType *resu, int n, Collective_Op op) const
{
#ifdef MPI_
  mp_collective_op_template<trustIdType, 2 /*long*/>(x, resu, n, op);
#endif
}
#endif

void Comm_Group_MPI::mp_collective_op(const double *x, double *resu, const Collective_Op *op, int n) const
{
#ifdef MPI_
  if (n <= 0)
    return;
  internal_collective(x, resu, n, op, n /* n different operations */, 0 /* recursion level */);
#endif
}

void Comm_Group_MPI::mp_collective_op(const float *x, float *resu, const Collective_Op *op, int n) const
{
#ifdef MPI_
  if (n <= 0)
    return;
  internal_collective(x, resu, n, op, n /* n different operations */, 0 /* recursion level */);
#endif
}

void Comm_Group_MPI::mp_collective_op(const int *x, int *resu, const Collective_Op *op, int n) const
{
#ifdef MPI_
  if (n <= 0)
    return;
  internal_collective(x, resu, n, op, n /* n different operations */, 0 /* recursion level */);
#endif
}

#if INT_is_64_ == 2
void Comm_Group_MPI::mp_collective_op(const trustIdType *x, trustIdType *resu, const Collective_Op *op, int n) const
{
#ifdef MPI_
  if (n <= 0)
    return;
  internal_collective(x, resu, n, op, n /* n different operations */, 0 /* recursion level */);
#endif
}
#endif


/*! @brief Synchronization point for all processors in the group (allows checking that everyone is present.
 *
 * ..). If check_enabled() is
 *  non-zero, the tag is used to verify that all processors
 *  are waiting on the same tag, otherwise it is a simple barrier.
 *  The tag must satisfy 0 <= tag < max_tag (i.e. 32).
 *
 */
void Comm_Group_MPI::barrier(int tag) const
{
#ifdef MPI_
  static const int max_tag = 32;
  statistics().begin_count(STD_COUNTERS::mpi_barrier);
  assert(tag >= 0 && tag < max_tag);
  if (check_enabled())
    {
      // We perform the barrier with mpmin and mpmax to verify
      // that the tag is the same on all processors:
      // WARNING: we need "int" and not "entier" !!!
      int tag_complet = get_new_tag() * max_tag + tag;
      int min_tag, amax_tag;
      mpi_error(MPI_Allreduce(& tag_complet, & min_tag, 1, MPI_ENTIER, MPI_MIN, mpi_comm_));
      mpi_error(MPI_Allreduce(& tag_complet, & amax_tag, 1, MPI_ENTIER, MPI_MAX, mpi_comm_));
      if (min_tag != tag_complet || amax_tag != tag_complet)
        {
          Cerr << "Error in Comm_Group_MPI::barrier(int tag)\n";
          Cerr << " the tag is not identical on all the processes.\n";
          Cerr << " (Loss of communications synchronisation)." << finl;
          Process::Journal() << "Comm_Group_MPI::barrier\n Error : tag = " << tag << finl;
          assert(0);
          exit();
        }
    }
  else
    {
      // Simple barrier without the tag:
      mpi_error(MPI_Barrier(mpi_comm_));
    }
  statistics().end_count(STD_COUNTERS::mpi_barrier);
#endif
}

/*! @brief Starts sending and receiving buffers.
 *
 * Buffers must remain valid until send_recv_finish() returns.
 *  The communication graph and buffer sizes must be correct!
 *
 *  send_list : list of processors (numbered within the current group) to send to
 *  send_size : size in bytes of each message
 *  send_buffers : address of data to send.
 *  recv_...  : same for data to receive.
 *
 *
 */
void Comm_Group_MPI::send_recv_start(const ArrOfInt& send_list,
                                     const ArrOfInt& send_size,
                                     const char * const * const send_buffers,
                                     const ArrOfInt& recv_list,
                                     const ArrOfInt& recv_size,
                                     char * const * const recv_buffers,
                                     TypeHint typehint) const
{
#ifdef MPI_
  statistics().begin_count(STD_COUNTERS::mpi_sendrecv);
  assert(mpi_nrequests_ < 0);

  const int tag = get_new_tag();
  int i, n;
  mpi_nrequests_ = 0;
  int msg_size = 0;

  int divisor = 0;
  MPI_Datatype datatype = MPI_CHAR;
  assert(sizeof(int) == sizeof(int)); // Otherwise MPI_ENTIER must be changed !!!
  switch(typehint)
    {
    case CHAR:
      divisor = 1;
      break;
    case INT:
      divisor = sizeof(int);
      datatype = MPI_ENTIER;
      break;
    case DOUBLE:
      divisor = sizeof(double);
      datatype = MPI_DOUBLE;
      break;
    case FLOAT:
      divisor = sizeof(float);
      datatype = MPI_FLOAT;
      break;
    default:
      Process::exit();
    }

  // Trick to maximize the chances of success: we declare
  // the reception first and the send afterwards.
  n = recv_list.size_array();
  for (i = 0; i < n; i++)
    {
      int source = recv_list[i];
      int sz   = recv_size[i];
      msg_size += sz;
      assert(source >= 0 && source < nproc());
      assert(mpi_nrequests_ < mpi_maxrequests_);
      assert(sz % divisor == 0);
      assert(mpi_requests_[mpi_nrequests_]==MPI_REQUEST_NULL);
      mpi_error(MPI_Irecv(recv_buffers[i], sz / divisor,
                          datatype,
                          source, tag, mpi_comm_,
                          & mpi_requests_[mpi_nrequests_]));
      mpi_nrequests_++;
    }

  n = send_list.size_array();
  for (i = 0; i < n; i++)
    {
      int dest = send_list[i];
      int sz   = send_size[i];
      msg_size += sz;
      assert(dest >= 0 && dest < nproc());
      assert(mpi_nrequests_ < mpi_maxrequests_);
      assert(sz % divisor == 0);
      mpi_error(MPI_Isend((char*) send_buffers[i], sz / divisor,
                          datatype,
                          dest, tag, mpi_comm_,
                          & mpi_requests_[mpi_nrequests_]));
      mpi_nrequests_++;
    }
  current_msg_size_ = msg_size;
#endif
}

/*! @brief Waits until all communications started by send_recv_start are finished.
 *
 */
void Comm_Group_MPI::send_recv_finish() const
{
#ifdef MPI_
  assert(mpi_nrequests_ >= 0);
  mpi_error(MPI_Waitall(mpi_nrequests_, mpi_requests_, mpi_status_));
  if (statistics().is_gpu_verbose_on() && Process::je_suis_maitre()) // Display
    {
      std::string clock(Process::is_parallel() ? "[clock]#" + std::to_string(Process::me()) : "[clock]  ");
      double ms = 0.001 * statistics().get_time_since_last_open(STD_COUNTERS::mpi_sendrecv) ;
      printf("%s %7.3f ms [MPI]   Comm_Group_MPI::exchange\n", clock.c_str(), ms);
      fflush(stdout);
    }
  statistics().end_count(STD_COUNTERS::mpi_sendrecv,mpi_nrequests_,current_msg_size_);
  /*
  for (int r=0;r<mpi_nrequests_;r++)
    {
      if ( mpi_requests_[r]!=MPI_REQUEST_NULL)
  {
    MPI_Request_free(&(mpi_requests_[r]));
  }
      mpi_requests_[r]=MPI_REQUEST_NULL;
      }
  */
  mpi_nrequests_ = -1;
#endif
}


/*! @brief Blocking send.
 *
 * To be sure the code is safe, we force
 *  a synchronous communication to enforce blocking in check mode (see check_enabled()).
 *  Otherwise, we use MPI_Send which is generally non-blocking for small messages
 *  (hence better performance).
 *
 */
void Comm_Group_MPI::send(int pe, const void *buffer, int size, int tag) const
{
#ifdef MPI_
  statistics().begin_count(STD_COUNTERS::mpi_send);
  assert(mpi_nrequests_ < 0);
  int dest = pe;
  assert(dest >= 0 && dest < nproc());
  // Problem: forced to cast (const void*) to (void*) because of
  // the MPI_Send prototype
  if (check_enabled())
    mpi_error(MPI_Ssend ((void*)buffer, size, MPI_CHAR, dest, tag, mpi_comm_));
  else
    mpi_error(MPI_Send ((void*)buffer, size, MPI_CHAR, dest, tag, mpi_comm_));
  statistics().end_count(STD_COUNTERS::mpi_send,1,size);
#endif
}

/*! @brief Blocking reception of a message.
 *
 */
void Comm_Group_MPI::recv(int pe, void *buffer, int size, int tag) const
{
#ifdef MPI_
  statistics().begin_count(STD_COUNTERS::mpi_recv);
  assert(mpi_nrequests_ < 0);
  MPI_Status status;
  int source = pe;
  assert(source >= 0 && source < nproc());
  mpi_error(MPI_Recv (buffer, size, MPI_CHAR, source, tag, mpi_comm_, & status));
  statistics().end_count(STD_COUNTERS::mpi_recv,1,size);
#endif
}

void Comm_Group_MPI::broadcast(void *buffer, int size, int pe_source) const
{
#ifdef MPI_
  statistics().begin_count(STD_COUNTERS::mpi_bcast);
  assert(mpi_nrequests_ < 0);
  mpi_error(MPI_Bcast (buffer, size, MPI_CHAR, pe_source, mpi_comm_));
  statistics().end_count(STD_COUNTERS::mpi_bcast,1,size);
#endif
}

void Comm_Group_MPI::all_to_all(const void *src_buffer, void *dest_buffer, int data_size) const
{
#ifdef MPI_
  statistics().begin_count(STD_COUNTERS::mpi_alltoall);
  assert(src_buffer != dest_buffer);
  void * ptr = (void *) src_buffer; // Cast needed because of the MPI_Alltoall interface
  mpi_error(MPI_Alltoall(ptr, data_size, MPI_CHAR, dest_buffer, data_size, MPI_CHAR, mpi_comm_));
  statistics().end_count(STD_COUNTERS::mpi_alltoall,1,data_size);
#endif
}

void Comm_Group_MPI::gather(const void *src_buffer, void *dest_buffer, int data_size, int root) const
{
#ifdef MPI_
  statistics().begin_count(STD_COUNTERS::mpi_gather);
  void * ptr = (void *) src_buffer; // Cast needed because of the MPI_Alltoall interface
  mpi_error(MPI_Gather(ptr, data_size, MPI_CHAR, dest_buffer, data_size, MPI_CHAR, root, mpi_comm_));
  statistics().end_count(STD_COUNTERS::mpi_gather,1,data_size);
#endif
}

void Comm_Group_MPI::all_gather(const void *src_buffer, void *dest_buffer, int data_size) const
{
#ifdef MPI_
  statistics().begin_count(STD_COUNTERS::mpi_allgather);
  void * ptr = (void *) src_buffer; // Cast needed because of the MPI_Alltoall interface
  mpi_error(MPI_Allgather(ptr, data_size, MPI_CHAR, dest_buffer, data_size, MPI_CHAR, mpi_comm_));
  statistics().end_count(STD_COUNTERS::mpi_allgather,1,data_size);
#endif
}

void Comm_Group_MPI::all_gatherv(const void *src_buffer, void *dest_buffer, int send_size, const int* recv_size, const int* displs) const
{
#ifdef MPI_
  statistics().begin_count(STD_COUNTERS::mpi_allgather);
  void * ptr = (void *) src_buffer; // Cast needed because of the MPI_Alltoall interface
  mpi_error(MPI_Allgatherv(ptr, send_size, MPI_CHAR, dest_buffer, recv_size, displs, MPI_CHAR, mpi_comm_));
  statistics().end_count(STD_COUNTERS::mpi_allgather,1,send_size);
#endif
}


#ifdef MPI_

/*! @brief Constructor for the "all" group. Warning: this constructor must only be called once.
 *
 *   The group is associated with trio_u_world_.
 *   If must_mpi_initialize_==false, it is assumed that MPI_Init has already been called.
 *   After calling init_group_trio, the group must be registered in PE_Groups.
 *   See PE_Groups::initialize().
 *
 */
void Comm_Group_MPI::init_group_trio()
{
  must_finalize_ = 0;

  if (mpi_status_ != 0)
    {
      Cerr << "Error : the construction of the global Comm_Group_MPI has already been done." << finl;
      exit();
    }

  if (must_mpi_initialize_)
    {
      if (trio_u_world_ != MPI_COMM_WORLD)
        {
          Cerr << "Error in Comm_Group_MPI::init_group_trio(...) : you cannot ask to initialize MPI\n"
               << " with something else than MPI_COMM_WORLD !" << finl;
          exit();
        }
      must_finalize_ = 1;
      int argc=0;
      char** argv=nullptr;
      int errcode = MPI_Init(&argc, &argv);
      //int errcode = MPI_Init(0,0); Error message on MPI Voltaire
      if (errcode != MPI_SUCCESS)
        {
          Cerr << "Error in Comm_Group_MPI::init_group_trio()\n"
               << " MPI_Init() failed (forget to run with mpirun ?)" << finl;
          mpi_error(errcode);
        }
    }

  int arank;
  int nbproc;

  mpi_error(MPI_Comm_size (trio_u_world_, & nbproc));
  mpi_error(MPI_Comm_rank (trio_u_world_, & arank));

  Comm_Group::init_group_trio(nbproc, arank);

  mpi_comm_ = trio_u_world_;
  MPI_Comm_group(mpi_comm_, &mpi_group_);

  // Initialization of the static member variables of the class
  // One send buffer and one receive buffer per processor,
  // hence the maximum...
  mpi_maxrequests_ = nbproc * 2;
  mpi_status_ = new MPI_Status[mpi_maxrequests_];
  mpi_requests_ = new MPI_Request[mpi_maxrequests_];

  for (int r=0; r<mpi_maxrequests_; r++)
    {
      mpi_requests_[r]=MPI_REQUEST_NULL;
    }
  if (arank == 0)
    {
      if (trio_u_world_ == MPI_COMM_WORLD)
        {
          Cerr << "Initialized MPI with MPI_COMM_WORLD (using all processors)" << finl;
        }
      else
        {
          Cerr << "Initialized MPI with communicator!=MPI_COMM_WORLD: using " << (int)nbproc << " processors" << finl;
        }
    }
}


// MPI_Group_free should be done before MPI_Finalize so not included into Comm_Group_MPI destructor
void Comm_Group_MPI::free()
{
  if (mpi_maxrequests_!=-1) // Group is created when mpi_maxrequests_>0 (avoid a crash with verifie_pere script)
    mpi_error(MPI_Group_free(& mpi_group_));
}


/*! @brief Free group and MPI communicator (to use for MPI subgroups only, MPI_COMM_WORLD can no be freed)
 *
 */
void Comm_Group_MPI::free_all()
{
  if (mpi_maxrequests_!=-1)
    {
      if (mpi_group_!=MPI_GROUP_NULL)
        mpi_error(MPI_Group_free(& mpi_group_));
      if (mpi_comm_!=MPI_COMM_NULL)
        mpi_error(MPI_Comm_free(&mpi_comm_));
    }
}


// Wrapper to MPI_Alltoallv. data type is MPI_CHAR
void Comm_Group_MPI::all_to_allv(const void *src_buffer, int *send_data_size, int *send_data_offset,
                                 void *dest_buffer, int *recv_data_size, int *recv_data_offset) const
{
  statistics().begin_count(STD_COUNTERS::mpi_alltoall);
  assert(src_buffer != dest_buffer);
  void * ptr = (void *) src_buffer; // Cast needed because of the MPI_Alltoall interface

  const int n = nproc();
  int size;

#ifdef INT_is_64_
  std::vector<int> send_data_size_int(n);
  std::vector<int> send_data_offset_int(n);
  std::vector<int> recv_data_size_int(n);
  std::vector<int> recv_data_offset_int(n);

  auto cast_func = [](int i) -> int { return static_cast<int>(i); };
  std::transform(send_data_size,   send_data_size + n,   send_data_size_int.begin(),   cast_func);
  std::transform(send_data_offset, send_data_offset + n, send_data_offset_int.begin(), cast_func);
  std::transform(recv_data_size,   recv_data_size + n,   recv_data_size_int.begin(),   cast_func);
  std::transform(recv_data_offset, recv_data_offset + n, recv_data_offset_int.begin(), cast_func);

  mpi_error(MPI_Alltoallv(ptr, send_data_size_int.data(), send_data_offset_int.data(), MPI_CHAR,
                          dest_buffer, recv_data_size_int.data(), recv_data_offset_int.data(), MPI_CHAR, mpi_comm_));
  size = send_data_offset_int[n-1] + send_data_size_int[n-1] + recv_data_size_int[n-1] + recv_data_offset_int[n-1];
#else
  mpi_error(MPI_Alltoallv(ptr, send_data_size, send_data_offset, MPI_CHAR,
                          dest_buffer, recv_data_size, recv_data_offset, MPI_CHAR, mpi_comm_));
  size = send_data_offset[n-1] + send_data_size[n-1] + recv_data_size[n-1] + recv_data_offset[n-1];

#endif
  statistics().end_count(STD_COUNTERS::mpi_alltoall,1,size);
}

/*! @brief So that trio_u uses only a subset of MPI_COMM_WORLD processors, a communicator must be provided before
 *
 *   calling init_group_trio.
 *
 */
void Comm_Group_MPI::set_trio_u_world(MPI_Comm world)
{
  if (mpi_status_ != 0)
    {
      Cerr << "Error : the construction of the global Comm_Group_MPI has already been done\n"
           << " set_trio_u_world call is forbidden" << finl;
      exit();
    }
#ifdef PETSCKSP_H
  PETSC_COMM_WORLD= world;
#endif
  trio_u_world_ = world;
}

MPI_Comm Comm_Group_MPI::get_trio_u_world()
{
  return trio_u_world_;
}



void Comm_Group_MPI::set_must_mpi_initialize(bool flag)
{
  if (mpi_status_ != 0)
    {
      Cerr << "Error : the construction of the global Comm_Group_MPI has already been done\n"
           << " set_must_mpi_initialize() call is forbidden." << finl;
      exit();
    }
  must_mpi_initialize_ = flag;
}

void Comm_Group_MPI::ptop_send_recv(const void * send_buf, int send_buf_size, int send_proc,
                                    void * recv_buf, int recv_buf_size, int recv_proc) const
{
  statistics().begin_count(STD_COUNTERS::mpi_sendrecv);
  assert(mpi_nrequests_ < 0);
  int dest = send_proc;
  int src = recv_proc;
  int tag = 1;
  MPI_Status status;
  if (send_proc < 0 && recv_proc < 0)
    {
      // do nothing
    }
  else if (send_proc < 0 && recv_proc >= 0)
    {
      mpi_error(MPI_Recv (recv_buf, recv_buf_size, MPI_CHAR, src, tag, mpi_comm_, &status));
    }
  else if (recv_proc < 0 && send_proc >= 0)
    {
      mpi_error(MPI_Send ((void*)send_buf, send_buf_size, MPI_CHAR, dest, tag, mpi_comm_));
    }
  else
    {
      assert(dest >= 0 && dest < nproc());
      assert(src >= 0 && src < nproc());
      // Problem: forced to cast (const void*) to (void*) because of
      // the MPI_Send prototype
      mpi_error(MPI_Sendrecv((void*)send_buf, send_buf_size, MPI_CHAR, dest, tag,
                             recv_buf, recv_buf_size, MPI_CHAR, src, tag, mpi_comm_,
                             &status));
    }
  statistics().end_count(STD_COUNTERS::mpi_sendrecv, 1, send_buf_size + recv_buf_size);
}

/*! @brief Builds the processor group from the list.
 *
 * See Comm_Group::init_group(const ArrOfInt &)
 *  Method called by PE_Groups::create_group()
 *
 */
void Comm_Group_MPI::init_group(const ArrOfInt& pe_list)
{
  must_finalize_ = 0;
  // The "all" group must exist
  assert(mpi_status_);

  Comm_Group::init_group(pe_list);

  const Comm_Group_MPI& cg = ref_cast(Comm_Group_MPI, PE_Groups::current_group());
  // Store a reference to the parent group: it is the current group at the time
  // of the call to init_group. The destructor must be called simultaneously
  // on all processors in the same group.
  groupe_pere_ = PE_Groups::current_group();
  // Build the MPI group
  const MPI_Group& current_mpi_group = cg.mpi_group_;
  const MPI_Comm& current_mpi_comm  = cg.mpi_comm_;
  // Copy pe_list in case int != int...
  const int nbproc = this->nproc();
  int *ranks = new int[nbproc];
  for (int i = 0; i < nbproc; i++)
    ranks[i] = pe_list[i];
  assert(mpi_group_==MPI_GROUP_NULL);
  mpi_error(MPI_Group_incl(current_mpi_group, nbproc, ranks, & mpi_group_));
  delete[] ranks;
  // Build the communicator
  // MPI_Comm_create returns MPI_COMM_NULL if the current processor
  // is not in the group.
  mpi_error(MPI_Comm_create(current_mpi_comm, mpi_group_, & mpi_comm_));
}



/*! @brief Building MPI communicator based on numa node (ie one communicator for each node)
 *
 */
void Comm_Group_MPI::init_comm_on_numa_node()
{
  must_finalize_ = 0;
  // The "all" group must exist
  assert(mpi_status_);
  assert(mpi_group_==MPI_GROUP_NULL);

  groupe_pere_ = PE_Groups::current_group();

  // Build the communicator
  const Comm_Group_MPI& cg = ref_cast(Comm_Group_MPI, PE_Groups::current_group());
  const MPI_Comm& current_mpi_comm  = cg.mpi_comm_;
  int current_rank = cg.rank();
  mpi_error(MPI_Comm_split_type(current_mpi_comm, MPI_COMM_TYPE_SHARED, current_rank, MPI_INFO_NULL, &mpi_comm_));
  mpi_error(MPI_Comm_group(mpi_comm_, &mpi_group_));

  int loc_rank;
  int nbproc;
  mpi_error(MPI_Comm_size(mpi_comm_, &nbproc));
  mpi_error(MPI_Comm_rank(mpi_comm_, &loc_rank));

  Comm_Group::init_group_node(nbproc, loc_rank, current_rank);

  // Getting rank of my node among all the other nodes:
  // we create a temporary communicator which gathers all masters of each node group
  // so that the rank of my node is the rank of my master inside this temporary communicator
  int master = loc_rank==0 ? 0 : MPI_UNDEFINED;
  MPI_Comm tmp;
  mpi_error(MPI_Comm_split(current_mpi_comm, master, current_rank, &tmp));
  if(tmp != MPI_COMM_NULL)
    {
      mpi_error(MPI_Comm_rank(tmp, &node_id_));
      mpi_error(MPI_Comm_size(tmp, &nb_nodes_));
    }
  // each master broadcasts id and size to their group
  mpi_error(MPI_Bcast(&node_id_, 1,  MPI_INT, 0, mpi_comm_));
  mpi_error(MPI_Bcast(&nb_nodes_, 1,  MPI_INT, 0, mpi_comm_));

  if (tmp!= MPI_COMM_NULL)
    mpi_error(MPI_Comm_free(&tmp));

}

/*! @brief Building MPI communicator containing only the master of my numa node (ie one different communicator for each node)
 *
 */
void Comm_Group_MPI::init_comm_on_node_master()
{
  must_finalize_ = 0;
  // The "all" group must exist
  assert(mpi_status_);
  assert(mpi_group_==MPI_GROUP_NULL);

  groupe_pere_ = PE_Groups::get_node_group();

  // Build the MPI communicator and group
  const Comm_Group_MPI& cg = ref_cast(Comm_Group_MPI, PE_Groups::get_node_group());
  const MPI_Comm& current_mpi_comm  = cg.mpi_comm_;
  const MPI_Group& current_mpi_group  = cg.mpi_group_;
  int master = 0;
  mpi_error(MPI_Group_incl(current_mpi_group, 1, &master, & mpi_group_));
  mpi_error(MPI_Comm_create(current_mpi_comm, mpi_group_, & mpi_comm_));

  int world_rank = ref_cast(Comm_Group_MPI, PE_Groups::current_group()).rank();
  int loc_rank = cg.rank() == 0 ? 0 : -1;
  Comm_Group::init_group_node(1, loc_rank, world_rank);
}

void Comm_Group_MPI::internal_collective(const int *x, int *resu, int nx, const Collective_Op *op, int nop, int level) const
{
  // For now, brute-force algorithm, to be optimized...
  for (int i = 0; i < nx; i++)
    {
      int j = (nop < 0) ? 0 : i;
      trustIdType xx = x[i], resu2 = -1;
      if (op[j] != COLL_PARTIAL_SUM)
        mp_collective_op(&xx, &resu2, 1, op[j]);
      else
        resu2 = mppartial_sum_impl(x[i]);
      assert(resu2 < std::numeric_limits<int>::max());
      resu[i] = static_cast<int>(resu2);
    }
}

#if INT_is_64_ == 2
void Comm_Group_MPI::internal_collective(const trustIdType *x, trustIdType *resu, int nx, const Collective_Op *op, int nop, int level) const
{
  // For now, brute-force algorithm, to be optimized...
  for (int i = 0; i < nx; i++)
    {
      int j = (nop < 0) ? 0 : i;
      if (op[j] != COLL_PARTIAL_SUM)
        mp_collective_op(x+i, resu+i, 1, op[j]);
      else
        resu[i] = mppartial_sum_impl(x[i]);
    }
}
#endif


void Comm_Group_MPI::internal_collective(const double *x, double *resu, int nx, const Collective_Op *op, int nop, int level) const
{
  // For now, brute-force algorithm, to be optimized...
  for (int i = 0; i < nx; i++)
    {
      int j = (nop < 0) ? 0 : i;
      if (op[j] != COLL_PARTIAL_SUM)
        mp_collective_op(x+i, resu+i, 1, op[j]);
      else
        {
          Cerr << "Error in Comm_Group_MPI: COLL_PARTIAL_SUM not coded for double" << finl;
          exit();
        }
    }
}

void Comm_Group_MPI::internal_collective(const float *x, float *resu, int nx, const Collective_Op *op, int nop, int level) const
{
  // For now, brute-force algorithm, to be optimized...
  for (int i = 0; i < nx; i++)
    {
      int j = (nop < 0) ? 0 : i;
      if (op[j] != COLL_PARTIAL_SUM)
        mp_collective_op(x+i, resu+i, 1, op[j]);
      else
        {
          Cerr << "Error in Comm_Group_MPI: COLL_PARTIAL_SUM not coded for float" << finl;
          exit();
        }
    }
}

/*! @brief Returns the sum of x over the preceding processors in the group (not including self).
 *
 * The result on the first processor of the group is therefore always 0.
 *  The result depends on the order in which the processors were
 *  provided in the constructor.
 *
 */
trustIdType Comm_Group_MPI::mppartial_sum_impl(trustIdType x) const
{
  statistics().begin_count(STD_COUNTERS::mpi_partialsum);
  trustIdType somme = 0;
  MPI_Status status;
  int tag = get_new_tag();
  int rang = rank();
  int np = nproc();

  if (rang > 0)
    {
      // Receives the partial sum from the previous processor
#ifndef INT_is_64_
      mpi_error(MPI_Recv(& somme, 1, MPI_INT, rang-1, tag, mpi_comm_, &status));
#else
      mpi_error(MPI_Recv(& somme, 1, MPI_LONG, rang-1, tag, mpi_comm_, &status));
#endif
    }
  if (rang+1 < np)
    {
      // Sends the partial sum to the next processor
      trustIdType s = somme + x;
#ifndef INT_is_64_
      mpi_error(MPI_Send(& s, 1, MPI_INT, rang+1, tag, mpi_comm_));
#else
      mpi_error(MPI_Send(& s, 1, MPI_LONG, rang+1, tag, mpi_comm_));
#endif
    }
  statistics().end_count(STD_COUNTERS::mpi_partialsum);
  return somme;
}

#endif
