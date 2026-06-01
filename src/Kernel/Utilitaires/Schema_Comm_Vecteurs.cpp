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
#include <Schema_Comm_Vecteurs.h>
#include <Comm_Group.h>
#include <communications.h>
#include <PE_Groups.h>
#include <sstream>
#include <comm_incl.h>

bool Schema_Comm_Vecteurs::buffer_locked_;
ArrOfDouble Schema_Comm_Vecteurs::tmp_area_double_;
ArrOfFloat Schema_Comm_Vecteurs::tmp_area_float_;
ArrOfInt Schema_Comm_Vecteurs::tmp_area_int_;
#if INT_is_64_ == 2
ArrOfTID Schema_Comm_Vecteurs::tmp_area_tid_;
#endif
Schema_Comm_Vecteurs_Static_Data Schema_Comm_Vecteurs::sdata_;
bool check_comm_vector = false;

Schema_Comm_Vecteurs_Static_Data::Schema_Comm_Vecteurs_Static_Data()
{
  buffer_base_ = 0;
  buffer_base_size_ = 0;
  buffer_base_device_size_ = 0;
  buf_pointers_ = 0;
  buf_pointers_size_ = 0;
}

void Schema_Comm_Vecteurs_Static_Data::init(int min_buf_size, bool bufferOnDevice)
{
  if (buf_pointers_size_ == 0)
    {
      int n = PE_Groups::groupe_TRUST().nproc();
      buf_pointers_size_ = n * 2; // Maximum required when exchanging with all procs in exchange()
      buf_pointers_ = new char*[n*2];
      for (int i = 0; i < n*2; i++)
        buf_pointers_[i] = 0;
    }
  // Is the global buffer large enough?
  if (buffer_base_size_ < min_buf_size)
    {
      if (buffer_base_device_size_ > 0)
        {
          deleteOnDevice(buffer_base_, buffer_base_device_size_);
          buffer_base_device_size_ = 0;
        }
      delete [] buffer_base_;
      buffer_base_ = new char[min_buf_size];
      // GF: zero-initialization added for mpiwrapper/valgrind but is it useful?
      for (int i = 0; i < min_buf_size; i++)
        buffer_base_[i] = 0;
      buffer_base_size_ = min_buf_size;
    }
  if (bufferOnDevice && buffer_base_device_size_ < min_buf_size)
    {
      // Allocate buffer_base_ on device:
      if (buffer_base_device_size_>0)
        deleteOnDevice(buffer_base_, buffer_base_device_size_);
      allocateOnDevice(buffer_base_, min_buf_size);
      buffer_base_device_size_ = min_buf_size;
    }
}

Schema_Comm_Vecteurs_Static_Data::~Schema_Comm_Vecteurs_Static_Data()
{
  /* ToDo OpenMP Fix crash when using AmgX: Failing in Thread:0
  call to cuInit returned error 4: Deinitialized
   deleteOnDevice(buffer_base_, buffer_base_size_);
  */
  delete[] buffer_base_;
  delete[] buf_pointers_;
}

Schema_Comm_Vecteurs::Schema_Comm_Vecteurs()
{
  status_ = RESET;
  const char* env_var = getenv("TRUST_USE_MPI_GPU_AWARE");
  use_gpu_aware_mpi_ = env_var != nullptr && std::stoi(env_var) == 1;
  if (use_gpu_aware_mpi_)
    {
#ifdef CRAY_MPICH_VER
      if (getenv("MPICH_GPU_SUPPORT_ENABLED") == nullptr)
        Process::exit("You try to enable GPU communications on Cray MPICH with TRUST_USE_MPI_GPU_AWARE=1 but forgot to set also MPICH_GPU_SUPPORT_ENABLED=1 !");
#endif
      std::cerr << "[MPI] Enabling GPU capability to communicate between devices." << std::endl;
      //Cerr << "[MPI] Warning! Only MPI calls with device pointers will benefit. Classic MPI calls with host pointers will be slower..." << finl;
    }
}

Schema_Comm_Vecteurs::~Schema_Comm_Vecteurs()
{
  assert (status_ == END_INIT || status_ == RESET);
}

/*! @brief Resets buffer sizes.
 *
 * Buffer sizes must then be defined with add_send/recv_area_...().
 *   This method must be called simultaneously on all processors in the group.
 *
 */
void Schema_Comm_Vecteurs::begin_init()
{
  assert(status_ == END_INIT || status_ == RESET);
  // Reset the sizes_ arrays
  const int np = Process::nproc();
  send_buf_sizes_.resize_array(np, RESIZE_OPTIONS::NOCOPY_NOINIT);
  send_buf_sizes_ = 0;
  recv_buf_sizes_.resize_array(np, RESIZE_OPTIONS::NOCOPY_NOINIT);
  recv_buf_sizes_ = 0;
  send_procs_.resize_array(0);
  recv_procs_.resize_array(0);
  sorted_ = 1;
  status_ = BEGIN_INIT;
  if (use_gpu_aware_mpi_)
    {
#if defined(TRUST_USE_CUDA) && !defined(MPIX_CUDA_AWARE_SUPPORT)
      Process::exit("MPI version is detected as not CUDA-Aware. You can't use TRUST_USE_MPI_GPU_AWARE=1");
#endif
    }
}

/*! @brief Once the data to exchange has been declared with add_send/recv_area_..(),
 *
 * initializes buffer offsets and allocates a global buffer of sufficient size.
 *   Must be called by all processors in the group.
 *
 */
void Schema_Comm_Vecteurs::end_init()
{
  assert(status_ == BEGIN_INIT);
  assert(Process::nproc() == send_buf_sizes_.size_array());
  // Verification of send and receive buffer sizes:
  if (Comm_Group::check_enabled())
    {
      const int n = send_buf_sizes_.size_array();
      ArrOfInt tmp(n);
      envoyer_all_to_all(send_buf_sizes_, tmp);
      int err = 0;
      for (int i = 0; i < n; i++)
        if (tmp[i] != recv_buf_sizes_[i])
          err++;
      if (Process::mp_sum(err))
        {
          Cerr << "Error in Schema_Comm_Vecteurs::end_init(): send_size_ and recv_size_ don't match, see log files" << finl;
          Process::Journal() << "Error in Schema_Comm_Vecteurs::end_init():\n"
                             << "send_sizes_ = " << send_buf_sizes_
                             << "\n recv_sizes_ = " << recv_buf_sizes_ << finl;
          Process::barrier();
          Process::exit();
        }
    }
  // Sort processor indices in ascending order
  if (!sorted_)
    {
      send_procs_.ordonne_array();
      recv_procs_.ordonne_array();
    }
  const int nsend = send_procs_.size_array();
  const int nrecv = recv_procs_.size_array();
  int offset = 0;
  int i;
  for (i = 0; i < nsend; i++)
    {
      int pe = send_procs_[i];
      const int size = (send_buf_sizes_[pe]+7)&(~7); // align size on 8 bytes
      offset += size;
      assert(pe >= i); // send_procs_ sorted in ascending order
      send_buf_sizes_[i] = size;
    }
  for (i = 0; i < nrecv; i++)
    {
      int pe = recv_procs_[i];
      const int size = (recv_buf_sizes_[pe]+7)&(~7); // align size on 8 bytes
      offset += size;
      assert(pe >= i); // recv_procs_ sorted in ascending order
      recv_buf_sizes_[i] = size;
    }
  min_buf_size_ = offset; // buffer size to allocate    recv_buf_offset_[i] = offset;

  send_buf_sizes_.resize_array(nsend);
  recv_buf_sizes_.resize_array(nrecv);
  status_ = END_INIT;
}

/*! @brief Starts a new data exchange (buffer sizes must have been initialized with begin_init()...end_init()).
 *
 * Sets sdata_.buf_pointers_ to the beginning of the buffers for each processor
 *    for which a "send" buffer was declared.
 *   After begin_comm(), the buffers must be filled using
 *    get_next_area_int() or get_next_area_double() in the same order
 *    as declared during the initialization phase, then exchange() must be called.
 *
 */
void Schema_Comm_Vecteurs::begin_comm(bool bufferOnDevice)
{
  assert(status_ == END_INIT);
  // Not an assert because the error is serious and presumably rare...
  if (buffer_locked_)
    {
      Cerr << "Internal error in Schema_Comm_Vecteurs::begin_comm(): buffers already locked by another communication" << finl;
      Process::exit();
    }
  buffer_locked_ = true;
  sdata_.init(min_buf_size_, bufferOnDevice);

  // Point the buffers to the beginning of the send_buffers
  char *ptr = sdata_.buffer_base_;

  const int nsend = send_procs_.size_array();
  for (int i = 0; i < nsend; i++)
    {
      const int pe = send_procs_[i];
      sdata_.buf_pointers_[pe] = ptr;
      ptr += send_buf_sizes_[i];
    }
  buffer_locked_ = true;
  status_ = BEGIN_COMM;
  bufferOnDevice_ = bufferOnDevice;
}

void Schema_Comm_Vecteurs::exchange(IsExchangeBlocking exchange_type, const std::string kernel_name)
{

  char * ptr = sdata_.buffer_base_;
  const Comm_Group& group = PE_Groups::current_group();
  const int nsend = send_procs_.size_array();
  const int nrecv = recv_procs_.size_array();

  if ((exchange_type == IsExchangeBlocking::DefaultBlocking)||(exchange_type == IsExchangeBlocking::NonBlockingStart))
    {

      // Copy buffer before MPI send
      if (bufferOnDevice_)
        {
          if (!use_gpu_aware_mpi_)
            copyFromDevice(sdata_.buffer_base_, min_buf_size_); // Copy buffer to host for MPI communication
          else
            {
              // Communication between devices. Use device buffer:
              ptr = addrOnDevice(sdata_.buffer_base_);
            }
        }

      assert(status_ == BEGIN_COMM);
      // Check that all buffers are full
      assert(check_buffers_full());
      // Exchange the data

      // Use the sdata_.buf_pointers_ array to store the addresses
      //  of buffers to pass to Comm_Group::send_recv_start()
      // (dimensioned to 2*nproc() so sufficient)
      assert(nsend + nrecv <= sdata_.buf_pointers_size_);
      char ** send_bufs = sdata_.buf_pointers_;
      char ** recv_bufs = sdata_.buf_pointers_ + nsend;
      for (int i = 0; i < nsend; i++)
        {
          send_bufs[i] = ptr;
          ptr += send_buf_sizes_[i];
        }
      for (int i = 0; i < nrecv; i++)
        {
          recv_bufs[i] = ptr;
          ptr += recv_buf_sizes_[i];
        }

      // We should be able to use int64 as type here because
      // the buffers are aligned on 8 bytes.

      if (exchange_type == IsExchangeBlocking::NonBlockingStart) start_gpu_timer(kernel_name);
      group.send_recv_start(send_procs_, send_buf_sizes_, send_bufs,
                            recv_procs_, recv_buf_sizes_, recv_bufs,
                            Comm_Group::INT);
    }

  if ((exchange_type == IsExchangeBlocking::DefaultBlocking)||(exchange_type == IsExchangeBlocking::NonBlockingFinish))
    {
      group.send_recv_finish();
      if (exchange_type == IsExchangeBlocking::NonBlockingFinish) end_gpu_timer(kernel_name);
      // Point the buffers to the received data
      char * recv_ptr = sdata_.buffer_base_;
      for (int i = 0; i < nsend; i++)
        recv_ptr += send_buf_sizes_[i];
      for (int i = 0; i < nrecv; i++)
        {
          const int pe = recv_procs_[i];
          sdata_.buf_pointers_[pe] = recv_ptr;
          recv_ptr += recv_buf_sizes_[i];
        }
      status_ = EXCHANGED;

      // Copy buffer to device after MPI recv if GPU-Aware MPI is not enabled:
      if (bufferOnDevice_ && !use_gpu_aware_mpi_) copyToDevice(sdata_.buffer_base_, min_buf_size_);
    }
}

void Schema_Comm_Vecteurs::end_comm()
{
  assert(status_ == EXCHANGED);
  // Check that all data has been read
  assert(check_buffers_full());
  status_ = END_INIT; // ready for a new begin_comm()
  buffer_locked_ = false;
  bufferOnDevice_ = false;
}

/*! @brief Depending on status_, verifies that all buffer pointers point to the end of the buffer allocated for each processor in send
 *
 *   or receive mode. Returns 0 on error (if a buffer has not been
 *   completely filled or emptied).
 *
 */
int Schema_Comm_Vecteurs::check_buffers_full() const
{
  char *ptr = sdata_.buffer_base_;
  const int nsend = send_procs_.size_array();
  int i;
  int ok = 1;
  if (status_ == BEGIN_COMM)
    {
      for (i = 0; i < nsend; i++)
        {
          ptr += send_buf_sizes_[i];
          const int pe = send_procs_[i];
          char *ptr2 = sdata_.buf_pointers_[pe];
          ALIGN_SIZE(ptr2, sizeof(double));
          if (ptr != ptr2)
            {
              Cerr << "Internal error in Schema_Comm_Vecteurs::check_buffers_full(): send buffer for processor "
                   << pe << " is not full" << finl;
              ok = 0;
            }
        }
    }
  else if (status_ == EXCHANGED)
    {
      for (i = 0; i < nsend; i++)
        ptr += send_buf_sizes_[i];
      const int nrecv = recv_procs_.size_array();
      for (i = 0; i < nrecv; i++)
        {
          ptr += recv_buf_sizes_[i];
          const int pe = recv_procs_[i];
          char *ptr2 = sdata_.buf_pointers_[pe];
          ALIGN_SIZE(ptr2, sizeof(double));
          if (ptr != ptr2)
            {
              Cerr << "Internal error in Schema_Comm_Vecteurs::check_buffers_full(): recv buffer for processor "
                   << pe << " has not been read entirely" << finl;
              ok = 0;
            }
        }
    }
  else
    {
      Cerr << "check_buffers_full: What ?" << finl;
      Process::exit();
    }
  return ok;
}

/*! @brief Verifies that there are at least byte_size bytes remaining in the buffer of processor pe.
 *
 */
int Schema_Comm_Vecteurs::check_next_area(int pe, int byte_size) const
{
  assert(byte_size >= 0);
  if (byte_size == 0)
    {
      return 1;
    }
  const ArrOfInt& procs = (status_ == BEGIN_COMM) ? send_procs_ : recv_procs_;
  const ArrOfInt& sizes = (status_ == BEGIN_COMM) ? send_buf_sizes_ : recv_buf_sizes_;
  const int n = procs.size_array();
  int i;
  char * ptr = sdata_.buffer_base_;
  // If in reception phase, the start of the receive buffers is located at
  // the end of the send buffers:
  if (status_ != BEGIN_COMM)
    {
      const int nsend = send_procs_.size_array();
      for (i = 0; i < nsend; i++)
        ptr += send_buf_sizes_[i];
    }
  for (i = 0; i < n; i++)
    {
      ptr += sizes[i]; // pointer to the end of this processor's buffer
      if (procs[i] == pe)
        return (sdata_.buf_pointers_[pe] + byte_size) <= ptr;
    }
  // no buffer declared for this processor
  return 0;
}
