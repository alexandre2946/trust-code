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
#include <Schema_Comm.h>
#include <PE_Groups.h>
#include <InOutCommBuffers.h>
#include <Comm_Group.h>
#include <communications.h>
#include <Comm_Group_MPI.h>

// ====================================================================
//                      class Schema_Comm
// ====================================================================

Schema_Comm::Static_Status Schema_Comm::status_ = UNINITIALIZED;
InOutCommBuffers   Schema_Comm::buffers_;
int                Schema_Comm::n_buffers_ = 0;

/*! @brief @brief Accessor to a member of the obuffers_ array (with verification).
 *
 * @param pe Processor index.
 * @return Reference to the OutputCommBuffer for processor pe.
 */
inline OutputCommBuffer& Schema_Comm::obuffer(int pe)
{
  assert(pe >= 0 && pe < n_buffers_);
  return buffers_.obuffers_[pe];
}

/*! @brief @brief Accessor to a member of the ebuffers_ array (with verification).
 *
 * @param pe Processor index.
 * @return Reference to the InputCommBuffer for processor pe.
 */
inline InputCommBuffer&   Schema_Comm::ebuffer(int pe)
{
  assert(pe >= 0 && pe < n_buffers_);
  return buffers_.ebuffers_[pe];
}

//static const int SET_GROUP_TAG = 0;
static const int BEGIN_COMM_TAG = 1;
static const int ECHANGE_MESSAGES_COMM_TAG = 2;
static const int END_COMM_TAG = 3;
static const int COPY_OPERATOR_TAG = 4;
//static const int CHECK_SEND_RCV_TAG = 5;

/*! @brief Constructs a new communication schema.
 *
 */
Schema_Comm::Schema_Comm()
{


  me_to_me_ = 0;
  use_all_to_allv_ = 0;
  // To verify later that we are still in the correct group,
  // we keep a reference to the current group.
  ref_group_ = PE_Groups::current_group();
  if (status_ == UNINITIALIZED)
    {
      n_buffers_ = Process::nproc();
      buffers_.obuffers_ = new OutputCommBuffer[n_buffers_];
      buffers_.ebuffers_ = new InputCommBuffer[n_buffers_];
      status_ = RESET;
    }
}

/*! @brief Destructs a communication schema.
 *
 */
Schema_Comm::~Schema_Comm()
{
}

/*! @brief Copy constructor (new schema placed in RESET mode).
 *
 * Warning: all members of the Comm_Group must execute this function simultaneously.
 *
 */
Schema_Comm::Schema_Comm(const Schema_Comm& schema)
{


  operator= (schema);
}

/*! @brief Copy operator: copies the list of communicating processors.
 *
 * The new schema is placed in the RESET state.
 *  Note: all members of the Comm_Group must execute this function simultaneously.
 *
 */
const Schema_Comm& Schema_Comm::operator=(const Schema_Comm& schema)
{
  assert(status_ == RESET);
  ref_group_ = schema.ref_group_;
  assert(&(ref_group_.valeur()) == &PE_Groups::current_group());
  send_pe_list_ = schema.send_pe_list_;
  recv_pe_list_ = schema.recv_pe_list_;
  me_to_me_     = schema.me_to_me_;
  use_all_to_allv_ = schema.use_all_to_allv_;
  const Comm_Group& group = ref_group_.valeur();
  if (group.check_enabled()) group.barrier(COPY_OPERATOR_TAG);
  return *this;
}

/*! @brief Obsolete method. The group associated with the schema is the current group at the time the schema is created.
 *
 * This method is only valid with the same group as the original group.
 *  The method does nothing.
 *
 */
void Schema_Comm::set_group(const Comm_Group& group)
{
  assert(&group == &(ref_group_.valeur()));
  assert(&group == &PE_Groups::current_group());
}

/*! @brief Returns the group associated with the schema.
 *
 */
const Comm_Group& Schema_Comm::get_group() const
{
  assert(ref_group_);
  return ref_group_.valeur();
}

/*! @brief Defines the list of processors to send data to and receive data from.
 *
 *  If me_to_me is non-zero, messages to oneself are allowed;
 *  otherwise not (optional argument: default me_to_me=0).
 *
 */
void Schema_Comm::set_send_recv_pe_list(const ArrOfInt& send_pe_list,
                                        const ArrOfInt& recv_pe_list,
                                        const int     me_to_me)
{
  assert(status_ == RESET);
  send_pe_list_ = send_pe_list;
  recv_pe_list_ = recv_pe_list;
  me_to_me_ = me_to_me;
  // Verification of the principle "you listen when I speak"
  const Comm_Group& group = ref_group_.valeur();
  assert(&group == &PE_Groups::current_group());
  if (group.check_enabled()) check_send_recv_pe_list();
}

/*! @brief Reserves communication buffers for a new communication.
 *
 * The schema transitions from RESET to WRITING; send_buffer() may now be called.
 * It is forbidden to call begin_comm() again on all communication objects
 * before finishing this communication with end_comm().
 *
 */
void Schema_Comm::begin_comm() const
{
  // Verify that no other communication is in progress.
  assert (status_ == RESET && ref_group_);
  status_ = WRITING;
  // Verify that all group members execute this.
  // If this crashes here, it means that not all declared members
  // are reaching the barrier.
  const Comm_Group& group = ref_group_.valeur();
  assert(&group == &PE_Groups::current_group());
  if (group.check_enabled()) group.barrier(BEGIN_COMM_TAG);
}

static void exchange_data(const ArrOfInt& send_list,
                          const ArrOfInt& send_size,
                          const char * const * const send_buffers,
                          const ArrOfInt& recv_list,
                          const ArrOfInt& recv_size,
                          char * const * const recv_buffers,
                          const Comm_Group& group,
                          bool use_all_to_all)
{
// GF added ifdef for builds without MPI where comm_group_mpi does not exist
#ifdef MPI_
  if (!use_all_to_all || !sub_type(Comm_Group_MPI, group))
#else
  //if (!use_all_to_all || !sub_type(Comm_Group_MPI, group))
  // do not use all_to_all as it does not exist without MPI
  if (1)
#endif
    {
      group.send_recv_start(send_list, send_size, send_buffers,
                            recv_list, recv_size, recv_buffers);
      group.send_recv_finish();
    }
  else
    {
#ifdef MPI_
      const int n = group.nproc();
      const int nsend = send_list.size_array();
      const int nrecv = recv_list.size_array();
      // Pack all send data in a single buffer:
      ArrOfInt a_send_size(n);
      ArrOfInt a_send_offset(n);
      ArrOfInt a_recv_size(n);
      ArrOfInt a_recv_offset(n);
      int i, offset;
      // Compute send_size array
      for (i = 0; i < nsend; i++)
        a_send_size[send_list[i]] = send_size[i];
      // Compute send_offset
      for (i = 0, offset = 0; i < n; i++)
        {
          a_send_offset[i] = offset;
          offset += a_send_size[i];
        }
      const int buf_size_send = offset;
      // Compute recv_size array
      for (i = 0; i < nrecv; i++)
        a_recv_size[recv_list[i]] = recv_size[i];
      // Compute recv_offset
      for (i = 0, offset = 0; i < n; i++)
        {
          a_recv_offset[i] = offset;
          offset += a_recv_size[i];
        }
      const int buf_size_recv = offset;
      // Allocate contiguous send and recv buffer:
      char *send_buffer = (char *) malloc(buf_size_send);
      char *recv_buffer = (char *) malloc(buf_size_recv);
      // Copy send data to send buffer
      for (i = 0; i < nsend; i++)
        memcpy(send_buffer + a_send_offset[send_list[i]], // dest
               send_buffers[i], // source
               send_size[i]); // size
      // Exchange data
      ref_cast(Comm_Group_MPI, group).all_to_allv(send_buffer, a_send_size.addr(), a_send_offset.addr(),
                                                  recv_buffer, a_recv_size.addr(), a_recv_offset.addr());
      // Copy recv data to recv_buffers
      for (i = 0; i < nrecv; i++)
        memcpy(recv_buffers[i], // dest
               recv_buffer + a_recv_offset[recv_list[i]], // source
               recv_size[i]); // size
      // [ABN] arrruhhhhhmmmmmm:
      free(send_buffer);
      free(recv_buffer);
#endif
    }
}

/*! @brief Transmits the size of messages to send to the processors that will receive them.
 *
 * The size is the number of bytes of the obuffers.
 *  send_pe_list and recv_pe_list must be initialized.
 *  The schema must be in the WRITING state.
 *
 */
void Schema_Comm::echange_taille(const ArrOfInt& send_size,
                                 ArrOfInt& recv_size) const
{
  static ArrOfInt send_sz;
  static ArrOfInt recv_sz;



  assert(status_ == WRITING && ref_group_);
  const Comm_Group& group = ref_group_.valeur();
  assert(&group == &PE_Groups::current_group());

  // Verify that all group members execute this.
  // If this crashes here, it means not all declared members are reaching the barrier.
  if (group.check_enabled()) group.barrier(ECHANGE_MESSAGES_COMM_TAG);

  const int n_send = send_pe_list_.size_array();
  const int n_recv = recv_pe_list_.size_array();

  recv_size.resize_array(n_recv);

  const char ** send_buffers = new const char* [n_send];
  char ** recv_buffers = new char* [n_recv];

  send_sz.resize_array(n_send);   // Size of an int (we exchange a size)
  send_sz = sizeof(int);
  int i;
  for (i = 0; i < n_send; i++)
    send_buffers[i] = (char*) (& send_size[i]);

  recv_sz.resize_array(n_recv);
  recv_sz = sizeof(int);
  recv_size.resize_array(n_recv);
  for (i = 0; i < n_recv; i++)
    recv_buffers[i] = (char*) (& recv_size[i]);

  exchange_data(send_pe_list_, send_sz, send_buffers,
                recv_pe_list_, recv_sz, recv_buffers,
                group,
                use_all_to_allv_);

  delete[] recv_buffers;
  delete[] send_buffers;
}

/*! @brief Launches the data exchange between all processors.
 *
 * The size of received messages must already be known.
 *  The schema transitions from WRITING to EXCHANGED.
 *
 */
void Schema_Comm::echange_messages(const ArrOfInt& send_size,
                                   const ArrOfInt& recv_size) const
{
  assert(status_ == WRITING && ref_group_);
  const Comm_Group& group = ref_group_.valeur();
  assert(&group == &PE_Groups::current_group());

  // Verify that all group members execute this.
  // If this crashes here, it means not all declared members are reaching the barrier.
  if (group.check_enabled()) group.barrier(ECHANGE_MESSAGES_COMM_TAG);

  const int n_send = send_pe_list_.size_array();
  const int n_recv = recv_pe_list_.size_array();
  const char ** send_buffers = new const char* [n_send];
  char ** recv_buffers = new char* [n_recv];

  int i;
  for (i = 0; i < n_recv; i++)
    {
      int pe = recv_pe_list_[i];
      int size = recv_size[i];
      InputCommBuffer& buf = ebuffer(pe);
      recv_buffers[i] = buf.reserve_buffer(size);
    }
  for (i = 0; i < n_send; i++)
    {
      int pe = send_pe_list_[i];
      OutputCommBuffer& buf = obuffer(pe);
      // Verify that the outgoing message size matches the size recorded in send_size_,
      // which guarantees that the receive size is also correct.
      assert(send_size[i] == buf.get_buffer_size());
      send_buffers[i] = buf.get_buffer();
    }

  exchange_data(send_pe_list_, send_size, send_buffers,
                recv_pe_list_, recv_size, recv_buffers,
                group,
                use_all_to_allv_);

  delete[] recv_buffers;
  delete[] send_buffers;

  // Create input streams from received buffers
  for (i = 0; i < n_recv; i++)
    {
      int pe = recv_pe_list_[i];
      InputCommBuffer& buf = ebuffer(pe);
      buf.create_stream();
    }

  // Special case of messages sent to oneself:
  if (me_to_me_)
    {
      int pe = Process::me();
      InputCommBuffer& buf = ebuffer(pe);
      OutputCommBuffer& obuf = obuffer(pe);
      buf.create_stream_from_output_stream(obuf);
    }

  status_ = EXCHANGED;
}

/*! @brief Launches the data exchange between all processors.
 *
 * The size of received messages does not need to be known a priori; it is transmitted.
 *  The schema transitions from WRITING to EXCHANGED.
 *
 */
void Schema_Comm::echange_taille_et_messages() const
{
  static ArrOfInt send_size;
  static ArrOfInt recv_size;
  const int n_send = send_pe_list_.size_array();
  send_size.resize_array(n_send);
  int i;
  for (i = 0; i < n_send; i++)
    {
      int pe = send_pe_list_[i];
      send_size[i] = obuffer(pe).get_buffer_size();
    }

  // Non-optimal method: in MPI one could use MPI_Probe but
  // this method seems unreliable.
  echange_taille(send_size, recv_size);
  echange_messages(send_size, recv_size);
}

/*! @brief Launches the data exchange.
 *
 * The size in bytes of received messages is provided in recv_size (array of the same size as recv_pe_list).
 *  In check_enabled mode, verifies that the size is correct.
 *
 */
void Schema_Comm::echange_messages(const ArrOfInt& recv_size) const
{
  const int n_send = send_pe_list_.size_array();
  ArrOfInt send_size(n_send);
  int i;
  for (i = 0; i < n_send; i++)
    {
      int pe = send_pe_list_[i];
      send_size[i] = obuffer(pe).get_buffer_size();
    }
  if (PE_Groups::current_group().check_enabled())
    {
      ArrOfInt check_recv_size;
      echange_taille(send_size, check_recv_size);
      if (!(check_recv_size == recv_size))
        {
          Cerr << "Error in Schema_Comm::echange_messages : bad recv_size" << finl;
          Process::exit();
        }
    }
  echange_messages(send_size, recv_size);
}

/*! @brief Clears the buffers and releases resources: reading of received data from buffers is complete.
 *
 *  The schema transitions from EXCHANGED to RESET.
 *
 */
void Schema_Comm::end_comm() const
{
  assert(status_ == EXCHANGED && ref_group_);
  const Comm_Group& group = ref_group_.valeur();
  assert(&group == &PE_Groups::current_group());

  // Check that all group members are executing this.
  // If it crashes here, it means not all declared members are at the barrier.
  if (group.check_enabled()) group.barrier(END_COMM_TAG);

  int i, n;
  n = send_pe_list_.size_array();
  for (i = 0; i < n; i++)
    {
      int pe = send_pe_list_[i];
      obuffer(pe).clear();
    }
  n = recv_pe_list_.size_array();
  for (i = 0; i < n; i++)
    {
      int pe = recv_pe_list_[i];
      ebuffer(pe).clear();
    }

  if (me_to_me_)
    {
      int pe = Process::me();
      obuffer(pe).clear();
      ebuffer(pe).clear();
    }

  status_ = RESET;
}

#ifndef NDEBUG
static int check_PE_in_list(int num_pe, const ArrOfInt& list)
{
  int i;
  int n = list.size_array();
  for (i = 0; i < n && list[i] != num_pe; i++);
  return (i < n);
}
#endif

/*! @brief Returns the buffer corresponding to processor num_PE to stack data to send.
 *
 * The schema must be in the WRITING state.
 *
 */
Sortie& Schema_Comm::send_buffer(int num_PE) const
{
  // If the following assert fails, it means we are trying to
  // put data in the buffer outside of the block
  //   begin_comm();
  //    ...
  //   echange_xxx();
  assert(status_ == WRITING && ref_group_);

  // Check that the requested PE is in the list of declared send PEs.
  assert((me_to_me_&&num_PE==Process::me()) || check_PE_in_list(num_PE, send_pe_list_));
  return obuffer(num_PE);
}

/*! @brief Returns the buffer corresponding to processor num_PE to read received data.
 *
 * The schema must be in the EXCHANGED state.
 *
 */
Entree& Schema_Comm::recv_buffer(int num_PE) const
{
  // If the following assert fails, it means we are trying to
  // read data from the buffers outside of the block
  //   echange_xxx();
  //    ...
  //   end_comm();
  assert(status_ == EXCHANGED && ref_group_);
  // Check that the requested PE is in the list of declared receive PEs.
  assert((me_to_me_&&num_PE==Process::me()) || check_PE_in_list(num_PE, recv_pe_list_));
  return ebuffer(num_PE);
}

const ArrOfInt& Schema_Comm::get_send_pe_list() const
{
  assert(ref_group_);
  return send_pe_list_;
}

const ArrOfInt& Schema_Comm::get_recv_pe_list() const
{
  assert(ref_group_);
  return recv_pe_list_;
}

/*! @brief Returns a reference to an array containing, for each processor in send_pe_list_, the size in bytes of the data
 *
 *   to send.
 *  TODO: TO FINISH !!!!
 *
 */
const ArrOfInt& Schema_Comm_statique::get_send_size() const
{
  assert(0);
  assert(status_ == EXCHANGED);
  return send_size_;
}

/*! @brief Returns a reference to an array containing, for each processor in send_pe_list_, the size in bytes of the data
 *
 *   received.
 *  TODO: TO FINISH !!!!
 *
 */
const ArrOfInt& Schema_Comm_statique::get_recv_size() const
{
  assert(0);
  assert(status_ == EXCHANGED);
  return recv_size_;
}

/*! @brief Verifies that send/recv_pe_list satisfy the property "you listen when I speak".
 *
 */
void Schema_Comm::check_send_recv_pe_list() const
{
  assert(status_ == RESET);
  // Check that processor indices are in the group
  int fail1 = 0;
  const int np = Process::nproc();
  const int n1 = send_pe_list_.size_array();
  int i;
  for (i = 0; i < n1; i++)
    if (send_pe_list_[i] < 0 || send_pe_list_[i] >= np)
      fail1 = 1;
  const int n2 = recv_pe_list_.size_array();
  for (i = 0; i < n2; i++)
    if (recv_pe_list_[i] < 0 || recv_pe_list_[i] >= np)
      fail1 = 1;
  int fail2 = 0;
  ArrOfInt recv_list;
  if (!fail1)
    {
      reverse_send_recv_pe_list(send_pe_list_, recv_list);
      // The array recv_pe_list_ is not necessarily sorted whereas recv_list is always sorted
      // Sort before comparing
      ArrOfInt copie(recv_pe_list_);
      copie.ordonne_array();
      fail2 = !(recv_list == copie);
    }
  if (Process::mp_sum(fail1+fail2))
    {
      if (Process::je_suis_maitre())
        Cerr << "Error in Schema_Comm::check_send_recv_pe_list(), see .log files" << finl;
      Process::Journal() << "Error in Schema_Comm::check_send_recv_pe_list() :\n"
                         << "send_list:\n" << send_pe_list_
                         << "recv_list:\n" << recv_pe_list_;
      if (fail1)
        Process::Journal() << "processor ranks not in current group: current group size = " << np;
      else if (fail2)
        Process::Journal() << "recv_list should be this one:\n" << recv_list << finl;
      else
        Process::Journal() << "OK on this processor" << finl;
      Process::barrier();
      Process::exit();
    }
}
