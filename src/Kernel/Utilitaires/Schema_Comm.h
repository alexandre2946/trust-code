/****************************************************************************
* Copyright (c) 2024, CEA
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

#ifndef Schema_Comm_included
#define Schema_Comm_included

#include <TRUSTArray.h>
#include <TRUST_Ref.h>

class Comm_Group;

// These objects store a communication graph between
// processors: each proc. has a list of processors to send to
// and a list of processors to receive from.
//
// The group is the set of processors on which it is GUARANTEED
// that the following methods will be called SIMULTANEOUSLY on all
// processors in the group:
// - begin_comm()
// - set_send_recv_pe_list(...)
// - echange_taille_et_messages()
// - end_comm()
//
// In particular: it is forbidden to use these methods inside
// a loop that is not executed the same number of times by all
// processors in the group, nor inside an "if() { }" whose
// execution is not identical on all processors in the group.
//
// It is forbidden to start a new communication while another is
// in progress (a communication always ends with "end_comm()")
// In particular, it is forbidden to use "envoyer", "recevoir" and
// shared disk files between the first call to schema.send_buffer()
// and "terminer()". Be very careful about all methods used between
// these two calls!
//
// On a given PE, an exchange sequence must be built as follows:
// schema.begin_comm()
// schema.send_buffer(pe1) << data_to_send;
// schema.send_buffer(pe2) << data_to_send;
// ...
// schema.echanger_taille_et_messages();
// schema.recv_buffer(pe2) >> data_to_recv;
// schema.recv_buffer(pe3) >> data_to_recv;
// ...
// schema.end_comm();
//
// Communication is not necessarily symmetric: a processor can send
// a message to one processor and receive from another.
// However, the user guarantees that the processor lists provided in
//  send_pe_list and recv_pe_list satisfy the principle "you listen when I speak!"
//  (i.e., processor A belongs to send_pe_list on processor B
//   if and only if processor B belongs to recv_pe_list on processor A).

// Modif BM 20/06/2013: adding set_all_to_allv_flag. If the flag is set, the communication
//  scheme uses MPI_alltoallv instead of ISend IRecv. To try to solve
//  problems encountered on supermuc in the file read/write routine
//  (schema where everyone writes to processor 0 => error allocating MPI_requests).

class OutputCommBuffer;
class InOutCommBuffers;
class InputCommBuffer;
class Comm_Group;
class Entree;
class Sortie;

class Schema_Comm
{
public:
  Schema_Comm();
  Schema_Comm(const Schema_Comm&);
  ~Schema_Comm();

  void set_group(const Comm_Group& group);  // Obsolete
  const Comm_Group& get_group() const;

  const Schema_Comm& operator= (const Schema_Comm&);
  void set_send_recv_pe_list(const ArrOfInt& send_pe_list, const ArrOfInt& recv_pe_list, const int me_to_me = 0);

  void begin_comm() const;                 // Status transitions to WRITING
  // Allowed when status_ == WRITING:
  Sortie& send_buffer(int num_PE) const;
  void echange_taille_et_messages() const; // Status transitions to EXCHANGED
  void echange_messages(const ArrOfInt& recv_size) const;  // Status transitions to EXCHANGED
  // Allowed when status_ == EXCHANGED:
  Entree& recv_buffer(int num_PE) const;
  void end_comm() const;                   // Status transitions to RESET
  // Accessors:
  const ArrOfInt& get_send_pe_list() const;
  const ArrOfInt& get_recv_pe_list() const;

  void set_all_to_allv_flag(int x) { use_all_to_allv_ = x; }

protected:
  void echange_taille(const ArrOfInt& send_size, ArrOfInt& recv_size) const;
  void echange_messages(const ArrOfInt& send_size, const ArrOfInt& recv_size) const;
  // Status transitions to EXCHANGED
  void check_send_recv_pe_list() const;

  // A single status for all exchanges: concurrent access to the class is not supported
  // because we want to limit the number of buffers and the number of "outstanding requests".
  // Therefore: it is forbidden to start a new communication if the buffers
  // are currently in use.
  enum Static_Status { UNINITIALIZED, RESET, WRITING, EXCHANGED };
  static Static_Status status_;
  static OutputCommBuffer& obuffer(int pe);
  static InputCommBuffer&   ebuffer(int pe);

  ArrOfInt send_pe_list_; // List of processors to send to
  ArrOfInt recv_pe_list_; // List of processors to receive from
  int   me_to_me_;     // Flag: is sending messages to oneself allowed?
  OBS_PTR(Comm_Group) ref_group_;// Group of processors that will communicate

  int use_all_to_allv_; // Flag, which type of communication should be used?
private:
  // Pointers are stored in a specific class (destructor of static members
  // called automatically at end of execution to free memory).
  static InOutCommBuffers buffers_;
  static int n_buffers_;
};

class Schema_Comm_statique : public Schema_Comm
{
public:
  //void echange_taille();
  //void echange_messages() const;           // Statut passe a EXCHANGED
  const ArrOfInt& get_send_size() const;
  const ArrOfInt& get_recv_size() const;
protected:
  ArrOfInt send_size_;    // Size of messages to send in bytes
  ArrOfInt recv_size_;    // Size of messages to receive in bytes
};

#endif
