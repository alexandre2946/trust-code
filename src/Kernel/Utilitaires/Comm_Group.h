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

#ifndef Comm_Group_included
#define Comm_Group_included

#include <TRUST_Deriv.h>
#include <TRUSTArray.h>
#include <assert.h>

/*! @brief : This class describes a group of processors on which
 *
 *   a portion of code executes simultaneously. It provides all the methods
 *   for exchanging data between the processors of the group (mpsum, send, recv, ...),
 *   and for synchronizing processors (barrier).
 *   It is specialized according to the network layer (MPI, PVM, ...).
 *   Note: these methods are reserved for low-level operations (TRUST kernel).
 *   In normal code, use the high-level communication class methods:
 *   (envoyer(), envoyer_broadcast(), class Schema_Comm, class Process, etc.)
 *   To create a new group and use it, see class PE_Groups.
 *   For the initialization procedure, see PE_Groups::Initialize().
 *
 */
class Comm_Group : public Objet_U
{
  Declare_base_sans_constructeur_ni_destructeur(Comm_Group);
public:
  Comm_Group();
  ~Comm_Group() override;
  virtual void   abort() const = 0;

  // COLL_SUM: sum over all procs
  // COLL_MIN: minimum
  // COLL_MAX: max
  // COLL_PARTIAL_SUM computes the partial sum of values over processors with rank
  // strictly less than me() (the result is always 0 on processor 0).
  enum Collective_Op { COLL_SUM, COLL_MIN, COLL_MAX, COLL_PARTIAL_SUM };
  virtual void mp_collective_op(const double *x, double *resu, int n, Collective_Op op) const = 0;
  virtual void mp_collective_op(const double *x, double *resu, const Collective_Op *op, int n) const = 0;
  virtual void mp_collective_op(const float *x, float *resu, int n, Collective_Op op) const = 0;
  virtual void mp_collective_op(const float *x, float *resu, const Collective_Op *op, int n) const = 0;
  virtual void mp_collective_op(const int *x, int *resu, int n, Collective_Op op) const = 0;
  virtual void mp_collective_op(const int *x, int *resu, const Collective_Op *op, int n) const = 0;
#if INT_is_64_ == 2
  virtual void mp_collective_op(const trustIdType *x, trustIdType *resu, int n, Collective_Op op) const = 0;
  virtual void mp_collective_op(const trustIdType *x, trustIdType *resu, const Collective_Op *op, int n) const = 0;
#endif

  virtual void barrier(int tag) const = 0;

  // Computes a new communication tag that allows identifying exchanges
  // uniquely across all groups.
  inline int get_new_tag() const;

  inline int rank() const;
  inline int nproc() const;

  inline int get_node_id() const;
  inline int get_number_of_nodes() const;


  // Do we want to perform additional checks on communications?
  // These checks imply extra communications, which modifies the program flow.
  // This is therefore a separate mechanism from "assert".
  inline static int check_enabled();

  enum TypeHint { CHAR, INT, DOUBLE, FLOAT };
  // Starts the exchange of buffers.
  // send_list / recv_list = list of PEs (ranks within the current group)
  // send_size / recv_size = size of messages in bytes
  // send_buffers / recv_buffers = address of buffers
  // Reception buffers must have sufficient size.
  // Note about const:
  //  send_buffers is completely const, nothing may be modified
  //  recv_buffers is const, recv_buffers[i] is const but *(recv_buffers[i])
  //               is not const because received data is stored there.
  virtual void send_recv_start(const ArrOfInt& send_list,
                               const ArrOfInt& send_size,
                               const char * const * const send_buffers,
                               const ArrOfInt& recv_list,
                               const ArrOfInt& recv_size,
                               char * const * const recv_buffers,
                               TypeHint typehint = CHAR) const = 0;
  // Waits until communications started by send_recv are finished.
  virtual void send_recv_finish() const = 0;

  // Blocking send/receive methods: each send must be matched
  // simultaneously by a recv on the destination processor.
  virtual void send(int pe, const void *buffer, int size, int tag) const = 0; // Blocking send
  virtual void recv(int pe, void *buffer, int size, int tag) const = 0; // Blocking receive

  // Broadcast methods: must be called on all processors simultaneously
  virtual void broadcast(void *buffer, int size, int pe_source) const = 0;

  // All-to-all methods
  virtual void all_to_all(const void *src_buffer, void *dest_buffer, int data_size) const = 0;
  virtual void all_gather(const void *src_buffer, void *dest_buffer, int data_size) const = 0;
  virtual void gather(const void *src_buffer, void *dest_buffer, int data_size, int root) const = 0;
  virtual void all_gatherv(const void *src_buffer, void *dest_buffer, int send_size, const int* recv_size, const int* displs) const = 0;

  static void set_check_enabled(int flag);
protected:
  Comm_Group(const Comm_Group&);  // forbidden!
  const Comm_Group& operator=(const Comm_Group&);   // forbidden!
  virtual void       init_group(const ArrOfInt& pe_list);
  void               init_group_node(int nproc, int loc_rank, int glob_rank);
  void               init_group_trio(int nproc, int rank);
  friend class PE_Groups;

  // ToDo gather that in a derived Comm_Group_MPI_Node class ?
  // id of my node among all the other nodes
  int node_id_  = -1;
  // total number of nodes
  int nb_nodes_ = -1;

private:
  static int check_enabled_;
  static int static_group_number_;

  // Rank of the local processor in the group, -1 if it is not in the group
  int rank_ = -1;
  // Number of processors in the group
  int nproc_ = -1;
  // For each PE in the full computation (array size = groupe_TRUST().nproc())
  //  index within the group if the PE is in it,
  //  -1 if the PE is not in the group
  ArrOfInt  local_ranks_;
  // List of processors in the group (indices of processors in groupe_TRUST())
  // (array size = nproc_)
  ArrOfInt  world_ranks_;

  // My group number (equal to static_group_number_ at the time the group was created).
  int group_number_ = -1;
  // The group_communication_tag_ is incremented by this amount at each
  // operation. It is a prime number, which allows different tags
  // for each group for a long time (until the tag number exceeds MAXINT...).
  int group_tag_increment_ = -1;
  // The tag is incremented at each operation, allowing verification that processes are properly synchronized.
  mutable int group_communication_tag_ = -1;
};

inline int Comm_Group::check_enabled()
{
  return check_enabled_;
}

/*! @brief Returns a new communication tag for the group.
 *
 * Side effect: increments the group_communication_tag_ member.
 *
 */
inline int Comm_Group::get_new_tag() const
{
  // B.M. This feature is ultimately of little practical use
  // and when the counter exceeds a limit, MPI crashes. Disabling:
  //group_communication_tag_ += group_tag_increment_;
  return group_communication_tag_;
}

/*! @brief Returns the rank of the local processor in the group *this.
 *
 * or -1 if this processor is not in the group.
 *
 */
inline int Comm_Group::rank() const
{
  return rank_;
}

/*! @brief Returns the number of processors in the group *this
 *
 */
inline int Comm_Group::nproc() const
{
  assert(nproc_ >= 0);
  return nproc_;
}

/*! @brief Retrieve ID of my numa node
 *
 */
inline int Comm_Group::get_node_id() const
{
  return node_id_;
}

inline int Comm_Group::get_number_of_nodes() const
{
  return nb_nodes_;
}


#endif
