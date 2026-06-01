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

#include <communications.h>
#include <Comm_Group.h>
#include <TRUST_Ref.h>
#include <PE_Groups.h>

Implemente_base_sans_constructeur_ni_destructeur(Comm_Group,"Comm_Group",Objet_U);

int Comm_Group::check_enabled_ = 0;
int Comm_Group::static_group_number_ = 0;

Sortie& Comm_Group::printOn(Sortie& os) const
{
  Process::exit();
  return os;
}

Entree& Comm_Group::readOn(Entree& is)
{
  Process::exit();
  return is;
}

Comm_Group::Comm_Group()
{
  static const int group_increment = 32;

  // Communication tags are different
  // for each group as long as there are no more than 32 groups.
  // Thus, each communication of each group will have a different tag.
  group_number_ = static_group_number_;
  group_tag_increment_ = group_increment;
  static_group_number_++;
  group_communication_tag_ = group_number_ % group_increment;
}

/*! @brief Copy constructor is forbidden!
 *
 */
Comm_Group::Comm_Group(const Comm_Group& a): Objet_U(a)
{
  Cerr << "Comm_Group::Comm_Group(const Comm_Group &) error" << finl;
  Process::exit();
}

/*! @brief Assignment is forbidden!
 *
 */
const Comm_Group& Comm_Group::operator=(const Comm_Group&)
{
  Process::exit();
  return *this;
}

/*! @brief Destructor (nothing to do for now).
 *
 */
Comm_Group::~Comm_Group()
{
}

/*! @brief This function must be called simultaneously by all PEs of the current_group with the same parameters.
 *
 *   The processors in pe_list are the ranks within current_group() of the processors of the new group.
 *   The master of the group is the first in the list. The rank of the current processor, if it is in
 *   the group, is determined by its position in the list. There must be no duplicates.
 *   This function is called by the init_group methods of derived classes.
 *
 */
void Comm_Group::init_group(const ArrOfInt& pe_list)
{
  // All processors in the current group must arrive here
  const Comm_Group& current = PE_Groups::current_group();
  current.barrier(0);

  nproc_ = pe_list.size_array();
  if (nproc_ == 0)
    {
      Cerr << "Comm_Group::set_group_properties : empty process list" << finl;
      Process::exit();
    }
  rank_ = -1;
  const Comm_Group& groupe_trio = PE_Groups::groupe_TRUST();
  local_ranks_.resize_array(groupe_trio.nproc());
  local_ranks_ = -1;
  world_ranks_.resize_array(nproc_);

  for (int i = 0; i < nproc_; i++)
    {
      // rank of the pe in current_group()
      const int pe = pe_list[i], me = Process::me();
      const bool in_group = std::find(pe_list.begin(), pe_list.end(), me) != pe_list.end();

      if (check_enabled() && nproc_ > 1 && in_group)
        {
          if(Process::me() == pe_list[0])
            {
              for(int j=1; j < nproc_; j++)
                envoyer(pe_list[i], Process::me(), pe_list[j], 1);
            }
          else
            {
              int rcv;
              recevoir(rcv, pe_list[0], Process::me(), 1);
              if (rcv != pe_list[i])
                {
                  Cerr << "Comm_Group::init_group : processes have different pe_lists" << finl;
                  Process::exit();
                }
            }
        }
      if (pe < 0 || pe >= current.nproc_)
        {
          Cerr << "Comm_Group::set_group_properties : process "
               << pe << " is not in current_group()" << finl;
          Process::exit();
        }
      // rank of the pe in groupe_TRUST
      const int world_rank = current.world_ranks_[pe];

      if (local_ranks_[world_rank] >= 0)
        {
          Cerr << "Comm_Group::set_group_properties : duplicate pe in pe_list" << finl;
          Process::exit();
        }

      world_ranks_[i] = world_rank;
      local_ranks_[world_rank] = i;
      if (pe == current.rank_)
        rank_ = i;
    }
}

/*! @brief Initializes groupe_TRUST().
 *
 * This method is called by init_group_trio() of derived classes.
 *
 */
void Comm_Group::init_group_trio(int nproc_tot, int arank)
{
  assert(arank >= 0 && arank < nproc_tot);
  rank_ = arank;
  nproc_ = nproc_tot;
  local_ranks_.resize_array(nproc_);
  world_ranks_.resize_array(nproc_);
  for (int i = 0; i < nproc_; i++)
    {
      local_ranks_[i] = i;
      world_ranks_[i] = i;
    }
}

/*! @brief Initialize all the information relative to world sizes and ranks for node communicator
 *
 * This method is called by derived classes when initializing communicator on node
 *
 */
void Comm_Group::init_group_node(int nproc, int loc_rank, int glob_rank)
{
  rank_ = loc_rank;
  nproc_ = nproc;

  const Comm_Group& groupe_trio = PE_Groups::groupe_TRUST();
  local_ranks_.resize_array(groupe_trio.nproc());
  local_ranks_ = -1;
  groupe_trio.all_gather(&loc_rank, local_ranks_.addr(), (int)sizeof(int));

  world_ranks_.resize_array(nproc_);
  world_ranks_ = glob_rank;
  if(nproc_ > 1) // the master of my node is of parallel type but only contains one proc, so can't perform MPI communications with it
    all_gather(&glob_rank, world_ranks_.data(), (int)sizeof(int));
}

void Comm_Group::set_check_enabled(int flag)
{
  if (flag)
    check_enabled_ = 1;
  else
    check_enabled_ = 0;
}
