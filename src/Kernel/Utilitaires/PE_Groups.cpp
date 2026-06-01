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

#include <Comm_Group.h>
#include <PE_Groups.h>
#include <TRUST_Ref.h>

// The following three variables save the group stack
// (see Comm_Group::enter_group(), Comm_Group::current_group(), Comm_Group::exit_group())
// The top of the stack is always groupe_TRUST(), provided to initialize().
// groups[0] points to groupe_trio.
static OBS_PTR(Comm_Group) groups[100];
static int ngroups = 0;
static int max_ngroups = 100;
const Comm_Group * PE_Groups::current_group_ = 0;
// node group is an isolated variable from all the other groups as it is only used for IO purposes
// and might be used throughout the code together with other groups
static OBS_PTR(Comm_Group) node_group;
static OBS_PTR(Comm_Group) node_master;
// For the user that defines his own group ! not done in the main, see the My_Comm_Group class !
static OBS_PTR(Comm_Group) user_defined_group;

int PE_Groups::check_current_group()
{
  assert(ngroups > 0);
  assert(current_group_ == &(groups[ngroups-1].valeur()));
  return 1;
}

/*! @brief Creates a new processor group (can be called anywhere in the code).
 *
 *   Must be called simultaneously on all processors of current_group() with the same
 *   liste_pe array. liste_pe is the list of ranks within the current group of processors
 *   to include in the new group. The first in the list will be the group master. The list
 *   must not contain duplicates and must include at least one processor.
 *   This method types and initializes the group object.
 *   enter_group() and exit_group() must then be called (as many times as desired).
 *
 * @param liste_pe List of PE ranks within current_group().
 * @param group The group object to initialize.
 * @param force_Comm_Group_NoParallel If non-zero, force a non-parallel group when possible.
 */
void PE_Groups::create_group(const ArrOfInt& liste_pe, OWN_PTR(Comm_Group) & group, int force_Comm_Group_NoParallel)
{
  if (liste_pe.size_array()==1 && force_Comm_Group_NoParallel)
    {
      // Create a non-parallel group if possible and if required
      group.typer("Comm_Group_NoParallel");
    }
  else
    {
      // Create a group of the same type as groupe_TRUST
      group.typer(groups[0]->que_suis_je());
    }
  group->init_group(liste_pe);
}

/*! @brief Initializes a new processor group that is already instantiated (can be called anywhere in the code).
 *
 *   Must be called simultaneously on all processors of current_group() with the same
 *   liste_pe array. liste_pe is the list of ranks within the current group of processors
 *   to include in the new group. The first in the list will be the group master. The list
 *   must not contain duplicates and must include at least one processor.
 *   enter_group() and exit_group() must then be called (as many times as desired).
 *
 * @param liste_pe List of PE ranks within current_group().
 * @param group The already-instantiated group object to initialize.
 */
void PE_Groups::init_group(const ArrOfInt& liste_pe, OWN_PTR(Comm_Group) & group)
{
  assert(group);
  group->init_group(liste_pe);
}

/*! @brief If the local processor belongs to the group, the current group for this processor becomes "group" and returns 1; otherwise returns 0.
 *
 *   A reference to the current group is saved and will be restored when exit_group() is called.
 *   This method must be called simultaneously on all processors of "group".
 *   Each enter_group() must be matched by a corresponding exit_group().
 *   However, it is allowed to enter the same group multiple times in a row, or to enter a
 *   larger group than the current one.
 *   Example: group1 and group2 form a partition of groupe_TRUST().
 *
 *    int sync_point(int x)
 *    {
 *      PE_Groups::enter_group(groupe_TRUST());
 *      int i = mp_sum(x);
 *      PE_Groups::exit_group();
 *      return i;
 *    }
 *    if (PE_Groups::enter_group(group1)) {
 *      s1 = mp_sum(x); // Sum over group 1
 *      // Sync point with the other group:
 *      s_all = sync_point(x);
 *    } else if (PE_Groups::enter_group(group2)) {
 *      s2 = mp_sum(x);
 *      // Sync point with the other group:
 *      s_all = sync_point(x);
 *    } else {
 *      Cerr << "Error: processor " << me() << " is not within a subgroup.";
 *      exit();
 *    }
 *    PE_Groups::exit_group(); // Exit the subgroup
 *
 * @param group The group to enter.
 * @return 1 if the local processor is in the group, 0 otherwise.
 */
int PE_Groups::enter_group(const Comm_Group& group)
{
  assert(&group != &current_group());
  if (ngroups >= max_ngroups-1)
    {
      Cerr << "Comm_Group::enter_group : fatal, too many groups" << finl;
      Process::exit();
    }
  int my_rank_in_group = rank_translate(current_group().rank(), current_group(), group);
  if (my_rank_in_group >= 0)
    {
      // Save the pointer to the current group and switch current_group():
      groups[ngroups] = group;
      current_group_ = &group;
      ngroups++;
      // Note: current_group() has changed!

      // Verify that all processors of the new group are present:
      if (Comm_Group::check_enabled())
        current_group().barrier(0);

      return 1;
    }
  else
    {
      return 0;
    }
}

/*! @brief Returns to the group that was active before the last successful enter_group() call (which returned 1).
 *
 *   This method must be called simultaneously on all processors of the current current_group() just before exit_group().
 *
 */
void PE_Groups::exit_group()
{
  if (Comm_Group::check_enabled())
    current_group().barrier(0);

  if (ngroups <= 1)
    {
      Cerr << "Comm_Group::exit_group() error : trying to exit from TRUST main group." << finl;
      Process::exit();
    }
  ngroups--;
  current_group_ = &(groups[ngroups-1].valeur());
}

/*! @brief Computes the rank in the current group of the processor with rank "rank" in "group".
 *
 * Requires 0 <= rank < group.nproc().
 *   Returns -1 if the processor is not in the current group.
 *
 * @param rank Rank in the source group.
 * @param group The source group.
 * @param dest_group The destination group.
 * @return Rank in the destination group, or -1 if not present.
 */
int PE_Groups::rank_translate(int rank, const Comm_Group& group,
                              const Comm_Group& dest_group)
{
  const int world_rank = group.world_ranks_[rank];
  const int local_rank = dest_group.local_ranks_[world_rank];
  return local_rank;
}

/*! @brief Returns a reference to the group containing all TRUST processors.
 *
 * @return Reference to the global TRUST Comm_Group.
 */
const Comm_Group& PE_Groups::groupe_TRUST()
{
  assert(ngroups > 0); // Initialized ?
  return groups[0].valeur();
}

/*! @brief Returns a reference to the node-level communicator group.
 *
 * @return Reference to the node Comm_Group.
 */
const Comm_Group& PE_Groups::get_node_group()
{
  assert(node_group);
  return node_group.valeur();
}

/*! @brief Returns the group containing the master of my node.
 *
 * @return Reference to the node-master Comm_Group.
 */
const Comm_Group& PE_Groups::get_node_master()
{
  assert(node_master);
  return node_master.valeur();
}

/*! @brief Returns a reference to the user-defined group.
 *
 * @return Reference to the user-defined Comm_Group.
 */
const Comm_Group& PE_Groups::get_user_defined_group()
{
  assert(user_defined_group);
  return user_defined_group.valeur();
}

/*! @brief Method to call at the beginning of execution (MAIN.cpp). Initializes current_group() with groupe_trio_u.
 *
 * @param groupe_trio_u The main TRUST communicator group to initialize with.
 */
void PE_Groups::initialize(const Comm_Group& groupe_trio_u)
{
  assert(ngroups == 0);
  ngroups = 1;
  groups[0] = groupe_trio_u;
  current_group_ = &groupe_trio_u;
}

/*! @brief Method to call after the initialization of trio_u_world and TRUST's statistical counters.
 *
 * @param ngrp The node-level communicator group.
 */
void PE_Groups::initialize_node(const Comm_Group& ngrp)
{
  assert(!node_group);
  node_group = ngrp;
}

void PE_Groups::initialize_user_defined_group(const Comm_Group& ngrp)
{
  assert(!user_defined_group);
  user_defined_group = ngrp;
}

bool PE_Groups::has_user_defined_group()
{
  return bool(user_defined_group);
}

/*! @brief Method to call after the initialization of trio_u_world, node_group, and TRUST's statistical counters.
 *
 * @param ngrp The node-master communicator group.
 */
void PE_Groups::initialize_node_master(const Comm_Group& ngrp)
{
  assert(!node_master);
  node_master = ngrp;
}

/*! @brief Method to call at the end of execution, once back in groupe_TRUST() and just before destroying the main Comm_Group.
 *
 */
void PE_Groups::finalize()
{
  assert(ngroups == 1);
  groups[0].reset();
  ngroups = 0;
  current_group_ = 0;
  node_group.reset();
  node_master.reset();
  user_defined_group.reset();
}

const int& PE_Groups::get_nb_groups()
{
  return ngroups ;
}
