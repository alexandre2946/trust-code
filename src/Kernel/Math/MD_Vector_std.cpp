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

#include <MD_Vector_std.h>
#include <Param.h>
#include <Array_tools.h>

Implemente_instanciable(MD_Vector_std,"MD_Vector_std",MD_Vector_mono);


/*! @brief Stupid ctor.
 *
 * Everything set to n.
 */
MD_Vector_std::MD_Vector_std(int n)
{
  nb_items_tot_ = n;
  nb_items_reels_ = n;
  nb_items_seq_tot_ = n;
  nb_items_seq_local_ = n;
  blocs_items_to_sum_.resize_array(2, RESIZE_OPTIONS::NOCOPY_NOINIT);
  blocs_items_to_sum_[0] = 0;
  blocs_items_to_sum_[1] = n;
  blocs_items_to_compute_.resize_array(2, RESIZE_OPTIONS::NOCOPY_NOINIT);
  blocs_items_to_compute_[0] = 0;
  blocs_items_to_compute_[1] = n;
}

/*! @brief Constructor.
 *
 * If nb_items_reels >= 0, items_to_compute will contain a block with the real items,
 *    otherwise, items_to_compute will cover all items (up to nb_items_tot).
 *
 * @param (nb_items_tot) value assigned to nb_items_tot_ (must be >= 0)
 * @param (nb_items_reels) value assigned to nb_items_reels_ (must be >= -1; -1 means there are no "real" items grouped at the start of the array)
 * @param (pe_voisins) list of neighbouring processors, sorted in ascending order (ArrsOfInt arrays must have the same size). This list is not taken as-is: processors for which all three item lists are empty are removed.
 * @param (items_to_send) for each pe_voisin, list of items to send (shared and virtual). count_items_to_send_to_items_ is computed via communication from items_to_recv.
 * @param (items_to_recv) list of individual items to receive from the different procs (often shared items but not necessarily, e.g. front-tracking)
 * @param (blocs_to_recv) list of blocks of items to receive (see MD_Vector_std::blocs_to_recv_) (often the virtual items...)
 */
MD_Vector_std::MD_Vector_std(int nb_items_tot, int nb_items_reels, const ArrOfInt& pe_voisins,
                             const ArrsOfInt& items_to_send, const ArrsOfInt& items_to_recv, const ArrsOfInt& blocs_to_recv)
{
  assert(nb_items_tot >= 0);
  assert(nb_items_reels >= -1 && nb_items_reels <= nb_items_tot);
  nb_items_tot_ = nb_items_tot;
  nb_items_reels_ = nb_items_reels;
  const int nb_voisins = pe_voisins.size_array();
  assert(items_to_send.size() == nb_voisins);
  assert(items_to_recv.size() == nb_voisins);
  assert(blocs_to_recv.size() == nb_voisins);

  // selection: list of indices of processors to keep in pe_voisins
  //  (procs with whom data is actually exchanged)
  ArrOfInt tmp, selection;


  int i;
  for (i = 0; i < nb_voisins; i++)
    if (items_to_send[i].size_array() > 0 || items_to_recv[i].size_array() > 0 || blocs_to_recv[i].size_array() > 0)
      {
        tmp.append_array(pe_voisins[i]);
        selection.append_array(i);
      }
  const int nb_voisins2 = tmp.size_array();
  // ** pe_voisins_ **
  pe_voisins_ = tmp;

  {
    ArrsOfInt tmpbis(nb_voisins2);

    for (i = 0; i < nb_voisins2; i++)
      tmpbis[i] = items_to_send[selection[i]];
    items_to_send_.set(tmpbis);

    for (i = 0; i < nb_voisins2; i++)
      tmpbis[i] = items_to_recv[selection[i]];
    items_to_recv_.set(tmpbis);

    for (i = 0; i < nb_voisins2; i++)
      tmpbis[i] = blocs_to_recv[selection[i]];
    blocs_to_recv_.set(tmpbis);
  }

  // Compute nb_items_to_items_
  {
    nb_items_to_items_.resize_array(nb_voisins2, RESIZE_OPTIONS::NOCOPY_NOINIT);
    tmp.resize_array(nproc(), RESIZE_OPTIONS::NOCOPY_NOINIT);
    tmp = 0;
    for (i = 0; i < nb_voisins2; i++)
      tmp[pe_voisins_[i]] = items_to_recv[selection[i]].size_array();
    ;
    envoyer_all_to_all(tmp, tmp);
    for (i = 0; i < nb_voisins2; i++)
      {
        int n = tmp[pe_voisins_[i]];
        assert(n <= items_to_send_.get_list_size(i));
        nb_items_to_items_[i] = n;
      }
  }
  // Compute blocs_items_count_
  {
    blocs_items_count_.resize_array(nb_voisins2, RESIZE_OPTIONS::NOCOPY_NOINIT);
    for (i = 0; i < nb_voisins2; i++)
      {
        const int nblocs = blocs_to_recv_.get_list_size(i) / 2;
        int count = 0;
        for (int j = 0; j < nblocs; j++)
          {
            int deb = blocs_to_recv_(i, j * 2);
            int fin = blocs_to_recv_(i, j * 2 + 1);
            assert(deb >= 0 && fin > deb && fin <= nb_items_tot_);
            count += fin - deb;
          }
        blocs_items_count_[i] = count;
      }
  }
  // Compute the blocks of sequential items (items not received from another proc)
  tmp.resize_array(nb_items_tot, RESIZE_OPTIONS::NOCOPY_NOINIT);
  {
    // Mark received items as zero
    tmp = 1;
    const ArrOfInt& items = items_to_recv_.get_data();
    int n = items.size_array();
    for (i = 0; i < n; i++)
      {
        assert(tmp[items[i]] == 1); // otherwise the same item is received twice
        tmp[items[i]] = 0;
      }
    const ArrOfInt& items2 = blocs_to_recv_.get_data();
    assert(items2.size_array() % 2 == 0);
    n = items2.size_array() / 2;
    for (i = 0; i < n; i++)
      {
        const int start = items2[i * 2];
        const int end = items2[i * 2 + 1];
        for (int j = start; j < end; j++)
          {
            assert(tmp[j] == 1); // otherwise, the same item is received twice
            tmp[j] = 0;
          }
      }

    // Build a list of sequential item blocks (those that remained 1)
    ArrOfInt blocs;

    int nb_seq = 0;
    for (i = 0; i < nb_items_tot; i++)
      {
        if (tmp[i])
          {
            if (nb_items_reels >= 0 && i >= nb_items_reels)
              {
                // On a un item non recu parmi les items virtuels !
                Cerr << "Error in MD_Vector_std::MD_Vector_std(...) [pe " << me() << "]\n "
                     "an item i(" << i << ") > nb_items_reels(" << nb_items_reels << ") is not received from any other processor" << finl;
                Process::exit();
              }
            append_item_to_blocs(blocs, i);
            nb_seq++;
          }
      }
    blocs_items_to_sum_ = blocs;
    nb_items_seq_local_ = nb_seq;
    nb_items_seq_tot_ = mp_sum(nb_seq);
  }

  // Bloc items to compute:
  // Default: a single block
  if (nb_items_reels_ >= 0)
    {
      // Array operators compute all real items
      blocs_items_to_compute_.resize(2, RESIZE_OPTIONS::NOCOPY_NOINIT);
      blocs_items_to_compute_[0] = 0;
      blocs_items_to_compute_[1] = nb_items_reels_;
    }
  else
    {
      // Array operators compute everything
      blocs_items_to_compute_.resize(2, RESIZE_OPTIONS::NOCOPY_NOINIT);
      blocs_items_to_compute_[0] = 0;
      blocs_items_to_compute_[1] = nb_items_tot_;
    }
}

/*! @brief method used to dump/restore a descriptor in a file Each process writes a different descriptor.
 *
 * See MD_Vector_tools::restore_vector_with_md()
 *
 */
Entree& MD_Vector_std::readOn(Entree& is)
{
  MD_Vector_mono::readOn(is);
  Param p(que_suis_je());
  ArrOfInt items_send_index, items_send_data, items_index, items_data, blocs_index, blocs_data;
  p.ajouter("pe_voisins", &pe_voisins_);
  p.ajouter("items_to_send_index", &items_send_index);
  p.ajouter("items_to_send_data", &items_send_data);
  p.ajouter("nb_items_to_items", &nb_items_to_items_);
  p.ajouter("items_to_recv_index", &items_index);
  p.ajouter("items_to_recv_data", &items_data);
  p.ajouter("blocs_to_recv_index", &blocs_index);
  p.ajouter("blocs_to_recv_data", &blocs_data);
  p.ajouter("blocs_items_count", &blocs_items_count_);
  p.lire_avec_accolades(is);
  if (items_send_data.size_array())
    items_to_send_.set_index_data(items_send_index, items_send_data);
  if (items_data.size_array())
    items_to_recv_.set_index_data(items_index, items_data);
  if (blocs_data.size_array())
    blocs_to_recv_.set_index_data(blocs_index, blocs_data);
  return is;
}

/*! @brief method used to dump/restore a descriptor in a file Each process writes a different descriptor.
 *
 * See MD_Vector_tools::dump_vector_with_md()
 *
 */
Sortie& MD_Vector_std::printOn(Sortie& os) const
{
  MD_Vector_mono::printOn(os);
  os << "{" << finl;
  os << "pe_voisins" << tspace << pe_voisins_ << finl;
  os << "items_to_send_index" << tspace << items_to_send_.get_index() << finl;
  os << "items_to_send_data" << tspace << items_to_send_.get_data() << finl;
  os << "nb_items_to_items" << tspace << nb_items_to_items_ << finl;
  os << "items_to_recv_index" << tspace << items_to_recv_.get_index();
  os << "items_to_recv_data" << tspace << items_to_recv_.get_data();
  os << "blocs_to_recv_index" << tspace << blocs_to_recv_.get_index() << finl;
  os << "blocs_to_recv_data" << tspace << blocs_to_recv_.get_data();
  os << "blocs_items_count" << tspace << blocs_items_count_ << finl;
  os << "}" << finl;
  return os;
}

namespace   // Anonymous namespace
{

/*! @brief returns i such that a[i]==x, or -1 if x is not found
 */
int find_in_array(const ArrOfInt& a, int x)
{
  int n;
  for (n = a.size_array() - 1; n >= 0; n--)
    if (a[n] == x)
      break;
  return n;
}

/*! @brief Returns the number of extracted items (not the number of blocks).
 */
int extract_blocs(const ArrOfInt src, const ArrOfInt& renum, ArrOfInt& dest)
{
  const int nblocs_src = src.size_array() / 2;

  dest.resize_array(0);
  int end_last_bloc = -1;
  int count = 0;
  for (int i = 0; i < nblocs_src; i++)
    {
      const int deb = src[i*2];
      const int fin = src[i*2+1];
      for (int j = deb; j < fin; j++)
        {
          const int rj = renum[j];
          if (rj >= 0)
            {
              count++;
              // Blocs construction assumes that renum is ordered
              assert(rj >= end_last_bloc);
              if (rj == end_last_bloc)
                {
                  // item is contiguous to last bloc in destination array, increase bloc size:
                  end_last_bloc++;
                  dest[dest.size_array()-1] = end_last_bloc;
                }
              else
                {
                  // non contiguous item, create new bloc
                  end_last_bloc = rj;
                  dest.append_array(end_last_bloc);
                  end_last_bloc++;
                  dest.append_array(end_last_bloc);
                }
            }
        }
    }
  return count;
}

} // end anonymous NS


/*! Append all items in "src" to "this", duplicating src items "multiplier" times
 *  and adding the "offset" to item indexes. Also merge pe_voisins_ lists.
 */
void MD_Vector_std::append_from_other_std(const MD_Vector_std& src, int offset, int multiplier)
{
  // List sizes for items_to_send_, items_to_recv_ and blocs_to_recv_ of new descriptor
  ArrOfInt x_sz, y_sz, z_sz;
  // Data of these lists:
  ArrOfInt x_data, y_data, z_data;
  ArrOfInt new_blocs_items_count, new_nb_items_to_items;
  ArrOfInt tmp;
  ArrOfInt pe_list(pe_voisins_);

  append_blocs(pe_list, src.pe_voisins_);
  array_trier_retirer_doublons(pe_list);
  const int np = pe_list.size_array();

  new_blocs_items_count.resize_array(np);
  new_nb_items_to_items.resize_array(np);

  for (int i = 0; i < np; i++)
    {
      // Number of items to recv, send and blocs to recv for this processor:
      int nx = 0, ny = 0, nz = 0;
      int blocs_count = 0;
      // Is processor pe_list[i] in initial "dest" descriptor ?
      int j;

      for (j = pe_voisins_.size_array()-1; j >= 0; j--)
        if (pe_voisins_[j] == pe_list[i])
          break;
      const int pe = pe_list[i];

      // Insert items to send and items to recv from dest
      j = ::find_in_array(pe_voisins_, pe);
      if (j >= 0)
        {
          const int count_items_single = nb_items_to_items_[j];
          nx += count_items_single;
          // take only single items to send
          tmp.resize_array(0);
          tmp.resize_array(count_items_single);
          for (int k = 0; k < count_items_single; k++)
            tmp[k] = items_to_send_(j, k);
          append_blocs(x_data, tmp);

          ny += items_to_recv_.get_list_size(j);
          items_to_recv_.copy_list_to_array(j, tmp);
          append_blocs(y_data, tmp);
        }

      // Insert items to send and items to recv from src, with multiplier
      j = find_in_array(src.pe_voisins_, pe);
      if (j >= 0)
        {
          const int count_items_single = src.nb_items_to_items_[j];
          nx += count_items_single * multiplier;
          // take only single items to send
          tmp.resize_array(0);
          tmp.resize_array(count_items_single);
          for (int k = 0; k < count_items_single; k++)
            tmp[k] = src.items_to_send_(j, k);
          append_items(x_data, tmp, offset, multiplier);

          ny += src.items_to_recv_.get_list_size(j) * multiplier;
          src.items_to_recv_.copy_list_to_array(j, tmp);
          append_items(y_data, tmp, offset, multiplier);
        }
      const int single_items_count = nx;
      // Insert blocs to send and blocs to recv from dest
      j = ::find_in_array(pe_voisins_, pe);
      if (j >= 0)
        {
          const int count_items_tot = items_to_send_.get_list_size(j);
          const int count_items_single = nb_items_to_items_[j];
          nx += count_items_tot - count_items_single;
          tmp.resize_array(0);
          tmp.resize_array(count_items_tot - count_items_single);
          for (int k = count_items_single; k < count_items_tot; k++)
            tmp[k - count_items_single] = items_to_send_(j, k);
          append_blocs(x_data, tmp);

          nz += blocs_to_recv_.get_list_size(j);
          blocs_to_recv_.copy_list_to_array(j, tmp);
          append_blocs(z_data, tmp);

          blocs_count += blocs_items_count_[j];
        }
      // Insert blocs to send and blocs to recv from src, with multiplier
      j = find_in_array(src.pe_voisins_, pe);
      if (j >= 0)
        {
          const int count_items_tot = src.items_to_send_.get_list_size(j);
          const int count_items_single = src.nb_items_to_items_[j];
          nx += (count_items_tot - count_items_single) * multiplier;
          tmp.resize_array(0);
          tmp.resize_array(count_items_tot - count_items_single);
          for (int k = count_items_single; k < count_items_tot; k++)
            tmp[k - count_items_single] = src.items_to_send_(j, k);
          append_items(x_data, tmp, offset, multiplier);

          nz += src.blocs_to_recv_.get_list_size(j);
          src.blocs_to_recv_.copy_list_to_array(j, tmp);
          append_blocs(z_data, tmp, offset, multiplier);

          blocs_count += src.blocs_items_count_[j] * multiplier;
        }

      x_sz.append_array(nx);
      y_sz.append_array(ny);
      z_sz.append_array(nz);
      new_blocs_items_count[i] = blocs_count;
      new_nb_items_to_items[i] = single_items_count;
    }

  items_to_send_.set_list_sizes(x_sz);
  items_to_send_.set_data(x_data);
  items_to_recv_.set_list_sizes(y_sz);
  items_to_recv_.set_data(y_data);
  blocs_to_recv_.set_list_sizes(z_sz);
  blocs_to_recv_.set_data(z_data);
  pe_voisins_ = pe_list;
  blocs_items_count_ = new_blocs_items_count;
  nb_items_to_items_ = new_nb_items_to_items;
}

void MD_Vector_std::fill_md_vect_renum(const IntVect& renum, MD_Vector& md_vect) const
{
  MD_Vector_std dest;

  // Compute nb_items_tot_ and nb_items_reels_
  // nb_items_tot_ is the number of items for which renum[i] >= 0
  // Real items are those that were real in the original array and are to be kept
  {
    const int src_size = get_nb_items_tot();
    const int src_size_r = get_nb_items_reels();
    dest.nb_items_tot_ = 0;
    dest.nb_items_reels_ = 0;
    for (int i = 0; i < src_size; i++)
      {
        if (renum[i] >= 0)
          {
            dest.nb_items_tot_++;
            if (i < src_size_r)
              dest.nb_items_reels_++;
          }
      }
    // If the nb_items_reels_ attribute is invalid in the source vector, it is also invalid in the
    // result vector:
    if (src_size_r < 0)
      dest.nb_items_reels_ = -1;
  }
  // Fill blocs_items_to_sum_ and blocs_items_to_compute_
  // These are the items for which renum[i] >= 0 and which are in the blocs_items_to_sum_
  // and blocs_items_to_compute_ of the source descriptor.
  ::extract_blocs(blocs_items_to_sum_, renum, dest.blocs_items_to_sum_);
  dest.nb_items_seq_local_ = ::extract_blocs(blocs_items_to_compute_, renum, dest.blocs_items_to_compute_);
  dest.nb_items_seq_tot_ = Process::mp_sum(dest.nb_items_seq_local_);

  // ********************************************************
  // Compute items to receive: these are items for which renum[i] >= 0 and which were to be
  // received in the source descriptor.
  // Each time a received item is found, the neighbouring processor that owns it must be notified
  // that this item needs to be sent (using schema_comm). To identify these items, their rank
  // in the array of items to send in the source descriptor is used (variable item_rank).
  //
  const int nb_pe_voisins = pe_voisins_.size_array();
  Schema_Comm schema_comm;
  schema_comm.set_send_recv_pe_list(pe_voisins_, pe_voisins_);
  schema_comm.begin_comm();
  // Compute the items_to_recv lists
  // index+data form the two arrays of a StaticIntLists.
  // We first build as many lists as there are neighbouring processors in the source descriptor.
  // Processors that do not exchange data will be removed at the end.
  ArrOfInt dest_items_recv_index;
  ArrOfInt dest_items_recv_data;
  ArrOfInt dest_blocs_recv_index;
  ArrOfInt dest_blocs_recv_data;

  ArrOfInt tmp;

  {
    dest_items_recv_index.resize_array(nb_pe_voisins + 1, RESIZE_OPTIONS::NOCOPY_NOINIT);
    dest_items_recv_index[0] = 0;
    dest_blocs_recv_index.resize_array(nb_pe_voisins + 1, RESIZE_OPTIONS::NOCOPY_NOINIT);
    dest_blocs_recv_index[0] = 0;
    dest.blocs_items_count_.resize_array(nb_pe_voisins, RESIZE_OPTIONS::NOCOPY_NOINIT);
    dest.nb_items_to_items_.resize_array(nb_pe_voisins, RESIZE_OPTIONS::NOCOPY_NOINIT);
    // Pre-allocate to maximum size
    dest_items_recv_data.resize_array(items_to_recv_.get_data().size_array(), RESIZE_OPTIONS::NOCOPY_NOINIT);

    dest_items_recv_data.resize_array(0);
    // The number of blocks cannot be predicted; there can be more than in the source

    for (int i_pe = 0; i_pe < nb_pe_voisins; i_pe++)
      {
        tmp.resize_array(0);
        int item_rank = 0; // Counter of received items in the source descriptor
        // 1) Extract individual items
        //    (fill dest_items_recv_index, dest_items_recv_data and tmp)
        const int n = items_to_recv_.get_list_size(i_pe);
        for (int i = 0; i < n; i++)
          {
            const int j = items_to_recv_(i_pe, i);
            const int rj = renum[j];
            if (rj >= 0)
              {
                // Add this item
                dest_items_recv_data.append_array(rj);
                // Request that the neighbouring processor sends this item
                tmp.append_array(item_rank);
              }
            item_rank++;
          }
        // Set the list size for items of this processor:
        dest_items_recv_index[i_pe+1] = dest_items_recv_data.size_array();

        // 2) Extract blocks of items received from this processor
        //    (fill dest_blocs_recv_index, dest_blocs_recv_data, dest.blocs_items_count_ and tmp)
        const int nblocs = blocs_to_recv_.get_list_size(i_pe) / 2;
        int received_count = 0; // number of items received from this processor
        for (int ibloc = 0; ibloc < nblocs; ibloc++)
          {
            const int jdeb = blocs_to_recv_(i_pe, ibloc * 2);
            const int jfin = blocs_to_recv_(i_pe, ibloc * 2 + 1);
            int last_added = -2; // can never equal rj-1
            // For each item in the block, if it is to be kept, report it in the blocks to receive:
            for (int j = jdeb; j < jfin; j++)
              {
                const int rj = renum[j];
                if (rj >= 0)
                  {
                    if (rj <= last_added)
                      {
                        // Items in the receive blocks must be sorted in ascending order
                        // (otherwise this algorithm does not work; a different one is needed to
                        //  build the blocks, and furthermore counting the real items may no longer make sense)
                        Cerr << "Error in creer_md_vect_renum: renum array is not sorted: cannot extract blocs to recv" << finl;
                        Process::exit();
                      }
                    // Is this item contiguous with the previous one?
                    if (last_added == rj-1)
                      {
                        // Yes, extend the last block
                        const int k = dest_blocs_recv_data.size_array();
                        // end of block = index of last + 1
                        dest_blocs_recv_data[k-1] = rj+1;
                      }
                    else
                      {
                        // No, create a new block
                        dest_blocs_recv_data.append_array(rj);
                        dest_blocs_recv_data.append_array(rj+1);
                      }
                    last_added = rj;
                    // Request that the neighbouring processor sends this item
                    tmp.append_array(item_rank);
                    received_count++;
                  }
                item_rank++;
              }
          }
        // Set the list size for blocks of this processor:
        dest_blocs_recv_index[i_pe+1] = dest_blocs_recv_data.size_array();
        // Number of items received from this processor
        dest.blocs_items_count_[i_pe] = received_count;
        // Send the tmp array to the neighbouring processor:
        schema_comm.send_buffer(pe_voisins_[i_pe]) << tmp;
      }
  }

  schema_comm.echange_taille_et_messages();

  // ********************************************************
  // Build dest.items_to_send with the received information
  //
  ArrOfInt dest_items_send_index;
  dest_items_send_index.resize_array(nb_pe_voisins + 1, RESIZE_OPTIONS::NOCOPY_NOINIT);
  ArrOfInt dest_items_send_data;
  ArrOfInt nb_items_to_items(nb_pe_voisins); // Initialized to zero

  // Allocate to upper-bound size:
  dest_items_send_data.resize_array(items_to_send_.get_data().size_array(), RESIZE_OPTIONS::NOCOPY_NOINIT);
  dest_items_send_index[0] = 0;
  {
    int count = 0;
    int error = 0;
    for (int i_pe = 0; i_pe < nb_pe_voisins; i_pe++)
      {
        // nb_item_received_single is the number of items in items_to_send_ that will be received
        //  by items_to_recv_. The subsequent ones are received in blocs_to_recv_.
        schema_comm.recv_buffer(pe_voisins_[i_pe]) >> tmp;
        const int n = tmp.size_array();
        for (int i = 0; i < n; i++)
          {
            // tmp contains indices pointing into items_to_send_
            int j = tmp[i];
            // k is the index in the old array of the item to send
            int k = items_to_send_(i_pe, j);
            int renum_k = renum[k];
            // A possible error here: a processor needs a virtual item and the
            // real item has been removed on the source processor. A different algorithm
            // must then be used to change the owner of the item...
            if (renum_k < 0)
              error = 1;
            dest_items_send_data[count++] = renum_k;
          }
        dest_items_send_index[i_pe+1] = count;
      }
    if (error)
      {
        Cerr << "Internal error in MD_Vector_tools::creer_md_vect_renum:\n"
             << " processor " << Process::me()
             << " received a request for distant items that are not in the new vector" << finl;
        Process::exit();
      }
    // Adjust to the final size

    dest_items_send_data.resize_array(count);
  }
  schema_comm.end_comm();

  /* new nb_items_to_items */
  schema_comm.begin_comm();
  for (int i_pe = 0; i_pe < nb_pe_voisins; i_pe++)
    schema_comm.send_buffer(pe_voisins_[i_pe]) << dest_items_recv_index(i_pe + 1) - dest_items_recv_index(i_pe);
  schema_comm.echange_taille_et_messages();
  for (int i_pe = 0; i_pe < nb_pe_voisins; i_pe++)
    schema_comm.recv_buffer(pe_voisins_[i_pe]) >> nb_items_to_items[i_pe];
  schema_comm.end_comm();

  // Build the list of neighbouring processors for dest
  //  (processors with whom data is exchanged)
  // And compress the StaticIntLists indexes to remove the
  //  deleted processors.
  {
    int pe_count = 0;
    dest.pe_voisins_.resize_array(nb_pe_voisins, RESIZE_OPTIONS::NOCOPY_NOINIT);
    for (int i_pe = 0; i_pe < nb_pe_voisins; i_pe++)
      {
        int flag = 0;
        if (dest_items_recv_index[i_pe+1] - dest_items_recv_index[i_pe] > 0)
          flag = 1;
        if (dest_blocs_recv_index[i_pe+1] - dest_blocs_recv_index[i_pe] > 0)
          flag = 1;
        if (dest_items_send_index[i_pe+1] - dest_items_send_index[i_pe] > 0)
          flag = 1;

        if (flag)
          {
            dest_items_recv_index[pe_count+1] = dest_items_recv_index[i_pe+1];
            dest_blocs_recv_index[pe_count+1] = dest_blocs_recv_index[i_pe+1];
            dest_items_send_index[pe_count+1] = dest_items_send_index[i_pe+1];
            dest.blocs_items_count_[pe_count] = dest.blocs_items_count_[i_pe];
            dest.pe_voisins_[pe_count] = pe_voisins_[i_pe];
            dest.nb_items_to_items_[pe_count] = nb_items_to_items[i_pe];
            pe_count++;
          }
      }

    dest_items_recv_index.resize_array(pe_count+1);

    dest_blocs_recv_index.resize_array(pe_count+1);

    dest_items_send_index.resize_array(pe_count+1);

    dest.pe_voisins_.resize_array(pe_count);

    dest.blocs_items_count_.resize_array(pe_count);

    dest.nb_items_to_items_.resize_array(pe_count);
  }
  dest.items_to_send_.set_index_data(dest_items_send_index, dest_items_send_data);
  dest.items_to_recv_.set_index_data(dest_items_recv_index, dest_items_recv_data);
  dest.blocs_to_recv_.set_index_data(dest_blocs_recv_index, dest_blocs_recv_data);

  md_vect.copy(dest);
}

