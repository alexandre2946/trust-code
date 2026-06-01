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

#include <Redistribute_Field.h>
#include <Linear_algebra_tools.h>
#include <Perf_counters.h>

void Redistribute_Field::initialize(const Domaine_IJK& input,
                                    const Domaine_IJK& output,
                                    const Domaine_IJK::Localisation loc)
{
  // empty map => ok for meshes of identical sizes
  VECT(IntTab) map;
  initialize(input, output, loc, map);
}

void Redistribute_Field::initialize(const Domaine_IJK& input,
                                    const Domaine_IJK& output,
                                    const Domaine_IJK::Localisation loc,
                                    const VECT(IntTab) & redistribute_maps)
{
  IntTab map;
  // Compute the send_blocs arrays
  for (int dir = 0; dir < 3; dir++)
    {
      if (redistribute_maps.size() == 0)
        {
          map.resize(1,3);
          map(0,0) = 0;
          map(0,1) = 0;
          map(0,2) = input.get_nb_items_global(loc, dir);
          if (output.get_nb_items_global(loc, dir) != map(0,2))
            {
              Cerr << "Error in Redistribute_Field::initialize: you must provide a redistribute map" << finl;
              Process::exit();
            }
        }
      else
        {
          map = redistribute_maps[dir];
        }
      compute_send_blocs(input, output, loc, dir, map, send_blocs_[dir]);
      // invert the map to compute the recv_blocs
      for (int i = 0; i < map.dimension(0); i++)
        {
          int x = map(i,0);
          map(i,0) = map(i,1);
          map(i,1) = x;
        }
      compute_send_blocs(output, input, loc, dir, map, recv_blocs_[dir]);
    }

  // Compute the communication schedule (list of neighbouring processors and data sizes)
  const FixedVector<IntTab, 3>& send_blocs = send_blocs_;
  const FixedVector<IntTab, 3>& recv_blocs = recv_blocs_;
  Schema_Comm_Vecteurs& schema = schema_comm_;
  schema.begin_init();
  size_buffer_for_me_= 0;
  int receive_size_buffer_for_me = 0;

  // Compute send blocs sizes and destination processors
  Int3 ibloc;
  for (ibloc[2] = 0; ibloc[2] < send_blocs[2].dimension(0); ibloc[2] ++)
    for (ibloc[1] = 0; ibloc[1] < send_blocs[1].dimension(0); ibloc[1] ++)
      for (ibloc[0] = 0; ibloc[0] < send_blocs[0].dimension(0); ibloc[0] ++)
        {
          const int dest_slice_i  = send_blocs[0](ibloc[0], 1);
          const int ni            = send_blocs[0](ibloc[0], 2);
          const int dest_slice_j  = send_blocs[1](ibloc[1], 1);
          const int nj            = send_blocs[1](ibloc[1], 2);
          const int dest_slice_k  = send_blocs[2](ibloc[2], 1);
          const int nk            = send_blocs[2](ibloc[2], 2);

          const int dest_pe = output.get_processor_by_ijk(dest_slice_i, dest_slice_j, dest_slice_k);
          if (dest_pe != Process::me())
            schema.add_send_area_template<double>(dest_pe, ni*nj*nk);
          else
            size_buffer_for_me_ += ni*nj*nk;
        }

  // Compute recv bloc sizes and source processors
  for (ibloc[2] = 0; ibloc[2] < recv_blocs[2].dimension(0); ibloc[2] ++)
    for (ibloc[1] = 0; ibloc[1] < recv_blocs[1].dimension(0); ibloc[1] ++)
      for (ibloc[0] = 0; ibloc[0] < recv_blocs[0].dimension(0); ibloc[0] ++)
        {
          const int src_slice_i   = recv_blocs[0](ibloc[0], 1);
          const int ni            = recv_blocs[0](ibloc[0], 2);
          const int src_slice_j   = recv_blocs[1](ibloc[1], 1);
          const int nj            = recv_blocs[1](ibloc[1], 2);
          const int src_slice_k   = recv_blocs[2](ibloc[2], 1);
          const int nk            = recv_blocs[2](ibloc[2], 2);

          const int src_pe = input.get_processor_by_ijk(src_slice_i, src_slice_j, src_slice_k);
          if (src_pe != Process::me())
            schema.add_recv_area_template<double>(src_pe, ni*nj*nk);
          else
            receive_size_buffer_for_me  += ni*nj*nk;
        }

  if (size_buffer_for_me_ != receive_size_buffer_for_me)
    {
      Cerr << "Internal error in Redistribute_Field::initialize (size_buffer_for_me)" << finl;
      Process::exit();
    }
  schema.end_init();
}

// Computes the intersection of segments [s1, s1 + n1 - 1] and [s2, s2 + n2 - 1],
// stores the result in s2 and n2.
// If the origin s2 is modified, updates s3 as:
//   new s3 = old s3 + new s2 - old s2
void Redistribute_Field::intersect(const int s1, const int n1, int& s2, int& n2, int& s3)
{
  // end of segment (plus one)
  int e1 = s1 + n1;
  int e2 = s2 + n2;
  int old_s2 = s2;
  int zero = 0;
  s2 = std::max(s1, s2);
  e2 = std::min(e1, e2);
  n2 = std::max(e2 - s2, zero);
  s3 += s2 - old_s2;
}

// Computes how to redistribute the local data of this processor from an input field
// into an output field (computation of the data to send to other processors)
// for direction "dir".
// global_index_mapping: array with three columns.
//   column 1: global index of the first segment element to copy from input
//   column 2: global index of the beginning of the destination segment (in output)
//   column 3: number of elements in the segment to copy
// To redistribute a simple field, pass a global_index_mapping array
//   with one row and three columns: { 0, 0, n } where n is the global element count of the field.
// To redistribute with periodisation (building the extended field for front-tracking),
//   pass an array with three rows:
//     row 0, the central zone:   { 0, extend_size, ni }  (ni is the non-extended field size)
//     row 1, the left zone:      { ni-extend_size, 0, extend_size }
//     row 2, the right zone:     { 0, ni + extend_size, extend_size }
//
// Fills the send_blocs array: list of data segments to send to other processors
//   (column 0: local index of the first element to send,
//    column 1: destination slice number,
//    column 2: number of consecutive elements to send)
void Redistribute_Field::compute_send_blocs(const Domaine_IJK& input,
                                            const Domaine_IJK& output,
                                            const Domaine_IJK::Localisation localisation,
                                            const int dir,
                                            const IntTab& global_index_mapping,
                                            IntTab& send_blocs)
{
  // Distribution of the output field over the different processor slices:
  ArrOfInt output_slice_offsets;
  ArrOfInt output_slice_size;
  output.get_slice_offsets(dir, output_slice_offsets);
  output.get_slice_size(dir, localisation, output_slice_size);

  // Which part of the input field is on this processor:
  const int input_slice_start = input.get_offset_local(dir);
  const int input_slice_size = input.get_nb_items_local(localisation, dir);
  send_blocs.resize(0,3);
  // Loop over the segments to copy (one per row of global_index_mapping)
  const int n_segments = global_index_mapping.dimension(0);
  for (int i_segment = 0; i_segment < n_segments; i_segment++)
    {
      int input_seg_start = global_index_mapping(i_segment, 0);
      int output_seg_start = global_index_mapping(i_segment, 1);
      int seg_length = global_index_mapping(i_segment, 2);

      // Intersection of the segment to copy with the local data of this processor:
      intersect(input_slice_start, input_slice_size, input_seg_start, seg_length, output_seg_start);

      // Loop over the output field slices:
      const int nb_output_slices = output_slice_offsets.size_array();
      for (int oslice = 0; oslice < nb_output_slices; oslice++)
        {
          int input_start = input_seg_start;
          int output_start = output_seg_start;
          int n = seg_length;
          // Does this slice oslice intersect with the current output_segment?
          intersect(output_slice_offsets[oslice], output_slice_size[oslice], output_start, n, input_start);
          if (n > 0)
            {
              const int input_start_local = input_start - input_slice_start;
              send_blocs.append_line(input_start_local, oslice, n);
            }
        }
    }
}

void Redistribute_Field::redistribute_(const IJK_Field_double& input_field,
                                       IJK_Field_double& output_field,
                                       bool add)
{
  statistics().create_custom_counter("Redistribute",2,"IJK");
  statistics().begin_count("Redistribute",statistics().get_last_opened_counter_level()+1);

  // For each direction, intersections between the input and output field:
  // Each line contain a bloc to copy from-to another processor
  // Column 0: local index in first field of the first element of the bloc
  // Column 1: slice number in the other field
  // Column 2: number of elements to copy
  const FixedVector<IntTab, 3>& send_blocs = send_blocs_;
  const FixedVector<IntTab, 3>& recv_blocs = recv_blocs_;
  Schema_Comm_Vecteurs& schema = schema_comm_;
  schema.begin_comm();

  ArrOfDouble buffer_for_me;
  ArrOfDouble tmp;
  buffer_for_me.resize_array(size_buffer_for_me_, RESIZE_OPTIONS::NOCOPY_NOINIT);
  int index_buffer_for_me = 0;
  // Prepare data to send:
  Int3 ibloc;
  for (ibloc[2] = 0; ibloc[2] < send_blocs[2].dimension(0); ibloc[2] ++)
    for (ibloc[1] = 0; ibloc[1] < send_blocs[1].dimension(0); ibloc[1] ++)
      for (ibloc[0] = 0; ibloc[0] < send_blocs[0].dimension(0); ibloc[0] ++)
        {
          const int input_i_start = send_blocs[0](ibloc[0], 0);
          const int dest_slice_i  = send_blocs[0](ibloc[0], 1);
          const int ni            = send_blocs[0](ibloc[0], 2);
          const int input_j_start = send_blocs[1](ibloc[1], 0);
          const int dest_slice_j  = send_blocs[1](ibloc[1], 1);
          const int nj            = send_blocs[1](ibloc[1], 2);
          const int input_k_start = send_blocs[2](ibloc[2], 0);
          const int dest_slice_k  = send_blocs[2](ibloc[2], 1);
          const int nk            = send_blocs[2](ibloc[2], 2);

          const int dest_pe = output_field.get_domaine().get_processor_by_ijk(dest_slice_i, dest_slice_j, dest_slice_k);
          // If the destination processor is this process, use a local array without going through MPI:
          ArrOfDouble *buf_ptr;
          if (dest_pe == Process::me())
            {
              // Point tmp to a zone of the local buffer of the correct size:
              tmp.ref_array(buffer_for_me, index_buffer_for_me /* start index */, ni*nj*nk /* size */);
              index_buffer_for_me += ni*nj*nk; // advance by the block size in the buffer
              buf_ptr = &tmp;
            }
          else
            {
              // The buffer is in the communication schedule:
              buf_ptr = &schema.get_next_area_template<double>(dest_pe, ni*nj*nk);
            }
          ArrOfDouble& buf = *buf_ptr;

          for (int k = 0; k < nk; k++)
            {
              for (int j = 0; j < nj; j++)
                {
                  for (int i = 0; i < ni; i++)
                    {
                      buf[k*nj*ni + j*ni + i] = input_field(i + input_i_start, j + input_j_start, k + input_k_start);
                    }
                }
            }
        }

  // Exchange data
  schema.exchange();
  index_buffer_for_me = 0;
  // Get data to receive:
  for (ibloc[2] = 0; ibloc[2] < recv_blocs[2].dimension(0); ibloc[2] ++)
    for (ibloc[1] = 0; ibloc[1] < recv_blocs[1].dimension(0); ibloc[1] ++)
      for (ibloc[0] = 0; ibloc[0] < recv_blocs[0].dimension(0); ibloc[0] ++)
        {
          const int output_i_start= recv_blocs[0](ibloc[0], 0);
          const int src_slice_i   = recv_blocs[0](ibloc[0], 1);
          const int ni            = recv_blocs[0](ibloc[0], 2);
          const int output_j_start= recv_blocs[1](ibloc[1], 0);
          const int src_slice_j   = recv_blocs[1](ibloc[1], 1);
          const int nj            = recv_blocs[1](ibloc[1], 2);
          const int output_k_start= recv_blocs[2](ibloc[2], 0);
          const int src_slice_k   = recv_blocs[2](ibloc[2], 1);
          const int nk            = recv_blocs[2](ibloc[2], 2);

          const int src_pe = input_field.get_domaine().get_processor_by_ijk(src_slice_i, src_slice_j, src_slice_k);
          ArrOfDouble *buf_ptr;
          if (src_pe == Process::me())
            {
              // Point tmp to a zone of the local buffer of the correct size:
              tmp.ref_array(buffer_for_me, index_buffer_for_me /* start index */, ni*nj*nk /* size */);
              index_buffer_for_me += ni*nj*nk; // advance by the block size in the buffer
              buf_ptr = &tmp;
            }
          else
            {
              // The buffer is in the communication schedule:
              buf_ptr = &schema.get_next_area_template<double>(src_pe, ni*nj*nk);
            }
          ArrOfDouble& buf = *buf_ptr;
          if (!add)
            for (int k = 0; k < nk; k++)
              {
                for (int j = 0; j < nj; j++)
                  {
                    for (int i = 0; i < ni; i++)
                      {
                        output_field(i + output_i_start, j + output_j_start, k + output_k_start) = buf[k*nj*ni + j*ni + i];
                      }
                  }
              }
          else
            for (int k = 0; k < nk; k++)
              {
                for (int j = 0; j < nj; j++)
                  {
                    for (int i = 0; i < ni; i++)
                      {
                        output_field(i + output_i_start, j + output_j_start, k + output_k_start) += buf[k*nj*ni + j*ni + i];
                      }
                  }
              }
        }
  schema.end_comm();
  statistics().end_count("Redistribute");

}

void Redistribute_Field::redistribute_(const IJK_Field_float& input_field,
                                       IJK_Field_float& output_field,
                                       bool add)
{
  Cerr << "Redistribute_Field::redistribute_(const IJK_Field_float & input_field,...) not implemented" << finl;
  Process::exit();
}
