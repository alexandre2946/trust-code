/****************************************************************************
* Copyright (c) 2022, CEA
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
#ifndef InputCommBuffer_included
#define InputCommBuffer_included
#include <Entree.h>
#include <sstream>
class OutputCommBuffer;
using std::istringstream;
/*! @brief : Helper class used exclusively by Schema_Comm.
 *
 * This is a class
 *   derived from Entree whose stream is of type istringstream (data
 *   read by operator>> are taken from a buffer).
 *   The class is used as follows:
 *   (1) reserve a buffer of a given size with
 *     char * buf = input_comm_buffer.reserve_buffer(size);
 *   (2) fill the buffer with data:
 *     for (i=0; i<size; i++)
 *        buf[i] = .....;
 *   (3) create a stream from the buffer:
 *     input_comm_buffer.create_stream();
 *   (4) data can then be read through operator>>:
 *     input_comm_buffer >> x >> y >> string >> ... ;
 *   (5) when done reading with operator>>, call
 *     input_comm_buffer.clear();
 *   and step (1) can be repeated.
 *
 */

class InputCommBuffer : public Entree
{
public:
  InputCommBuffer();
  ~InputCommBuffer() override;
  // Specifies the buffer size and returns its address.
  // The user must then fill it with data (the entire buffer
  // must be filled).
  char * reserve_buffer(int bufsize);
  void   create_stream_from_output_stream(OutputCommBuffer&);
  // Create the read stream from the buffer.
  void create_stream();
  void clear();
private:
  istringstream * stream_;
  char * buffer_;
  int size_;
  int memorysize_; // memorysize_ >= size_ (allocated capacity)
};
#endif
