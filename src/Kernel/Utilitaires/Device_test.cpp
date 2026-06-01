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

#include <TRUSTTrav.h>
#include <TRUSTTab_parts.h>
#include <Device.h>
#include <string>
#include <sstream>
#include <comm_incl.h>
#include <DeviceMemory.h>

// Device_test: intensively tests the methods of the Device.h interface:
bool self_tested_ = false;
void self_test()
{
  if (self_tested_)
    return;
  else
    self_tested_ = true;
  {
    // local_operations_vect_bis_generic
    // Critical unit test case:
    {
      DoubleTab a(2);
      a(0) = 1;
      a(1) = 2;
      mapToDevice(a);
      double mp = mp_norme_vect(a);
      double sol = sqrt(Process::nproc() * 5);
      if (mp != sol)
        {
          Cerr << "What! " << mp << " != " << sol << finl;
          Cerr
              << "Weird bug seen during a sum reduction with Kokkos::parallel_reduce for TRUST production build on Cuda 11.6. Fixed with Cuda 12.1 (or before...)."
              << finl;
          Cerr << "Update Cuda please !" << finl;
          Process::exit();
        }
    }
  }
#ifndef NDEBUG
  // Test mapToDevice(arr)
  // Status
  // Status before		Status after		Copy ?
  // DataLocation::HostOnly	        DataLocation::HostDevice	        Yes
  // Host		        DataLocation::HostDevice	        Yes
  // DataLocation::HostDevice	        DataLocation::HostDevice	        No
  // Device		        Device		        No
  {
    DoubleTab a(10);
    assert(a.get_data_location() == DataLocation::HostOnly);
    mapToDevice(a);
    assert(a.get_data_location() == DataLocation::HostDevice);
    a[1]=0;
    assert(a.get_data_location() == DataLocation::Host);
    mapToDevice(a);
    assert(a.get_data_location() == DataLocation::HostDevice);
  }

  {
    DoubleTab a(10);
    a = 23;
    assert(a.get_data_location() == DataLocation::HostOnly);
    {
      DoubleTab b;
      b.ref(a);
      mapToDevice(b);
      assert(b.get_data_location() == DataLocation::HostDevice);
      assert(a.get_data_location() == DataLocation::HostDevice);
    }
  }

  int N = 10;
  // Test access methods on the device:
  DoubleTab inco(N);
  inco = 1;
  mapToDevice(inco); // copy
  assert(inco.get_data_location() == DataLocation::HostDevice);
  assert(inco.ref_count() == 1);
  {
    // Example: 1st operator
    DoubleTab a;
    a.ref(inco); // Must take the state of inco
    assert(a.get_data_location() == DataLocation::HostDevice);
    assert(a.ref_count() == 2);
    assert(inco.ref_count() == 2);
    DoubleTab b(N);

    CDoubleArrView a_v = static_cast<const ArrOfDouble&>(a).view_ro();
    DoubleArrView b_v = static_cast<ArrOfDouble&>(b).view_wo();
    Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int i)
    {
      b_v[i] = a_v[i];
    });
    Kokkos::fence();

    const DoubleTab& const_b = b;
    const DoubleTab& const_a = a;
    assert(const_b[5] == const_a[5]);
    assert(const_b[5] == 1);
    //assert(b[5] == a[5]); // Argh double& TRUSTArray<double>::operator[](int i) called for a and thus falls back to host
    // How to detect whether operator[](int i) is used for read or write? Possible? No, unless const is used.
    assert(a.get_data_location() == DataLocation::HostDevice);
    assert(b.get_data_location() == DataLocation::HostDevice);
    assert(inco.get_data_location() == DataLocation::HostDevice);
  }
  assert(inco.get_data_location() == DataLocation::HostDevice);
  assert(inco.ref_count() == 1);
  {
    // Example: 2nd operator
    DoubleTab a;
    a.ref(inco); // Must take the state of inco which is still DataLocation::HostDevice
    assert(a.get_data_location() == DataLocation::HostDevice);
    assert(a.ref_count() == 2);
    assert(inco.ref_count() == 2);

    DoubleTab b(N);
    CDoubleArrView a_v = static_cast<const ArrOfDouble&>(a).view_ro();
    DoubleArrView b_v = static_cast<ArrOfDouble&>(b).view_wo();
    Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int i)
    {
      b_v[i] = a_v[i];
    });
    Kokkos::fence();

    const DoubleTab& const_b = b;
    const DoubleTab& const_a = a;
    assert(const_b[5] == const_a[5]);
    assert(const_b[5] == 1);
    assert(a.get_data_location() == DataLocation::HostDevice);
    assert(b.get_data_location() == DataLocation::HostDevice);
  }
  assert(inco.ref_count() == 1);

  // Update of the unknown on the device:
  inco += 1;
  assert(inco.get_data_location() == DataLocation::Device);
  assert(inco.ref_count() == 1);
  {
    // Next time step, new operator:
    DoubleTab a;
    a.ref(inco); // Must take the state of inco
    assert(a.get_data_location() == DataLocation::Device);
    assert(a.ref_count() == 2);
    assert(inco.ref_count() == 2);

    DoubleTab b(N);
    CDoubleArrView a_v = static_cast<const ArrOfDouble&>(a).view_ro();
    DoubleArrView b_v = static_cast<ArrOfDouble&>(b).view_wo();
    Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int i)
    {
      b_v[i] = a_v[i];
    });
    Kokkos::fence();
    const DoubleTab& const_b = b;
    const DoubleTab& const_a = a;
    assert(const_b[5] == const_a[5]);
    assert(const_b[5] == 2);
    assert(a.get_data_location() == DataLocation::HostDevice);
    assert(b.get_data_location() == DataLocation::HostDevice);
    assert(inco.get_data_location() == DataLocation::HostDevice); // Because a refs inco
  }
  assert(inco.get_data_location() == DataLocation::HostDevice);
  // Update of the unknown on the device
  {
    inco += 1;
    assert(inco.get_data_location() == DataLocation::Device);

    DoubleTab a;
    a.ref(inco); // Must take the state of inco
    assert(a.get_data_location() == DataLocation::Device);
    assert(a.ref_count() == 2);
    assert(inco.ref_count() == 2);

    DoubleTab b(N);
    CDoubleArrView a_v = static_cast<const ArrOfDouble&>(a).view_ro();
    DoubleArrView b_v = static_cast<ArrOfDouble&>(b).view_wo();
    Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int i)
    {
      b_v[i] = a_v[i];
    });
    Kokkos::fence();
    const DoubleTab& const_b = b;
    const DoubleTab& const_a = a;
    assert(const_b[5] == const_a[5]);
    assert(const_b[5] == 3);
    assert(a.get_data_location() == DataLocation::HostDevice);
    assert(b.get_data_location() == DataLocation::HostDevice);
    assert(inco.get_data_location() == DataLocation::HostDevice);
  }
  assert(inco.get_data_location() == DataLocation::HostDevice);

  // Test ArrOfDouble operations on GPU
  {
    ArrOfDouble a(10), b(10);
    a=1;
    b=2;
    mapToDevice(a);
    mapToDevice(b);
    b+=a; // TRUSTArray& operator+=(const TRUSTArray& y) on the device
    b+=3; // TRUSTArray& operator+=(const _TYPE_ dy)
    b-=2; // TRUSTArray& operator-=(const _TYPE_ dy)
    b-=a; // TRUSTArray& operator-=(const TRUSTArray& y)
    // ToDo fix: Multiple definition of 'nvkernel__ZN10TRUST
    //b*=10; // TRUSTArray& operator*= (const _TYPE_ dy)
    //b/=2;  // TRUSTArray& operator/= (const _TYPE_ dy)
    assert(a.get_data_location() == DataLocation::HostDevice);
    assert(b.get_data_location() == DataLocation::Device);
    const ArrOfDouble& const_b = b;
    // Operations on the device:
    // Return to host to verify the result
    copyFromDevice(b);
    assert(const_b[0] == 3);
  }

  // Constructors on device:
  {
    DoubleTab a(10);
    allocateOnDevice(a);
    assert(a.get_data_location() == DataLocation::Device);
    a = 1; // On the device
    assert(a.get_data_location() == DataLocation::Device);
    DoubleTab b(a); // b must also be allocated on the device and the copy done on the device
    assert(b.get_data_location() == DataLocation::Device);
    copyFromDevice(b);
    const ArrOfDouble& const_b = b;
    assert(const_b[0] == 1);
    assert(const_b[b.size() - 1] == 1);
  }
  {
    // Array copies:
    DoubleTab a(10);
    a=1;
    mapToDevice(a); // a on the device
    DoubleTab b;
    b = a; // b must also be allocated/filled on the device by copying a:
    assert(b.get_data_location() == DataLocation::Device);
    const ArrOfDouble& const_b = b;
    copyFromDevice(b);
    assert(const_b[0] == 1);
    assert(const_b[b.size() - 1] == 1);
  }
  // operator_vect_vect_generic for DoubleTab::operator+-*/
  {
    DoubleTab a(10), b(10);
    const ArrOfDouble& const_a = a;
    const ArrOfDouble& const_b = b;
    a=1;
    b=10;
    mapToDevice(a);
    assert(a.get_data_location() == DataLocation::HostDevice);
    mapToDevice(b);
    assert(b.get_data_location() == DataLocation::HostDevice);
    a=b;  // TRUSTArray<_TYPE_>::inject_array(v) done on the device (a=10)
    assert(a.get_data_location() == DataLocation::Device);
    a+=b; // operator_vect_vect_generic(ADD_) done on the device (a=20)
    assert(a.get_data_location() == DataLocation::Device);
    b-=a; // operator_vect_vect_generic(SUB_) done on the device (b=-10)
    assert(b.get_data_location() == DataLocation::Device);
    // Return to host to verify results
    copyFromDevice(a);
    copyFromDevice(b);
    assert(a.get_data_location() == DataLocation::HostDevice);
    assert(b.get_data_location() == DataLocation::HostDevice);
    assert(const_a[0] == 20);
    assert(const_a[a.size()-1] == 20);
    assert(const_b[0] == -10);
    assert(const_b[b.size()-1] == -10);
  }
  // DoubleTrav
  {
    DoubleTrav a(10*N);
    a = 1;
    mapToDevice(a); // copy
    assert(a.get_data_location() == DataLocation::HostDevice);
    assert(a.ref_count() == 1);
    DeviceMemory::printMemoryMap();
  }
  // Second DoubleTrav
  {
    DeviceMemory::printMemoryMap();
    DoubleTrav a(10*N); // a is initialized to 0 on the device because the previous DoubleTrav was HostDevice!
    assert(a.get_data_location() == DataLocation::Device);
    const ArrOfDouble& const_a = a;
    assert(const_a[0] == 0);
  }
  // Copy constructor DoubleTab
  {
    DoubleTab a(N);
    a=-1;
    mapToDevice(a); // On the device
    DoubleTrav b(a); // b must be on the device (=0)
    assert(b.get_data_location()==DataLocation::Device);
    b+=1; // Operation must be done on the device (=1)
    assert(b.get_data_location()==DataLocation::Device);
    copyFromDevice(b);
    const ArrOfDouble& const_b = b;
    assert(const_b[0] == 1);
    assert(const_b[b.size()-1] == 1);
  }
  // max/min methods
  {
    DoubleTab a(3);
    a(0)=1;
    a(1)=3;
    a(2)=-10;
    mapToDevice(a);
    // Change on the host for testing:
    a.data()[0]=0;
    a.data()[1]=0;
    a.data()[2]=0;
    a.set_data_location(DataLocation::Device);
    assert(local_max_vect(a)==3);
    assert(local_min_vect(a)==-10);
    assert(local_max_abs_vect(a)==10);
    assert(local_min_abs_vect(a)==1);
    //assert(local_imax_vect(a)==1);
    //assert(local_imin_vect(a)==2);
    assert(a.isDataOnDevice());
    // Test on host the two methods imin,imax not ported to GPU:
    copyFromDevice(a);
    assert(a.get_data_location()==DataLocation::HostDevice);
    assert(local_imax_vect(a)==1);
    assert(local_imin_vect(a)==2);
  }
  // ref_array
  {
    DoubleTab a(N);
    a=1;
    mapToDevice(a); // On the device
    assert(a.get_data_location()==DataLocation::HostDevice);
    DoubleTab b;
    b.ref_array(a);
    assert(b.get_data_location()==DataLocation::HostDevice); // b must be on the device
  }
  // ref_tab/ref_array on a chunk array
  {
    DoubleTab a(2*N);
    a=1;
    DoubleTab b;
    b.ref_tab(a, 0, N); // partial reference on a
    mapToDevice(b); // On the device
    assert(b.get_data_location()==DataLocation::HostDevice); // b must be on the device
    assert(a.get_data_location()==DataLocation::HostDevice); // a must be on the device

    DoubleArrView a_v = static_cast<ArrOfDouble&>(a).view_wo();
    Kokkos::parallel_for(2*N, KOKKOS_LAMBDA(const int i)
    {
      a_v[i] = 2;
    });
    Kokkos::fence();
    // Return to the device and verify that a was completely on the device:
    assert(a.get_data_location()==DataLocation::Device);
    assert(a(0)==2);
    assert(a.get_data_location()==DataLocation::Host);
    assert(a(2*N-1)==2);
  }
  double * ptr_host;
  {
    DoubleTab a(N);
    a=1;
    {
      DoubleTab b;
      b.ref_array(a);
      mapToDevice(b); // On the device
      assert(a.data()==b.data());
      assert(b.get_data_location() == DataLocation::HostDevice);
      assert(a.get_data_location() == DataLocation::HostDevice); // a is also considered on the device
    }
    assert(a.get_data_location() == DataLocation::HostDevice);
    ptr_host = a.data();
    assert(isAllocatedOnDevice(ptr_host)); // Check that the array has a memory area on the device
  }
  assert(!isAllocatedOnDevice(ptr_host)); // Check that the array no longer has a memory area on the device
  // Test of resize_array_ in TRUSTArray<_TYPE_, _SIZE_>::resize_array_
  {
    DoubleTab a(1);
    mapToDevice(a);
    a=1;
    a.resize(2);
    a+=1;
    assert(a.isDataOnDevice());
    assert(a(0)==2); // Check that after resize data on device is keep
    assert(a(1)==1); // Check that resize initialize new elements to 0
  }
  // ToDo: How to handle DoubleTab_Parts? Not easy, so for now
  // the copy constructor does a copyFromDevice on the DoubleTab...
  /*
    {
        DoubleTab pression;
        mapToDevice(pression);
        DoubleTab_parts P0P1(pression); // P0P1 Device
        DoubleTab& P0 = P0P1[0];
        DoubleTab& P1 = P0P1[1];
        double moyenne_K = mp_moyenne_vect(P0); // P0 DataLocation::HostDevice but not P1!
        P1 -= moyenne_K; // Big problem as P1 is always on Device!
    }
    */
  //if (Process::me()==0) std::cerr << ptr_host << std::endl;

  // DoubleTrav copy Constructor:
  {
    DoubleTrav b(100);
    b = 123;
    mapToDevice(b);
    DoubleVect b2(b);
    assert(b2.get_mem_storage() == STORAGE::TEMP_STORAGE);
    assert(b2.get_data_location() == DataLocation::Device);

    DoubleVect b3(b2);
    assert(b3.get_mem_storage() == STORAGE::TEMP_STORAGE);
    assert(b3.get_data_location() == DataLocation::Device);
  }

  {
    // Provisional: reproduce the deadlock?
    //DoubleTab a(10);
    //allocateOnDevice(a);
    //allocateOnDevice(a);
  }
#endif
}


