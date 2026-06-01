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

#include <TRUSTVect.h>
#include <Memoire.h>
#include <Nom.h>

static int compteur=0;
Memoire* Memoire::_instance = 0;
int Memoire::prems=0;
int Memoire::step=4096;
//struct mallinfo minfo;
static int max_sz_mem=0;
static int min_sz_mem=0;

/*! @brief Returns a pointer to the memory instance. Creates a new memory object if no instance has been created yet.
 *
 * @return (Memoire*) pointer to the memory instance
 */
Memoire& Memoire::Instance()
{
  if (_instance == 0) _instance = new Memoire;

  return *_instance;
}

/*! @brief Constructor. Initializes a working area for Objet_U, "double", and "int" objects.
 *
 */
Memoire::Memoire() :
  size(step), data(new Memoire_ptr[step])
{
  for(int i=prems; i<size; i++)
    data[i].next=i+1;
}

/*! @brief Adds an Objet_U to the TRUST memory.
 *
 * @param (Objet_U* obj) pointer to the Objet_U to add
 * @return (int) the index assigned to the object in memory
 */
int Memoire::add(Objet_U* obj)
{
  assert(_instance != 0);
  compteur++;
  if(prems<size)
    {
      const int num=prems;
      assert(data[num].libre());
      data[num].set(obj);
      prems=data[num].next;
      data[num].next = -1;
      assert(prems>=0);
      //      minfo = mallinfo();
      //      int sz_mem = minfo.arena;
      //      if(sz_mem > max_sz_mem) max_sz_mem = sz_mem;
      return num;
    }
  int old_size=size;
  int i;
  size*=2;
  // Cerr<<"size "<<size<<" "<<size/step<<finl;
  Memoire_ptr* newdata=new Memoire_ptr[size];
  if(!newdata)
    {
      Cerr << "Not enough memory " << finl;
      Process::exit();
    }
  Memoire_ptr  newptr(obj);
  for(i=0; i<old_size; i++)
    {
      newdata[i]=data[i];
      data[i].set(0);
    }
  newdata[old_size]=newptr;
  for(i=old_size; i<size; i++)
    newdata[i].next=i+1;
  if(data)
    delete[] data;
  data=newdata;
  prems=old_size+1;
  //   minfo = mallinfo();
  //   int sz_mem = minfo.arena;
  //   if(sz_mem > max_sz_mem) max_sz_mem = sz_mem;
  return old_size;
}

/*! @brief Removes from memory the Objet_U with index num. The Objet_U itself is not deleted; only its pointer in memory is removed.
 *
 * @param (int num) the index of the Objet_U to remove
 * @return (int) return code, always returns 1
 */
int Memoire::suppr(int num)
{
  //   Cerr << "Suppression de " << num << finl;
  assert(!data[num].libre());
  data[num].next=prems;
  prems=num;
  data[num].set(0);
  compteur--;
#ifndef _COMPILE_AVEC_PGCC
  if((size>4*compteur)&&(size>step))
    //if((size-compteur)>step)
    {
      compacte(); // Obscure crash with the NVidia compiler at end of computation
    }
#endif
  /*
    static int deb=0;
    if (compteur >500) deb=1;
    if ((deb)&&(compteur<3))
    {
    Cerr<<"here "<<compteur<<finl;

    int i;
    for(i=0; i<size; i++)
    if(! data[i].libre())
    {
    Cerr << finl;
    Cerr << i << " ";
    const Objet_U& obj=objet_u(i);
    Cerr << " :: ";
    Cerr << "TYPE :" << obj.le_type()<<finl;
    if(sub_type(Nom,obj)) Cerr << " NAME : " << obj.le_nom() <<finl;

    }
    }
  */

  //GF when there are no more objet_u, destroy everything

  if (compteur==0)
    {
      //Cerr<<"delete data"<<finl;
      if (data)
        {
          delete[] data;
          data=0;
          delete _instance;
          _instance=0;

        }
    }
  return 1;
}


/*! @brief Returns the index in memory of the object with the given type and name.
 *
 * @param (const Nom& type) the type of the object
 * @param (const Nom& nom) the name of the object
 * @return (int) the index of the object if found in memory, -1 otherwise
 */
int Memoire::rang(const Nom& type, const Nom& nom) const
{
  Memoire_ptr* x=data;
  for(int i=0; i<size; i++)
    {
      if( !( x->libre() ) )
        if ( ((x->obj()).que_suis_je()==type) &&
             ((x->obj()).le_nom()==nom) )
          return i;
      x++;
    }
  return -1;
}


/*! @brief Returns the index in memory of the object with the given name.
 *
 * @param (const Nom& nom) the name of the object
 * @return (int) the index of the object if found in memory, -1 otherwise
 */
int Memoire::rang(const Nom& nom) const
{
  Memoire_ptr* x=data;
  for(int i=0; i<size; i++)
    if( !( x->libre() ) )
      {
        if (((x->obj()).le_nom())==nom)
          return i;
        x++;
      }
  return -1;
}


/*! @brief Returns a reference to the Objet_U at index num in memory.
 *
 * @param (int num) the index of the object in memory
 * @return (Objet_U&) reference to the found Objet_U
 * @throws Exits with an error if the memory has an error at index num
 */
Objet_U& Memoire::objet_u(int num)
{
  assert(num >=0);
  assert(num < size);
  assert(!data[num].libre());
  Objet_U& objet = data[num].obj();
  const int num_obj = objet.numero();
  if (num_obj != num)
    {
      Cerr << "Error in Objet_U & Memoire::objet_u(int num) " << finl;
      Cerr << " num                       = " << num << finl;
      std::cerr << " data[num].o_ptr           = " << &objet << std::endl;
      Cerr << " data[num].o_ptr->numero() = " << num_obj << finl;
      assert(0);
      Process::exit();
    }
  return objet;
}


/*! @brief Returns a const reference to the Objet_U at index num in memory.
 *
 * @param (int num) the index of the object in memory
 * @return (const Objet_U&) const reference to the found Objet_U
 * @throws Exits with an error if the memory has an error at index num
 */
const Objet_U& Memoire::objet_u(int num) const
{
  assert(num >=0);
  assert(num < size);
  assert(!data[num].libre());
  Objet_U& objet = data[num].obj();
  const int num_obj = objet.numero();
  if (num_obj != num)
    {
      Cerr << "Error in const Objet_U & Memoire::objet_u(int num) " << finl;
      Cerr << " num                       = " << num << finl;
      std::cerr << " data[num].o_ptr           = " << &objet << std::endl;
      Cerr << " data[num].o_ptr->numero() = " << num_obj << finl;
      assert(0);
      Process::exit();
    }
  return objet;
}


/*! @brief Returns a pointer to the Objet_U at index num in memory.
 *
 * @param (int num) the index of the object in memory
 * @return (Objet_U*) pointer to the found Objet_U
 * @throws Exits with an error if the memory has an error at index num
 */
Objet_U* Memoire::objet_u_ptr(int num)
{
  assert(num >=0);
  assert(num < size);
  assert(!data[num].libre());
  Objet_U& objet = data[num].obj();
  const int num_obj = objet.numero();
  if (num_obj != num)
    {
      Cerr << "Error in Objet_U * Memoire::objet_u(int num) " << finl;
      Cerr << " num                       = " << num << finl;
      std::cerr << " data[num].o_ptr           = " << &objet << std::endl;
      Cerr << " data[num].o_ptr->numero() = " << num_obj << finl;
      assert(0);
      Process::exit();
    }
  return &objet;
}


/*! @brief Returns a const pointer to the Objet_U at index num in memory.
 *
 * @param (int num) the index of the object in memory
 * @return (const Objet_U*) const pointer to the found Objet_U
 * @throws Exits with an error if the memory has an error at index num
 */
const Objet_U* Memoire::objet_u_ptr(int num)const
{
  assert(num >=0);
  assert(num < size);
  assert(!data[num].libre());
  Objet_U& objet = data[num].obj();
  const int num_obj = objet.numero();
  if (num_obj != num)
    {
      Cerr << "Error in const Objet_U * Memoire::objet_u(int num) " << finl;
      Cerr << " num                       = " << num << finl;
      std::cerr << " data[num].o_ptr           = " << &objet << std::endl;
      Cerr << " data[num].o_ptr->numero() = " << num_obj << finl;
      assert(0);
      Process::exit();
    }
  return &objet;
}

/*! @brief Compacts memory. This compaction is performed automatically when it becomes necessary.
 *
 */
void Memoire::compacte()
{
  int i, compte=0;
  int* newnum=new int[size];
  if(!newnum)
    {
      Cerr << "Not enough memory " << finl;
      Process::exit();
    }
  for(i=0; i<size; i++)
    if(data[i].libre())
      newnum[i]=-1;
    else
      newnum[i]=compte++;
  for(i=0; i<size; i++)
    if(newnum[i]>-1)
      {
        objet_u(i).change_num(newnum);
      }
  Memoire_ptr* newdata=new Memoire_ptr[compte];
  if(!newdata)
    {
      Cerr << "Not enough memory " << finl;
      Process::exit();
    }
  for(i=0; i< size; i++)
    if(newnum[i]>-1)
      newdata[newnum[i]]=data[i];
  delete[] data;
  delete[] newnum;
  data=newdata;
  size=prems=compte;
  for(i=0; i<size; i++)
    data[i].next=i+1;
  verifie();
}

/*! @brief Prints a memory state summary to the error output stream.
 *
 * @return (int) return code; always returns 1
 */
int Memoire::imprime() const
{
  assert(prems <=size);
  int i, tmp=0;
  for(i=0; i<size; i++)
    if(! data[i].libre()) tmp++;
  Cerr << "\n There are " << size << " memory squares";
  Cerr << "\n of which " << tmp << " ==(const char*)" << compteur << " are used" << finl;
  for(i=0; i<size; i++)
    {
      Cerr << i << " :: " << flush;
      if(! (data[i].libre()) )
        {
          const Objet_U& obj=data[i].obj();
          if (obj.numero() != i)
            Cerr << "error";
          else
            {
              Cerr << "TYPE :" << obj.le_type()<< flush;
              if (sub_type(Objet_U_ptr, obj))
                {
                  const Objet_U_ptr& x=ref_cast(Objet_U_ptr, obj);
                  Cerr << " key: " << x << flush;
                }
              else
                Cerr << " NAME : " << obj.le_nom() << flush;
              Cerr << " address : " << (long)(&(obj))<< flush;
              Cerr << " SIZE : " << (int)(obj.taille_memoire())<< " bytes"<< flush;
            }
        }
      else
        Cerr << "free ";
      Cerr << finl;
    }
  return 1;
}

/*! @brief Verifies the content of all memory slots.
 *
 * @return (int) return code; always returns 1
 */
int Memoire::verifie() const
{
  assert(prems <=size);
  for(int i=0; i<size; i++)
    {
      if(! (data[i].libre()) )
        {
          if ((data[i].obj()).numero() != i)
            {
              Cerr << "ERROR at the square " << i << finl;
              imprime();
              Process::exit();
            }
        }
    }
  return 1;
}

/*! @brief Output operator that prints the state of memory mem to output stream os.
 *
 * @param (Sortie& os) the output stream to use
 * @param (const Memoire& mem) the memory to examine
 * @return (Sortie&) the modified output stream
 */
Sortie& operator << (Sortie& os, const Memoire& mem)
{
  int i;
  int total=0;
  int tmp=0;
  int tmp1=0;
  for(i=0; i<mem.size; i++)
    if(! mem.data[i].libre()) tmp++;
  os << "\n There are " << mem.size << " memory slots";
  os << "\n of which " << tmp << " ==(const char*)" << compteur << " are used" << finl;
  for(i=0; i<mem.size; i++)
    if(! mem.data[i].libre())
      {
        os << i << " ";
        const Objet_U& obj=mem.objet_u(i);
        os << " :: ";
        os << "TYPE :" << obj.le_type();
        if(sub_type(Nom,obj)) os << " NAME : " << obj.le_nom() ;
        os << " address : " << (long)(&(obj));

        const ArrOfInt* intA = dynamic_cast<const ArrOfInt*>(&obj);
        const ArrOfDouble* intAD = dynamic_cast<const ArrOfDouble*>(&obj);
        if(intA)
          {
            const ArrOfInt& toto = *intA;
            tmp = obj.taille_memoire() + (int) ((toto.size_array()*sizeof(int))/toto.ref_count());
            os << " SIZE ArrOfInt:" << tmp<< " bytes";

            const IntVect* intV = dynamic_cast<const IntVect*>(&obj);
            if(intV)
              {
                const IntVect& titi = *intV;
                tmp1 = obj.taille_memoire() + (int) ((titi.size()*sizeof(int))/titi.ref_count());
                os << " of which:" << (tmp-tmp1) << " virtual bytes ";
              }
            os << "\n";
          }
        else if (intAD)
          {
            const ArrOfDouble& toto = *intAD;
            tmp = obj.taille_memoire() + (int) ((toto.size_array()*sizeof(double))/toto.ref_count());
            os << " SIZE ArrOfDouble:" << tmp<< " bytes";
            /* PL: Plante en P1Bulle donc je commente
             * const DoubleVect* intVD = dynamic_cast<const DoubleVect*>(&obj);
               if(intVD){
               const DoubleVect& titi = *intVD;
               tmp1 = obj.taille_memoire() + (int) ((titi.size()*sizeof(double))/titi.ref_count());
               os << " of which:" << (tmp-tmp1) << " virtual bytes ";
               } */
            os << "\n";
          }
        else
          {
            const ArrOfFloat* intAF = dynamic_cast<const ArrOfFloat*>(&obj);
            if(intAF)
              {
                const ArrOfFloat& toto = *intAF;
                tmp = obj.taille_memoire() + (int) ((toto.size_array()*sizeof(float))/toto.ref_count());
                os << " SIZE ArrOfFloat:" << tmp<< " bytes";
                os << "\n";
              }
            else
              {
                os << " SIZE : " << (tmp=obj.taille_memoire())<< " bytes\n ";
              }
          }
        total += tmp;
        os.flush();
      }
  os << "Max memory size (MB): " <<  max_sz_mem << finl;
  os << "Min memory size (MB): " <<  min_sz_mem << finl;
  os << "Used memory size (MB): " <<  (max_sz_mem-min_sz_mem)/1024/1024 << finl;

  return os << "Total occupied memory size (MB): " << total/1024/1024 << finl;
}

