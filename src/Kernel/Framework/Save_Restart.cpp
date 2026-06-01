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

#include <EcritureLectureSpecial.h>
#include <Entree_Fichier_base.h>
#include <Entree_Fichier_base.h>
#include <Discretisation_base.h>
#include <EcrFicCollecteBin.h>
#include <LecFicDiffuseBin.h>
#include <communications.h>
#include <FichierHDFPar.h>
#include <Probleme_base.h>
#include <Sortie_Nulle.h>
#include <Save_Restart.h>
#include <TRUST_2_PDI.h>
#include <Ecrire_YAML.h>
#include <sys/stat.h>
#include <Avanc.h>
#include <Perf_counters.h>
#define CHECK_ALLOCATE 0
#ifdef CHECK_ALLOCATE
#include <unistd.h> // For access to int close(int fd); with PGI
#include <fcntl.h>
#include <errno.h>
#endif

// Returns the version of the save format
// 151 to say that it is the version initiated at version 1.5.1 of TRUST
inline int version_format_sauvegarde() { return 184; }

// Version using PDI library
inline int version_format_PDI() { return 196; }

/*! @brief Initialization of file_size, bad_allocate, nb_pb_total, num_pb
 *
 */
long int Save_Restart::File_size_=0;        // file_size is the disk space in bytes necessary to write the XYZ files
int Save_Restart::Bad_allocate_=1;        // bad_allocate is an int that tells us if the allocation has already taken place
int Save_Restart::Nb_pb_total_=0;        // nb_pb_total is the total number of problems
int Save_Restart::Num_pb_=1;                // num_pb is the number of the current problem

void Save_Restart::assoscier_pb_base(const Probleme_base& pb)
{
  pb_base_ = pb;
}

void Save_Restart::allocation() const
{
  if(pb_base_->schema_temps().file_allocation() && EcritureLectureSpecial::Active)        // Permet de tester l'allocation d'espace disque
    {
      if (Bad_allocate_==1)                                        // Si l'allocation n'a pas eut lieu
        if (Process::je_suis_maitre())                                // Qu'avec le proc maitre
          {
            if (Num_pb_==1)                                                // Si le probleme est le premier
              if (!allocate_file_size(File_size_))                        // je tente une allocation d'espace disque de taille 2*file_size
                Bad_allocate_=0;                                        // Si cela echoue, j'indique au code que l'allocation a deja eut lieu et n'a pas fonctionner
              else
                Num_pb_=Nb_pb_total_;                                        // Si OK, je modifie num_pb pour que les autres pb ne tentent pas d'allocation
            else
              Num_pb_-=1;                                                // Si le probleme n'est pas le premier, je decremente le numero de probleme
          }
      const int canal = 2007;
      if (Process::je_suis_maitre())                                // le processeur maitre envoi bad_allocate a tout le monde
        for (int p=1; p<Process::nproc(); p++)
          envoyer(Bad_allocate_,p,canal);
      else
        recevoir(Bad_allocate_,0,canal);

      if (Bad_allocate_==0)                                        // Si l'allocation a echoue
        {
          sauver_xyz(1);
          if (Num_pb_==Nb_pb_total_)                                        // Si le numero de probleme correspond au nombre total de probleme
            {
              if (Process::je_suis_maitre())
                {
                  Cerr << finl;                                                // j'arrete le code de facon claire
                  Cerr << "***Error*** " << error_ << finl;                // et je sort l'erreur du code
                  Cerr << "A xyz backup was made because you do not have enough disk space" << finl;
                  Cerr << "to continue the current calculation. Free up disk space and" << finl;
                  Cerr << "restart the calculation thanks to the backup just made." << finl;
                  Cerr << finl;
                }
              Process::barrier();
              Process::exit();
            }
          Num_pb_+=1;                                                // I increment the problem number
        }
    }
}

/*! @brief Verifies that the necessary space exists on the hard disk.
 *
 * @param the required disk space
 * @return (int) returns 1 if disk space is sufficient, 0 otherwise
 */
int Save_Restart::allocate_file_size(long int& size) const
{
#ifndef MICROSOFT
#ifndef __APPLE__
#ifndef RS6000
#ifdef CHECK_ALLOCATE
  Nom Fichier_File_size(Objet_U::nom_du_cas());
  Fichier_File_size+="_File_size";
  const char *file = Fichier_File_size;                        // Allocation file
  //  if (size<1048576)                                        // If size is too small we set it to 1 MB
  //     size=1048576;
  off_t taille = off_t(size+size);                        // Conversion of file size 2*size

  int fichier = open(file, O_WRONLY | O_CREAT, 0666);        // Opening of File_size file
  if (fichier == -1)                                        // Opening error
    {
      error_="Open of ";
      error_+=file;
      error_+=" : ";
      error_+=strerror(errno);                                // Error on opening
      close(fichier);                                        // file closing
      remove(file);                                        // Destruction of File_size file
      return 0;                                                // Allocation failure because file not opened
    }

  if (posix_fallocate(fichier, 0, taille) != 0)                // Disk space allocation error
    {
      error_="Allocation of ";
      error_+=file;
      error_+=" : ";
      error_+=strerror(errno);                                // Error on allocation
      close(fichier);                                        // file closing
      remove(file);                                        // Destruction of File_size file
      return 0;                                                // Allocation failure because not enough space
    }
  close(fichier);                                        // file closing
  remove(file);                                                // Destruction of File_size file
#endif
#endif
#endif
#endif
  return 1;
}

void Save_Restart::preparer_calcul()
{
#ifndef RS6000
  if (pb_base_->schema_temps().file_allocation() && EcritureLectureSpecial::Active)
    {
      Nom nom_fich_xyz(".xyz");
      sauver_xyz(0);
      if (Process::je_suis_maitre())
        {
          ifstream fichier(nom_fich_xyz); // Calculation of disk space taken by XYZ file of current problem
          fichier.seekg(0, std::ios_base::end);
          File_size_ += fichier.tellg(); // Increments the disk space already necessary
          fichier.close();
          remove(nom_fich_xyz);
        }
      Nb_pb_total_ += 1; // Allows knowing the total number of problems at the end of preparer_calcul
    }
#endif

  for(int i=0; i<pb_base_->nombre_d_equations(); i++)
    pb_base_->equation(i).init_save_file();
}

void Save_Restart::setTinitFromLastTime(double last_time)
{
  // Set the time to restart the calculation
  pb_base_->schema_temps().set_temps_courant() = last_time;
  // Initialize tinit and current time according last_time
  if (pb_base_->schema_temps().temps_init() > -DMAXFLOAT)
    {
      Cerr << "tinit was defined in .data file to " << pb_base_->schema_temps().temps_init() << ". The value is fixed to " << last_time << " accroding to resume_last_time_option" << finl;
    }
  pb_base_->schema_temps().set_temps_init() = last_time;
  pb_base_->schema_temps().set_temps_precedent() = last_time;
  Cerr << "==================================================================================================" << finl;
  Cerr << "In the backup file, we find the last time: " << last_time << " and read the fields." << finl;
}

void Save_Restart::checkVersion(const Nom& nomfic)
{
  if (Process::mp_min(restart_version_) != Process::mp_max(restart_version_))
    {
      Cerr << "The version of the format backup/resumption is not the same in the resumption files " << nomfic << finl;
      Process::exit();
    }
  if (restart_version_ > version_format_sauvegarde())
    {
      Cerr << "The format " << restart_version_ << " of the resumption file " << nomfic << " is posterior" << finl;
      Cerr << "to the format " << version_format_sauvegarde() << " recognized by this version of TRUST." << finl;
      Cerr << "Please use a more recent version." << finl;
      Process::exit();
    }

  // Writing of restart format
  Cerr << "The version of the resumption format of file " << nomfic << " is " << restart_version_ << finl;
}

void Save_Restart::prepare_PDI_restart(int resume_last_time)
{
  TRUST_2_PDI::set_PDI_restart(1);
  TRUST_2_PDI pdi_interface;

  int last_iteration = -1;
  double tinit = -1.;

  // Restart from the last time
  if (resume_last_time)
    {
      // Look for the last time saved in checkpoint file to init current computation
      pdi_interface.prepareRestart(restartComm_, last_iteration, tinit, 1 /*resume_last_time */);

      // set last time found in checkpoint file to tinit if tinit not set
      setTinitFromLastTime(tinit);
    }
  else // resume from the requested time
    {
      // looking for tinit in backup file
      tinit = pb_base_->schema_temps().temps_init();
      pdi_interface.prepareRestart(restartComm_, last_iteration, tinit, 0 /* reprise */);
    }
}

void Save_Restart::sauver_xyz(int verbose) const
{
  statistics().begin_count(STD_COUNTERS::backup_file,statistics().get_last_opened_counter_level()+1);
  Nom nom_fich_xyz("");
  if (verbose)
    {
      nom_fich_xyz += Objet_U::nom_du_cas();
      nom_fich_xyz += "_";
      nom_fich_xyz += pb_base_->le_nom();
      nom_fich_xyz += ".xyz";
      Cerr << "Creation of " << nom_fich_xyz << " (" << EcritureLectureSpecial::get_Output() << ") for resumption of a calculation with a different number of processors." << finl;
    }
  else
    {
      nom_fich_xyz = ".xyz";
    }
  // Create the XYZ file for the current problem
  ficsauv_.typer(EcritureLectureSpecial::get_Output());
  ficsauv_->ouvrir(nom_fich_xyz);
  // New for xyz since version 155: the backup format is written in the header
  if (Process::je_suis_maitre())
    ficsauv_.valeur() << "format_sauvegarde:" << finl << version_format_sauvegarde() << finl;

  EcritureLectureSpecial::mode_ecr = 1;
  int bytes = pb_base_->sauvegarder(ficsauv_.valeur());
  EcritureLectureSpecial::mode_ecr = -1;

  if (Process::je_suis_maitre())
    ficsauv_.valeur() << Nom("fin");
  (ficsauv_.valeur()).flush();
  (ficsauv_.valeur()).syncfile();
  ficsauv_.detach();
  Cout << "[IO] " << statistics().get_time_since_last_open(STD_COUNTERS::backup_file) << " s to write xyz file." << finl;
  statistics().end_count(STD_COUNTERS::backup_file,1,bytes);
}


void Save_Restart::lire_pdi_sauvegarde_reprise(Entree& is, Motcle& motlu, Nom& restart_file_name, Nom& yaml_fname)
{
  Nom nom;
  is >> nom;
  motlu = nom;
  if(motlu==Motcle("{"))
    {
      Motcles compris(3);
      compris[0]="}";
      compris[1]="checkpoint_fname";
      compris[2]="yaml_fname";
      int ind = -1;
      while (ind!=0)
        {
          is >> motlu;
          ind = compris.rang(motlu);
          if (ind==1)
            is >> restart_file_name;
          else if (ind==2)
            {
              Cerr << "[Save_Restart] lire_pdi_sauvegarde_reprise :: You have provided your own yaml file to initialize PDI ! " << finl;
              is >> yaml_fname;

              // Check to see if the file exists
              LecFicDiffuse test;
              if (!test.ouvrir(yaml_fname))
                {
                  Cerr << "[Save_Restart] lire_pdi_sauvegarde_reprise :: Error! The provided file " << yaml_fname << " does not exist " << finl;
                  Process::exit();
                }
            }
          else if (ind==-1)
            {
              Cerr << "[Save_Restart] lire_pdi_sauvegarde_reprise :: " << motlu << " is not understood. Keywords are:" << finl;
              Cerr << compris << finl;
              Process::exit();
            }
        }
    }
  else
    {
      Cerr << "[Save_Restart] lire_pdi_sauvegarde_reprise :: " << motlu << " is not understood. Expected { :" << finl;
      Process::exit();
    }
}


/////////////////////////////////////////////
// Reading restart options for a computation
/////////////////////////////////////////////
void Save_Restart::lire_reprise(Entree& is, Motcle& motlu)
{
  int resume_last_time = (motlu == "resume_last_time" ? 1 : 0);
  // reset to zero to allow a standard restart after an xyz restart
  EcritureLectureSpecial::mode_lec = 0;
  Motcle format_rep;
  is >> format_rep;
  if ((format_rep != "formatte") && (format_rep != "binaire") && (format_rep != "xyz") && (format_rep != "single_hdf") && (format_rep != "pdi") && (format_rep != "pdi_expert"))
    {
      Cerr << "Restarting calculation... : keyword " << format_rep << " not understood. Waiting for:" << finl << motlu << " formatte|binaire|xyz|single_hdf|pdi|pdi_expert Filename" << finl;
      Process::exit();
    }

  // XXX Elie Saikali : for PolyMAC_CDO => only .sauv files are possible
  if (pb_base_->discretisation().is_PolyMAC_MPFA() && format_rep != "binaire")
    {
      Cerr << "Error in Save_Restart::" << __func__ << " !! " << finl;
      Cerr << "Only the binary format is currently supported to resume a simulation with the discretization " << pb_base_->discretisation().que_suis_je() << " ! " << finl;
      Cerr << "Please update your data file and use a .sauv file !" << finl;
      Process::exit();
    }

  // Read the filename:
  Nom nom_yaml;
  if( format_rep == "pdi_expert" )
    {
      lire_pdi_sauvegarde_reprise(is, motlu, restart_filename_, nom_yaml);
      format_rep = "pdi";
    }
  else
    is >> restart_filename_;
  // Force hdf restart beyond a certain number of MPI ranks:
  if (format_rep != "xyz" && Process::force_single_file(Process::nproc(), restart_filename_))
    format_rep = "pdi";

  if(format_rep == "pdi")
    {
      std::string yaml_fname = nom_yaml.getString();
      if(yaml_fname == "??")
        {
          Ecrire_YAML yaml_file;
          yaml_file.add_pb_base(pb_base_, restart_filename_);
          yaml_fname = "restart_" + pb_base_->le_nom().getString() + ".yml";
          yaml_file.write_restart_file(yaml_fname);
        }
      TRUST_2_PDI::init(yaml_fname);

      // Prepare restart
      prepare_PDI_restart(resume_last_time);

      Entree useless;
      pb_base_->reprendre(useless);

      TRUST_2_PDI::finalize();
    }
  else if(format_rep == "single_hdf")
    {
      // !! DEPRECATED HDF5 FILE !!
      Cerr << "==============================================================================" << finl;
      Cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
      Cerr << "WARNING::you are using a deprecated backup file format. Please switch to PDI." << finl;
      Cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
      Cerr << "==============================================================================" << finl;
      LecFicDiffuse test;
      if (!test.ouvrir(restart_filename_))
        {
          Cerr << "Error! " << restart_filename_ << " file not found ! " << finl;
          Process::exit();
        }
      FichierHDFPar fic_hdf;
      Entree_Brute input_data;
      fic_hdf.open(restart_filename_, true);
      fic_hdf.read_dataset("/sauv", Process::me(),input_data);

      if(resume_last_time)
        {
          double last_time = -1;
          last_time = get_last_time(input_data);
          setTinitFromLastTime(last_time);
          fic_hdf.read_dataset("/sauv", Process::me(), input_data);
        }

      input_data >> motlu;
      if (motlu=="format_sauvegarde:")
        {
          input_data >> restart_version_;
          checkVersion(restart_filename_);
        }
      else
        {
          Cerr<<"This .sauv file is too old and the format is not supported anymore."<<finl;
          Process::exit();
        }
      fic_hdf.close();

      // Restart computation from checkpoint file
      pb_base_->reprendre(input_data);
    }
  else
    {
      OWN_PTR(Entree_Fichier_base) fic;
      if (format_rep == "formatte")
        fic.typer("LecFicDistribue");
      else if (format_rep == "binaire")
        fic.typer("LecFicDistribueBin");
      else if (format_rep == "xyz")
        {
          EcritureLectureSpecial::mode_lec = 1;
          fic.typer(EcritureLectureSpecial::Input);
        }
      fic->ouvrir(restart_filename_);
      if (fic->fail())
        {
          Cerr << "Error during the opening of the restart file : " << restart_filename_ << finl;
          Process::exit();
        }

      // Restart from the last time
      if (resume_last_time)
        {
          // Look for the last time and set it to tinit if tinit not set
          double last_time = -1.;
          last_time = get_last_time(fic);
          setTinitFromLastTime(last_time);

          fic->close();
          fic->ouvrir(restart_filename_);
        }

      // Read the backup format version if this is a standard restart
      // Since 1.5.1, the backup format is marked at the header of backup files
      // to allow easier evolution of the format in the future.
      // Moreover with 1.5.1, faces are numbered differently, so restarting
      // from an older backup file is incorrect; this is a way to warn users:
      // they must perform an xyz restart to continue a computation started with an older version.
      // Since 1.5.5, there is no format version for xyz
      fic.valeur() >> motlu;
      if (motlu != "FORMAT_SAUVEGARDE:")
        {
          if (format_rep == "xyz")
            {
              // We close and re-open the file:
              fic->close();
              fic->ouvrir(restart_filename_);
              restart_version_ = 151;
            }
          else
            {
              Cerr << "-------------------------------------------------------------------------------------" << finl;
              Cerr << "The resumption file " << restart_filename_ << " can not be read by this version of TRUST" << finl;
              Cerr << "which is a later version than 1.5. Indeed, the numbering of the faces have changed" << finl;
              Cerr << "and it would produce an erroneous resumption. If you want to use this version," << finl;
              Cerr << "you must do a resumption of the file .xyz saved during the previous calculation" << finl;
              Cerr << "because this file is independent of the numbering of the faces." << finl;
              Cerr << "The next backup will be made in a format compatible with the new" << finl;
              Cerr << "numbering of the faces and you can then redo classical resumptions." << finl;
              Cerr << "-------------------------------------------------------------------------------------" << finl;
              Process::exit();
            }
        }
      else
        {
          fic.valeur() >> restart_version_;
          checkVersion(restart_filename_);
        }

      // Restart computation from checkpoint file
      pb_base_->reprendre(fic.valeur());
    }

  restart_done_ = true;
  restart_in_progress_ = true;
}

////////////////////////////////////////////////
// Reading save options for a computation
////////////////////////////////////////////////
void Save_Restart::lire_sauvegarde(Entree& is, Motcle& motlu)
{
  // restart_file=1: the file is overwritten at each backup (and thus contains only one instant)
  if (motlu == "sauvegarde_simple")
    simple_restart_ = true;
  is >> checkpoint_format_;
  if ((Motcle(checkpoint_format_) != "binaire") && (Motcle(checkpoint_format_) != "formatte") &&
      (Motcle(checkpoint_format_) != "xyz") && (Motcle(checkpoint_format_) != "single_hdf") &&
      (Motcle(checkpoint_format_) != "pdi") && (Motcle(checkpoint_format_) != "pdi_expert")  )
    {
      checkpoint_filename_ = checkpoint_format_;
      checkpoint_format_ = "binaire";
    }
  else
    {
      if( Motcle(checkpoint_format_) == "pdi_expert" )
        {
          lire_pdi_sauvegarde_reprise(is, motlu, checkpoint_filename_, yaml_fname_);
          checkpoint_format_ = "pdi";
        }
      else
        {
          if( Motcle(checkpoint_format_) == "single_hdf" )
            {
              Cerr << "==============================================================================" << finl;
              Cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
              Cerr << "WARNING::you are using a deprecated backup file format. Please switch to PDI." << finl;
              Cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
              Cerr << "==============================================================================" << finl;
            }
          is >> checkpoint_filename_;
        }
    }
}

void Save_Restart::lire_sauvegarde_reprise(Entree& is, Motcle& motlu)
{
  // XXX Elie Saikali : for PolyMAC_MPFA => No xyz for the moment
  if (pb_base_->discretisation().is_PolyMAC_MPFA())
    {
      Cerr << "Problem "  << pb_base_->le_nom() << " with the discretization "
           << pb_base_->discretisation().que_suis_je() <<  " => EcritureLectureSpecial = 0 !" << finl;
      EcritureLectureSpecial::Active = 0;
    }
  checkpoint_format_ = "binaire";
  checkpoint_filename_ = Objet_U::nom_du_cas();
  checkpoint_filename_ += "_";
  checkpoint_filename_ += pb_base_->le_nom();
  checkpoint_filename_ += ".sauv";

  while (1)
    {
      if ((motlu == "reprise") || (motlu == "resume_last_time"))
        lire_reprise(is, motlu);
      else if (motlu == "sauvegarde" || motlu == "sauvegarde_simple")
        lire_sauvegarde(is, motlu);
      else if (motlu == "}")
        break;
      else
        {
          Cerr << "Error in Save_Restart::lire_sauvegarde_reprise" << finl;
          Cerr << "We expected } instead of " << motlu << " to mark the end of the data set" << finl;
          Process::exit();
        }
      is >> motlu;
    }

  ficsauv_.detach();
  // Force hdf backup beyond a certain number of MPI ranks:
  if (checkpoint_format_ != "xyz" && Process::force_single_file(Process::nproc(), checkpoint_filename_))
    checkpoint_format_ = "pdi";

  if ((Motcle(checkpoint_format_) != "binaire") && (Motcle(checkpoint_format_) != "formatte") &&
      (Motcle(checkpoint_format_) != "xyz") && (Motcle(checkpoint_format_) != "pdi") &&
      (Motcle(checkpoint_format_) != "single_hdf"))
    {
      Cerr << "Error of backup format ! We expected formatte, binaire, xyz, or pdi/pdi_expert (which replace single_hdf). single_hdf format is still available but is deprecated, please use pdi instead." << finl;
      Process::exit();
    }

  if (pb_base_->schema_temps().temps_init() <= -DMAXFLOAT)
    {
      pb_base_->schema_temps().set_temps_init() = 0;
      pb_base_->schema_temps().set_temps_courant() = 0;
    }

  if (reprise_effectuee())
    {
      // check if dt_ev exists, otherwise set reprise to 2
      // we will recreate the header in dt_ev otherwise the header is wrong on pb_couple restart
      Nom fichier(Objet_U::nom_du_cas());
      fichier += ".dt_ev";
      struct stat f;
      if (stat(fichier, &f))
        reprise_effectuee() = 2;
    }
}

/*! @brief Writes to file for restart (backup).
 *
 */
int Save_Restart::sauver() const
{
  int pdi_format = Motcle(checkpoint_format_) == "pdi";
  if(pdi_format)
    {
      if(!TRUST_2_PDI::is_PDI_initialized())
        {
          std::string yaml_fname = yaml_fname_.getString();
          if(yaml_fname == "??")
            {
              Ecrire_YAML yaml_file;
              yaml_fname = "save_" + pb_base_->le_nom().getString() + ".yml";
              yaml_file.add_pb_base(pb_base_, checkpoint_filename_);
              yaml_file.write_checkpoint_file(yaml_fname);
            }
          TRUST_2_PDI::init(yaml_fname);
        }

      if(!config_file_created_)
        {
          int nb_proc = Process::nproc();
          IntTab nodeRanks(nb_proc);
          nodeRanks = PE_Groups::get_node_group().get_node_id();
          envoyer_gather(nodeRanks, nodeRanks, 0);
          // Creating and filling the configuration checkpoint file
          // which contains all the information about the nodes partition used for checkpoint
          if (Process::je_suis_maitre())
            {
              TRUST_2_PDI pdi_interface;
              pdi_interface.write("nb_proc", &nb_proc);
              pdi_interface.trigger("InitConfig");

              int nb_nodes = PE_Groups::get_node_group().get_number_of_nodes();
              pdi_interface.write("nb_nodes", &nb_nodes);
              pdi_interface.trigger("WriteConfig");
              pdi_interface.TRUST_start_sharing("nodeRanks", nodeRanks.data());
              pdi_interface.trigger("WriteNodeRanks");
              pdi_interface.stop_sharing_last_variable();
            }
          config_file_created_ = true;
        }

      // if we are dealing with a coupled problem, the initialization might have been done twice
      // in which case we don't want to overwrite the file
      if(Process::node_master() && !ficsauv_created_)
        {
          TRUST_2_PDI pdi_interface;
          // if a file with the same name already exists, delete it and create a new one
          int non_const_sr = simple_restart_;
          pdi_interface.TRUST_start_sharing("simple_sauvegarde", &non_const_sr);
          std::string event = "init_" + pb_base_->le_nom().getString();
          pdi_interface.trigger(event);
          pdi_interface.stop_sharing_last_variable();

          // format information
          int version = version_format_PDI();
          pdi_interface.write("version", &version);
          ficsauv_created_ = true;
        }
    }
  else if (!ficsauv_ && !osauv_hdf_)
    {
      // If the backup file has not been opened yet, create the backup file:
      if (Motcle(checkpoint_format_) == "formatte")
        {
          ficsauv_.typer("EcrFicCollecte");
          //
          // Even in 64b, a save/restart SAUV file never actually requires 64b indices since all the information
          // saved is per proc. So we might as well save some space (not so much actually, since most of the data
          // saved are double values).
          //
          ficsauv_->set_64b(false);
          ficsauv_->ouvrir(checkpoint_filename_);
          ficsauv_->setf(ios::scientific);
        }
      else if (Motcle(checkpoint_format_) == "binaire")
        {
          ficsauv_.typer("EcrFicCollecteBin");
          ficsauv_->set_64b(false); // see comment above!
          ficsauv_->ouvrir(checkpoint_filename_);
        }
      else if (Motcle(checkpoint_format_) == "xyz")
        {
          ficsauv_.typer(EcritureLectureSpecial::get_Output());
          ficsauv_->ouvrir(checkpoint_filename_);
        }
      else if (Motcle(checkpoint_format_) == "single_hdf")
        osauv_hdf_ = new Sortie_Brute;
      else
        {
          Cerr << "Error in Save_Restart::sauver() " << finl;
          Cerr << "The format for the backup file must be either binary or formatted or pdi/pdi_expert (which replace single_hdf). single_hdf is still available but is deprecated, please use pdi instead." << finl;
          Cerr << "But it is :" << checkpoint_format_ << finl;
          Process::exit();
        }
      // If this is the first backup, write the backup format in the header
      if (Motcle(checkpoint_format_) == "xyz")
        {
          if (Process::je_suis_maitre())
            ficsauv_.valeur() << "format_sauvegarde:" << finl << version_format_sauvegarde() << finl;
        }
      else if ((Motcle(checkpoint_format_) == "single_hdf"))
        *osauv_hdf_ << "format_sauvegarde:" << finl << version_format_sauvegarde() << finl;
      else
        ficsauv_.valeur() << "format_sauvegarde:" << finl << version_format_sauvegarde() << finl;
    }

  // Perform the backup writing
  int bytes;
  EcritureLectureSpecial::mode_ecr = (Motcle(checkpoint_format_) == "xyz");
  TRUST_2_PDI::set_PDI_checkpoint(pdi_format);
  if(pdi_format)
    {
      TRUST_2_PDI pdi_interface;
      int tmp = 1;
      pdi_interface.TRUST_start_sharing("TYPES", &tmp);

      Sortie_Nulle useless;
      bytes = pb_base_->sauvegarder(useless);

      // backup of the unknown fields (which are local to each processor so it will involve a parallel writing)
      std::string f_event = "local_backup_" + pb_base_->le_nom().getString();
      pdi_interface.trigger(f_event);
      if(Process::node_master())
        {
          // backup of the data that are global to everyone so we just need one proc (the master of the node) to write it
          std::string s_event = "global_backup_" + pb_base_->le_nom().getString();
          pdi_interface.trigger(s_event);

          // we save the types of fields we want to save (not all of them actually, only the ones that are necessary for restart)
          std::string t_event = pb_base_->le_nom().getString() + "_get_types";
          pdi_interface.trigger(t_event);
        }
      pdi_interface.stop_sharing();
    }
  else if (Motcle(checkpoint_format_) == "single_hdf")
    bytes = pb_base_->sauvegarder(*osauv_hdf_);
  else
    bytes = pb_base_->sauvegarder(ficsauv_.valeur());
  EcritureLectureSpecial::mode_ecr = -1;
  TRUST_2_PDI::set_PDI_checkpoint(0);

  // If this is a simple backup, close the file immediately and properly
  if (simple_restart_ && !pdi_format)
    {
      if (Motcle(checkpoint_format_) == "xyz")
        {
          if (Process::je_suis_maitre())
            ficsauv_.valeur() << Nom("fin");
          (ficsauv_.valeur()).flush();
          (ficsauv_.valeur()).syncfile();
        }
      else if (Motcle(checkpoint_format_) == "single_hdf")
        {
          *osauv_hdf_ << Nom("fin");
          FichierHDFPar fic_hdf;
          fic_hdf.create(checkpoint_filename_);
          fic_hdf.create_and_fill_dataset_MW("/sauv", *osauv_hdf_);
          fic_hdf.close();
          delete osauv_hdf_;
          osauv_hdf_ = 0;
        }
      else
        {
          ficsauv_.valeur() << Nom("fin");
          (ficsauv_.valeur()).flush();
        }
      ficsauv_.detach();
    }
  return bytes;
}

void Save_Restart::finir()
{
  // Close the backup file properly
  // If it is a simple backup, fin was written at each call to ::sauver()
  if(Motcle(checkpoint_format_) == "pdi" && TRUST_2_PDI::is_PDI_initialized())
    TRUST_2_PDI::finalize();
  else  if (!simple_restart_ && (ficsauv_ || osauv_hdf_) )
    {
      if (Motcle(checkpoint_format_) == "xyz")
        {
          if (Process::je_suis_maitre())
            ficsauv_.valeur() << Nom("fin");
          (ficsauv_.valeur()).flush();
          (ficsauv_.valeur()).syncfile();
        }
      else if (Motcle(checkpoint_format_) == "single_hdf")
        {
          *osauv_hdf_ << Nom("fin");
          FichierHDFPar fic_hdf;
          fic_hdf.create(checkpoint_filename_);
          fic_hdf.create_and_fill_dataset_MW("/sauv", *osauv_hdf_);
          fic_hdf.close();
          delete osauv_hdf_;
          osauv_hdf_ = 0;
        }
      else
        {
          ficsauv_.valeur() << Nom("fin");
          (ficsauv_.valeur()).flush();
        }

      ficsauv_.detach();
    }
  // If the backup is a standard one and the user has not disabled the final xyz backup,
  // then perform the final xyz backup
  if (Motcle(checkpoint_format_) != "xyz")
    {
      if (EcritureLectureSpecial::Active)
        sauver_xyz(1);
      else
        Cerr << "As saving .xyz file disabled since 1.9.7, add into your datafile \"EcritureLectureSpecial 1\" to enable it again if wanted." << finl;
    }

  for(int i=0; i<pb_base_->nombre_d_equations(); i++)
    pb_base_->equation(i).close_save_file();
}
