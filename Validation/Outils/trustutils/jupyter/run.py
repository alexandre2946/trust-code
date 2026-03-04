"""
Victor Banon Garcia / Adrien Bruneton / Teo Boutin
CEA Saclay - DM2S/STMF/LGLS
03/2021

We provide here a Python package that can be used to run validation files from jupyterlab.

IMPORTANT: we assume that the methods of this module are invoked from the root directory of a
validation form (i.e. a directory containing at least a 'src' subdirectory).

"""

import os
import subprocess
import argparse

from string import Template

from time import time, sleep, strftime 

def saveFormOutput():
    """ Dummy method to indicate that the output of the notebook should be saved.
    This method does nothing here, but its invokation is detected by the pre_save hook
    registered in the TRUST Jupyter configuration.
    By default, outputs of the validation forms are not saved.
    """
    pass

def getLastLines_(err_pth):
    """ Get last 20 lines of an error file ...
    """
    if os.path.exists(err_pth):
        with open(err_pth) as f:
            s = "(....)"
            txt = f.readlines()
            err = "".join(txt[-20:])
    else:
        err = "Error file %s can not be opened!" % err_pth
    return err


def displayMD(txt):
    """ Display text as markdown in Jupyter notebook.
    """
    try:
        from IPython.display import display, Markdown

        ok = True
    except:
        ok = False
    if ok:
        display(Markdown(txt))


def useMEDCoupling():
    """ Load MEDCoupling environment in the PYTHONPATH so that 'import medcoupling' can work.
    """
    import sys

    mcr = os.environ["TRUST_MEDCOUPLING_ROOT"]
    arch = os.environ["TRUST_ARCH"]
    sub = "%s_opt/lib/python%d.%d/site-packages" % (arch, sys.version_info.major, sys.version_info.minor)
    sys.path.append(os.path.join(mcr, sub))
    sys.path.append(os.path.join(mcr, "bin"))
    try:
        import medcoupling
    except:
        raise Exception("Could not load MEDCoupling environment!")
        
def useLataTools():
    """ Load LataTools environment in the PYTHONPATH so that 'import LataLoader' can work.
    """
    import sys
    useMEDCoupling()
    trust_root = os.environ["TRUST_ROOT"]
    sys.path.append(f"{trust_root}/exec/lata_tools/lib")
    # for testing. no need to try/catch, that would just shadow the initial exception
    import LataLoader
    
######## PRIVATE STUFF #########

def _print(*args, also_to_nb=False):
    """
    Private. When jupyter is run non-interactively, print to terminal rather than to the notebook.
    Useful for debugging in parallel runs.
    Usage is the same as python builtin print.
    Also takes kwarg also_to_nb=False/True, to also print to the notebook cell.
    
    Printing only to console and not to jupyter cell output is pretty hard.
    
    This is because of ipythonkernel PR 630 (https://github.com/ipython/ipykernel/pull/630)
    where the goal was to capture outputs of calls such as os.system(...).
    
    This was partially solved in PR https://github.com/ipython/ipykernel/pull/752
    but requires adding options to kernel used by jupyter, I do not know how to do that dynamically
    (requires editing some user config files afaik)
    
    Solution used here is drastic, as it completely disables the feature of PR 630
    It comes from https://github.com/ipython/ipykernel/issues/795
    
    To change when a better solution is added to jupyter (or add it yourself if you are motivated ;)
    
    """
    import sys
    
    if isExtractingNR():
        return
    
    for std, __std__ in [
        (sys.stdout, sys.__stdout__),
        (sys.stderr, sys.__stderr__),
    ]:
        if getattr(std, "_original_stdstream_copy", None) is not None:
            # redirect captured pipe back to original FD
            os.dup2(std._original_stdstream_copy, __std__.fileno())
            std._original_stdstream_copy = None
            
            
    terminal = sys.__stdout__

        
    print(*args, file=terminal, flush=True)

    if also_to_nb:
        print(*args)

def _initBuildDir():
    """
    Private. Initialize build directory: this can be overriden with the "-dest" option in env variable JUPYTER_OPTIONS
    Change the build directory if -dest option provided. This is mostly for NR test and validation
    """
    ret = os.path.join(ORIGIN_DIRECTORY,"build")
    # Compute correct build directory if overriden:
    opt = os.environ.get("JUPYTER_RUN_OPTIONS", "")
    a = opt.split(" ")
    if '' in a:
        a=a.remove('') or []
    if "-dest" in a:
        idx = a.index("-dest")
        if idx >= 0 and len(a) >= idx + 2:
            ret = os.path.join(a[idx + 1], "build")
    
    
    # Choose run option
    # can be parallel, sequential or not_run
    # default is sequential
    _set_run_sequential()
    
    # if Sserver detected, use it
    if ("-parallel_sjob" in opt) or _detectSserver():
        _set_use_sserver()
        _set_run_parallel()
    
    if ("-not_run" in opt):
        _set_not_run()
        
    if ("-parallel_run" in opt):
        _set_run_parallel()
        

        

    # not used yet
    # ~ parser=argparse.ArgumentParser()
    # ~ parser.add_argument("-dest")
    # ~ parser.add_argument("-parallel_sjob", action="store_true")
    # ~ parser.add_argument("-parallel_run", action="store_true")
    # ~ parser.add_argument("-export_pdf", action="store_true")
    # ~ parser.add_argument("-timeout")  # TODO for timeout values
    # ~ parser.print_help()
    # ~ if a and len(a)>0:
        # ~ JUPYTER_RUN_OPTIONS=parser.parse_args(a)
        
    return ret


def _runCommand(cmd, verbose):
    """ Private method to run a command, print the logs if verbose or if it fails, and throw an exception
    if it fails.
    """
    # Run by redirecting stderr to stdout
    # add universal_newlines as in https://stackoverflow.com/questions/41171791/how-to-suppress-or-capture-the-output-of-subprocess-run
    complProc = subprocess.run(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.STDOUT,universal_newlines = True)
    if verbose:
        print(cmd)
        print(complProc.stdout)
    # Throw if return code non-zero:
    if complProc.returncode != 0: # and complProc.stdout:
        # Display message through a custom exception so that jupyter-nbconvert also shows it properly in the console:
        msg = "\nExecution of following command failed!!\n"
        msg += "  " + cmd
        msg += "\nwith return code %d\n" % complProc.returncode
        msg += "and with following output:\n\n"
        msg += (complProc.stdout if complProc.stdout else "")
        _print(msg, also_to_nb=True)
        raise RuntimeError(msg)

def _detectSserver():
    """ Detect whether the Sserver is running on the machine
    """
    squeue = os.path.join(os.environ["TRUST_ROOT"], "bin", "Sjob", "Squeue")
    try:
        subprocess.check_call([squeue], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False

######## End PRIVATE STUFF #########

class TRUSTCase(object):
    """ 
    Class which allows the user to load and execute a validation case
    """

    _UNIQ_ID_START = -1

    def __init__(self, directory, datasetName, nbProcs=1, execOptions="", excluNR=False, pre_run=None, post_run=None):
        """ 
        Initialisation of the class

        Parameters
        ---------

        directory: str 
            Path of the file
        datasetName: str
            Path of the case we want to run, relatively to the 'src' folder.
        nbProcs : int
            Number of processors to use to run the case.
        execOptions : str
            TRUST Options to add at the execution of the test case 
        excluNR : bool 
            to remove of the non regression test cases

        Returns
        ------

        """
        TRUSTCase._UNIQ_ID_START += 1
        self.id_ = TRUSTCase._UNIQ_ID_START  # Unique ID
        self.dir_ = os.path.normpath(directory)
        self.name_ = datasetName.split(".data")[0]  # always w/o the trailing .data
        self.dataFileName_ = self.name_ #Change in case of parallel case
        self.nbProcs_ = nbProcs
        self.last_run_ok_ = -255  # exit status of the last run of the case
        self.last_run_err_ = ""  # error message returned when last running the case
        self.execOptions = execOptions
        self.excluNR = excluNR
        
        self.has_python_pre_run=False
        if pre_run is not None:
            if not callable(pre_run):
                raise ValueError("pre_run should be a function")
            self.has_python_pre_run=True
            self.pre_run=pre_run
        
        self.has_python_post_run=False
        if post_run is not None:
            if not callable(post_run):
                raise ValueError("post_run should be a function")
            self.has_python_post_run=True
            self.post_run=post_run

    def _fullDir(self):
        """
        full directory of the test case
        """
        return os.path.normpath(os.path.join(BUILD_DIRECTORY, self.dir_))

    def _fullPath(self):
        """
        full path of the test case in the build directory
        """
        fullPath = os.path.join(self._fullDir(), self.dataFileName_)
        return fullPath + ".data"

    def _relPath(self):
        """
        relative path of the test case in the build directory
        """
        fullPath = os.path.join(self._fullDir(), self.dataFileName_)
        relPath = os.path.relpath(fullPath, start = BUILD_DIRECTORY)
        return relPath + ".data"

    def _fullPath_OutFile(self):
        """
        full path of the test case in the build directory
        """
        fullPath = os.path.join(self._fullDir(), self.dataFileName_)
        return fullPath + ".out"

    def _fullPath_ErrFile(self):
        """
        full path of the test case in the build directory
        """
        fullPath = os.path.join(self._fullDir(), self.dataFileName_)
        return fullPath + ".err"

    def substitute(self, find, replace):
        """ 
        Substitute (in place and in the build directory) a part of the dataset

        Parameters
        ---------

        find: str 
            Text we want to substitute.
        replace: str
            Text to insert in replacement.

        """
        path = self._fullPath()
        text_out = ""
        with open(path, "r") as f:
            text_in = f.read()
            text_out += text_in.replace(find, replace)

        with open(path, "w") as f:
            f.write(text_out)

    def substitute_template(self, subs_dict):
        """ 
        Substitute (in place and in the build directory) a part of the dataset

        Parameters
        --------

        subs_dict: dict 
            Text we want to substitute with python Template formalism (character identified with $ in datafile)
        """
        path = self._fullPath()
        with open(path, "r") as file: filedata = Template(file.read())
        try:
            result = filedata.substitute(subs_dict)
        except KeyError as e:
            _print(f"TRUSTCase: in template substitution of case {self.dir_}/{self.name_}.data", also_to_nb = True)
            _print("Missing a keyword in dict given:", e, also_to_nb = True)
            pattern="$"+str(e).replace("'", "")
            index = filedata.template.find(pattern)
            ln = 1+filedata.template[0:index].count("\n")
            _print("At line", ln, also_to_nb = True)
            lt = filedata.template.split("\n")[ln-1]
            raise Exception(f"Missing keyword {e} in dict given to case template {self.dir_}/{self.name_}.data.\nOccurs in dataset at line {ln}:\n{lt}")
        except ValueError as e:
            raise Exception(f"Error in template substitution for case {self.dir_}/{self.name_}.data:\n{e}")
            
        with open(path, "w") as file: file.write(result)

    def copy(self, targetName, targetDirectory=None, nbProcs=1, execOptions="", excluNR=False, pre_run=None, post_run=None):
        """ 
            Copy a TRUST Case. 
            Only copies the datafile and not other optionnal atributes, such as excluNR, nbProcs, pre and post_run...
            These must be given again to the copy method
            

        Parameters
        ---------

        targetName: str
            New name
        targetDirectory: str
            New directory
        nbProcs: int
            number of procs
        execOptions : str
            TRUST Options to add at the execution of the test case 
        excluNR : bool 
            to remove of the non regression test cases

        """
        if targetDirectory is None:
            targetDirectory = self.dir_
        # Create the directory if not there:
        fullDir2 = os.path.join(BUILD_DIRECTORY, targetDirectory)
        if not os.path.exists(fullDir2):
            os.makedirs(fullDir2, exist_ok=True)
        pthTgt = os.path.join(BUILD_DIRECTORY, targetDirectory, targetName)
        # And copy the .data file:TRUSTSuite
        from shutil import copyfile

        copyfile(self._fullPath(), pthTgt)

        return TRUSTCase(targetDirectory, targetName, nbProcs=nbProcs, execOptions=execOptions, excluNR=False, pre_run=pre_run, post_run=post_run)

    def dumpDataset(self, user_keywords=[]):
        """ 
        Print out the .data file. TODO handle upper / lower case

        Parameters
        --------
        
        list_Trust_user_words: list(str) 
            List of word the user wants to color in red.
        """
        ###############################################################
        # Voici un lien qui explique comment on manipule les couleurs #
        # https://stackabuse.com/how-to-print-colored-text-in-python/ #
        ###############################################################

        ### Class that I found interesting to change theirs colors. Victor ###
        list_Trust_classe = [
            "fields",
            "Partition",
            "Read",
            "Conduction",
            "solveur implicite",
            "solveur gmres",
            "Post_processing",
            "diffusion",
            "initial_conditions",
            "boundary_conditions",
            "Probes",
        ]

        ### Import and change de layout of TRUST Keyword ###
        keywords = os.path.join(os.getenv("TRUST_ROOT"), "doc", "TRUST", "Keywords.txt")
        f = open(keywords, "r")
        tmp = []
        for i in f.readlines():
            i = i.replace("|", " ")
            i = i.replace("\n", " ")
            if len(i) == 2:
                i = i.replace(i, "")
            else:
                tmp.append(i)
        list_Trust_keywords = tmp
        f.close()

        ### Import and Underline important words of the data file  ###
        f = open(self._fullPath(), "r")
        tmp = f.readlines()
        test = []
        comment_flag = False
        for j in tmp:
            flag = 0
            if j.count("#") == 1 and  not comment_flag:
                comment_flag = True
                flag = 1
                j = "\033[38;5;88m" + j + "\033[0;0m"
            elif j.count("#") == 1 and  comment_flag :
                comment_flag = False
                j = "\033[38;5;88m" + j + "\033[0;0m"
                flag = 1
            elif j.count("#") == 2 or comment_flag:
                j = "\033[38;5;88m" + j + "\033[0;0m"
                flag = 1


            if "end" in j.lower():
                flag = 1
                j = "\033[38;5;88m" + j + "\033[0;0m"
            j = " " + j
            j = j.replace("\t", "\t ")
            if flag == 0:
                j = j.replace("{", "\033[196;5;30m{\033[0;0m") # en gras
                j = j.replace("}", "\033[196;5;30m}\033[0;0m") # en gras

                for i in list_Trust_keywords:
                    j = j.replace(i, "\033[38;5;28m" + i + "\033[0;0m")

                for i in list_Trust_classe:
                    j = j.replace(i, "\033[38;5;4m" + i + "\033[0;0m")

                for i in user_keywords:
                    # tmp=tmp.replace(i, "\033[196;5;30m"+i+"\033[0;0m ") # en gras
                    j = j.replace(i, "\033[38;5;196m" + i + "\033[0;0m")  # en rouge
            test.append(j)
        f.close()

        ### Print the coloured .data file ###
        print("".join(test))

    def dumpDatasetMD(self):
        """
        Print out the .data file in a basic MarkDown format.
        This will not highlight TRUST keywords.
        """

        f = open(self._fullPath(), "r").read()
        displayMD( "```\n" + f + "\n```" )

    def _runScript(self, scriptName, verbose=False):
        """ Internal. Run a shell script if it exists.
        """
        pth = "./%s" % scriptName
        if os.path.exists(pth):
            subprocess.check_output("chmod u+x %s" % pth, shell=True)
            cmd = "%s %s" % (pth, self.name_)
            _runCommand(cmd, verbose)

    def _preRun(self, verbose):
        if self.has_python_pre_run:
            case_path=os.path.join(BUILD_DIRECTORY, self.dir_)
            if os.path.exists("pre_run"):
                pre_run_path=os.path.join(case_path,'pre_run')
                raise Exception(f"Error in case {self.name_}: pre_run script exists when pre_run function was given to the test case. You should delete {pre_run_path}")
            self.pre_run(case_path, self.name_)
        else:
            self._runScript("pre_run", verbose)

    def _postRun(self, verbose):
        if self.has_python_post_run:
            if os.path.exists("post_run"):
                post_run_path=os.path.join(case_path,'post_run')
                raise Exception(f"Error in case {self.name_}: post_run script exists when pre_run function was given to the test case. You should delete {post_run_path}")
            case_path=os.path.join(BUILD_DIRECTORY, self.dir_)
            self.post_run(case_path, self.name_)
        else:
            self._runScript("post_run", verbose)

    def _generateExecScript(self):
        """ Generate a shell script doing the
            - pre_run
            - launching the case
            - and doing the post_run
        """
        uniq_id = "{:04d}".format(self.id_)
        
        # script names has additionnal details about the report and case names
        # this is useful to check the state of a running validation via Sserver
        scriptName=".cmds_%s_%s.sh" % (uniq_id, self.name_)
        scriptFl = os.path.join(BUILD_DIRECTORY, scriptName)

        logName = scriptName.replace(".sh", ".log")
        n = self.dataFileName_
        d = self.dir_
        e = self.execOptions
        out_file = self._fullPath_OutFile()
        err_file = self._fullPath_ErrFile()
        fullD, fullL = os.path.join(BUILD_DIRECTORY, d), os.path.join(BUILD_DIRECTORY, logName)
        para = ""
        if self.nbProcs_ != 1:
            para = str(self.nbProcs_)
        with open(scriptFl, "w") as f:
            s = "#!/bin/bash\n"
            s += "( echo;\n"
            s += '  echo "-> Running the calculation of the %s data file in the %s directory ...";\n' % (n, d)
            s += "  cd %s ; \n" % fullD
            
            # Running and checking pre_run was OK:
            if not self.has_python_pre_run:
                s += "  if [ -f pre_run ]; then\n"
                s += "     chmod +x pre_run\n"
                s +=f'     echo "-> Running the pre_run script in the {d} directory ..."\n'
                s +=f"     (./pre_run {n} || (echo '  FAILED pre_run!' && exit -1)) || exit -1 \n"
                s += "fi\n"

            # Running case
            s += f"  trust {n} {para} {e} 1>{out_file} 2>{err_file};\n"
            s += "  if [ ! $? -eq 0 ]; then exit 1; fi; \n"

            # Running and checking post_run was OK:
            if not self.has_python_post_run:
                s += "  if [ -f post_run ]; then\n"
                s += "     chmod +x post_run\n"
                s +=f'     echo "-> Running the post_run script in the {d} directory ..."\n'
                s +=f"     (./post_run {n} || (echo '  FAILED post_run!' && exit -1)) || exit -1\n"
                s += "fi\n"

            s += "  exit 0;"
            s += ") 1>%s 2>&1 \n" % fullL
            f.write(s)
        os.chmod(scriptFl, 0o755)
        
        baseName = os.path.join(d, n)

        saveFileAccumulator(f"{baseName}.dt_ev")
        saveFileAccumulator(f"{baseName}.newton_evol")
        saveFileAccumulator(f"{baseName}.data")
        saveFileAccumulator(f"{baseName}.out")
        saveFileAccumulator(f"{baseName}_*.out")
        saveFileAccumulator(f"{baseName}.err")
        saveFileAccumulator(f"{baseName}_*.son")
        saveFileAccumulator(f"{baseName}.TU")
        saveFileAccumulator(f"{baseName}_csv.TU")
        
        return scriptFl, fullL

    def partition(self, verbose=False, overwritePartition=True):
        """ 
        Apply partitioning of specified test case with trust -partition if nbProcs>1 only
        This avoids to modify trust -partition to be able to call it with 1 cpu

        Parameters
        ---------
        verbose: bool
        overwritePartition: bool
            Whether to overwrite the partition from the dataset when running trust -partition
            If True, uses 
                trust -partition xxx.data nbProcs
            which will overwrite (or create) the partition in the dataset, but will cut the domain in a single direction.
            This is the old and default behavior.
            If False, uses the partition from the dataset by using
                trust -partition xxx.data
            without specifying nbProcs
            Warning: this will not crash if there is no commented partition block in the dataset
        
        """
        ok = True
        path = os.getcwd()
        os.chdir(self._fullDir())
        self.dataFileName_ = "PAR_"+self.name_

        # when getting list of test cases, no need for partition. But we need to get inside this partition()
        # method to have the correct dataset name (see above). so only skip from here
        if isExtractingNR_ListOnly():
            return
        opt = os.environ.get("JUPYTER_RUN_OPTIONS", "")
        if "-not_run" in opt:
            return

        if self.nbProcs_ == 1:
            err_msg = "Not allowed to call case.partition() on case = addCase(%s,..., nbProcs=1) \n" % (self.name_)
            err_msg += "You cannot run parallel computation on 1 proc!"
            raise ValueError(err_msg)

        err_file = self.name_ + "_partition.err"
        out_file = self.name_ + "_partition.out"
        if overwritePartition:
            # old default behavior
            cmd = "trust -partition %s %s 2>%s 1>%s" % (self.name_, str(self.nbProcs_), err_file, out_file)
        else:
            # possibility to keep user defined partition
            # Warning: no error if there is no partition in the dataset
            cmd = "trust -partition %s 2>%s 1>%s" % (self.name_, err_file, out_file)
        output = subprocess.run(cmd, shell=True, executable="/bin/bash", stderr=subprocess.STDOUT)
        if verbose:
            print(cmd)
        if output.returncode != 0:
            ok = False
            err = getLastLines_(err_file)

        os.chdir(path)

        baseName = os.path.join(self.dir_, self.dataFileName_)

        saveFileAccumulator(f"{baseName}.dt_ev")
        saveFileAccumulator(f"{baseName}.newton_evol")
        saveFileAccumulator(f"{baseName}.data")
        saveFileAccumulator(f"{baseName}.out")
        saveFileAccumulator(f"{baseName}_*.out")
        saveFileAccumulator(f"{baseName}.err")
        saveFileAccumulator(f"{baseName}_*.son")
        saveFileAccumulator(f"{baseName}.TU")
        saveFileAccumulator(f"{baseName}_csv.TU")

        return ok
        
    def run(self, verbose=False):
        """ 
        Execute the current test case.
        
        In parallel mode (default):
            - runs the pre_run, which may also run some TRUSTCases.
            - if that happens, add the case to a waiting list. It will be started by TRUSTSuite.runCases() when pre_run is finished.
            - otherwise, directly submit the job via subprocess.Popen(...) and store it in _RUNNING_CASES.
            - post_run is called when the case is finished by TRUSTSuite.runCases()
        
        In sequential mode (when preventConcurrent=True is passed to TRUSTSuite.runCases()):
            - runs the pre_run
            - wait for potential cases started by pre_run to finish
            - run the case and wait for it to finish
            - run the post_run and wait for potential cases started by post_run to finish
            
        Sserver usage is independent of parallel/sequential mode.
        Jobs are submitted to SSserver if it is detected, so that preventConcurrent does not overloads cpus during complete validation
        
        Parameters
        ---------
        verbose: bool
        
        The results of the run are stored in members self.last_run_ok_ and self.last_run_err_
        """
        
        # Very specific to the validation process. Sometimes we want the core
        # method 'runCases()' not to do anything ... see script 'archive_resultat' for example.
        if _NOT_RUN:
            return    
        if isExtractingNR_ListOnly():
            return
        
        (script, logFile,) = self._generateExecScript()  # Generate the shell script doing pre_run, case and post_run
        cmdLst = ["bash", script, "&"]
        
        if _USE_SSERVER:
            # Invoke Salloc to schedule test case execution:
            salloc = os.path.join(os.environ["TRUST_ROOT"], "bin", "Sjob", "Salloc")
            cmdLst = [salloc, "-n", str(self.nbProcs_), script]
            
        case_popen=lambda: subprocess.Popen(cmdLst, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        case_info={"case":self, "script":script, "logFile":logFile, "callback":lambda: None, "callbackDone": False, "popen": case_popen}
        
        if self.has_python_post_run:
            case_info["callback"]=lambda: self._postRun(verbose)
        
        if _RUN_PARALLEL:
            
            nj=len(_RUNNING_CASES)
            
            if self.has_python_pre_run:
                self._preRun(verbose)
                
            if len(_RUNNING_CASES)>nj:
                # that means pre_run started some jobs. We store the case in a waiting list with newest jobs as dependency
                lj=[a for a in range(nj, len(_RUNNING_CASES))] # the list of dependencies
                
                _print(f"pre_run of case {self._relPath()} started {len(lj)} jobs. Adding the case to waiting list.", also_to_nb=verbose)
                _print(f"Deps of case {self._relPath()} are {lj}.")
                
                _RUNNING_CASES.append({"process":None, **case_info, "depends":lj})
                
            else:
                # Launch the case. 
                
                # Wait for cpu availability if Sserver is not managing the jobs
                if not(_USE_SSERVER):
                    _wait_for_available_procs(self.nbProcs_)
                
                # Launch and store in _RUNNING_CASES
                p = case_popen()
                _print(f"Starting case {self._relPath()}", also_to_nb=verbose)
                _RUNNING_CASES.append({"process":p, **case_info, "depends":[]})
                
        elif _RUN_SEQUENTIAL:
            # case of sequential run
            # we add to the _RUNNING_CASES global list because to facilitate keeping track of failed runs
            if self.has_python_pre_run:
                self._preRun(verbose)
                
            # wait till all runs from pre_run are done
            for r in _RUNNING_CASES:
                r["process"].wait()
                
            p = case_popen()
            _print(f"Starting case {self._relPath()}", also_to_nb=verbose)
            _RUNNING_CASES.append({"process":p, **case_info})
            
            # wait till I am done
            p.wait()
            
            # do the post_run myself, otherwise it would only be done after all case ran
            _RUNNING_CASES[-1]["callbackDone"] = True
            if self.has_python_post_run:
                self._postRun(verbose)
                
            # wait till all runs from post_run are done
            for r in _RUNNING_CASES:
                r["process"].wait()
                
            
        else:
            raise Exception("Invalid run option, this should not be reachable.")
                
            
    
            
    def _addPerfToTable(self, zeTable):
        """ Extract performances for this case and add it to the global table
            passed in parameter.
        """
        os.chdir(self._fullDir())

        opt = os.environ.get("JUPYTER_RUN_OPTIONS", None)
        # Very specific to the validation process - we need to keep build there:
        if not opt is None and "-not_run" in opt:
            cmd = '.'
        else:
            cmd = os.environ["TRUST_ROOT"] + "/Validation/Outils/Genere_courbe/scripts/extract_perf " + self.dataFileName_
            _runCommand(cmd, False)

        f = open(self.dataFileName_ + ".perf", "r")
        row = f.readlines()[0].replace("\n", "").split(" ")[1:]
        f.close()

        if len(row) == 6:
            row_arrange = row[0:5]
            row_arrange[4] = str(row[4])
            if (row[5] != row[4]):
                row_arrange[4] += "-" + str(row[5])
            row = row_arrange

        zeTable.addLine([row[:5]], self.dir_ + "/" + self.dataFileName_)
        os.chdir(ORIGIN_DIRECTORY)

        ## Save the file
        saveFileAccumulator(self.dir_ + "/" + self.dataFileName_ + ".perf")


class TRUSTSuite(object):
    """ A set of TRUST cases to be run for the validation form.
    """
    def __init__(self, runPrepare=True):
        self.cases_ = []
        os.chdir(ORIGIN_DIRECTORY)
        self.copySrc()
        if runPrepare:
            executeScript("prepare")

    def copySrc(self):
        """ Copy content of src directory into build directory.
        WARNING: this is where we expect to be at the root of the validation form.
        """
        if not os.path.exists("src"):
            raise Exception("Not a coherent validation form directory: 'src' subdirectory not found.")
        # Mimick what is done in 'prepare_gen' TRUST script:
        subprocess.run("mkdir -p %s && cp -a src/* %s 2>/dev/null" % (BUILD_DIRECTORY, BUILD_DIRECTORY), shell=True)  # Note the '*' !!
        # Add image file that can be in src file for archiving
        os.chdir(BUILD_DIRECTORY)
        import pathlib
        desktop = pathlib.Path("./")
        image_files = list(desktop.rglob("*.png"))
        image_files.append(list(desktop.rglob("*.jpg")))
        image_files.append(list(desktop.rglob("*.jpeg")))
        for image in image_files:
            saveFileAccumulator(str(image))
        

    def addCase(self, case):
        self.cases_.append(case)

        baseName = os.path.join(case.dir_, case.name_)

        saveFileAccumulator(f"{baseName}.dt_ev")
        saveFileAccumulator(f"{baseName}.newton_evol")
        saveFileAccumulator(f"{baseName}.data")
        saveFileAccumulator(f"{baseName}.out")
        saveFileAccumulator(f"{baseName}_*.out")
        saveFileAccumulator(f"{baseName}.err")
        saveFileAccumulator(f"{baseName}_*.son")
        saveFileAccumulator(f"{baseName}.TU")
        saveFileAccumulator(f"{baseName}_csv.TU")

    def getCases(self):
        return self.cases_



    def runCases(self, verbose=False, preventConcurrent=False):
        """ Launch all cases for the current suite.
        
        Parameters
        ---------
        verbose: bool 
            whether to print the console output of the run.
        preventConcurrent: bool 
            run the cases in the order they were provided, even if the Sserver is up and running, and the -parallel_sjob option was passed.
        """
            
        # Very specific to the validation process. Sometimes we want the core
        # method 'runCases()' not to do anything ... see script 'archive_resultat' for example.
        if _NOT_RUN:
            return
            
        if preventConcurrent:
            _set_run_sequential()
            _print("runCases: Running sequential because preventConcurrent = True (Run_fiche -parallel_run will be ignored)", also_to_nb=verbose)
        
        if _USE_SSERVER and not(_detectSserver()):
            raise Exception("Sserver not detected when runCases was asked to use it.")

        if _RUN_PARALLEL and not _USE_SSERVER:
            _print("runCases: Running parallel without Sserver", also_to_nb=verbose)
            
        if _RUN_PARALLEL and _USE_SSERVER:
            _print("runCases: Running parallel with Sserver", also_to_nb=verbose)

            

        t0 = time()
        stream = os.popen("echo Returned output")
        if verbose:
            print(stream.read())


        allOK = True

        
        allOK = _wait_for_prepare(verbose)

                
        # start the cases
        if allOK:
            lstC = self.getCases()
            for case in lstC:
                if allOK:
                    case.run(verbose)
                    allOK = _handle_error_in_RUNNING_CASES()
        

        if allOK:
            allOK = wait_run(verbose)
            



        t1 = time()
        if allOK:
            dt = t1 - t0
            # Validation process - we do not output variable information:
            if os.environ.get("JUPYTER_IS_VALIDATION", False):
                dt = 0
                # Extract total time for all the runs of datasets
                tr = os.environ.get("TRUST_ROOT", "")
                gtt = os.path.join(tr, "Validation", "Outils", "Genere_courbe", "scripts", "get_total_time")
                l = [c._fullPath() for c in lstC]
                cmd = gtt + " " + " ".join(l)
                os.chdir(BUILD_DIRECTORY)
                _runCommand(cmd, False)
                os.chdir(ORIGIN_DIRECTORY)

            _print("  => A total of %d cases were (successfully) run in %.1fs." % (len(lstC), dt), also_to_nb=True)
        else:
            raise ValueError ("at least one case has failed ! See previous logs to get more information")
        
    def tablePerf(self):
        """ Prints the table of performance
        """
        from . import plot
        columns = ["host", "system", "Total CPU Time", "CPU time/step", "number of cells"]
        zeTable = plot.Table(columns)
        for case in self.getCases():
            try:
                case._addPerfToTable(zeTable)
            except Exception as e:
                raise e
        zeTable.sum("Total CPU Time")
        return zeTable.df

    def extractNRCases(self, verbose=False):
        """
        Prints out the list of cases in a suitable format for processing by validation and lance_test tools.
        
        verbose does nothing, just in case we forgot to remove it when actually extracting the cases 
        (we do a complicated sed on the notebook at this step)
        
        Teo Boutin: IMO, this shouldn't exist. 
        concatenating the python parts of the report, doing various replacements and executing that is just a bad idea.
        test case extraction/lauch should be controlled by options given ro Run_fiche
        passed through the env var JUPYTER_RUN_OPTIONS which is parser in this file (should use argparse btw)
        
        with args such as 
        -get_list_cas_nr
        -get_nb_cas_nr
        -run_single_nr_case <test_name.data>
        
        for now I will change this method to take care of executing all pre_runs of registered tests
        
        WARNING:
        do not modify this without looking at scripts get_list_cas_nr and get_nb_cas_nr in Validation/Outils/Genere_Courbe/scripts
        """
        import numpy as np
        
        # have to wait for cases launched in the prepare to complete (probably)
        # at least to avoid errors when deleting directories
        allOK = _wait_for_prepare()
        if not allOK: 
            raise RuntimeError("Case failed in prepare")
        # do the pre_run of the cases if they are in python and not bash scripts
        
        lstC = self.getCases()
        for case in lstC:
            if allOK and case.has_python_pre_run:
                case._preRun(verbose)
                allOK = _handle_error_in_RUNNING_CASES()

        if allOK: 
            allOK = wait_run()

        if not allOK: 
            raise RuntimeError("Case failed in preRun")
        
        
        list_exclu_nr = []
        if os.path.exists("src/liste_cas_exclu_nr"):
            list_cases = np.loadtxt("src/liste_cas_exclu_nr", dtype=str)
            if len(list_cases) == 0:
                raise RuntimeError("Empty file 'liste_cas_exclu_nr' in your src directory. Please remove it.")
            list_exclu_nr = list(map(lambda a: os.path.normpath(a), list_cases))

        for c in self.getCases():
            if c.dir_ != ".":
                t = os.path.join(c.dir_, c.name_ + ".data")
            else:
                t = c.name_ + ".data"
            t = os.path.normpath(t)
            
            if c.excluNR: continue
            if t in list_exclu_nr: continue
            
            print("@@@CAS_NR_JY@@@ " + t)

    def printCases(self):
        """
        display testCases
        """
        text = "### Test cases \n"
        for c in self.getCases():
            text += "* " + c.dir_ + "/" + c.dataFileName_ + ".data "
            if (int(c.nbProcs_) > 1):
                text += "with " + str(c.nbProcs_)  + " procs"
            text += "\n"
        displayMD(text)

def saveFileAccumulator(data):
    """ Method for saving files for the Non Regression test

    Parameters
    ---------
    data: str
        name of the file we want to save.
    """
    from .filelist import FileAccumulator

    path = os.getcwd()
    os.chdir(BUILD_DIRECTORY)

    from pathlib import Path

    data=sanitizePathToBUILD_DIRECTORY(data, relative=True)

    FileAccumulator.active = True
    new = FileAccumulator.Append(data)
    if new==True:
        FileAccumulator.WriteToFile("used_files")

    os.chdir(path)


def sanitizePathToBUILD_DIRECTORY(input_path, relative=False):
    """ Method for checking that a path given to trustutils is correctly located inside BUILD_DIRECTORY and for getting a sanitized version


    Parameters
    ---------
    input_path: str
        path that must be either absolute and located inside BUILD_DIRECTORY or relative (will check that it does not leave BUILD_DIRECTORY)
    relative: bool
        if False (default), path is returned as absolute
        if True, path is returned as relative to BUILD_DIRECTORY (for use in saveFileAccumulator mostly)
    """
    from pathlib import Path
    from os.path import normpath
    path=Path(input_path)

    if path.is_absolute():
        output_path=Path(normpath(input_path))

        try:
            p=output_path.relative_to(Path(BUILD_DIRECTORY), walk_up=False)
        except ValueError as e:
            raise ValueError(f"Invalid path given to trustutils:'{input_path}' \npath must be absolute and located inside BUILD_DIRECTORY '{BUILD_DIRECTORY}' or relative (will then be searched inside BUILD_DIRECTORY).")

    else:
        output_path=Path(normpath(os.path.join(BUILD_DIRECTORY, input_path)))

    if relative:
        return output_path.relative_to(Path(BUILD_DIRECTORY), walk_up=False).as_posix()
    else:
        return output_path.as_posix()

def introduction(auteur, creationDate=None):
    """ Function that creates an introduction cell Mardown

    Parameters
    ---------

    author: str 
        Name of the author of the test case.
    creationDate: date 
        format dd/mm/YYYY 
    """
    from datetime import datetime

    format = "%d/%m/%Y"
    today = datetime.today()
    # Validation process - we do not output variable information:
    if os.environ.get("JUPYTER_IS_VALIDATION", False): dat = ""
    else:                                              dat = today.strftime(format)

    displayMD("## Introduction \n Validation made by : " + auteur + "\n")
    if creationDate:
        try:
            datetime.strptime(creationDate, format)
        except ValueError:
            print("This is the incorrect date string format. It should be DD/MM/YYYY")
        displayMD("\n Report created : " + creationDate + "\n")
    displayMD("\n Report generated " + dat)


def TRUST_parameters(version="", param=[]):
    """ Function that creates a cell Mardown giving TRUST parameters (version, binary)
    
    Parameters
    ----------

    version: str 
        version of TRUST - if void TRUST_VERSION read
    param: list(str) 
        List of Parameter used in this test case
    """
    # Validation process - we do not output variable information:
    builtOn = BUILD_DIRECTORY
    if os.environ.get("JUPYTER_IS_VALIDATION", False):
        binary, version, builtOn = "BINARY", "VERSION", "BUILDDIR"
    elif version == "":
        binary, version = os.environ.get("exec", "UNKNOWN"), os.environ.get("TRUST_VERSION", "UNKNOWN")
    else:
        binary = os.environ.get("exec", "UNKNOWN")

    text = "### TRUST parameters \n * Version TRUST: " + version + "\n * Binary used: " + binary + " (built in directory " + builtOn + ")"
    for i in param:
        text = text + "\n" + i
    displayMD(text)


def dumpDataset(fiche, list_Trust_user_words=[]):
    """ 
    Print out the .data file.

    Parameters
    ---------
    fiche: str 
        Path of the file
    list_Trust_user_words: list(str) 
        List of word the user wants to color in red.
    """
    c = TRUSTCase(directory=BUILD_DIRECTORY, datasetName=fiche)
    c.dumpDataset(list_Trust_user_words)


def dumpDatasetMD(data):
    """
    Print the .data file in MD format.

    Parameters
    ---------
    datafile: str
        relative path (to build directory) + name of datafile
    """

    c = TRUSTCase(directory=BUILD_DIRECTORY, datasetName=data)
    c.dumpDatasetMD()


def dumpText(fiche, list_keywords=[]):
    """ Print out the file.

    Parameters
    --------
    fiche: str 
        Path of the file
    """
    ## Save the file
    name = os.path.join(BUILD_DIRECTORY, fiche)
    saveFileAccumulator(fiche)

    f = open(name, "r")
    tmp = f.readlines()
    test = []
    for j in tmp:
        flag = 0
        if j[0] == "#":
            flag = 1
            j = "\033[38;5;88m" + j + "\033[0;0m"
        j = " " + j
        j = j.replace("\t", "\t ")
        if flag == 0:
            for i in list_keywords:
                # tmp=tmp.replace(i, "\033[196;5;30m"+i+"\033[0;0m ") # en gras
                j = j.replace(i, "\033[38;5;28m" + i + "\033[0;0m")  # en rouge
        test.append(j)
    f.close()

    print("".join(test))

def addCaseFromTemplate(templateData, targetDirectory, dic, nbProcs=1, targetData=None, execOptions="",excluNR=False, pre_run=None, post_run=None):
    """ Add a case to run to the list of globally recorded cases.
    
    Parameters
    ----------
    targetDirectory: str 
        targetDirectory where the case is stored (relative to build/)
    templateData: str
        Name of the template datafile
    dic: dictionary
        substitution of term (key identified with $ in datafile) by value in a new data file
    targetData : str
        if provided, templateData will be copied as targetData + dictionary applied + added to TRUSTSuite()
    nbProcs : int 
        Number of processors
    excluNR : bool 
        to remove of the non regression test cases

    Returns
    -------
    a new TRUSTcase instance we want to launch.
    """
    global defaultSuite_
    if defaultSuite_ is None:
        # When called for the first time, will copy src to build, and execute 'prepare' script if any
        defaultSuite_ = TRUSTSuite()

    if targetData is None:
        targetData = templateData
    else:
        if targetData[-5:] != ".data":
            raise ValueError("targetData should ends with .data in addCaseFromTemplate!!!")

    fullDir = os.path.join(BUILD_DIRECTORY, templateData)
    fullDir2 = os.path.join(BUILD_DIRECTORY, targetDirectory)
    if not os.path.exists(fullDir2):
        os.makedirs(fullDir2, exist_ok=True)
    pthTgt = os.path.join(BUILD_DIRECTORY, targetDirectory, targetData)
    # And copy the .data file:
    from shutil import copyfile

    copyfile(fullDir, pthTgt)
    tc = addCase(targetDirectory, targetData, nbProcs, execOptions, excluNR, pre_run=pre_run, post_run=post_run)
    tc.substitute_template(dic)
    return tc


def addCase(directoryOrTRUSTCase, datasetName="", nbProcs=1, execOptions="", excluNR=False, pre_run=None, post_run=None):
    """ 
    Add a case to run to the list of globally recorded cases.

    Parameters
    ---------

    directoryOrTRUSTCase: str 
        directory where the case is stored (relative to build/)
    datasetName: str 
        Name of the case we want to run.
    nbProcs : int 
        Number of processors
    execOptions : str
        TRUST Options to add at the execution of the test case 
    excluNR : bool 
        to remove of the non regression test cases

    Returns
    -------
    TRUSTcase instance we want to launch.
    """
    global defaultSuite_
    if isinstance(directoryOrTRUSTCase, TRUSTCase):
        if datasetName != "":
            raise ValueError("addCase() method can either be called with a single argument (a TRUSTCase object) or with at least 2 arguments (directory and case name)")
        tc = directoryOrTRUSTCase
        initCaseSuite()
        defaultSuite_.addCase(tc)
        return
    elif isinstance(directoryOrTRUSTCase, str):
        if datasetName == "":
            raise ValueError("addCase() method can either be called with a single argument (a TRUSTCase object) or with at least 2 arguments (directory and case name)")
        tc = TRUSTCase(directoryOrTRUSTCase, datasetName, nbProcs,execOptions=execOptions,excluNR=excluNR, pre_run=pre_run, post_run=post_run)
        initCaseSuite()
        defaultSuite_.addCase(tc)
        return tc

def initCaseSuite():
    """ Instantiate a default suite of cases """
    global defaultSuite_
    if defaultSuite_ is None:
        # When called for the first time, will copy src to build, and execute 'prepare' script if any
        defaultSuite_ = TRUSTSuite()

def reset(keepBuildDir=True,keepSubdirectories=False):
    """ 
    Wipe out build directory completly and reset default suite.
    """
    opt = os.environ.get("JUPYTER_RUN_OPTIONS", None)
    # Very specific to the validation process - we need to keep build there:
    if not opt is None and "-not_run" in opt:
        return

    import shutil


    global defaultSuite_
    defaultSuite_ = None
    if os.path.exists(BUILD_DIRECTORY):
        # only delete content. to avoid killing your terminal if you are inside BUILD_DIRECTORY
        # can even choose to keep empty dirs with the optional variable
        for root, dirs, files in os.walk(BUILD_DIRECTORY, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            if not keepSubdirectories:
                for name in dirs:
                    os.rmdir(os.path.join(root, name))
        if not keepBuildDir:
            shutil.rmtree(BUILD_DIRECTORY)


def getCases():
    global defaultSuite_
    if defaultSuite_ is None:
        return []
    return defaultSuite_.getCases()

def executeScript(scriptName, verbose=False, nonRegression=False):
    """ 
    Execute a script shell in the BUILD_DIRECTORY

    Parameters
    ----------

    scriptName: str
        Name of the exec script
    verbose: bool
    nonRegression: bool
        by default executeScript is inactive when option -not_run is applied (for non regression test)
    """
    opt = os.environ.get("JUPYTER_RUN_OPTIONS", "")
    # Very specific to the validation process. Sometimes we want the core
    # method 'runCases()' not to do anything ... see script 'archive_resultat' for example.
    if "-not_run" in opt and not nonRegression :
        return

    os.chdir(BUILD_DIRECTORY)
    pth = "./" + scriptName
    if os.path.exists(pth):
        cmd = pth
        subprocess.check_output("chmod u+x %s" % pth, shell=True)
        _runCommand(cmd, verbose)
    os.chdir(ORIGIN_DIRECTORY)

def executeCommand(cmd, verbose=False, nonRegression=False):
    """ Execute a bash command in the BUILD_DIRECTORY
    Parameters
    ----------

    cmd: command to execute
    verbose: bool
    nonRegression: bool
    """
    opt = os.environ.get("JUPYTER_RUN_OPTIONS", "")
    # Very specific to the validation process. Sometimes we want the core
    # method 'runCases()' not to do anything ... see script 'archive_resultat' for example.
    if "-not_run" in opt and not nonRegression :
        return
    os.chdir(BUILD_DIRECTORY)
    _runCommand(cmd, verbose)
    os.chdir(ORIGIN_DIRECTORY)

def extraitCoupe(data_file, probe_name, time=None, verbose=False, directory=None):
    """ execute script extrait_coupe of $TRUST_ROOT/bin 
    Parameters
    ----------
    directory: str
    data_file: str
    probe_name: str
    verbose: bool
    time: value or 'all' or nothing
    """
    cmd = ""
    if directory:
        cmd = "cd " + directory + " && "
    cmd += " extrait_coupe " + data_file + " " + probe_name
    if time:
        cmd+=" " + str(time)
    executeCommand(cmd, verbose)



def printCases():
    """
    display testCases
    """
    global defaultSuite_
    return defaultSuite_.printCases()


def extractNRCases(verbose=False):
    """
    Prints out the list of cases in a suitable format for processing by validation and lance_test tools.
    
    verbose does nothing, just in case we forgot to remove it when actually extracting the cases 
    (we do a complicated sed on the notebook at this step)
    
    WARNING:
    do not modify this without looking at scripts get_list_cas_nr and get_nb_cas_nr in Validation/Outils/Genere_Courbe/scripts
    """
    global defaultSuite_
    return defaultSuite_.extractNRCases()

def isExtractingNR():
    """
    Check environnement variable IS_EXTRACTING_NR to know if extracting_nr is active
    """

    if (os.getenv("IS_EXTRACTING_NR") == '1'): return True
    return False

def isExtractingNR_ListOnly():
    """
    Check global variable IS_EXTRACTING_NR_LIST_ONLY to know if we want only the list of nr cases
    Should only be used internally.
    This global variable is read from os.environ at the begining (import) of this module
    """

    return IS_EXTRACTING_NR_LIST_ONLY

def is_NOT_RUN():
    """
    Check global variable _NOT_RUN to know if option -not_run was set.

    Useful if lata files are read by the user (eg. with lataloader).
    And probably a lot of other use cases
    """

    return _NOT_RUN

def needPrepareForNRCaseList():
    """
    Sets IS_EXTRACTING_NR_LIST_ONLY to False.

    To use when the trust runs launched manually before run.runCases are necessary in order
    to generate the list of non regression test cases.
    """
    global IS_EXTRACTING_NR_LIST_ONLY
    IS_EXTRACTING_NR_LIST_ONLY=False
    os.environ['IS_EXTRACTING_NR_LIST_ONLY'] = '0' # just in case, also reset the env variable


def runCases(verbose=False, preventConcurrent=False):
    """ Launch all TRUST cases for the current validation form.

    Parameters
    ---------
    verbose: bool 
        whether to print the console output of the run.
    preventConcurrent: bool
        run the cases in the order they were provided, even if the Sserver is up and running, and the -parallel_sjob option was passed.
    """
    global defaultSuite_
    if defaultSuite_ is None:
        raise Exception("No test cases currently recorded! Call 'addCase' first ...")

    defaultSuite_.runCases(verbose, preventConcurrent)


def tablePerf():
    """ Prints the table of performance
    """
    global defaultSuite_
    return defaultSuite_.tablePerf()

def initBuildDirectory():
    """ triggers build directory creation and copy src stuff into it
    """
    global defaultSuite_
    defaultSuite_ = TRUSTSuite()



# Handling of parallelism of trust cases
# at global level. Other solution would be to handle this at the trust suite level. but then rogue cases in eg pre_run couldn't be handled
# default is parallel. Sserver usage is independent of parallelism
_NOT_RUN=False
_RUN_SEQUENTIAL=False
_RUN_PARALLEL=True
_USE_SSERVER=False

def _set_not_run():
    global _NOT_RUN
    global _RUN_SEQUENTIAL
    global _RUN_PARALLEL
    _NOT_RUN=True
    _RUN_SEQUENTIAL=False
    _RUN_PARALLEL=False
    _print("Not running trust cases")
    
def _set_run_sequential():
    global _NOT_RUN
    global _RUN_SEQUENTIAL
    global _RUN_PARALLEL
    _NOT_RUN=False
    _RUN_SEQUENTIAL=True
    _RUN_PARALLEL=False
    _print("Set sequential run")
    
def _set_run_parallel():
    global _NOT_RUN
    global _RUN_SEQUENTIAL
    global _RUN_PARALLEL
    _NOT_RUN=False
    _RUN_SEQUENTIAL=False
    _RUN_PARALLEL=True
    _print("Set parallel run")
    
def _set_use_sserver(v=True):
    global _USE_SSERVER
    _USE_SSERVER=v
    if v:
        _print("Using Sserver")
    else:
        _print("Disabling Sserver")
        
    
"""
_RUNNING_CASES global variable:

For storing informations about the subprocesses handling the trust cases

See usage in TRUSTCase.run method
Also in wait_run and _wait_for_prepare
Will store dict objects with the following structure:
{
    "case": a ref to the TRUSTCase object 
    "script": path to the launch script (.cmdsxxx)
    "logFile": path to the resulting log file  
    "callback": a lambda function that either do nothing or calls the post_run function 
    "callbackDone": to know whether the callback has already been called 
    "popen": a lambda function that can start the trust simulation with a subprocess.popen
    "process": Initially None, then the popen object returned by the 'popen' lambda function of this dict (see line above)
    "depends": a list of indices of cases (for this array _RUNNING_CASES) on which this specific cases depends. 
                For example, cases that were launched during the pre_run of this function
}

This object is stores both the cases that are running, but also those that are waiting for their pre_run, or waiting for available procs. Thus, the name may not be the most pertinent. Feel free to change it if you find a better one.

See the functions _count_running, _count_running_deps, _has_waiting_cases, _count_waiting_cases
to understand conditions that correspond to the different states

"""
_RUNNING_CASES=[]

def _count_running():
    # A case is considered over when the callback (post_run) is done
    c=0
    for r in _RUNNING_CASES:
        if r["process"] != None and not(r["callbackDone"]):
            c+=1
    return c
def _count_running_deps(deps):
    # A case is considered over when the callback (post_run) is done
    c=0
    for index in deps:
        r=_RUNNING_CASES[index]
        if r["process"] == None or not(r["callbackDone"]):
            c+=1
    return c

def _has_waiting_cases():
    for r in _RUNNING_CASES:
        if r["process"] == None:
            return True
    return False

def _count_waiting_cases():
    c=0
    for r in _RUNNING_CASES:
        if r["process"] == None:
            c+=1
    return c

def _count_procs_usage():
    return sum([(r["case"].nbProcs_ if (r["process"] and r["process"].poll() == None) else 0) for r in _RUNNING_CASES])

def _wait_for_available_procs(n_procs):
    """ 
    When using the parallel_run mode (without Sserver), we must manually manage proc usage from here
    This function returns when enough procs are available.

    To allow for occasional runs that ask for more procs than the total, 
    this also stops when all procs become available.

    TODO: a potentially better approach may be to launch all cases that fit into the procs available, rather than waiting on each individual cases. Right now, cases where we have 3 cases that uses 2, 4, 2 procs, with 4 available will lauch the first, wait till it finishes to lauch the second, then wait again and lauch the third, when we could have done this in 2 batches. 

    Still, current handling is sufficient for use cases of this option (and it is also way better than it was in the past...)
    
    
    """
    max_procs = int(os.environ.get("TRUST_NB_PROCS","")) # if not set, then what ? should not happen anyway...
    used = _count_procs_usage()
    free = max_procs - used
    if n_procs > max_procs:
        _print("\nAsking for more procs than TRUST_NB_PROCS. Case will run when all procs are free.\n")
    
    # this loop ends if either there are at least <n_procs> which are free, 
    # or if all procs are free (to allow the case n_cpu > max_used_procs)
    while used > 0 and free < n_procs:
        _print(f"\nWaiting for {n_procs} free procs: {used}/{max_procs} used.\n")
        sleep(1)
        used = _count_procs_usage()
        free = max_procs - used

# p should be a Popen object or None
def _is_process_running(p):
    # p is not None means the case was started (we ran Popen and stored the resulting Popen Object)
    # p.poll() is not None when the subprocess has finished
    return p is not None and p.poll() is None

# p should be a Popen object or None
def _is_process_finished(p):
    # p is not None means the case was started (we ran Popen and stored the resulting Popen Object)
    # p.poll() is not None when the subprocess has finished
    return p is not None and p.poll() is not None

def _handle_error_in_RUNNING_CASES():

    sep = "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    err_msg = sep
    err_msg += "\nCase '%s' FAILED !! Below you will find the last 20 lines of the trust cerr and the log of the wrapper:"
    err_msg += "\n(If you don't see anything suspicious, also check pre/post_run scripts!!)\n"
    err_msg += sep

    allOK=True
    for i,r in enumerate(_RUNNING_CASES):
        p=r["process"]
        
        if _is_process_finished(p):
            case = r["case"]

            # if it failed, print error and abort
            if p.returncode != 0:
                allOK = False
                _print(err_msg % (case._relPath()), also_to_nb=True)
                err_file=case._fullPath_ErrFile()
                if os.path.isfile(err_file):
                    print(sep)
                    _print("errfile path:" ,err_file, "\n", also_to_nb=True)
                    _print(getLastLines_(err_file), also_to_nb=True)
                logFile=r["logFile"]
                if os.path.isfile(logFile):
                    print(sep)
                    _print("logFile path:" ,logFile, also_to_nb=True)
                    print(sep)
                    _print(getLastLines_(logFile), also_to_nb=True)
                # only report about the first failed case.
                # otherwise, too much clutter
                # user will learn about each failed case after fixing previous one
                # debatable. 
                break 
    return allOK

def wait_run(verbose=False):
    """ 
    Wait until all TRUST Cases have run, including cases launched from their post_run functions.
    This is called in TRUSTSuite methods runCases and extractNRCases.
    
    This function waits till all cases in the global list _RUNNING_CASES are completed.
    post_run of such cases may append to the list, thus the need for this function.

    In case of failure of a single case, it should correctly abort all others cases and display a clear error message.
    
    Parameters
    ---------
    verbose: bool
    
    """


    _print("Starting wait_run")
    allOK = True
    
    
    run_count=_count_running()
    waiting=_count_waiting_cases()
    reps=0
    while allOK and (waiting>0 or (run_count>0)):
        reps+=1
        # Wait for a bit to avoid polling too frequently
        if reps>1 and not(isExtractingNR_ListOnly()):
            sleep(1)

        if os.getenv("TRUSTUTILS_VERBOSE") == '1':
            tf=strftime('%H:%M:%S')
            _print(f"\n[{tf}]")
            _print("Running:", run_count, "| Waiting for pre_run:", waiting, "| Finished:", len(_RUNNING_CASES) - run_count - waiting, "| Total:", len(_RUNNING_CASES), "\n")


        allOK=_handle_error_in_RUNNING_CASES()
        
        if allOK:
            for i,r in enumerate(_RUNNING_CASES):
                p=r["process"]
                
                if _is_process_finished(p):
                    case = r["case"]

                    # failed cases handled in _handle_error_in_RUNNING_CASES
                    # if we get here, raise an error
                    if p.returncode != 0:
                        raise RuntimeError("A failed case was not handled correctly. Contact TRUST Team")

                    if allOK and not(r["callbackDone"]):
                        _print("Finished case", case._relPath(), also_to_nb=verbose)
                        r["callback"]()
                        r["callbackDone"]=True
                    
        if allOK: # handle deps of running cases
            for i,r in enumerate(_RUNNING_CASES):
                p=r["process"]
                
                if p is None:
                    case = r["case"]
                    deps = r["depends"]
                    deps_running=_count_running_deps(deps)
                    if deps_running == 0:
                        # Wait for cpu availability if Sserver is not managing the jobs
                        if not(_USE_SSERVER):
                            _wait_for_available_procs(case.nbProcs_)
                            
                        _print("Starting case", case._relPath(), "from waiting list", also_to_nb=verbose)
                        
                        r["process"] = r["popen"]()
                
        
        if not allOK:
            _print("ABORTING all running cases, because a case failed")
            for r in _RUNNING_CASES:
                if r["process"] is not None and p.poll() is None: r["process"].terminate()
        
        # IMPORTANT
        # update variables used in while condition
        run_count=_count_running()
        waiting=_count_waiting_cases()
        
        
    _print("Total runs:", len(_RUNNING_CASES))
    
    return allOK
    
    
def _wait_for_prepare(verbose=False):
    """ 
    Wait till all runs that may have been started manually before runCases are done 
    (the part that we call 'prepare' as it replaces the old prepare from prm reports)
    pre and post_run at this point are not handled. Maybe done later.

    Parameters
    ---------
    verbose: bool
    
    """
    
    _print("Starting _wait_for_prepare")
    allOK=True
    
    # err_msg = "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    # err_msg += "Case '%s/%s.data' FAILED !! Here are the last 20 lines of the log file:\n"
    # err_msg += "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    
    if _has_waiting_cases():
        raise Exception("Cases with pre_run that launches other cases are not allowed in prepare.")
    _print(f"\nrunCases: Waiting for cases in prepare\n")
    for r in _RUNNING_CASES:
        allOK = _handle_error_in_RUNNING_CASES()
        if not allOK:
            break
        p=r["process"]
        case = r["case"]
        p.wait()
        
        if not(r["callbackDone"]):
            _print("Finished prepare case", case._relPath(), also_to_nb=verbose)
            # post run is called, but post_run are not allowed to start more TRUSTCases
            nj= len(_RUNNING_CASES)
            r["callback"]()
            r["callbackDone"]=True
            if len(_RUNNING_CASES) > nj or _has_waiting_cases():
                raise Exception(f"Cases with post_run that launches other cases are not allowed in prepare.\n From case {case._relPath()}")

    if not allOK:
        _print("ABORTING, a case failed in prepare")
        for r in _RUNNING_CASES:
            if r["process"]: r["process"].terminate()
        
    return allOK
    
ORIGIN_DIRECTORY = os.getcwd()

JUPYTER_RUN_OPTIONS=None

IS_EXTRACTING_NR_LIST_ONLY=(os.getenv("IS_EXTRACTING_NR_LIST_ONLY") == '1')

defaultSuite_ = None  # a TRUSTSuite instance



import sys
# add imports directory to path and save contained scripts
import_path=os.path.join(ORIGIN_DIRECTORY, "src", "python_modules")
if os.path.exists(import_path):
    sys.path.append(import_path)
# imports = list(pathlib.Path(import_path).rglob("*.py"))
# for f in imports:
#     saveFileAccumulator(str(f))


# I think it is IMPORTANT to call this at the very end. Not sure though, but be careful when expanding this file
BUILD_DIRECTORY = _initBuildDir()

