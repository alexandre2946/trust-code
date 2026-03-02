import sys
import os
import argparse

def usage():
   print("Extract list of test cases used in a Jupyter Notebook.")
   print("This is done by extracting first all Python cells from the notebook and then replacing calls to 'runCases'")
   print("by a specific instruction outputing the test case list into a special file (./list_cases_jy)")
   print("   python get_nr_cases_jup.py <notebook.ipynb>")
   print("")
   print("No input notebook specified!")

def extract_cases(args):
   import json
   import re

   with open(args.notebook, "r", encoding='utf-8') as f:
     root = json.loads(f.read())
     # Now parse tree to retrieve all Python up to the invocation of "runCases()":
     cells = root.get("cells", [])
     
     done = False
     s = "os.environ['IS_EXTRACTING_NR'] = '1'\n"

     # another flag, for the case where we want to only get the list of cases
     # because copie_cas_test will also use this, but needs to run cases that are lauched before run.runCases
     # we use os.environ to pass this message to the same script who later uses os.environ.get, which may seem absurd
     # but this is necessary as the script may also be executed from the jupyter notebook, and we need to use a method that provides a default
     if args.listonly:
       s+="os.environ['IS_EXTRACTING_NR_LIST_ONLY'] = '1'\n"
     else:
       s+="os.environ['IS_EXTRACTING_NR_LIST_ONLY'] = '0'\n"

     for c in cells:
       s += "\n"
       if done: break
       if c["cell_type"] == "code":
         if "source" in c:
           for loc in c["source"]:
             # for now, skip shell code lines:
             if loc.strip().startswith("!") or loc.strip().startswith("#"):
               continue
             # skip also display of plot: 
             if "plot(" in loc: 
               continue
             if ".runCases(" in loc:
               done = True
               s += re.sub(r"^(.*).runCases\([^\)]*\)(.*)$", "\\1.extractNRCases()\\2", loc)
               break
             s += loc
   # Now execute this code - it will perform whatever action the user wants to have as 'prepare' steps
   # and it will finally print out the list of test cases.

   print("=========== Here is the part extracted from the notebook ===========")
   print(s)
   print("======== Here is the output when executing previous python code ==========")

   exec(s, globals())

if __name__ == "__main__":

   parser=argparse.ArgumentParser(prog="get_nr_cases_jup")
   parser.add_argument("notebook", help="Path to the notebook") 
   parser.add_argument("-l", "--listonly", action='store_true', help="Indicate if we only want the list of nr cases. In that case, an env variable will be set to indicate to trustutils that we only want the list.") 
   args = parser.parse_args()

   extract_cases(args)

