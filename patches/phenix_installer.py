#
# TOP-LEVEL INSTALLER SCRIPT FOR MAC AND LINUX
#
# Invoked by indirectly by make_dist() in libtbx\auto_build\create_installer.py
# when generating zipped distribution on Windows
#
# This subclasses the base installer class included in CCTBX in this file:
#   libtbx/auto_build/install_distribution.py
#

from __future__ import division
from __future__ import print_function
import glob
import shutil
import os.path
op = os.path
import sys
libtbx_path = os.path.join(
  os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "lib")
if (not libtbx_path in sys.path) :
  sys.path.append(libtbx_path)
from libtbx.auto_build import install_distribution
from libtbx.auto_build.installer_utils import *
from libtbx.auto_build import bootstrap

def get_modules(builder):
  modules = list()
  for class_var in ['CODEBASES', 'CODEBASES_EXTRA', 'HOT', 'HOT_EXTRA', 'LIBTBX',
                    'LIBTBX_EXTRA']:
    modules += getattr(builder, class_var, [])
  modules = set(modules)
  # msgpack has a different directory name
  modules.discard('msgpack')
  modules.add('msgpack-3.1.1')
  modules.add('qrefine')
  # remove voyager_regression for Windows
  if sys.platform == 'win32':
    modules.discard('voyager_regression')
    modules.discard('phasertng')
  return modules

class installer (install_distribution.installer) :
  product_name = "Phenix"
  dest_dir_prefix = "phenix"
  make_apps = ["phenix"]
  include_gui_packages = True
  base_package_options = ["--all"]
  installer_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
  base_modules = []
  modules = get_modules(bootstrap.PhenixReleaseBuilder)
  configure_modules = modules

  if sys.platform == "win32": # exclude some modules on Windows if necessary
    configure_modules = list( set(configure_modules) - set(bootstrap.windows_remove_list) )
    modules = list( set(modules) - set(bootstrap.windows_remove_list) )

  def add_product_specific_options (self, parser) :
    parser.add_option("--rosetta", default=False, action="store_true",
      help="Also compile Rosetta with Phenix extensions")
    parser.add_option("--without-pymol", default=False, action="store_true",
      help="Disable installation of PyMOL executable")
    parser.add_option("--without-dials", default=False, action="store_true",
      help="Disable installation of DIALS-related executables (dials.*, dxtbx.*, xia2.*, iota.*")

  def product_specific_setup_before_compile (self, log) :
    phaser_version_cc = op.join(self.modules_dir, "phaser", "Version.cc")
    if not op.isfile(phaser_version_cc) :
      open(phaser_version_cc, "w").write("")

  def product_specific_source_install (self, log) :
    out = self.out
    if (self.options.rosetta) :
      print("Attempting to compile Rosetta with Python embedded...", file=out)
      rosetta_build_bin = op.join(self.build_dir, "bin",
        "rosetta.build_phenix_interface")
      try :
        log = open(op.join(self.tmp_dir, "compile_rosetta.log"), "w")
        call([rosetta_build_bin, "nproc=%s" % self.options.nproc], log=log)
      except RuntimeError :
        print("", file=out)
        print("Error: Rosetta build failed.  Will continue with normal", file=out)
        print("Phenix installation.  You can view the command output", file=out)
        print("here:", file=out)
        print(op.join(self.tmp_dir, "compile_rosetta.log"), file=out)
        print("", file=out)
      else :
        print("  Rosetta compiled successfully.", file=out)
        print("", file=out)
    #if (not self.options.without_pymol) :
    #  self.install_pymol(out=out)

  def product_specific_dispatcher_epilogue (self) :
    if sys.platform == "win32":
      return [] # QuantumBio is only available on Linux
    return [
      # XXX QuantumBio environment variables
      """if [ ! -z "$QB_PYTHONPATH" ]; then""",
      """  export PYTHONPATH=$PYTHONPATH:$QB_PYTHONPATH""",
      """fi""",
      """if [ ! -z "$QB_LD_LIBRARY_PATH" ]; then""",
      """  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$QB_LD_LIBRARY_PATH""",
      """fi""",
      """if [ ! -z "$QB_DYLD_LIBRARY_PATH" ]; then""",
      """  export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$QB_DYLD_LIBRARY_PATH""",
      """fi""",
      # XXX end of QBio stuff
      """if [ "$PHENIX_MTYPE" != "mac-ppc-osx" ] && \\""",
      """   [ "$PHENIX_MTYPE" != "mac-intel-osx" ] && \\""",
      """   [ "$PHENIX_MTYPE" != "mac-intel-osx-x86_64" ]; then""",
      """  export PYMOL_PATH=$PHENIX/pymol""",
      """fi""",
    ]

  def product_specific_finalize_install (self, log) :
    out = self.out

    # solve_resolve license
    license_file = os.path.join(self.modules_dir, "solve_resolve", "ext_ref_files", "solve2.access")
    f = open(license_file, "w")
    print("Generating SOLVE/RESOLVE license file...", file=out)
    f.write("Access code for solve for PHENIX users\n")
    f.write("31087535G5\n")
    f.write("\n")
    f.close()

    # rebuild docs
    print("Rebuilding HTML documentation...", file=out)
    try:
      call([os.path.join(self.build_dir, 'bin', 'phenix_html.rebuild_docs')], log=log)
    except RuntimeError:
      print("Error rebuilding documentation!", file=out)

    # generate rotarama pickle files
    rotarama_dir = op.join(self.modules_dir, "chem_data", "rotarama_data")
    os.chdir(rotarama_dir)
    print("Creating rotamer/Ramachandran database files...", end=' ', file=out)
    call([os.path.join(self.build_dir, 'bin', 'mmtbx.rebuild_rotarama_cache')], log=log)
    no_update_file = op.join(rotarama_dir, "NO_UPDATE")
    open(no_update_file, "w").write("1")
    print("ok", file=out)
    # verify monomer library
    print("Checking for monomer library files...", end=' ', file=out)
    try :
      call([os.path.join(self.build_dir, 'bin', 'mmtbx.verify_mon_lib_data')], log=log)
    except RuntimeError :
      raise InstallerError("Incomplete monomer library installation! "+
        "This usually indicates a corrupted installer; please "+
        "try downloading again and re-installing.")
    print("ok", file=out)

    if (not self.options.without_pymol) :
      self.setup_pymol(out=out)

    if self.options.without_dials:
      print("Removing DIALS-related executables...", file=out)
      cwd = os.getcwd()
      os.chdir(os.path.join(self.build_dir, 'bin'))
      for module in ['dials', 'dxtbx', 'iota', 'xia2']:
        filenames = glob.iglob(module + '.*')
        for filename in filenames:
          print("  Removing %s ... " % filename, end=' ', file=out)
          try:
            os.remove(filename)
            print("done", file=out)
          except Exception:
            print("failed", file=out)
      print('',file=out)
      os.chdir(cwd)

  def product_specific_reduce_installation_size (self, out) :
    # remove CABLAM raw text files (unnecessary since we have pickles)
    cablam_dir = op.join(self.modules_dir, "chem_data", "cablam_data")
    for file_name in os.listdir(cablam_dir) :
      if file_name.endswith(".stat") :
        full_path = op.join(cablam_dir, file_name)
        os.remove(full_path)
    # top-level directories that we don't need
    for dir_name in ["phenix_regression", "phenix_html"] :
      full_path = op.join(self.modules_dir, dir_name)
      if op.isdir(full_path) :
        shutil.rmtree(full_path)
    # chem_data stuff
    lig_lib_dir = op.join(self.modules_dir, "chem_data", "ligand_lib")
    for file_name in lig_lib_dir :
      if file_name.startswith("PDB.") :
        full_path = op.join(lig_lib_dir, file_name)
        os.remove(full_path)
    rotarama_dir = op.join(self.modules_dir, "chem_data", "rotarama_data")
    for file_name in os.listdir(rotarama_dir) :
      if ((file_name.endswith(".data")) and
          (not file_name.startswith("rama500"))) :
        full_path = op.join(rotarama_dir, file_name)
        pkl_file = full_path.replace(".data", ".pickle")
        if op.isfile(pkl_file) :
          os.remove(full_path)

  def install_pymol (self, out) :
    base_bin_dir = op.join(self.base_dir, "bin")
    pymol_dir = op.join(self.installer_dir, "pymol")
    if (not op.isdir(pymol_dir)) :
      print("PyMOL not available!", file=out)
      return False
    print("Installing PyMOL...", end=' ', file=out)
    try :
      log = open(op.join(self.tmp_dir, "pymol.log"), "w")
      os.chdir(self.dest_dir)
      if (sys.platform == "darwin") :
        pymol_pkg = op.join(pymol_dir, "MacPyMOL.app.tar.gz")
        untar(pymol_pkg, log=log, check_output_path=False)
      else :
        pymol_pkg = op.join(pymol_dir,
          "pymol-0_99rc6-bin-linux-x86-glibc23.tar.gz")
        untar(pymol_pkg, log=log, check_output_path=False)
    except Exception as e :
      print("failed!", file=out)
      print("Error message:", file=out)
      print("  %s" % str(e), file=out)
      return False
    print("ok", file=out)
    return True

  # configure an already in-place PyMOL (Linux only)
  def setup_pymol (self, out) :
    base_bin_dir = op.join(self.base_dir, "bin")
    pymol_dir = op.join(self.dest_dir, "pymol")
    if (not op.isdir(pymol_dir)) :
      print("PyMOL not available!", file=out)
      return False
    print("Creating PyMOL launcher...", end=' ', file=out)
    os.chdir(pymol_dir)
    args = [ "./setup.sh" ]
    log = open(op.join(self.tmp_dir, "pymol_setup.log"), "w")
    call(args=args, log=log)
    shutil.copyfile(op.join(os.getcwd(), "pymol"),
                    op.join(base_bin_dir, "pymol"))
    print("ok", file=out)
    return True

  def display_final_message (self) :
    if sys.platform == "win32":# not applicable for the setup.exe installer
      return
    print("""
==========================================================================

                      PHENIX installation complete
                      ----------------------------

 You can begin using PHENIX now by setting your environment with the
 'source' command:

   csh users:
     source %s/phenix_env.csh

   bash users:
     source %s/phenix_env.sh

 To use PHENIX, go to a work directory and type:

   phenix

 You may wish to put the source statement in your .cshrc or .bashrc
 file.

"""%(self.dest_dir, self.dest_dir), file=self.out)

if (__name__ == "__main__") :
  installer(sys.argv[1:]).install()
