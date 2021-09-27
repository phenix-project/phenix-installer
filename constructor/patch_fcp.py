"""
Script for patching fcp.py in constructor to avoid checking if a directory is writable
"""
import argparse
import os
import sys

def patch_fcp(site_packages=None):
  '''
  Comment out line that checks if the directory is writable
  '''
  fcp_file = os.path.join(site_packages, 'constructor', 'fcp.py')
  assert os.path.isfile(fcp_file)
  with open(fcp_file, 'r') as f:
    lines = f.readlines()
  line_to_patch = 'assert pc.is_writable, download_dir'
  print(fcp_file)
  with open(fcp_file, 'w') as f:
    for line in lines:
      if line_to_patch in line:
        print('Removed line:', line)
        continue
      else:
        f.write(line)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)

  parser.add_argument('--site-packages', default=None, type=str,
    help='The site-packages directory, usually ${CONDA_PREFIX}/lib/python${PY_VER}/site-packages')

  # show help if no arguments are provided
  if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

  namespace = parser.parse_args()

  # patch
  patch_fcp(site_packages=namespace.site_packages)
