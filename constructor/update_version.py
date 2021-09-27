"""
Script for updating the version in construct.yaml
"""
import argparse
import os
import sys

def update_version(filename, version):
  assert os.path.isfile(filename)
  with open(filename, 'r') as f:
    lines = f.readlines()
  with open(filename, 'w') as f:
    for line in lines:
      if line.startswith('version:'):
        f.write(f'version: {version}')
      else:
        f.write(line)

if __name__ == '__main__':

  parser = argparse.ArgumentParser(description=__doc__)

  parser.add_argument('--filename', default='construct.yaml', type=str,
    help='The file to be patched')
  parser.add_argument('--version', default=None, type=str,
    help='The version')

  # show help if no arguments are provided
  if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

  namespace = parser.parse_args()

  # patch
  update_version(namespace.filename, namespace.version)
