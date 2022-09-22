"""
Script for updating the version in construct.yaml
"""
import argparse
import datetime
import os
import sys

def update_version(filename, pkg_name, version):
  assert os.path.isfile(filename)
  with open(filename, 'r') as f:
    lines = f.readlines()
  with open(filename, 'w') as f:
    for line in lines:
      if line.startswith('pkg_name'):
        f.write(f'pkg_name: {pkg_name}-{version}')
      elif line.startswith('version:'):
        f.write(f'version: {version}')
      else:
        f.write(line)

if __name__ == '__main__':

  now = datetime.datetime.now()
  default_version = '.'.join([now.year, now.month, now.day])

  parser = argparse.ArgumentParser(description=__doc__)

  parser.add_argument('--filename', default='construct.yaml', type=str,
    help='The file to be patched')
  parser.add_argument('--pkg-name', default='phenix', type=str,
    help='The package name (the version will be appended to this)')
  parser.add_argument('--version', default=default_version, type=str,
    help='The version')

  # show help if no arguments are provided
  if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

  namespace = parser.parse_args()

  # patch
  update_version(namespace.filename, namespace.pkg_name, namespace.version)
