"""
Script for updating the version in construct.yaml
"""
import argparse
import datetime
import os
import sys

def get_version_from_file(filename):
  version = os.path.basename(filename).split('-')[1]  # dev.1234 or 1.23.4567
  split_version = version.split('.')
  build_number = split_version[-1]  # 1234
  if len(split_version) > 2:
    version = '.'.join(split_version[:2])  # 1.23
  else:
    version = split_version[0]  # dev
  version = '-'.join([version, build_number])  # dev-1234 or 1.23-4567
  return version

def update_version(filename, pkg_name, version):
  # check if version is a filename (e.g. conda package with version)
  # version follows the Phenix format (e.g. dev-1234 or 1.23-4567)
  # conda package naming does not allow hyphens (e.g. dev.1234 or 1.23.4567)
  if os.path.isfile(version):
    version = get_version_from_file(version)
  elif version == 'windows':  # special handling for Windows
    search_dir = 'D:\\a\\1\\win-64'
    contents = os.listdir(search_dir)
    for f in contents:
      f = os.path.join(search_dir, f)
      if os.path.isfile(f) and (f.endswith('conda') or f.endswith('bz2')):
        version = get_version_from_file(f)
        break
  assert os.path.isfile(filename)
  with open(filename, 'r') as f:
    lines = f.readlines()
  with open(filename, 'w') as f:
    for line in lines:
      if line.startswith('pkg_name') and pkg_name is not None:
        f.write(f'pkg_name: {pkg_name}-{version}\n')
      elif 'REPLACEME' in line:
        f.write(line.replace('REPLACEME', version))
      else:
        f.write(line)

  return version

if __name__ == '__main__':

  now = datetime.datetime.now()
  default_version = '.'.join([str(now.year), str(now.month), str(now.day)])

  parser = argparse.ArgumentParser(description=__doc__)

  parser.add_argument('--filename', default='construct.yaml', type=str,
    help='The file to be patched')
  parser.add_argument('--pkg-name', default=None, type=str,
    help='The package name (the version will be appended to this)')
  parser.add_argument('--version', default=default_version, type=str,
    help='The version')

  # show help if no arguments are provided
  if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

  namespace = parser.parse_args()

  # patch
  version = update_version(namespace.filename, namespace.pkg_name, namespace.version)

  print(version)
