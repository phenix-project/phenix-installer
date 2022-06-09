'''
Script for uploading the installer to Google Drive
'''

import sys
from argparse import ArgumentParser
from xfel.command_line.upload_mtz import pydrive2_interface

if __name__ == '__main__':

  parser = ArgumentParser(description=__doc__)
  parser.add_argument(
    '--credentials', help='Your credentials JSON from Google'
  )
  parser.add_argument(
    '--folder_id', help='The folder id on the Google Drive'
  )
  parser.add_argument(
    '--folder_list', help='The folder_list argument for pydrive2_interface'
  )
  parser.add_argument(
    '--files', help='The files argument for pydrive2_interface'
  )

  if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

  namespace = parser.parse_args(sys.argv[1:])

  pi = pydrive2_interface(
    cred_file=namespace.credentials,
    folder_id=namespace.folder_id
  )
  pi.upload(folder_list=[namespace.folder_list], files=[namespace.files])
  # pi = pydrive2_interface(
  #   cred_file='/Users/bkpoon/Downloads/phenix-lbl-4f3f6ca7f212.json',
  #   folder_id='1DVRkUc_nlS4i19zyI7OVdUcHcTI50YTz')
  # pi.upload(folder_list=['test'], files=['cctbx.xfel-VERSION-MacOSX-arm64.sh'])
