'''
Script for uploading the installer to Google Drive

conda create -n gd python=3.10
conda activate gd
pip install google-api-python-client google-auth-httplib2 google-auth-oauthlib oauth2client

https://developers.google.com/resources/api-libraries/documentation/drive/v3/python/latest/index.html
'''

import argparse
import os
import sys
import time

from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
from googleapiclient.http import MediaFileUpload

from oauth2client.service_account import ServiceAccountCredentials

# =============================================================================
# "My Drive" for service account
my_drive_id = '0AG16Yc-WJrbLUk9PVA'

# "Installers" folder in "Phenix Installers" drive
installers_id = '1m2p6k6UdF798vlH1D6eemvTqngj5Ei-q'

# =============================================================================
def get_drive_id(name=None, credentials=None):
  try:
    service = build('drive', 'v3', credentials=credentials)
    results = service.drives().list().execute()
    items = results.get('drives', [])
    if len(items) == 0:
      return None
    for item in items:
      if name == item['name']:
        return item['id']
    return None
  except HttpError as error:
    print(f'An error occurred: {error}')

# -----------------------------------------------------------------------------
def get_folder_id(name=None, parent=None, driveId=None, pageSize=50, credentials=None):
  if os.path.isfile(name):
    name = os.path.basename(name)
  try:
    service = build('drive', 'v3', credentials=credentials)
    results = service.files().list(
      pageSize=pageSize,
      q=f"'{parent}' in parents",
      fields='nextPageToken, files(id, name, parents)',
      includeItemsFromAllDrives=True,
      supportsAllDrives=True,
      corpora='drive',
      driveId=driveId
    ).execute()
    items = results.get('files', [])
    if len(items) == 0:
      return None
    for item in items:
      if name == item['name']:
        return item['id']
    return None
  except HttpError as error:
    print(f'An error occurred: {error}')

# -----------------------------------------------------------------------------
def upload_file(name=None, parent=None, driveId=None, chunksize=-1, credentials=None, retries=20):
  if not os.path.exists(name):
    raise IOError(f'The "{name}" file does not exist.')

  name = os.path.basename(name)

  file_metadata = {'name': name,
                   'parents': parent,
                   'driveID': driveId}

  service = build('drive', 'v3', credentials=credentials)
  new_file_id = None
  for retry in range(retries):
    try:
      media = MediaFileUpload(
        name,
        mimetype='application/octet-stream',
        chunksize=chunksize,
        resumable=True)
      results = service.files().create(
        body=file_metadata,
        media_body=media,
        supportsAllDrives=True,
      ).execute()
      new_file_id = results['id']

    except Exception as error:
      print(f'An error occurred (upload): {error}')

    # exit for loop if file has been uploaded
    if new_file_id is not None:
      break

    # sleep for retry_attempt * 2 * 60 s
    print()
    print('='*79)
    print(f'Retry {retry + 1} will start after {(retry + 1)*120} seconds.')
    print(time.asctime())
    print('='*79)
    print()
    time.sleep((retry + 1)*120)

  # stop if file was not uploaded
  if new_file_id is None:
    return None

  # otherwise, move file to final location
  for retry in range(retries):
    try:
      # parent defaults to "My Drive"
      results = service.files().update(
        fileId=new_file_id,
        removeParents=my_drive_id,
        addParents=parent,
        fields='name, id, parents',
        supportsAllDrives=True,
      ).execute()

      return new_file_id

    except Exception as error:
      print(f'An error occurred (move): {error}')

    # sleep for retry_attempt * 2 * 60 s
    print()
    print('='*79)
    print(f'Retry {retry + 1} will start after {(retry + 1)*120} seconds.')
    print(time.asctime())
    print('='*79)
    print()
    time.sleep((retry + 1)*120)

# -----------------------------------------------------------------------------
def main():

  parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('--credentials', help='JSON credentials for service account', required=True)
  parser.add_argument('--drive', help='Google Drive name', required=True)
  parser.add_argument('--installers-id', help='Google Drive id of parent folder', default=installers_id)
  parser.add_argument('--folder', help='Folder name for version')
  parser.add_argument('--subfolder', help='Subfolder in version folder')
  parser.add_argument('--file', help='File to upload')
  parser.add_argument('--retries', type=int, default=20, help='Number of retries to upload')
  parser.add_argument('--cleanup', action='store_true', help='Clean up "My Drive"')

  if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

  namespace = parser.parse_args(sys.argv[1:])

  # load credentials
  credentials = None
  if os.path.exists(namespace.credentials):
    credentials = ServiceAccountCredentials.from_json_keyfile_name(namespace.credentials)
  else:
    raise IOError(f'The credentials file, {namespace.credentials}, cannot be read.')

  # check drive
  drive_id = get_drive_id(name=namespace.drive, credentials=credentials)
  if drive_id is None:
    raise RuntimeError(f'The "{namespace.drive}" drive could not be found.')

  # clean up "My Drive" of service account
  if namespace.cleanup:
    service = build('drive', 'v3', credentials=credentials)
    results = service.files().list(
      pageSize=50,
      q=f"'{my_drive_id}' in parents",
      fields='nextPageToken, files(id, name, parents)',
      includeItemsFromAllDrives=True,
      supportsAllDrives=True,
      corpora='drive',
      driveId=drive_id
    ).execute()
    items = results.get('files', [])
    print(f'Cleaning up {len(items)} items')
    for item in items:
      print(item)
      service.files().delete(fileId=item['id']).execute()
    print()
    return 0

  # check folder (version)
  version_folder_id = get_folder_id(
    name=namespace.folder,
    parent=namespace.installers_id,
    driveId=drive_id,
    credentials=credentials)
  if version_folder_id is None:
    raise RuntimeError(f'The "{namespace.folder}" folder could not be found.')

  # check subfolder
  if namespace.subfolder is not None:
    subfolder_id = get_folder_id(
      name=namespace.subfolder,
      parent=version_folder_id,
      driveId=drive_id,
      credentials=credentials)
    if subfolder_id is None:
      raise RuntimeError(f'The "{namespace.subfolder}" subfolder could not be found.')
  else:
    subfolder_id = version_folder_id

  # check file
  new_file_id = None
  if namespace.file is not None:
    new_file_id = get_folder_id(
      name=namespace.file,
      parent=subfolder_id,
      driveId=drive_id,
      credentials=credentials)

  # id summary
  print('ID Summary')
  print('==========')
  for name, id in zip([namespace.drive, namespace.folder, namespace.subfolder],
                      [drive_id, version_folder_id, subfolder_id]):
    if name is not None and id is not None:
      print(f'{name:>20}: {id}')
  print()

  # upload file to subfolder if it does not already exist
  if new_file_id is None:
    new_file_id = upload_file(
      name=namespace.file,
      parent=subfolder_id,
      driveId=drive_id,
      credentials=credentials,
      retries=namespace.retries)
  else:
    print(f'{namespace.file} was not uploaded because it already exists.')
    print()

  # verify
  service = build('drive', 'v3', credentials=credentials)
  metadata = service.files().get(
    fileId=new_file_id,
    fields='name, id, parents, size, md5Checksum',
    supportsAllDrives=True,
  ).execute()

  # final summary
  print('Uploaded file')
  print('=============')
  for key in ['name', 'id', 'parents', 'size', 'md5Checksum']:
    print(f'{key:>20}: {metadata[key]}')
  print()

# =============================================================================
if __name__ == '__main__':
  sys.exit(main())
