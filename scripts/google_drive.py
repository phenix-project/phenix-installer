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
  try:
    service = build('drive', 'v3', credentials=credentials)
    results = service.files().list(
      pageSize=pageSize,
      q="'{}' in parents".format(parent),
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
def upload_file(name=None, parent=None, driveId=None, chunksize=-1, credentials=None, retries=5):
  if not os.path.exists(name):
    raise IOError('The "{}" file does not exist.'.format(name))

  name = os.path.basename(name)

  file_metadata = {'name': name,
                   'parents': parent,
                   'driveID': driveId}

  for retry in range(retries):
    try:
      media = MediaFileUpload(
        name,
        mimetype='application/octet-stream',
        chunksize=chunksize,
        resumable=True)
      service = build('drive', 'v3', credentials=credentials)
      results = service.files().create(
        body=file_metadata,
        media_body=media,
        supportsAllDrives=True,
      ).execute()
      new_file_id = results['id']

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
      print(f'An error occurred: {error}')

    # sleep for retry_attempt * 2 * 60 s
    time.sleep(retry*120)

# -----------------------------------------------------------------------------
def main():

  parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('--credentials', help='JSON credentials for service account')
  parser.add_argument('--drive', help='Google Drive name')
  parser.add_argument('--folder', help='Folder name for version')
  parser.add_argument('--subfolder', help='Subfolder in version folder')
  parser.add_argument('--file', help='File to upload')
  parser.add_argument('--retries', type=int, default=5, help='Number of retries to upload')

  if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

  namespace = parser.parse_args(sys.argv[1:])

  # load credentials
  credentials = None
  if os.path.exists(namespace.credentials):
    credentials = ServiceAccountCredentials.from_json_keyfile_name(namespace.credentials)
  else:
    raise IOError('The credentials file, {}, cannot be read.'.format(namespace.credentials))

  # check drive
  drive_id = get_drive_id(name=namespace.drive, credentials=credentials)
  if drive_id is None:
    raise RuntimeError('The "{}" drive could not be found.'.format(namespace.drive))

  # check folder (version)
  version_folder_id = get_folder_id(
    name=namespace.folder,
    parent=installers_id,
    driveId=drive_id,
    credentials=credentials)
  if version_folder_id is None:
    raise RuntimeError('The "{}" folder could not be found.'.format(namespace.folder))

  # check subfolder
  subfolder_id = get_folder_id(
    name=namespace.subfolder,
    parent=version_folder_id,
    driveId=drive_id,
    credentials=credentials)
  if subfolder_id is None:
    raise RuntimeError('The "{}" subfolder could not be found.'.format(namespace.subfolder))

  # id summary
  print('ID Summary')
  print('==========')
  for name, id in zip([namespace.drive, namespace.folder, namespace.subfolder],
                      [drive_id, version_folder_id, subfolder_id]):
    print('{:>20}: {}'.format(name, id))
  print()

  # upload file to subfolder
  new_file_id = upload_file(
    name=namespace.file,
    parent=subfolder_id,
    driveId=drive_id,
    credentials=credentials,
    retries=namespace.retries)

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
    print('{:>20}: {}'.format(key, metadata[key]))
  print()

  # clean up "My Drive" of service account
  results = service.files().list(
    pageSize=50,
    q="'{}' in parents".format(my_drive_id),
    fields='nextPageToken, files(id, name, parents)',
    includeItemsFromAllDrives=True,
    supportsAllDrives=True,
    corpora='drive',
    driveId=drive_id
  ).execute()
  items = results.get('files', [])
  for item in items:
    print(item)
    # service.files().delete(fileId=item['id']).execute()

# =============================================================================
if __name__ == '__main__':
  sys.exit(main())
