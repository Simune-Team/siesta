#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script for releasing SIESTA on GitLab.
Other uses: release deletion and release description replacement.

For more information, execute
    $ python3 release.py -h
"""

__author__ = "J.M. Escartin"
__email__ = "jm.escartin@icn2.cat"
__version__ = "0.8"
__date__ = "10th February 2021"
__license__ = "GPLv3"

# TODO:
# a) Check that signature files uploaded to generic packages are correct.
# b) Make strict test of contents of list vs. dict.
# c) Check that signatures are valid.

################################################################################

import os
import sys
import requests

################################################################################

# Check returned status code.
def check_status_code(response, name, check=200):

    if not response.status_code == check:
        if r.content:
            print(r.json())
        sys.exit('Request "%s" returned status code %s.\nAborting.' \
            % (name, response.status_code))
    elif args.verbose:
        print('Request "%s" returned status code %s.' \
            % (name, response.status_code))


# Compute hashes for large files.
# http://stackoverflow.com/questions/1131220/get-md5-hash-of-a-files-without-open-it-in-python
def file_hash(fname, alg, block_size=2**20):

    import hashlib

    if alg.lower() == 'md5':
        h = hashlib.md5()
    elif alg.lower() == 'sha256':
        h = hashlib.sha256()
    elif alg.lower() == 'sha512':
        h = hashlib.sha512()
    else:
        sys.exit('Cannot find algorithm "%s".\nAborting.' % alg)

    with open(fname, 'rb') as f:
        while True:
            data = f.read(block_size)
            if not data:
                break
            h.update(data)
    return h.hexdigest()


# Check checksums file using local paths.
def check_checksums(cfname='checksums.txt'):

    import re

    with open(cfname, 'r') as cf:
        # Cycle over the lines of the checksums file.
        for l in cf.readlines():
            # Read info from that line.
            alg, f, rh = re.split('(|)= ', l)
            # Compute corresponding hash.
            ch = file_hash(f, alg)
            # Compare both digests.
            if fh != ch:
                sys.exit('Error checking checksums in file %s' % cfname
                    + '%s checksum failed for file %s:\n' % (alg, f)
                    + 'Expected digest: %s.\n'          % rh
                    + 'Computed digest: %s.\nAborting.' % ch)

    # No error found when processing the complete file.
    print('All checksums in file %s verified!' % cfname)

################################################################################

###### 0. Parse input ######

import argparse

### Create parser.
input_parser = argparse.ArgumentParser( \
    description=(
"""
Script for releasing SIESTA on GitLab.

Main assumptions:
  1) Tag <version> already exists in the repository.
  2) The following release files already exist in the local directory:
      - siesta tarball,
      - manual files (2 pdfs),
      - checksums.txt,
      - their signatures (4 .asc).
  3) The tarball with all the signatures does not exist in the local directory.
  4) There are no other "packages" for this version yet.

The script may also be used to delete a release from GitLab
or to replace the description (text) of an existing release
(in both cases, related generic packages or uploads will remain unchanged).
"""), \
    formatter_class=argparse.RawTextHelpFormatter)

### Add to the parser all the input arguments:

# Token.
input_parser.add_argument('token', type=str, action='store', \
    help="A GitLab API-enabled token.")

# Repository.
input_parser.add_argument('repo', type=str, action='store', \
    help='Repository including namespace\n(e.g., "siesta-project/siesta").')

# Version (= release name).
input_parser.add_argument('version', type=str, action='store', \
    help='Release version.\nThe repository must contain tag <version>.')

# Release description or deletion option (they are mutually exclusive).
text_group = input_parser.add_mutually_exclusive_group(required=True)
text_group.add_argument('-f', type=str, \
    action='store', dest='file', \
    help=('File that contains the description of the release.\n'
    'Markdown allowed.'))
text_group.add_argument('-m', type=str, \
    action='store', dest='message',
    help='Description of the release. Markdown allowed.')
text_group.add_argument('--delete', '-d', \
    action='store_true',
    help=('Delete an existing release. Unlike the web UI,\n'
    'this preserves the git tag associated to the release.'))

# Method for storing files on Gitlab: generic packages or uploads.
input_parser.add_argument('--no-generic-packages', \
    action='store_true', dest='no_generic',
    help=('Do not store files as generic packages\n'
    '(use unmanageable uploads instead).\n'
    'Generic packages are the default,\n'
    'but they require "X.Y.Z" version names\n'
    '(3 fully numeric components separated by points).'))

# Fix release description option.
input_parser.add_argument('--fix-description', \
    action='store_true', dest='fix_description', \
    help=('Fix an existing release by uploading a new description.\n'
    'All the release links and the JSON release evidence\n'
    'will be preserved.\n'
    'This option is incompatible with the --delete option.\n'))

# Verbose output.
input_parser.add_argument('--verbose', '-v', action='store_true', \
    help="Enable verbose output")

### Parse arguments
args = input_parser.parse_args()
if args.verbose:
    print(args)

### Process parsed arguments.

# My GitLab API token.
token = args.token

# Project name:
project_name = args.repo

# Version (= release name)
version = args.version
# and version with leading 'v' removed (if that is the case).
if version[0] == 'v' and version[1].isdigit():
    nov_version = version[1:]
else:
    nov_version = version

# Delete release?
if args.delete:
    delete_release = True
else:
    delete_release = False
    # Release will be created, so text needed:
    if args.file:
        with open(args.file, 'r') as f:
            text = f.read()
    else:
        # If you try to be smart and use escape sequences in the message
        # you also need to remember to use bash ANSI-C quoting ($'string')
        # so that the escape sequences are interpreted by the shell, see
        # https://www.gnu.org/software/bash/manual/html_node/ANSI_002dC-Quoting.html
        # Remember also about line breaks in Markdown - $'Foo\nBar' is futile.
        text = args.message
    if args.verbose:
        print(text)

# Generic packages?
generic = not args.no_generic

# Fix description?
fix_description = args.fix_description
if fix_description:
    if args.no_generic:
        print('Warning: --no-generic-packages option not used\n'
            'since --fix-description option provided.')
    if args.delete:
        sys.exit('Error: options --fix-description and --delete'
            'are incompatible.\nAborting.')

###### 1. Obtain projectID ######

# Header with only the token.
token_header = { 'PRIVATE-TOKEN': token }
# Headers for JSON calls.
json_headers = { 'PRIVATE-TOKEN': token,
    'Content-Type': 'application/json' }

# URL-encode the project name.
import urllib
safe_project_name = urllib.parse.quote_plus(project_name)

# Base URL.
base_base_url = 'https://gitlab.com/api/v4/projects/%s'
base_url = base_base_url % safe_project_name

# Obtain project information.
r = requests.get( \
    url = base_url, \
    headers = token_header)

# Check status code.
check_status_code(r, 'project')

# Extract project ID.
projectID = str(r.json()['id'])
if args.verbose:
    print('projectID: %s' % projectID)

# Redefine base url to use project ID (safer than name).
base_url = base_base_url % projectID

###### 2. Simpler options ######

  #### 2.1 Delete release (if requested) ####

if delete_release:
    r = requests.delete( \
        url = base_url + '/releases/' + version, \
        headers = token_header)

    # Check status code.
    check_status_code(r, 'delete release')

    # Exit successfully.
    print('The script thinks this is a successful release deletion.')
    print('Check the following link:')
    print('https://gitlab.com/%s/-/releases/' % project_name )
    sys.exit()

  #### 2.2 Fix description (if requested) ####

if fix_description:

    # Create JSON payload dictionary with only one key.
    payload = { 'description' : text }

    if args.verbose:
        import json
        print(json.dumps(payload, indent=4))

    # Update the release.
    # API: https://docs.gitlab.com/ee/api/releases/#update-a-release
    r = requests.put( \
        url = base_url + '/releases/' + version, \
        headers = token_header, \
        json = payload)

    if args.verbose:
        print(r.json())

    # Check returned status code.
    check_status_code(r, 'update release PUT with json payload')

    # Exit successfully.
    print('The script thinks this is a successful release edition.')
    print('Check the following link:')
    print('https://gitlab.com/%s/-/releases/' % project_name )
    sys.exit()


###### 3. Check and prepare local files ######

# Files:
# (Note: release links will have the order of release_files).
siesta_tarball = 'siesta-%s.tar.gz' % nov_version
manuals = [ 'siesta.pdf', 'tbtrans.pdf']
#    'siesta-screen.pdf', 'tbtrans-screen.pdf' ]
checksums = 'checksums.txt'
release_files = [ siesta_tarball ] + manuals + [ checksums ]
signatures = [ g + '.asc' for g in release_files ]

# Manual names
manual_name = { \
    'siesta.pdf'         : 'Siesta manual (pdf)', \
    'tbtrans.pdf'        : 'Tbtrans manual (pdf)' }
#    'siesta-screen.pdf'  : 'Siesta manual (screen-friendly pdf)', \
#    'tbtrans-screen.pdf' : 'Tbtrans manual (screen-friendly pdf)' }

# Check that all files and their signatures exist in cwd.
for f in release_files + signatures:
    if not os.path.isfile(f):
        sys.exit('Cannot find file "%s".\nAborting.' % f)

# Check checksums file,
check_checksums

# and compute its checksum.
checksums_checksum = file_hash(checksums, 'sha512')

# Check that a tarball with all signatures does not exist,
signatures_tarball = 'signatures.tar.gz'
if os.path.exists(signatures_tarball):
    sys.exit('File "%s" already exists.\nAborting' % signatures_tarball)

# create it,
import tarfile
with tarfile.open(signatures_tarball, 'w:gz') as tgz:
    for f in signatures:
        tgz.add(f)

# add it to the list of release files,
release_files.append(signatures_tarball)

# and compute its checksum.
signatures_checksum = file_hash(signatures_tarball, 'sha512')

###### 4. Upload files ######

# Create dictionary of files to be linked from the release.
file_urls = {}

  #### 4.1. Generic packages ####

if generic:

    ## 4.1.1. Upload files as generic packages ##

    for f in reversed(release_files):

        # Determine generic "package":
        if f == siesta_tarball:
            package = 'tarball'
        elif f in manuals:
            package = 'manuals'
        elif f == checksums:
            package = 'checksums'
        else: # signatures tarball
            package = 'signatures_tarball'

        # Loop over the actual file and the signature file:
        for ext in [ '', '.asc' ]:
            f1 = f + ext

            # The signatures tarball is not signed.
            if f1 == signatures_tarball + '.asc':
                continue

            # Upload the file.
            # As of 2021-01-31, only X.Y.Z versions are supported, see
            # https://docs.gitlab.com/ee/user/packages/generic_packages/#publish-a-package-file
            # This restriction may be lifted in the future, see GitLab issue
            # https://gitlab.com/gitlab-org/gitlab/-/issues/255234 .
            with open(f1, 'rb') as data:
                r = requests.put( \
                    url = base_url + \
                        '/packages/generic/%s/%s/%s' % (package, nov_version, f1), \
                    headers = token_header, \
                    data = data)

            # Check status code.
            check_status_code(r, 'upload file %s' % f1, 201) # 201: Created

    ## 4.1.2. Obtain links of relevant package files ##

    # Obtain list of "packages" in this version.
    r = requests.get( \
        url = base_url + '/packages', \
        headers = token_header)

    # Check status code.
    check_status_code(r, 'list of packages')

    # Filter list of dicts: we are only interested in the packages of this version.
    list_of_packages = [ d for d in r.json() if d['version'] == nov_version ]

    # Populate dictionary of links by looping over files in these packages.
    for d in list_of_packages:

        # Determine package ID.
        packageID = d['id']

        # Request list of files in the package.
        r = requests.get( \
            url = base_url + '/packages/%s/package_files' % packageID,
            headers = token_header)

        # Check returned status code.
        check_status_code(r, 'list of files of %s package %s (version %s)' \
            % (d['package_type'], d['name'], d['version']))

        # Add files to dictionary with their download links.
        for f in r.json():
            filename = f['file_name']
            if filename not in signatures:

                # Check for multiple versions of the same file within
                # the packages for this version.
                if filename in file_urls:
                    sys.exit('Packages contains multiple versions of file\n'
                        + '"%s" for version %s.\nAborting.' % (filename, version))

                # Build link to package file.
                myurl = 'https://gitlab.com/%s/-/package_files/%s/download' \
                    % (project_name, f['id'])

                # Add to dictionary.
                file_urls[filename] = myurl

else:

  #### 4.2. Unmanageable uploads ####

    for f in release_files:

        if f in file_urls:
            sys.exit('File "%s" already uploaded.\nAborting.' % f)

        # Upload the file.
        # API: https://docs.gitlab.com/ee/api/projects.html#upload-a-file
        with open(f, 'rb') as myfile:
            r = requests.post( \
                url = base_url + '/uploads', \
                headers = token_header, \
                files = dict(file=myfile))

        # Check status code.
        check_status_code(r, 'upload file %s' % f, 201) # 201: Created

        # Build link to uploaded file.
        myurl = 'https://gitlab.com' + r.json()['full_path']

        # Add files to dictionary with their download links.
        file_urls[f] = myurl

###### 5. In both cases: consistency check ######

# Soft consistency check:     # TODO: strict check
if len(file_urls) != len(release_files):
    print(release_files)
    print(file_urls)
    sys.exit('Consistency check failed.\nAborting.')

###### 6. Build list of links for release ######

# It is an ordered list:
list_of_links = []

for f in reversed(release_files): # Based on observation

    # Obtain download link.
    link_url = file_urls[f]

    # Select whether this file is the "package" or "other", and
    # build its file name.
    if f == siesta_tarball:
        link_type = 'package'
        link_name = f          # Keeps the tarball name
    elif f in manuals:
        link_type = 'other'
        link_name = manual_name[f]
    elif f == checksums:
        link_type = 'other'
        link_name = 'Checksums'
    else: # Signatures tarball
        link_type = 'other'
        link_name = 'Signatures tarball'

    # Add link to list.
    list_of_links.append( { \
        'name': link_name, \
        'url': link_url, \
        'link_type': link_type, \
        'filepath': '/' + f } )

###### 7. Create the release ######

# Build JSON payload.
# Note no 'ref': version tag <version> must already exist!
payload = { 'id': projectID, \
    'tag_name': version, \
    'assets': { 'links': list_of_links } }
if text:
    payload['description'] = text

if args.verbose:
    import json
    print(json.dumps(payload, indent=4))

# Headers for JSON calls.
json_headers = { 'PRIVATE-TOKEN': token,
    'Content-Type': 'application/json' }

# Call to create release.
# API: https://docs.gitlab.com/ee/api/releases/#create-a-release
r = requests.post(base_url + '/releases', \
    headers = json_headers, \
    json = payload)

if args.verbose:
    print(r.json())

# Check returned status code.
check_status_code(r, 'release post with json payload', 201) # 201: Created

###### 8. Download files and check checksums ######

import time
import pathlib
import cgi

# Create auxiliary directory named after current Unix time.
auxdir = pathlib.Path.cwd() / str(int(time.time()))
pathlib.Path.mkdir(auxdir)

# Change to this auxiliary directory.
os.chdir(auxdir)

# Loop over the returned link-type assets:
for link in r.json()['assets']['links']:
    link_url = link['url']

    # Download the link.
    rd = requests.get(link_url)
    check_status_code(rd, 'download %s' % link_url)

    # Determine the name of the downloaded file according to GitLab's headers.
    # (This should work, especially if the file name only contains
    # ASCII characters).
    filename = cgi.parse_header( \
        rd.headers.get('content-disposition'))[1]['filename']

    # Write the downloaded data to the file.
    with open(filename, "wb") as f:
        f.write(rd.content)

# Check that the checksums file has not changed.
downloaded_checksums_checksum = file_hash(checksums, 'sha512')
if downloaded_checksums_checksum != checksums_checksum:
    sys.exit( \
        'Error checking checksums of downloaded "%s" file:\n' % checksums
        + 'Expected digest: %s\n' % checksums_checksum
        + 'Computed digest: %s\n' % downloaded_checksums_checksum)

# Check that the signatures tarball also has not changed.
downloaded_signatures_checksum = file_hash(signatures_tarball, 'sha512')
if downloaded_signatures_checksum != signatures_checksum:
    sys.exit( \
        'Error checking checksums of downloaded "%s" file:\n' % signatures_tarball
        + 'Expected digest: %s\n' % signatures_checksum
        + 'Computed digest: %s\n' % downloaded_signatures_checksum)

# Check the contents of the checksums file (covers the rest of the files).
check_checksums

###### 9. Final messages ######

print('The script thinks this is a successful release.')
print('Check the following link(s):')
print('https://gitlab.com/%s/-/releases/%s' % (project_name, version) )
if generic:
    print('https://gitlab.com/%s/-/packages/' % project_name)

