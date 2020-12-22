#!/usr/bin/env python3

from ftplib import FTP
from termcolor import colored
import ftplib
import urllib.parse
import sys, os, stat
import hashlib
import csv
import re
import progressbar
import datetime

# =================================================================================================
#     Settings
# =================================================================================================

# CSV table listing all the FTP server/user combinations that we want to get files from.
# The table needs to have at least the following columns:
#   * active: "true" if that dataset shall be downloaded, or "false" if not (skip it)
#   * date: The date when this was downloaded. Needed to keep track.
#   * target: The local target directory to write files to. Recommended to start with the date as well.
#   * url: The url to download from. Can be just a host, or contain protocal and/or path as well.
#   * username, password: If needed to connect to the host. Empty strings are okay if not needed.
run_table = "runs.csv"

# If the main directory of the FTP server where we land only has a single subdirectory which
# contains all the data, there is no need to keep that subdirectory name in our local files as well.
# Let's just pretend in that case that this single subdirectory is our actual main directory.
descend_into_single_dir = True

# Log file where all FTP download attempts are logged
ftp_download_log_file = "ftp-download.log"

# File name search pattern (regex) for finding a file with md5 hashes of the files on the server.
md5_file_re = "(.*/)?md5(sum)?\.txt"

# Summary of all processed files
summary = {}

# =================================================================================================
#     Structures
# =================================================================================================

# Plain data class of information about files we are downloading.
# `status` is supposed to be a single character, which we use as follows:
#  - 'N': file does not exists locally
#  - 'S': file exists and has the correct size and md5 hash, so we can skip it
#  - 'R': file exists, but has a different size or md5 hash than on the server, so we need to replace/re-download it
#  - 'D': new download of non existing file
#  - 'E': some error occurred
class FileInfo:

    # Init function that expects the minimum data that we need to get the file.
    def __init__(self, url, user, remote_path, local_path):
        # Where the file is downloaded from, and remote file information
        self.url = url
        self.user = user
        self.remote_path = remote_path
        self.remote_size = None
        self.remote_md5_hash = None

        # Local file information, and overall status
        self.local_path = local_path
        self.local_size = None
        self.local_md5_hash = None
        self.status = None

# Translate our status codes into meaningful text.
def file_status(status):
    stati = {
        'N': "None Existing ",
        'S': "Skipped       ",
        'R': "Replaced      ",
        'D': "Downloaded    ",
        'E': "Errored       ",
    }
    return stati.get(status, "Invalid")

# =================================================================================================
#     General Helpers
# =================================================================================================

# We log each FTP download and its properties, to be sure we don't miss anything.
# Expects a FileInfo object.
def write_ftp_download_log(fileinfo):
    with open( ftp_download_log_file, "a") as logfile:
        now = datetime.datetime.now().strftime("%Y-%m-%d\t%H:%M:%S")
        logfile.write(
            now +
            "\t" + str(fileinfo.url) +
            "\t" + str(fileinfo.user) +
            "\t" + str(fileinfo.status) +
            "\t" + str(fileinfo.local_md5_hash) +
            "\t" + str(fileinfo.local_size) +
            "\t" + str(fileinfo.local_path) + "\n"
        )

# Given a file as produced by the unix `md5sum` command, return a dict from file names to
# their hashes. The file is typically named `md5.txt`, and its expected file format consists
# of rows of the format `<hash>  <filename>`.
def get_md5_hash_dict(md5_file):
    hashdict = {}
    with open(md5_file) as fp:
        for line in fp:
            sl = [item for item in re.split("\s+", line) if item]
            if len(sl) == 0:
                continue
            elif len(sl) > 2:
                raise Exception("md5 file " + md5_file + " has a line with more than 2 columns.")
            assert len(sl) == 2
            if len(sl[0]) != 32:
                raise Exception("md5 file " + md5_file + " has a line with an invalid md5 hash.")

            if sl[1] in hashdict:
                raise Exception("md5 file " + md5_file + " has multiple entries for file " + sl[1])
            else:
                hashdict[sl[1]] = sl[0]
    return hashdict

# Compute the md5 hash of a given local file, efficiently (hopefully?! not all to sure about large
# file handling in python...) by using blocks of data instead of reading the (potentially huge)
# files all at once in to memory.
def get_file_md5(filename, blocksize=65536):
    if not os.path.isfile(filename):
        raise Exception("Cannot compute md5 hash for path \"" + filename + "\"")
    md5_hash = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(blocksize), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()

# Check local file properties against remote: file size and md5 hash have to match.
# Expects to be given a FileInfo object.
# Return True if they match (file is good), or False if either is wrong (file needs to be
# downloaded [again]). Also, fill in the values in the FileInfo while doing so.
def get_and_check_file_properties( fileinfo ):
    if not os.path.exists(fileinfo.local_path):
        raise Exception("Local path \"" + fileinfo.local_path + "\" does not exists.")
    if not os.path.isfile(fileinfo.local_path):
        raise Exception("Local path \"" + fileinfo.local_path + "\" exists, but is not a file.")

    # Get and check file size.
    fileinfo.local_size = os.stat(fileinfo.local_path).st_size
    if fileinfo.local_size != fileinfo.remote_size:
        print(colored(
            "Local file \"" + fileinfo.local_path + "\" exists, but has size " +
            str(fileinfo.local_size) + " instead of remote file size " + str(fileinfo.remote_size) + ".",
            "yellow"
        ))
        return False

    # Compute local file md5 hash.
    if fileinfo.remote_md5_hash:
        fileinfo.local_md5_hash = get_file_md5( fileinfo.local_path )
    else:
        fileinfo.local_md5_hash = None

    # Check that we got the correct md5 hash.
    if fileinfo.local_md5_hash != fileinfo.remote_md5_hash:
        print(colored(
            "Local file \"" + fileinfo.local_path + "\" exists, " +
            "but its md5 hash does not match the remote md5 hash.",
            "yellow"
        ))
        return False

    # If we are here, everything is good.
    assert fileinfo.local_size == fileinfo.remote_size
    assert fileinfo.local_md5_hash == fileinfo.remote_md5_hash
    return True

# =================================================================================================
#     FTP Helpers
# =================================================================================================

# Test whether a name is a file or not - so, probably a directory? FTP is messy...
# But that's the best we can do with the limitations of that protocol.
def ftp_is_file(ftp, name):
    try:
        fs = ftp.size(name)
    except:
        return False
    return fs is not None

# Get all names (files and dirs) in the given or the current working directory.
def ftp_get_list(ftp, dir=None):
    try:
        if dir:
            names = ftp.nlst(dir)
        else:
            names = ftp.nlst()
    except ftplib.error_perm as resp:
        # No files found / Can't check for file existence
        if str(resp).startswith( "550" ):
            return []
        else:
            raise
    return names

# Get all file names in the current working directory.
def ftp_get_files(ftp, dir=None):
    names = ftp_get_list(ftp, dir)
    return [ n for n in names if ftp_is_file(ftp, n) ]

# Get all directory names in the current working directory. Hopefully.
def ftp_get_dirs(ftp, dir=None):
    names = ftp_get_list(ftp, dir)
    return [ n for n in names if not ftp_is_file(ftp, n) ]

# =================================================================================================
#     FTP Download File
# =================================================================================================

# Download a specific file, and fill in the respective FileInfo data.
# This is the inner function, that we use to keep the control flow simple.
# See ftp_download_file() for the actual function to be called.
def ftp_download_file_inner(ftp, fileinfo):
    # Init the (expected) remote file size.
    fileinfo.remote_size = ftp.size(fileinfo.remote_path)
    if fileinfo.remote_size is None:
        raise Exception(
            "Cannot work with a server that does not support to retreive file sizes. " +
            "Feel free however to refactor this script accordingly."
        )

    # Check that we do not overwrite files accidentally, that is, only download again if the file
    # does not match its expectations. If all is good, we can skip the file.
    if os.path.exists(fileinfo.local_path):
        if get_and_check_file_properties(fileinfo):
            print(colored(
                "Local file \"" + fileinfo.local_path + "\" exists and is good. Skipping.", "green"
            ))
            fileinfo.status='S'
            return
        else:
            print( "Will download the file again." )
            fileinfo.status='R'

    # Make the target dir if necessary.
    if not os.path.exists(os.path.dirname( fileinfo.local_path )):
        os.makedirs(os.path.dirname( fileinfo.local_path ))

    # Report progress while downloading. We have gigabytes of data, so that is important.
    print(colored("\nDownloading \"" + fileinfo.local_path + "\"...", "blue"), flush=True)
    pbar = progressbar.ProgressBar( max_value = (
        fileinfo.remote_size if fileinfo.remote_size is not None else progressbar.UnknownLength
    ))
    pbar.start()

    # Open the file locally, and define a callback that writes to that file while reporting progress.
    filehandle = open(fileinfo.local_path, 'wb')
    def file_write_callback(data):
        filehandle.write(data)
        nonlocal pbar
        pbar += len(data)

    # Go go gadget!
    try:
        ftp.retrbinary("RETR " + fileinfo.remote_path, file_write_callback)
    except Exception as ex:
        print(colored("Error downloading file: " + str(ex), "red"))
        fileinfo.status='E'
    pbar.finish()
    filehandle.close()

    # Check that we got the correct size, and the correct md5 hash,
    # and if so, make it read-only, and return.
    if get_and_check_file_properties(fileinfo):
        print(colored("Done. File passed checks.", "green"))
        if fileinfo.status != 'R':
            fileinfo.status='D'
        os.chmod( fileinfo.local_path, stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH )
    else:
        print(colored("Error downloading file!", "red"))
        fileinfo.status='E'

# Download a specific file, fill in the respective FileInfo data,
# and add the file to the log and summary.
def ftp_download_file( ftp, fileinfo ):
    ftp_download_file_inner( ftp, fileinfo )
    write_ftp_download_log(fileinfo)

    # Summary of all downloads. Cumbersome, because Python...
    global summary
    if fileinfo.status in summary:
        summary[fileinfo.status] += 1
    else:
        summary[fileinfo.status] = 1

# =================================================================================================
#     FTP Download All
# =================================================================================================

# Download all files from an FTP server into a target directory.
def ftp_download_all(url, user, passwd, target_dir):
    # Get host name components. We allow with or without protocol, for simplicity.
    # But that means that we have to do a bit more parsing here...
    parsed_url = urllib.parse.urlparse(url)
    host = parsed_url.netloc
    path = parsed_url.path
    if not host:
        host = url.split('/', 1)[0]
        path = url.split('/', 1)[1] if len(url.split('/', 1)) > 1 else ""

    # User output.
    print(colored(
        "===================================================================================",
        "blue"
    ))
    print(colored("Connecting to host " + host + ( " as user " + user if user else "" ), "blue"))
    print()

    # Check target.
    if os.path.exists(target_dir):
        if not os.path.isdir(target_dir):
            raise Exception("Local path", target_dir, "exists, but is not a directory.")
    else:
        os.mkdir(target_dir)

    # Connect to FTP server. If the remote host is not available, for example becaue the sequencing
    # center already deleted the data, we simply skip it with a warning.
    try:
        ftp = FTP( host )
        if user or passwd:
            ftp.login( user=user, passwd=passwd )
        else:
            ftp.login()
        if path:
            ftp.cwd( path )
    except:
        print(colored("Cannot connect to host, skipping.", "red"))
        return

    # We work through all directories on the server, and store them in a queue
    # that we process dir by dir, pushing new (sub)dirs as we discover them.
    # Initialize with the current dir (after login) of the server.
    queue = ftp_get_dirs(ftp)

    # If there is only a single dir on the server, we might want to descend into that, and use
    # that as our new main directory, to keep our local file structure a bit easier.
    # That only works if there is just that directory, and no files.
    # In that case, change to that dir, and load its subdirectories again.
    if descend_into_single_dir and len(queue) == 1 and len(ftp_get_files(ftp)) == 0:
        maindir = queue.pop(0)
        print("Descending into directory " + maindir)
        print()
        ftp.cwd(maindir)
        queue = ftp_get_dirs(ftp)

    # We of course also want to download files from the current directory (either the main, or
    # the one we descended into). In fact, make this the first directory to process.
    if not "." in queue:
        queue.insert(0, ".")

    while len(queue) > 0:
        remote_dir = queue.pop(0)
        if remote_dir == ".." or remote_dir.startswith("./") or remote_dir.endswith("/.") or remote_dir.endswith("/.."):
            continue
        print(colored(
            "-----------------------------------------------------------------------------------",
            "blue"
        ))
        print(colored("Processing directory " + remote_dir, "blue"))

        # Add all subdirs of the current one to the queue.
        for f in ftp_get_dirs( ftp, remote_dir ):
            queue.append( f )

        # Get list of all files in the dir.
        files = ftp_get_files( ftp, remote_dir )
        # print("Files:", files)

        # If there is an md5 hash file in that directory for the files in there, get that first,
        # so that we can check hashes for each downloaded file.
        md5_hashes = {}
        if md5_file_re is not None:
            # See if there is a file in the list that fits our regular expression.
            md5_regex = re.compile(md5_file_re)
            md5_match_list = list(filter(md5_regex.match, files))
            if len(md5_match_list) > 1:
                raise Exception("Multiple md5 hash files found. Refine your regex to find the file.")
            elif len(md5_match_list) == 1:
                # Get the md5 txt file, as produced by the unix `md5sum` command.
                # First, prepare is properties.
                md5_remote_file = md5_match_list[0]
                print("Using md5 hash check file", md5_remote_file)
                md5_local_file = os.path.join( target_dir, md5_remote_file )
                md5_fileinfo = FileInfo( url, user, md5_remote_file, md5_local_file )

                # Now, download it, and remove it from the file list, so that we don't download
                # it again. Then, extract a dict of all hashes for the files.
                ftp_download_file( ftp, md5_fileinfo )
                files.remove(md5_remote_file)
                md5_hashes = get_md5_hash_dict(md5_local_file)

        # Download them all!
        print()
        for f in files:
            # Initialize a FileInfo where we capture all info as we process that file.
            fileinfo = FileInfo( url, user, f, os.path.join( target_dir, f ))

            # See if there is an md5 hash that we can use to check the file contents.
            fbn = os.path.basename(fileinfo.remote_path)
            if fbn in md5_hashes:
                fileinfo.remote_md5_hash = md5_hashes[fbn]
            else:
                fileinfo.remote_md5_hash = None

            # Download the file, do all checks, and write a log line about it.
            ftp_download_file( ftp, fileinfo )
        print()

    # We are polite, and close the connection respectfully. Bye, host.
    ftp.quit()

# =================================================================================================
#     Main function to process the download table
# =================================================================================================

# Some nice error reporting on a common problem...
if progressbar.__author__ != "Rick van Hattem (Wolph)":
    raise Exception("You are using `progressbar`, instead of `progressbar2`. Please update.")

if __name__ == "__main__":
    with open( run_table ) as csvfile:
        runreader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
        exp_cols = ['active', 'date', 'target', 'url', 'username', 'password']
        if(not all(x in runreader.fieldnames for x in exp_cols)):
            raise Exception(
                "Run table " + run_table + " does not contain all needed fields [" +
                ", ".join(exp_cols) + "], but instead contains [" + ", ".join(runreader.fieldnames) + "]"
            )

        for row in runreader:
            if row["active"] == "true":
                ftp_download_all(
                    row["url"], row["username"], row["password"], row["target"]
                )
            else:
                print("Skipping inactive run", row["url"], "from", row["date"])

    print(colored(
        "===================================================================================\n",
        "blue"
    ))
    print("Summary:")
    for key, val in summary.items():
        print(file_status(key) + ": " + str(val))
    if 'E' in summary:
        print(colored("There were errors in the downloads! Please check!", "red"))
