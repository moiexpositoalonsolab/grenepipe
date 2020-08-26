#!/usr/bin/env python3

import sys, os, stat
import csv
from ftplib import FTP
import progressbar
import datetime

# =================================================================================================
#     Settings
# =================================================================================================

# CSV table listing all the FTP server/user combinations that we want to get files from.
# The table needs to have at least the following columns: host,username,password, as well as a column
# that specifies the target directory to which the files from that server/user are downloaded to.
# The name of that column can be set with run_table_target_col
run_table = "runs.csv"
run_table_target_col = "submission_date"

# Log file where all FTP download attempts are logged
ftp_download_log_file = "ftp-download.log"

# Summary of all processed files
summary = {}

# =================================================================================================
#     General Helpers
# =================================================================================================

# We log each FTP download, to be sure we don't miss anything.
# Status is supposed to be a single character, which we use as follows:
#  - 'S': file exists and has the correct size, so we skip it
#  - 'R': file exists, but has a different size than on the server, so we replace it
#  - 'D': new download of non existing file
#  - 'E': some error occurred
def write_ftp_download_log(host, user, status, size, file):
    with open( ftp_download_log_file, "a") as logfile:
        now = datetime.datetime.now().strftime("%Y-%m-%d\t%H:%M:%S")
        logfile.write(
            now + "\t" + str(host) + "\t" + str(user) + "\t" +
            str(status) + "\t" + str(size) + "\t" + str(file) + "\n"
        )

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
    except( ftplib.error_perm, resp ):
        if str(resp) == "550 No files found":
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
#     FTP Download
# =================================================================================================

# Download a specific file. We return a status code and the remove file size. This is quick and dirty.
# If we want to handle these information properly, and want more logging etc in the future,
# refactor this to use a proper data structure that captures all details of a file download.
def ftp_download_file(ftp, remote_file, local_file):
    # Init our return values, which is status "Downloaded" (all good), and the file size.
    status='D'
    rsize = ftp.size(remote_file)
    if rsize is None:
        raise Exception("Cannot work with a server that does not support to retreive file sizes.")

    # Check that we do not overwrite files accidentally.
    if os.path.exists(local_file):
        if not os.path.isfile(local_file):
            raise Exception("Local path \"" + local_file + "\" exists, but is not a file.")

        # If the file sizes are identical, we can skip. If not, we download again.
        lsize = os.stat(local_file).st_size
        if rsize != lsize:
            print(
                "Local file \"" + local_file + "\" exists, but has size", str(lsize),
                "instead of remote file size", str(rsize) + ".", flush=True
            )
            status='R'
        else:
            print("Local file \"" + local_file + "\" exists. Skipping.")
            status='S'
            return [status, rsize]

    # Make the target dir if necessary.
    if not os.path.exists(os.path.dirname( local_file )):
        os.mkdir(os.path.dirname( local_file ))

    # Report progress while downloading. We have gigabytes of data, so that is important.
    print("\nDownloading \"" + remote_file + "\"", flush=True)
    pbar = progressbar.ProgressBar(max_value=(rsize if rsize is not None else progressbar.UnknownLength))
    pbar.start()

    # Open the file locally, and define a callback that writes to that file while reporting progress.
    file = open(local_file, 'wb')
    def file_write(data):
        file.write(data)
        nonlocal pbar
        pbar += len(data)

    # Go go gadget
    try:
        ftp.retrbinary("RETR " + remote_file, file_write)
    except ex:
        print("Error downloading file:", str(ex))
        status='E'
    pbar.finish()
    file.close()

    # Check that we got the correct size, if so, make it read-only, and return.
    lsize = os.stat(local_file).st_size
    if rsize != lsize:
        print("Error downloading file: Downloaded size", str(lsize), "does not match remote size", str(rsize))
        status='E'
    else:
        os.chmod( local_file, stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH )
    return [status, rsize]

# Download all files from an FTP server into a target directory.
def ftp_download_all(host, user, passwd, target_dir):
    # Check target.
    if os.path.exists(target_dir):
        if not os.path.isdir(target_dir):
            raise Exception("Local path", target_dir, "exists, but is not a directory.")
    else:
        os.mkdir(target_dir)

    # Connect to FTP server.
    ftp = FTP( host )
    ftp.login( user=user, passwd=passwd )

    # We work ourselves through all directories on the server, and store them in a queue
    # that we process dir by dir, pushing new (sub)dirs as we go.
    # Initialize with the current dir (after login) of the server.
    queue = ftp_get_dirs(ftp)
    while len(queue) > 0:
        dir = queue.pop(0)
        print("-----------------------------------------------------------------------------------")
        print("Processing", dir)
        print()

        # Add all subdirs of the current one to the queue.
        for f in ftp_get_dirs( ftp, dir ):
            queue.append( f )

        # Get all files in the dir, and download them.
        files = ftp_get_files( ftp, dir )
        # print("Files:", files)
        for f in files:
            result = ftp_download_file( ftp, f, os.path.join( target_dir, f ))
            write_ftp_download_log(host, user, result[0], result[1], f)

            # Summary of all downloads. Cumbersome, because Python...
            if result[0] in summary:
                summary[result[0]] += 1
            else:
                summary[result[0]] = 1

        print()

    # We are polite, and close the connection respectfully. Bye, host.
    ftp.quit()

# =================================================================================================
#     Table of Sequencing Runs
# =================================================================================================

with open( run_table ) as csvfile:
    runreader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    for row in runreader:
        print("===================================================================================")
        print("Connecting to " + row["host"] + " as " + row["username"])
        print()

        ftp_download_all( row["host"], row["username"], row["password"], row[run_table_target_col] )

print("Summary:")
for key, val in summary.items():
    print(key + ": " + str(val))
