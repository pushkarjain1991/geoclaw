#import remove_file
import os
import shutil

def removefile(files):
    """Removes one or more files or directories"""
    if isinstance(files,str): #is files a string
        files = [files]
    if not isinstance(files,list):
        "Error"
    for file in files:
        if os.path.isdir(file):
            shutil.rmtree(file)
        elif os.path.isfile(file):
            os.remove(file)



def take(dirname):
    """Delete the directory *dirname* if exists, and create a new directory"""
    if os.access(dirname,os.F_OK):
        removefile(dirname)
    os.mkdir(dirname, 0755)


