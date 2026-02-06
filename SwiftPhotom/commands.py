import os
import subprocess
import sys


def run(_command, verbose=False):
    pid = subprocess.Popen(
        _command, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    output, _ = pid.communicate()
    returncode = pid.returncode

    if verbose or returncode != 0:
        print('CMD:', _command)
        if output:
            print(output.decode('utf-8', errors='replace'))
        if returncode != 0:
            print('Return code:', returncode)

    if returncode != 0:
        print('ERROR while running ' + _command)
        sys.exit()

    return output


def uvotimsum(_in, _out, _exclude='none', traceback=False):
    _in = os.path.abspath(_in)
    _out = os.path.abspath(_out)
    comm = 'uvotimsum %s %s exclude=%s' % (_in, _out, _exclude)
    # if traceback:
    #     comm += ' %tb'
    run(comm, verbose=True)  # TODO: set verbose=False after debugging


def uvotmaghist(_in, _reg, _bgreg, _out, _gif):
    comm = 'uvotmaghist %s srcreg=%s bkgreg=%s  outfile=%s plotfile=%s ' \
           'coinfile=caldb zerofile=caldb exclude=none chatter=0 clobber=yes ' \
           'logtime=no psffile=caldb apercorr=curveofgrowth ' \
           'skipbad=yes' % (_in, _reg, _bgreg, _out, _gif)

    run(comm)


def fappend(_in, _out):
    # Use absolute paths so HEASoft can find files regardless of working directory
    _in = os.path.abspath(_in)
    _out = os.path.abspath(_out)
    comm = 'fappend "' + _in + '" "' + _out + '"'
    run(comm)


def fcopy(_in, _out):
    # Use absolute paths so HEASoft can find files regardless of working directory
    _in = os.path.abspath(_in)
    _out = os.path.abspath(_out)
    comm = 'fcopy "' + _in + '" "' + _out + '"'
    run(comm)


# def find_attfile(_infile):
#    obsID=_infile.split('/')[0]
#    aux=glob.glob(obsID+'/auxil/sw'+obsID+'sat.fits*')
#    if os.path.isfile(aux[0]):
#        return aux[0]
#    else:
#        return 0
#
# def uvotskycorr(_what,_infile,_corrfile,_attfile,_outfile):
#    comm = 'uvotskycorr what=%s skyfile=%s corrfile=%s attfile=%s outfile=%s ' \
#        'starid="mag.err=5 rot.error=60" catspec=usnob1.spec options=NONE ' \
#        'flagfile=NONE' % (_what, _infile, _corrfile, _attfile, _outfile)
#    run(comm)
#
# def aspect_correction(_infile):
#    attfile=find_attfile(_infile)
#    if not attfile:
#        print('attfile for '+_infile+'not found. Cannot perform aspect correction.')
#        return
#
#    corrfile='reduction/auxiliary/'+_infile.split('/')[0]
#    uvotskycorr('ID',_infile,'NONE',attfile,corrfile)
#
#    uvotskycorr('ID',_infile,corrfile,attfile,'NONE')
