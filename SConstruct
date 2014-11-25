import os
import SCons

#get the mode flag from the command line
#default to 'release' if the user didn't specify
mymode = ARGUMENTS.get('mode', 'release')   #holds current mode
Export('mymode')

# only 'debug' or 'release' allowed
if not (mymode in ['debug', 'release']):
   print "Error: expected 'debug' or 'release', found: " + mymode
   Exit(1)
   
#tell the user what we're doing
print '**** Compiling in ' + mymode + ' mode...'

debugcflags = ['-g']   #extra compile flags for debug
releasecflags = ['-O2', '-march=native', '-flto', '-DNDEBUG']         #extra compile flags for release

EIGEN_PATH='/usr/include/eigen3'
flags = ["-std=c++11"]
if mymode == 'debug':
	BUILD_DIR = 'dbg'
	flags += debugcflags
else:
	BUILD_DIR = 'rel'
	flags += releasecflags
Export('BUILD_DIR')
	
env = Environment(CPPPATH = ['#', EIGEN_PATH], CXXFLAGS=flags)
Export('env')

GTEST_LIBS = ['gtest', 'pthread']
Export('GTEST_LIBS')
OTHER_LIBS = ['fftw3', 'boost_system', 'lapack', 'gsl', 'gslcblas']
Export('OTHER_LIBS')
PREFIX = '/usr/local'

def call(subdir):
        return SConscript(os.path.join(subdir, 'SConscript'), variant_dir = os.path.join(subdir, BUILD_DIR), duplicate = 0)

## Build libraries
math = call('math')
Export('math')
qsd = call('qsd')
Export('qsd')
protein_chain = call('protein_chain')
Export('protein_chain')
cfa = call('cfa')
Export('cfa')
pseudomode = call('pseudomode')
Export('pseudomode')

## Build tests
call('math-test')
call('qsd-test')
call('protein_chain-test')

## Build calculation programs
cfa_transport = call('cfa_transport')
strunz_transport = call('strunz-transport')
pseudomode_absorption = call('pseudomode-absorption')

# define the custom function for "chmod"
from SCons.Script.SConscript import SConsEnvironment
SConsEnvironment.Chmod = SCons.Action.ActionFactory(os.chmod, lambda dest, mode: 'Chmod("%s", 0%o)' % (dest, mode))
        
def InstallPerm(env, dest, files, perm):
    obj = env.Install(dest, files)
    for i in obj:
        env.AddPostAction(i, env.Chmod(str(i), perm))
    return dest

# put this function "in" scons
SConsEnvironment.InstallPerm = InstallPerm


env.Alias('install-lib', env.InstallPerm(os.path.join(PREFIX, 'lib'), [math, qsd, protein_chain, cfa], 0644))
env.Alias('install-bin', env.InstallPerm(os.path.join(PREFIX, 'bin'), [pseudomode_absorption, cfa_transport, strunz_transport], 0755))
