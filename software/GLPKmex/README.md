# GLPKMEX

GLPKMEX is a Matlab MEX Interface for the GLPK library developed by
Andrew Makhorin. GLPKMEX is developed by Nicolo' Giorgetti, email
giorgetti  at  ieee.org.
GLPK is currently being mantained by Niels Klitgord, email
niels  at  bu.edu.

This version is maintained by Benoît Legat, email
benoit.legat  at  gmail.com.

To install mps2mat.py script, please see INSTALL_mps2mat file.

## Precompiled binaries

You don't need to compile GLPKMEX by yourself nor install GLPK.
They are both compiled into `glpkcc.mex*` files.

* `a64` is for Linux 64 bits
* `glx` is for Linux 32 bits
* `maci64` is for Mac OS 64 bits
* `w32` is for Windows 32 bits (in the win32 subdir)
* `w64` is for Windows 64 bits (in the win64 subdir)

On Linux 64 bits and Windows 64 bits, GLPK v4.62 is used, but Linux 32 bits and Mac OS versions use GLPK v4.48 (see [this issue](https://github.com/blegat/glpkmex/issues/3)).

Extension | OS      | Architecture | GLPK version |
--------- | ------- | ------------ | ------------ |
`a64`     | Linux   | 64 bits      | v4.62        |
`glx`     | Linux   | 32 bits      | v4.48        |
`maci64`  | Mac OS  | 64 bits      | v4.48        |
`w32`     | Windows | 32 bits      | v4.60        |
`w64`     | Windows | 64 bits      | v4.62        |

## Quick installation instructions

Open MATLAB (for GNU/Linux user, you will need roots privileges so open a terminal (`CTRL+ALT+T`) and enter `sudo matlab` or `sudo matlab -glnx86` if you have problems with architecture (see [this wiki](https://help.ubuntu.com/community/MATLAB))).
Then click on _File/Set Path..._ and add the path to where you
put __glpkmex__. Then hit save.
For GNU/Linux users, you can now quit MATLAB, don't run it with `sudo` again except for _Set Path..._.

To check everything, try the two commands below in MATLAB

    >> glpktest1
    >> glpktest2
which should output a file `SimpleLP.mps` without error.

## Instructions for compiling from source

### Standard installation procedure (suitable for most 32-bit linux/windows users)

1. Download and install GLPK version 4.60 or higher:
       http://ftp.gnu.org/gnu/glpk/

2. Start Matlab and run makeglpkmex.m. Specify correct path to both GLPK
   directory and eventually to the GLPK include and library directories.

3. Test the interface on the examples included (`glpktest1.m`, `glpktest2.m`, `glpksparse.m`). Everything should works fine.

### Linux 64bit machine
Niels's Install Notes

1. if you plan on compling glpk with gmp, recompile gmp with `CFLAGS+=-fPIC`. Then recompile `glpk-4.60` (or newer) with `CFLAGS+=-fPIC`.
    This is do deal with a weird issue with ld in 64bit format

2. Directly compile glpkmex with (default):

        mex -I<path to>/include glpkcc.cpp <path to>/libglpk.a

(with gmp compiled into glpk and large arrays (64bit stuff...)):

        mex -largeArrayDims -I<path to>/include glpkcc.cpp <path to>/libglpk.a <path to>/libgmp.a

**note -largArrayDims is optional, but allows matlab to use large arrays with glpk on 64bit machines**

this particular pipeline has worked for me on a number of different 64bit linux boxes:

1. in the linux comandline cd to glpk dir
2. run (**note: update the dir install dir if you want, but be sure to do the same below**):

        make clean
        ./configure
        make CFLAGS+=-fPIC
        make check
        make  prefix=/usr/local install

3. cd glpkmex dir
4. run from the comand line:

        mex -largeArrayDims -I/usr/local/include glpkcc.cpp /usr/local/lib/libglpk.a

5. start matlab, add glpkmex dir to path and run glpktest1 and glpktest2 in matlab environment.  if these work you should be good to go.

### Install GLPK on 64bit Linux with 32bit MATLAB
First, you will need to install GLPK.
Get the source [here](http://ftp.gnu.org/gnu/glpk/), unpack it with the ''Archive Manager'' and then open a terminal and go in the folder where you un packed it (for example, if it is in download, run `cd ~/Downloads/glpk-4.60`)
Then run

    make clean
    ./configure --build=i686-pc-linux-gnu "CFLAGS=-m32" "CXXFLAGS=-m32" "LDFLAGS=-m32"
    make CFLAGS+=-fPIC
    make check
    make prefix=/usr/local install
The `CFLAGS+=-fPIC` is needed if you plan to install it in 64 bits with gmp (if you're not sure, just use it).
The ` --build=i686-pc-linux-gnu "CFLAGS=-m32" "CXXFLAGS=-m32" "LDFLAGS=-m32"` is needed if you are running MATLAB 32 bits on a 64 bits machine.

Open MATLAB and navigate to the folder of GLPKMEX and run

    mex CXXFLAGS='$CXXFLAGS -m32' -I/usr/local/include glpkcc.cpp /usr/local/lib/libglpk.a

Only use `CXXFLAGS='$CXXFLAGS -m32'` if you have a 32 bits MATLAB.
I don't understand why it needs this option since `mex` is run inside MATLAB so it should know that it is in 32 bits however it is needed for me.

### Installing on a PC
1. if you used winglpk, you still need to compile the `glpk.lib` file.
   to do this go to w32 or w64 folder in winglpk and run `Build_GLPK_with_VC9.bat`.
   see the `readme.txt` file in w32 or w64 for more details on how to do this.
2. note windows uses `.lib`, not `.a` files so change `libglpk.a` to `glpk.lib` (simular for gmp if used)
3. full paths to directions (like include) often need to be in quotes.
   ex:

        mex -I'C:\\some path\include' glpkcc.cpp 'C:\\some other path\glpk.lib'

### Installing with cygwin (alternate method on a pc):
Requirements:
A) The last version of cygwin, downloadable from http://www.cygwin.com
   (or download mingw)
B) Gnumex: http://gnumex.sourceforge.net/

1. Download GLPK and install it by specifying CFLAGS="-O3 -mno-cygwin" as
   argument of make, i.e.:

        ./configure;
        make CFLAGS="-O3 -mno-cygwin"
        make install
   This avoid to build a glpk library which needs of cygwin1.dll to run.

2. Start Matlab and run gnumex. Make a mexopts.bat with option 'cygwin-mingw' if
   you are using cygwin or 'mingw' if you are using mingw.

3. run makeglpkmex.m. Specify correct path to GLPK directory and eventually to
   the GLPK include and library directories. Moreover, you have to specify the
   location of mexopts.bat file generated at step 2).

4. Test the interface on the examples included (glpktest1.m, glpktest2.m,
   glpksparse.m). Everything should works fine.


### Installing on a MAC
On mac, it works out of the box with the file `glpkmex.mexmaci64`.
It was compiled with glpk-4.60.

However, if you want to compile by yourself for whatever reason.
I'd like to mention that everything worked as expected for me except for mex for which
I had to change the compiler from gcc-4.2 to gcc and I had to remove `-syslibroot /Developer/SDKs/MacOSX10.6.sdk`
So here is the command I had to enter in MATLAB:

        mex -largeArrayDims -I/usr/local/include glpkcc.cpp /usr/local/lib/libglpk.a -v CC='gcc' CXX='g++' CFLAGS='-fno-common -no-cpp-precomp -arch x86_64 -mmacosx-version-min=10.5  -fexceptions' CXXFLAGS='-fno-common -no-cpp-precomp -fexceptions -arch x86_64 -mmacosx-version-min=10.5' LD='gcc' LDFLAGS='-Wl,-twolevel_namespace -undefined error -arch x86_64 -Wl -mmacosx-version-min=10.5 -bundle -Wl,-exported_symbols_list,/Applications/MATLAB_R2012a_Student.app/extern/lib/maci64/mexFunction.map'

Another special notes for installing on a Mac,( worked for Leapard )
add GLPK lib install directory to your DYLD_LIBRARY_PATH shell variable (before you start matlab)
example: 'export DYLD_LIBRARY_PATH=/usr/local/lib:$DYLD_LIBRARY_PATH' in .bashrc
continue as above for linux install.
