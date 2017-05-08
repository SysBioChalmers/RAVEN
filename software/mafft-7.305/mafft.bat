@echo off

setlocal

set ROOTDIR=%~d0%~p0
set PATH=/usr/bin/:%PATH%
set MAFFT_BINARIES=/usr/lib/mafft
set TMPDIR=%ROOTDIR%/tmp

"%ROOTDIR%\usr\bin\bash" "%ROOTDIR%\usr\bin\mafft" %*

:EOF
