@echo off
set ignored_warnings= -wd4201 -wd4244 -wd4365 -wd4514 -wd4710 -wd5045 -wd4711
cl -nologo -FC -Wall %ignored_warnings% -O2 main.cpp /link /out:rtr.exe
set ignored_warnings=
del *.obj