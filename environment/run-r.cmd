@echo off
rem wrapper: clear LC_* vars so R 4.5.1 can read non-ASCII paths, then run script
set "LC_ALL="
set "LANG="
set "LC_CTYPE="
set "LC_COLLATE="
set "LC_TIME="
set "LC_MONETARY="
"C:\Program Files\R\R-4.5.1\bin\x64\Rscript.exe" --vanilla %*
exit /b %ERRORLEVEL%
