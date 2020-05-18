@echo off

set /p file="SRA IDs list file: "
set /p outfile="Ourput directory: "

echo
echo "Start procesing ... "
echo "You can take a nap, if you want to download many files zzZ"
echo

E:\SRA_download\sratoolkit.2.9.6\bin\prefetch.exe --option-file %file% --output-directory ./public/%outfile% --verbose

echo
echo "Congratulations! All is finished! Please check the md5 values."
echo
pause

goto :eof
