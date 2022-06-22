## Test environments
* local Windows 11 install, R 4.2.0
* R Hub Windows Server 2022, R-devel, 64 bit
* R Hub Ubuntu Linux 20.04.1 LTS, R-release, GCC
* R Hub Fedora Linux, R-devel, clang, gfortran

## R CMD check results
0 errors, 0 warnings, 0 notes on local Windows.

On R Hub Windows R-devel, I get: note on checking Rd cross-references - Package unavailable to check Rd xrefs: 'mice', and a note - checking for detritus in the temp directory - Found the following files/directories: 'lastMiKTeXException'. Not sure why mice package would be unavailable.
