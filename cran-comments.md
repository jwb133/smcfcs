## Test environments
* local Windows 11 install, R 4.2.2
* R Hub Windows Server 2022, R-devel, 64 bit
* R Hub Ubuntu Linux 20.04.1 LTS, R-release, GCC
* R Hub Fedora Linux, R-devel, clang, gfortran

## R CMD check results
On Windows local machine, 0 errors, 0 warnings, 1 note about being unable to verify current time, which I think is because http://worldclockapi.com/ is down.

On R Hub I get 4 notes: 1) checking CRAN incoming feasibility ... [12s] NOTE because my maintainer email address has changed. 2) Note on checking Rd cross-references - Package unavailable to check Rd xrefs: 'mice'. 3) Note - checking for detritus in the temp directory - Found the following files/directories: 'lastMiKTeXException'. 4) Checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found (on Linux only).
