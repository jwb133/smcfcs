## Test environments
* local Mac install, R 4.1.2
* local Windows 11, R 4.1.2
* R hub builder Windows Server 2022, R-devel, 64 bit
* R hub Ubuntu Linux 20.04.1 LTS, R-release, GCC
* R hub Fedora Linux, R-devel, clang, gfortran

## R CMD check results
No notes or warnings on Mac. On R hub I get a 'checking Rd cross-references note about the mice package being unavailable to check Rd xrefs, but I'm not sure why (the package is on CRAN). On the R hub Windows build I get a note about lastMiKTeXException being found when checking for detritus. On the R hub Windows I also get a qpdf warning, but from my previous submissions I don't think this is a problem for CRAN.
