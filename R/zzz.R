# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Package initialization
#

#
# Package attaching
#
.onAttach <- function(libname, pkgname){
    packageStartupMessage("Thank you for using 'spfrontier'")
}