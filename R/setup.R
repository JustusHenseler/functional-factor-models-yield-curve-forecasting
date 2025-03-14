
packages <- c("zoo",
              "remotes", 
              "vars", 
              "tsbox", 
              "gsheet", 
              "rugarch", 
              "rmgarch",
              "xts",
              "gridExtra",
              "ggplot2",
              "plotly",
              "magrittr",
              "reshape2")

packages_to_install <- packages[!(packages %in% installed.packages()[,"Package"])]

if(length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

if(!("dffm" %in% installed.packages()[,"Package"])) {
  install_github("ottosven/dffm")
}

