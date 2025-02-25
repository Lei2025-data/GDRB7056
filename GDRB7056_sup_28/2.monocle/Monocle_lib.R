source("site_lib/lib_load.R", chdir = T)

lib_conf <- yaml::yaml.load_file("conf/Monocle_lib.yaml")

attach_lib(lib_conf)

#backup_lib(lib_conf, outfile = "test.R")

