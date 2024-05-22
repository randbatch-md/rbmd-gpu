message ("src")

set(VAR1 "CACH1 VAR1" CACHE STRING "CACHE DOC")
set(VAR1 "CACH1 VAR2" CACHE STRING "CACHE DOC" FORCE)
message("var1 = ${VAR1}")

option(opt1 "opt1 doc" OFF)
message("opt1 = ${opt1}")
option(opt1 "opt1 doc" ON)
message("opt1 = ${opt1}")
