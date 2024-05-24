


set_property(DIRECTORY . PROPERTY DIR1 "DIR001")
get_property(dir DIRECTORY . PROPERTY DIR1)
message("dir = ${dir}")

message("use double = $CACHE{ENABLE_DOUBLE}")
message("use 64 BIT = $CACHE{ENABLE_64BIT_IDS}")