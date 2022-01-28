file(GLOB file_SRC CONFIGURE_DEPENDS "*.f" "*.f90")
add_library(Utili ${file_SRC})
