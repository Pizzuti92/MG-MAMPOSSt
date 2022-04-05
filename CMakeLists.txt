cmake_minimum_required(VERSION "3.17.1")
project("gomamposstopt")

enable_language (Fortran)
if (gfortran)
	set (CMAKE_Fortran_FLAGS "-fallow-invalid-boz")
#	set(dialect "-fallow-invalid-boz")
#    add_compile_options("-fallow-invalid-boz")
endif()


add_subdirectory(GamI)
add_subdirectory(Powell)
add_subdirectory(Utili)
add_subdirectory(JJin)
add_subdirectory(Newuoa)

file(GLOB file_SRC CONFIGURE_DEPENDS "*.f")
add_library(API ${file_SRC})
target_link_libraries(API
	GamI
	Powell
	Utili
	JJin
	Newuoa
	m
)



# add the executable
add_executable(${PROJECT_NAME} gomamposstoptS.f)
target_link_libraries(${PROJECT_NAME}
	GamI
	Powell
	Utili
	JJin
	Newuoa
	API
	m
)
	

find_package (Python3 COMPONENTS Interpreter REQUIRED)
install(TARGETS "${PROJECT_NAME}")

