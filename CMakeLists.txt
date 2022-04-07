cmake_minimum_required(VERSION "3.17.1")
project("gomamposstopt")
include(CheckFortranCompilerFlag)
enable_language (Fortran)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    CHECK_Fortran_COMPILER_FLAG("-fallow-invalid-boz" prova)

    if(${prova})
		set (CMAKE_Fortran_FLAGS "-fallow-invalid-boz")
#		set(dialect "-fallow-invalid-boz")
#		add_compile_options("-fallow-invalid-boz")
    endif()
endif()


add_subdirectory(GamI)
add_subdirectory(Powell)
add_subdirectory(Utili)
add_subdirectory(JJin)
add_subdirectory(Newuoa)

file(GLOB file_SRC CONFIGURE_DEPENDS "*.f")
add_library(MAM ${file_SRC})
target_link_libraries(MAM
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
	MAM
	m
)
	

find_package (Python3 COMPONENTS Interpreter REQUIRED)
install(TARGETS "${PROJECT_NAME}")

