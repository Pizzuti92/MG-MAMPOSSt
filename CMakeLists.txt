cmake_minimum_required(VERSION "3.17.1")
project("gomamposstopt")

enable_language (Fortran)


add_subdirectory(GamI)
add_subdirectory(Powell)
add_subdirectory(Utili)
add_subdirectory(JJin)
add_subdirectory(Newuoa)

# add the executable
add_executable(${PROJECT_NAME} gomamposstoptS.f)
target_link_libraries(${PROJECT_NAME}
	GamI
	Powell
	Utili
	JJin
	Newuoa
	m
)
	

find_package (Python3 COMPONENTS Interpreter REQUIRED)
install(TARGETS "${PROJECT_NAME}")

