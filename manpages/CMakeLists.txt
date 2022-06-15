set(A2X_ROOT ${PROJECT_SOURCE_DIR}/manpages)

find_program(A2X_EXECUTABLE NAMES a2x a2x.py)
set(A2X_OPTS
  -D ${PROJECT_BINARY_DIR}/manpages
  -d manpage
  -f manpage
)

file(GLOB a2xsource ${A2X_ROOT}/*.1.txt)

foreach (gs ${a2xsource})
    get_filename_component(mpname ${gs} NAME_WE)
    ADD_CUSTOM_COMMAND(
        COMMAND  ${A2X_EXECUTABLE} ${A2X_OPTS} ${gs}
        OUTPUT   ${PROJECT_BINARY_DIR}/manpages/${mpname}.1
        DEPENDS  ${gs}
    )
    list(APPEND MAN_FILES ${PROJECT_BINARY_DIR}/manpages/${mpname}.1)
    ADD_CUSTOM_COMMAND(
        COMMAND  ${A2X_EXECUTABLE} -d manpage -f pdf -D ${PROJECT_BINARY_DIR}/manpages ${gs}
        OUTPUT   ${PROJECT_BINARY_DIR}/manpages/${mpname}.pdf
        DEPENDS  ${gs}
    )
    list(APPEND PDF_FILES ${PROJECT_BINARY_DIR}/manpages/${mpname}.pdf)
endforeach(gs)

add_custom_target(groffdocs ALL DEPENDS ${MAN_FILES})
add_custom_target(pdfdocs ALL DEPENDS ${PDF_FILES})

install(FILES ${MAN_FILES} DESTINATION ${MANDIR})
