add_executable(LogicFilter LogicFilter.cu)

set_target_properties(LogicFilter PROPERTIES
  CUDA_STANDARD 17
  CUDA_STANDARD_REQUIRED ON
)

target_include_directories(LogicFilter PUBLIC 
  ${Boost_INCLUDE_DIRS}
  ${FairMQ_INCDIR}
  ${HIREDIS_HEADER}
  ${REDIS_PLUS_PLUS_HEADER}
  ${UHBOOK_HEADER}
  ${CMAKE_SOURCE_DIR} # nestdaq-user-impl root directory
  ${CMAKE_CURRENT_SOURCE_DIR} # this subdirectory
  ${fmt_INCDIR}
)

target_link_directories(LogicFilter PUBLIC
  ${Boost_LIBRARY_DIRS}
  ${FairMQ_LIBDIR}
)

# branching for FairMQ version compatibility
if ((FairMQ_VERSION VERSION_GREATER_EQUAL 1.8.0) AND (FairMQ_VERSION VERSION_LESS_EQUAL 1.8.9))
  target_link_libraries(LogicFilter PUBLIC 
    ${Boost_LIBRARIES}
    ${HIREDIS_LIB}
    ${REDIS_PLUS_PLUS_LIB}
    ${UHBOOK_LIB}
    fairmq
    FairLogger
    Utility
    stdc++fs
    AmQTdcModule
    sqlite3
    ChannelMap::ChannelMap
    ROOT::RIO
    ROOT::Tree
  )
elseif (((FairMQ_VERSION VERSION_GREATER_EQUAL 1.4.55) AND (FairMQ_VERSION VERSION_LESS_EQUAL 1.4.56)) OR (FairMQ_VERSION VERSION_EQUAL 1.5.0))
  target_link_libraries(LogicFilter PUBLIC 
    ${Boost_LIBRARIES}
    ${HIREDIS_LIB}
    ${REDIS_PLUS_PLUS_LIB}
    ${UHBOOK_LIB}
    FairMQ
    FairMQStateMachine
    FairLogger
    Utility
    stdc++fs
    AmQTdcModule
    sqlite3
    ChannelMap::ChannelMap
    ROOT::RIO
    ROOT::Tree
  )
elseif ((FairMQ_VERSION VERSION_EQUAL 1.10.0))
  target_link_libraries(LogicFilter PUBLIC 
    ${Boost_LIBRARIES}
    ${HIREDIS_LIB}
    ${REDIS_PLUS_PLUS_LIB}
    ${UHBOOK_LIB}
    fairmq
    FairLogger
    Utility
    stdc++fs
    AmQTdcModule
    sqlite3
    ChannelMap::ChannelMap
    ROOT::RIO
    ROOT::Tree
  )
else()
  message(FATAL_ERROR "Unsupported FairMQ version ${FairMQ_VERSION}")
endif()

install(TARGETS LogicFilter
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)