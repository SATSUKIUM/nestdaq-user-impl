#==============================================================================

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

message(STATUS "CMAKE_CUDA_ARCHITECTURES: ${CMAKE_CUDA_ARCHITECTURES}")

set(INSTALL_BINDIR ${CMAKE_SOURCE_DIR}/bin)
set(INSTALL_LIBDIR ${CMAKE_SOURCE_DIR}/lib)

# Find cudamap library
find_library(CUDAMAP_LIBRARY
    NAMES cudamap
    PATHS "$ENV{NESTDAQ}/lib"
    REQUIRED
)


add_executable(LogicFilter LogicFilter.cu;)
target_compile_definitions(LogicFilter PRIVATE
  LOGICFILTER_ELAPSED_DIR="${CMAKE_SOURCE_DIR}/test"
)
target_sources(LogicFilter PRIVATE
    SignalParser.cxx
)
set_target_properties(LogicFilter PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_include_directories(LogicFilter PUBLIC 
  ${Boost_INCLUDE_DIRS};  
  ${FairMQ_INCDIR};
  ${FairMQ_INCDIR}/fairmq;
  ${FairLogger_INCDIR};
  ${CMAKE_SOURCE_DIR};
  ${HIREDIS_HEADER};
  ${REDIS_PLUS_PLUS_HEADER};
  ${UHBOOK_HEADER};
  "$ENV{NESTDAQ}/include"
  ./;
)

target_link_directories(LogicFilter PUBLIC
  ${Boost_LIBRARY_DIRS}; 
  ${FairMQ_LIBDIR};
)

if ((FairMQ_VERSION VERSION_GREATER_EQUAL 1.8.0)
  AND (FairMQ_VERSION VERSION_LESS_EQUAL 1.8.9))
  target_link_libraries(LogicFilter PUBLIC 
    ${Boost_LIBRARIES};  
    fairmq;
    FairLogger;
    ${HIREDIS_LIB};
    ${REDIS_PLUS_PLUS_LIB};
    pthread;
    "${CUDAMAP_LIBRARY}"
  )
elseif ((FairMQ_VERSION VERSION_GREATER_EQUAL 1.4.55)
  AND (FairMQ_VERSION VERSION_LESS_EQUAL 1.4.56))
  target_link_libraries(LogicFilter PUBLIC 
    ${Boost_LIBRARIES};  
    FairMQ;
    FairMQStateMachine;
    FairLogger;
    ${HIREDIS_LIB};
    ${REDIS_PLUS_PLUS_LIB};
    pthread;
    "${CUDAMAP_LIBRARY}"
  )
elseif((FairMQ_VERSION VERSION_EQUAL 1.10.0))
  target_link_libraries(LogicFilter  PUBLIC
    fairmq;
    FairLogger;
    ${Boost_LIBRARIES};
    Utility;
    sqlite3;
    ${HIREDIS_LIB};
    ${REDIS_PLUS_PLUS_LIB};
    "${CUDAMAP_LIBRARY}"
  )
else()
  message(FATAL_ERROR "Unsupported FairMQ version ${FairMQ_VERSION}")
endif()

install(TARGETS LogicFilter
  #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib64
  ARCHIVE DESTINATION lib64
)
