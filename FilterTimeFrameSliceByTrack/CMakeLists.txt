# Tracking filter and its detector/track reconstruction implementation.
add_executable(FilterTimeFrameSliceByTrack
  FilterTimeFrameSliceByTrack.cxx
  DCHit.cxx
  DCPairHitCluster.cxx
  DCLTrackHit.cxx
  DCMathTools.cxx
  DCLocalTrack.cxx
)

target_include_directories(FilterTimeFrameSliceByTrack PRIVATE
  ${Boost_INCLUDE_DIRS}
  ${FairMQ_INCDIR}
  ${HIREDIS_HEADER}
  ${REDIS_PLUS_PLUS_HEADER}
  ${UHBOOK_HEADER}
  ${PROJECT_SOURCE_DIR}
  ${fmt_INCDIR}
)

target_link_directories(FilterTimeFrameSliceByTrack PRIVATE
  ${Boost_LIBRARY_DIRS}
  ${FairMQ_LIBDIR}
)

target_link_libraries(FilterTimeFrameSliceByTrack PRIVATE
  ChannelMap::ChannelMap
  ROOT::RIO
  ROOT::Tree
  ${Boost_LIBRARIES}
  ${HIREDIS_LIB}
  ${REDIS_PLUS_PLUS_LIB}
  ${UHBOOK_LIB}
  Utility
  stdc++fs
  AmQTdcModule
  sqlite3
)

if ((FairMQ_VERSION VERSION_GREATER_EQUAL 1.8.0 AND
     FairMQ_VERSION VERSION_LESS_EQUAL 1.8.9) OR
    FairMQ_VERSION VERSION_EQUAL 1.10.0)
  target_link_libraries(FilterTimeFrameSliceByTrack PRIVATE fairmq FairLogger)
elseif ((FairMQ_VERSION VERSION_GREATER_EQUAL 1.4.55 AND
         FairMQ_VERSION VERSION_LESS_EQUAL 1.4.56) OR
        FairMQ_VERSION VERSION_EQUAL 1.5.0)
  target_link_libraries(FilterTimeFrameSliceByTrack PRIVATE
    FairMQ FairMQStateMachine FairLogger)
else()
  message(FATAL_ERROR "Unsupported FairMQ version ${FairMQ_VERSION}")
endif()

install(TARGETS FilterTimeFrameSliceByTrack
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
