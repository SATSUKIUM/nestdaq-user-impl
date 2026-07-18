/*
  DCTdcCalibMan.hh

  2018/12  K.Shirotori
*/

#ifndef DCTdcCalibMan_h
#define DCTdcCalibMan_h 1

#include <map>
#include <string>
#include <vector>

typedef std::vector <int> IntVec;

struct DCTdcCalMap;

class DCTdcCalibMan
{
public:
  DCTdcCalibMan( const std::string & filename );
  ~DCTdcCalibMan();

private:
  DCTdcCalibMan( const DCTdcCalibMan & );
  DCTdcCalibMan & operator = ( const DCTdcCalibMan & );

private:
  std::string MapFileName_;

  mutable std::map <unsigned int, DCTdcCalMap *> Cont_;

public:
  bool Initialize( void );
  bool GetTime( int PlaneId, double WireId, int tdc, double & time ) const;
  bool GetTdc( int PlaneId, double WireId, double time, int & tdc ) const;

private:
  DCTdcCalMap * getMap( int PlaneId, double WireId ) const;
  void clearElements( void );
};

#endif 
