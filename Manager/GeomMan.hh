/*
  GeomMan.hh

  2024/04  K.Shirotori
*/

#ifndef GeomMan_h
#define GeomMan_h 1

#include "ThreeVector.hh"
#include <string>
#include <vector>
#include <map>

class GeomRecord;

class GeomMan 
{
private:
  GeomMan();
public:
  ~GeomMan();

public:
  void SetFileName( const char *filename ) { filename_=filename; }
  void SetFileName( const std::string &filename ) { filename_=filename; }

  bool Initialize( void );
  bool Initialize( const char *filename )
      { filename_=filename; Initialize(); return 0; }
   bool Initialize( const std::string &filename )
      { filename_=filename; Initialize(); return 0; }

  static GeomMan & GetInstance( void );
  double GetLocalZ( int lnum ) const;
  double GetResolution( int lnum ) const;

  double GetTiltAngle( int lnum ) const;
  double GetRotAngle1( int lnum ) const;
  double GetRotAngle2( int lnum ) const;
  const ThreeVector & GetGlobalPosition( int lnum ) const;
  ThreeVector NormalVector( int lnum ) const;
  ThreeVector UnitVector( int lnum ) const;
  const GeomRecord *GetRecord( int lnum ) const;

  ThreeVector Local2GlobalPos( int lnum, const ThreeVector &in ) const;
  ThreeVector Global2LocalPos( int lnum, const ThreeVector &in ) const;
  ThreeVector Local2GlobalDir( int lnum, const ThreeVector &in ) const;
  ThreeVector Global2LocalDir( int lnum, const ThreeVector &in ) const;

  double calcWirePosition( int lnum, double wire ) const;
  int calcWireNumber( int lnum, double position ) const;

  std::vector<int> GetDetectorIDList( void ) const;
  int GetLTofId( void ) const { return LTOFid_; }
  int GetDetectorId( const std::string &detName ) const;
  
private:
  static GeomMan *geomMan_;
  std::string filename_;
  mutable std::map <int, GeomRecord *> geomRecord_;
  int LTOFid_;

private:
  void clearElements( void );
};


#endif
