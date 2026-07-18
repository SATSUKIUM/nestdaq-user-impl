/*
  DCDriftParamMan.hh

  2018/12  K.Shirotori
*/

#ifndef DCDriftParamMan_h 
#define DCDriftParamMan_h 1

#include <map>
#include <string>

struct DCDriftParamRecord;

class DCDriftParamMan
{
public:
  DCDriftParamMan( const std::string & filename );
  ~DCDriftParamMan();

private:
  DCDriftParamMan( const DCDriftParamMan & );
  DCDriftParamMan & operator = ( const DCDriftParamMan & );

private:
  std::string ParameterFileName_;

  mutable std::map <unsigned int, DCDriftParamRecord *> Cont_;

public:
  bool Initialize( void );

  bool calcDrift( int PlaneId, double WireId, double ctime, 
		  double & dt, double & dl ) const;

private:
  DCDriftParamRecord * getParameter( int PlaneId, double WireId ) const;
  void clearElements( void );

  static double DriftLength1( double dt, double vel );
  static double DriftLength2( double dt, double p1, double p2, double p3,
			      double st, double p5, double p6 );
  static double DriftLength3( double dt, double p1, double p2, int PlaneId);
  static double DriftLength4( double dt, double p1, double p2, double p3);
  static double DriftLength5( double dt, double p1, double p2, double p3, double p4, double p5);
  static double DriftLength6(int PlaneId, double dt, double p1, double p2 ,double p3, double p4, double p5);

};

#endif
