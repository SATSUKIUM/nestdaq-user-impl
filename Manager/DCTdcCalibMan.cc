/*
  DCTdcCalibMan.cc

  2018/12  K.Shirotori
*/

#include "DCTdcCalibMan.hh"

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstdlib>

const int BufSize = 200;

struct DCTdcCalMap {
  DCTdcCalMap( double q0, double q1 )
    : p0(q0), p1(q1)
  {}
  double p0, p1;
};

DCTdcCalibMan::DCTdcCalibMan( const std::string & filename )
  : MapFileName_(filename)
{}

DCTdcCalibMan::~DCTdcCalibMan()
{
  clearElements();
}

void DCTdcCalibMan::clearElements( void )
{
  std::map <unsigned int, DCTdcCalMap *>::iterator itr;
  for( itr=Cont_.begin(); itr!=Cont_.end(); ++itr ){
    delete itr->second;
  }
  Cont_.clear();
}

inline unsigned int SortKey( int plane, double wire )
{
  return (plane<<10) | int(wire);
}


bool DCTdcCalibMan::Initialize( void )
{
  static const std::string funcname = "[DCTdcCalibMan::Initialize]";

  int pid, wid;
  double p0, p1;
  FILE *fp;
  char buf[BufSize];

  if((fp=fopen(MapFileName_.c_str(),"r"))==0){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }

  while(fgets(buf,BufSize,fp)){
    if(buf[0]!='#'){
      if(sscanf(buf,"%d %d %lf %lf",&pid,&wid,&p1,&p0)==4){
	unsigned int key=SortKey(pid,wid);
	DCTdcCalMap *p=new DCTdcCalMap(p0,p1);
	if(p){
	  if(Cont_[key]) delete Cont_[key];
	  Cont_[key]=p;
	}
	else{
	  std::cerr << funcname << ": new fail. "
		    << " Plane=" << std::setw(3) << std::dec << pid 
		    << " Wire=" << std::setw(3) << std::dec << wid 
		    << std::endl;
	}
      }
      else{
	std::cerr << funcname << ": Bad format => "
		  << std::string(buf) << std::endl;
      }
    }
  }
  fclose(fp);

  std::cout << funcname << ": Initialization finished" << std::endl;
  return true;
}

DCTdcCalMap * DCTdcCalibMan::getMap( int PlaneId, double WireId ) const
{
  unsigned int key=SortKey(PlaneId,WireId);
  return Cont_[key];
}  

bool DCTdcCalibMan::GetTime( int PlaneId, double WireId,
			     int tdc, double & time ) const
{
  static const std::string funcname = "[DCTdcaCalibMan::GetTime]";
  DCTdcCalMap *p=getMap(PlaneId,WireId);
  if(p){ 
    time=(p->p0)+(p->p1)*tdc;
    return true;
  }
  else{
    //    std::cerr << funcname << ": No record. "
    //	      << " PlaneId=" << std::setw(3) << std::dec << PlaneId
    //	      << " WireId=" << std::setw(3) << std::dec << WireId 
    //	      << std::endl;
    return false;
  }
}

bool DCTdcCalibMan::GetTdc( int PlaneId, double WireId,
			    double time, int &tdc ) const
{
  static const std::string funcname = "[DCTdcaCalibMan::GetTdc]";
  DCTdcCalMap *p=getMap(PlaneId,WireId);
  if(p){
    tdc=int((time-(p->p0))/(p->p1));
    return true;
  }
  else{
    //    std::cerr << funcname << ": No record. "
    //	      << " PlaneId=" << std::setw(3) << std::dec << PlaneId
    //	      << " WireId=" << std::setw(3) << std::dec << WireId 
    //	      << std::endl;
    return false;
  }
}
