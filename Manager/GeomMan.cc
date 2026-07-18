/*
  GeomMan.cc

  2024/04  K.Shirotori
*/

#include "GeomMan.hh"
#include "GeomRecord.hh"

#include <string>
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <cstdlib>

const int MaxChar = 200;

GeomMan *GeomMan::geomMan_=0;

GeomMan::GeomMan()
  : LTOFid_(21) 
{}

GeomMan::~GeomMan()
{}


GeomMan & GeomMan::GetInstance( void )
{
  if( !geomMan_ ){
    geomMan_ = new GeomMan();
  }
  return *geomMan_;
}

double GeomMan::GetLocalZ( int lnum ) const
{
   static const std::string funcname = "[GeomMan::GetLocalZ(int)]"; 
   GeomRecord *pGeo = geomRecord_[lnum];
   if( pGeo ) return pGeo->length_;
   else{
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
}

double GeomMan::GetResolution( int lnum ) const
{
   static const std::string funcname = "[GeomMan::GetResolution(int)]"; 
   GeomRecord *pGeo = geomRecord_[lnum];
   if( pGeo ) return pGeo->resol_;
   else{
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
}

double GeomMan::GetTiltAngle( int lnum ) const
{
   static const std::string funcname = "[GeomMan::GetTiltAngle(int)]";
   GeomRecord *pGeo = geomRecord_[lnum];
   if( pGeo ) return pGeo->tiltAngle_;
   else{
      std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
  }
}

double GeomMan::GetRotAngle1( int lnum ) const
{
   static const std::string funcname = "[GeomMan::GetRotAngle1(int)]";
   GeomRecord *pGeo = geomRecord_[lnum];
   if( pGeo ) return pGeo->rotAngle1_;
   else{
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
}

double GeomMan::GetRotAngle2( int lnum ) const
{
   static const std::string funcname = "[GeomMan::GetRotAngle2(int)]";
   GeomRecord *pGeo = geomRecord_[lnum];
   if( pGeo ) return pGeo->rotAngle2_;
   else{
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
}

const ThreeVector & GeomMan::GetGlobalPosition( int lnum ) const
{
   static const std::string funcname = "[GeomMan::GetGlobalPosition(int)]";
   GeomRecord *pGeo = geomRecord_[lnum];
   if( pGeo ) return pGeo->pos_;
   else{
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
}

ThreeVector GeomMan::NormalVector( int lnum ) const
{
   static const std::string funcname = "[GeomMan::NormalVector(int)]";
   GeomRecord *pGeo = geomRecord_[lnum];
   if( pGeo ) return pGeo->NormalVector();
   else{
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
}

ThreeVector GeomMan::UnitVector( int lnum ) const
{
   static const std::string funcname = "[GeomMan::UnitVector(int)]";
   GeomRecord *pGeo = geomRecord_[lnum];
   if( pGeo ) return pGeo->UnitVector();
   else{
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
}

const GeomRecord * GeomMan::GetRecord( int lnum ) const
{
   static const std::string funcname = "[GeomMan::GetRecord(int)]";
   GeomRecord *pGeo = geomRecord_[lnum];
   if( pGeo ) return pGeo;
   else{
      std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
}


double GeomMan::calcWirePosition( int lnum, double wire ) const
{
   static const std::string funcname = "[GeomMan::calcWirePosition()]";
   GeomRecord *pGeo = geomRecord_[lnum];
   if( pGeo ){
      return pGeo->WirePos(wire);
   }
   else{
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
}

int GeomMan::calcWireNumber( int lnum, double pos ) const
{
   static const std::string funcname = "[GeomMan::calcWireNumber()]";
   GeomRecord *pGeo = geomRecord_[lnum];
   if( pGeo ){
      return pGeo->WireNumber(pos);
   }
   else{
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
}


void GeomMan::clearElements( void )
{
   //  for_each( geomRecord_.begin(), geomRecord_.end(), DeleteObject() );
   std::map <int, GeomRecord *>::iterator itr;
   for( itr=geomRecord_.begin(); itr!=geomRecord_.end(); ++itr ){
      delete itr->second;
   }
   geomRecord_.clear();
   LTOFid_=51;
}


bool GeomMan::Initialize( void )
{
   static const std::string funcname = "[GeomMan::Initialize]";
   char str[MaxChar];
   char cname[MaxChar];
   int id;
   double xs, ys, zs, ta, ra1, ra2, l, res, w0, dd, ofs;
   double X_offset = 0.0;
   double Y_offset = 0.0;
   double Z_offset = 0.0;
   
   FILE *fp;
   
   if( ( fp = fopen( filename_.c_str(), "r" ) ) == 0 ){
    throw std::invalid_argument(funcname+": file open fail");
   }
   
   clearElements();
   
   while( fgets( str, MaxChar, fp ) != 0 ){
      if( str[0]!='#' ){
         if( sscanf( str, "%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                     &id, cname, &xs, &ys, &zs, &ta, &ra1, &ra2, &l, &res,
                     &w0, &dd, &ofs )
             == 13 ){
            GeomRecord *pRec = 
               new GeomRecord( id, cname, xs, ys, zs, 
                               ta, ra1, ra2, l, res,
                               w0, dd, ofs );
            GeomRecord *pOld = geomRecord_[id];
            geomRecord_[id] = pRec;
            
            if(strcmp(cname,"LTOF")==0) {
               LTOFid_=id;
            }
            
            if( pOld ){
               std::cerr << funcname << ": duplicated id number. "
                         << " following record is deleted." << std::endl;
               std::cerr << "Id=" << pOld->id_ << " " << pOld->pos_
                         << " ) ... " << std::endl;
               delete pOld;
            }
         }
         else {
            std::string strtemp=str;
            std::cerr << funcname << ": Invalid format " << strtemp << std::endl;
         }
      }
   }
   
   fclose(fp);
   
   return true;
}


std::vector <int> GeomMan::GetDetectorIDList( void ) const
{
   std::vector<int> vlist;
   vlist.reserve(geomRecord_.size());
   std::map <int, GeomRecord *>::const_iterator 
      itr=geomRecord_.begin(), end=geomRecord_.end();
   
   for(; itr!=end; ++itr ){
      vlist.push_back( itr->first );
   }
   
   return vlist;
}

ThreeVector GeomMan::Local2GlobalPos( int lnum, 
                                      const ThreeVector &in ) const
{
   static const std::string funcname = 
      "[GeomMan::Local2GlobalPos(ThreeVecor &)]";
   
   GeomRecord *pGeo = geomRecord_[lnum];
   if( !pGeo ) 
      throw std::out_of_range(funcname+": No record" );
   
   double x = pGeo->dxds_*in.x() + pGeo->dxdt_*in.y()
      + pGeo->dxdu_*in.z() + pGeo->pos_.x();
   double y = pGeo->dyds_*in.x() + pGeo->dydt_*in.y()
      + pGeo->dydu_*in.z() + pGeo->pos_.y();
   double z = pGeo->dzds_*in.x() + pGeo->dzdt_*in.y()
      + pGeo->dzdu_*in.z() + pGeo->pos_.z();
   
   return ThreeVector( x, y, z );
}

ThreeVector GeomMan::Global2LocalPos( int lnum,
                                      const ThreeVector &in ) const
{
   static const std::string funcname = 
      "[GeomMan::Global2LocalPos(ThreeVecor &)]";
   
   GeomRecord *pGeo = geomRecord_[lnum];
   if( !pGeo ){
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
   
   double x 
      = pGeo->dsdx_*(in.x()-pGeo->pos_.x())
      + pGeo->dsdy_*(in.y()-pGeo->pos_.y())
      + pGeo->dsdz_*(in.z()-pGeo->pos_.z());
   double y 
      = pGeo->dtdx_*(in.x()-pGeo->pos_.x())
      + pGeo->dtdy_*(in.y()-pGeo->pos_.y())
      + pGeo->dtdz_*(in.z()-pGeo->pos_.z());
   double z 
      = pGeo->dudx_*(in.x()-pGeo->pos_.x())
      + pGeo->dudy_*(in.y()-pGeo->pos_.y())
      + pGeo->dudz_*(in.z()-pGeo->pos_.z());
   
   return ThreeVector( x, y, z );
}

ThreeVector GeomMan::Local2GlobalDir( int lnum,
                                      const ThreeVector &in ) const
{
   static const std::string funcname = 
      "[GeomMan::Local2GlobalDir(ThreeVecor &)]";
   
   GeomRecord *pGeo = geomRecord_[lnum];
   if( !pGeo ){
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
   
   double x = pGeo->dxds_*in.x() + pGeo->dxdt_*in.y()
      + pGeo->dxdu_*in.z();
   double y = pGeo->dyds_*in.x() + pGeo->dydt_*in.y()
      + pGeo->dydu_*in.z();
   double z = pGeo->dzds_*in.x() + pGeo->dzdt_*in.y()
      + pGeo->dzdu_*in.z();
   
   return ThreeVector( x, y, z );
}

ThreeVector GeomMan::Global2LocalDir( int lnum,
                                      const ThreeVector &in ) const
{
   static const std::string funcname = 
      "[GeomMan::Global2LocalDir(ThreeVecor &)]";
   
   GeomRecord *pGeo = geomRecord_[lnum];
   if( !pGeo ){
      std::cerr << funcname << ": No record. Layer#=" 
                << lnum << std::endl;
      throw std::out_of_range(funcname+": No record" );
   }
   
   double x = pGeo->dsdx_*in.x() + pGeo->dsdy_*in.y()+ pGeo->dsdz_*in.z();
   double y = pGeo->dtdx_*in.x() + pGeo->dtdy_*in.y()+ pGeo->dtdz_*in.z();
   double z = pGeo->dudx_*in.x() + pGeo->dudy_*in.y()+ pGeo->dudz_*in.z();
   
   return ThreeVector( x, y, z );
}

int GeomMan::GetDetectorId( const std::string &detName ) const
{
   const std::string funcName = "GeomMan::GetDetectorId";
   
   std::map <int, GeomRecord *>::const_iterator 
      itr=geomRecord_.begin(), end=geomRecord_.end();
   
   for(; itr!=end; ++itr ){
      if (itr->second->name_ == detName)
         return itr->second->id_;
   }
   
   std::cerr << funcName << " : No such detector " << detName << std::endl;
   exit(-1);
}
