#include "DCLocalTrack.h"
#include "DCLTrackHit.h"
#include "DCMathConstants.h"
#include "DCMathTools.h"
#include "DCConstants.h"
#include <string>

using nestdaq::DCLocalTrack;
using nestdaq::DCLTrackHit;

DCLTrackHit* DCLocalTrack::GetHit(std::size_t iHit) const {
    if(iHit < dclthits.size()){
        return dclthits[iHit];
    }
    else{
        return nullptr;
    } 
} // DCLTrackHit* nestdaq::DCLocalTrack::GetHit(std::size_t iHit) const

bool DCLocalTrack::DoFit(){
// True/False condition
// 1. Number of hits >= DCLocalMinNHits(⇦comes from DCConstants.h)
// 2. Gauss elimination is successful
    const std::string_view funcname = "[nestdaq::DCLocalTrack::DoFit] ";

    std::size_t nHits = dclthits.size();
    if(nHits < DCConstants::DCLocalMinNHits){
        status = false;
        return status;
    }

    std::vector<double> z, res, w, ct, st;
    z.reserve(nHits); // position along z-axis
    res.reserve(nHits); // 1/sigma^2, where sigma is the resolution of the hit
    w.reserve(nHits); // position along measurement axis
    ct.reserve(nHits); // cos of wire angle
    st.reserve(nHits); // sin of wire angle

    // fill observables
    for(std::size_t iHit=0; iHit<nHits; ++iHit){
        DCLTrackHit* dclhit = dclthits[iHit];
        if(dclhit != nullptr){
            z.push_back(dclhit->GetGlobalZ());
            res.push_back(1.0 / (dclhit->GetResolution() * dclhit->GetResolution()));
            w.push_back(dclhit->GetWirePosition() + dclhit->GetLeftRight() * dclhit->GetDriftLength());
            ct.push_back(cos(DCMath::Deg2Rad * dclhit->GetWireAngle()));
            st.push_back(sin(DCMath::Deg2Rad * dclhit->GetWireAngle()));
        } // if(dclhit != nullptr)
    } // for(std::size_t iHit=0; iHit<nHits; ++iHit)

    size_t nn = z.size();

    // prepare matrix and vector for linear equation
    double matrx[16], *mtp[4], fitp[4];
    mtp[0]=&matrx[0]; mtp[1]=&matrx[4]; mtp[2]=&matrx[8]; mtp[3]=&matrx[12];
    for( int i=0; i<4; ++i ){
        fitp[i]=0.0;
        for( int j=0; j<4; ++j ){
            mtp[i][j]=0.0;
        }
    }

    for( std::size_t i=0; i<nn; ++i ){
        double res_=res[i], z_=z[i], w_hat=w[i], ctt=ct[i], stt=st[i];
        mtp[0][0] += res_*ctt*ctt;
        mtp[0][1] += res_*z_*ctt*ctt;
        mtp[0][2] += res_*ctt*stt;
        mtp[0][3] += res_*z_*ctt*stt;
        mtp[1][1] += res_*z_*z_*ctt*ctt;
        mtp[1][2] += res_*z_*ctt*stt;
        mtp[1][3] += res_*z_*z_*ctt*stt;
        mtp[2][2] += res_*stt*stt;
        mtp[2][3] += res_*z_*stt*stt;
        mtp[3][3] += res_*z_*z_*stt*stt;

        fitp[0] += res_*w_hat*ctt;
        fitp[1] += res_*z_*w_hat*ctt;
        fitp[2] += res_*w_hat*stt;
        fitp[3] += res_*z_*w_hat*stt;
    }
    mtp[1][0]=mtp[0][1]; mtp[2][0]=mtp[0][2]; mtp[3][0]=mtp[0][3];
    mtp[2][1]=mtp[1][2]; mtp[3][1]=mtp[1][3]; mtp[3][2]=mtp[2][3];

    std::vector<int> indxc(nn), indxd(nn), ipiv(nn);

#if 0
  std::cout<<"             Vector:  Q_i               "<<std::endl;
  std::cout<<"   A1=  "<<fitp[0]<<"     A2=  "<<fitp[1]<<"     A3=  "<<fitp[2]
	   <<"     A4=  "<<fitp[3]<<std::endl;

  std::cout<<"             original matrix: M_ij        "<<std::endl;
  std::cout<<"    A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2]
      	   <<"    A14="<<mtp[0][3]<<std::endl;
  std::cout<<"    A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2]
	   <<"    A24="<<mtp[1][3]<<std::endl;
  std::cout<<"    A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2]
	   <<"    A34="<<mtp[2][3]<<std::endl;
  std::cout<<"    A41="<<mtp[3][0]<<"    A42="<<mtp[3][1]<<"    A43="<<mtp[3][2]
	   <<"    A44="<<mtp[3][3]<<std::endl;
#endif

    double Org[4][4]={0},Red[4][4]={0},Final[4][4]={0};
    double Org_vec[4]={0}, Solution_vec[4]={0};
    for(int l=0; l<4;l++){
        for(int m=0; m<4; m++){
        Org[l][m]=mtp[l][m];
        }
    }

    for (int i=0; i<4; i++){
        Org_vec[i]=fitp[i];
    }

    if( DCMathTools::GaussJordan(mtp,4,fitp,&indxc[0],
                    &indxd[0],&ipiv[0])==false ){
        std::cerr << funcname << ": Fitting fails" << std::endl;
        return status_=false;
    }
    x0=fitp[0]; y0=fitp[2]; u0=fitp[1]; v0=fitp[3];

#if 0
  std::cout<<"             reduced matrix        "<<std::endl;
  std::cout<<"     A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2]
	   <<"     A14="<<mtp[0][3]<<std::endl;
  std::cout<<"     A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2]
	   <<"     A24="<<mtp[1][3]<<std::endl;
  std::cout<<"     A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2]
	   <<"     A34="<<mtp[2][3]<<std::endl;
  std::cout<<"     A41="<<mtp[3][0]<<"    A42="<<mtp[3][1]<<"    A43="<<mtp[3][2]
	   <<"     A44="<<mtp[3][3]<<std::endl;
#endif

#if 0
  std::cout<<"             Solution from Gauss-Jordan             "<<std::endl;
  std::cout<<"     A1="<<fitp[0]<<"    A2="<<fitp[1]<<"    A3="<<fitp[2]
	   <<"     A4="<<fitp[3]<<std::endl;
#endif


    for(int l=0; l<4;l++){
        for(int m=0; m<4; m++){
        Red[l][m]=mtp[l][m];
        }
    }

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
        Solution_vec[i] +=Red[i][j]*Org_vec[j];
        }
    }

#if 0
  std::cout<<"             Solution from explicite calculation       "<<std::endl;
  std::cout<<"    A1="<<Solution_vec[0]<<"    A2="<<Solution_vec[1]<<"    A3="<<Solution_vec[2]
	   <<"    A4="<<Solution_vec[3]<<std::endl;
#endif

    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
        for(int k=0; k<4; k++){
        Final[i][j] += Red[i][k]*Org[k][j];
        }
        }
    }
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
        if(Final[i][j]<1.0e-10) Final[i][j]=0.0;
        }
    }

#if 0
  std::cout<< ""<<std::endl;
  std::cout<<"            final matrix        "<<std::endl;
  std::cout<<"       A11="<<std::setw(10)<<Final[0][0]<<"       A12="<<std::setw(10)<<Final[0][1]
	   <<"       A13="<<std::setw(10)<<Final[0][2]<<"       A14="<<std::setw(10)<<Final[0][3]
	   <<std::endl;
  std::cout<<"       A21="<<std::setw(10)<<Final[1][0]<<"       A22="<<std::setw(10)<<Final[1][1]
	   <<"       A23="<<std::setw(10)<<Final[1][2]<<"       A24="<<std::setw(10)<<Final[1][3]
	   <<std::endl;
  std::cout<<"       A31="<<std::setw(10)<<Final[2][0]<<"       A32="<<std::setw(10)<<Final[2][1]
	   <<"       A33="<<std::setw(10)<<Final[2][2]<<"       A34="<<std::setw(10)<<Final[2][3]
	   <<std::endl;
  std::cout<<"       A41="<<std::setw(10)<<Final[3][0]<<"       A42="<<std::setw(10)<<Final[3][1]
	   <<"       A43="<<std::setw(10)<<Final[3][2]<<"       A44="<<std::setw(10)<<Final[3][3]
	   <<std::endl;

#endif

    double csr=0.0;
    for( std::size_t i=0; i<nn; ++i ){
        double w_hat=w[i], z_=z[i];
        double scal=CalcX(z_)*ct[i]+CalcY(z_)*st[i];
        csr += w_hat*(s[i]-scal)*(s[i]-scal);

#if 0
    if(1){
      //      std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
      //     std::cout<<"x coordinate = "<<(x0+u0*z_)<<std::endl;
      //     std::cout<<"y coordinate = "<<(y0+v0*z_)<<std::endl;
      std::cout<<std::setw(10)<<"layer = "<<i<<
	std::setw(10)<<"scal = "<<scal<<
	std::setw(10)<<"sdata = "<<s[i]<<std::endl;
      std::cout<<std::setw(10)<<"Res = "<<s[i]-scal<<std::endl;
      std::cout<<std::setw(10)<<"X = "<< CalcX(z_)<<" Y = "<< CalcY(z_)<<std::endl;
      std::cout<<std::setw(10)<<"chisqr = "<<chisqr<<std::endl;
    }
#endif
    } // for( std::size_t i=0; i<nn; ++i )
    csr /= nn-4.; // data points - number of parameters
  /*  
  if(chisqr<2){
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
    std::cout << "chisqr = " << chisqr << " nn-4 = " << nn-4 << std::endl;
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
  }
  */
    chisqr = csr;
    for( std::size_t i=0; i<nn; ++i ){
        DCLTrackHit *dclthit = dclthits[i];
        if(dclthit != nullptr){
            dclthit->SetCalPosition(CalcX(dclthit->GetGlobalZ()), CalcY(dclthit->GetGlobalZ()));
        }
    } // for( std::size_t i=0; i<nn; ++i )

    return status_=true;
} // bool nestdaq::DCLocalTrack::DoFit()

bool DCLocalTrack::AngleCorrection( void )
{
// True/False condition
// 1. Number of hits >= DCLocalMinNHits(⇦comes from DCConstants.h)
// 2. Gauss elimination is successful
    const std::string funcname = "[DCLocalTrack::AngleCorrection()] ";

    std::size_t nHits = dclthits.size();

    const double u0 = this->GetU0();
    const double v0 = this->GetV0();

    std::vector <double> z, res, w, ct, st, scaleFactor;
    // z: global z position of the hit layer
    // res: 1/sigma^2, where sigma is the resolution of the hit
    // w: position along measurement axis
    // ct: cos(tilt angle)
    // st: sin(tilt angle)
    // scaleFactor: for angle correction, used in the way that drift length = drift length * scaleFactor
    z.reserve(nHits); res.reserve(nHits); w.reserve(nHits);
    ct.reserve(nHits); st.reserve(nHits);

    // ================================
    // Prepare the data for fitting
    // ================================
    for( std::size_t i=0; i<nHits; ++i ){
        DCLTrackHit *hitp = dclthits[i];
        if( hitp ){
        double resolution = hitp->GetResolution();
        double zz = hitp->GetGlobalZ();
        double aa = hitp->GetWireAngle()*DCMath::Deg2Rad;
        double cc = cos(aa);
        double ss = sin(aa);
        double tan = u0*cc + v0*ss;
        double scale = sqrt(1.0 + tan*tan);

        z.push_back( zz ); res.push_back( 1./(resolution*resolution*scale*scale) ); 
        // Correct drift-length contribution using the current track angle.
        w.push_back( hitp->GetWirePosition() + static_cast<double>(hitp->GetLeftRight()) * hitp->GetDriftLength() * scale );
        ct.push_back( cc ); st.push_back( ss );
        scaleFactor.push_back( scale );

    #if 0
        std::cout << std::setw(10) << "layer = " << lnum 
            << std::setw(10) << "wire  = " << hitp->GetWire() << " "
        << std::setw(20) << "WirePosition = "<<hitp->GetWirePosition() << " "
        << std::setw(20) << "DriftLength = "<<hitp->GetDriftLength() << " "
            << std::endl;
    #endif
        }
    } // for( std::size_t i=0; i<nHits; ++i )

    std::size_t nn = z.size();


    ///std::cout << "nn = " << nn << std::endl;


    double matrx[16], *mtp[4], fitp[4];
    mtp[0]=&matrx[0]; mtp[1]=&matrx[4]; mtp[2]=&matrx[8]; mtp[3]=&matrx[12];

    for( int i=0; i<4; ++i ){
        fitp[i]=0.0;
        for( int j=0; j<4; ++j ){
        mtp[i][j]=0.0;
        }
    }

    for( std::size_t i=0; i<nn; ++i ){
        double res_ = res[i], zz=z[i], w_hat=w[i], ctt=ct[i], stt=st[i];
        mtp[0][0] += res_*ctt*ctt;
        mtp[0][1] += res_*zz*ctt*ctt;
        mtp[0][2] += res_*ctt*stt;
        mtp[0][3] += res_*zz*ctt*stt;
        mtp[1][1] += res_*zz*zz*ctt*ctt;
        mtp[1][2] += res_*zz*ctt*stt;
        mtp[1][3] += res_*zz*zz*ctt*stt;
        mtp[2][2] += res_*stt*stt;
        mtp[2][3] += res_*zz*stt*stt;
        mtp[3][3] += res_*zz*zz*stt*stt;

        fitp[0] += res_*w_hat*ctt;
        fitp[1] += res_*zz*w_hat*ctt;
        fitp[2] += res_*w_hat*stt;
        fitp[3] += res_*zz*w_hat*stt;
    }
    mtp[1][0]=mtp[0][1]; mtp[2][0]=mtp[0][2]; mtp[3][0]=mtp[0][3];
    mtp[2][1]=mtp[1][2]; mtp[3][1]=mtp[1][3]; mtp[3][2]=mtp[2][3];

    std::vector<int> indxc(nn), indxd(nn), ipiv(nn);

    #if 0
    std::cout<<"             Vector:  Q_i               "<<std::endl;
    std::cout<<"   A1=  "<<fitp[0]<<"     A2=  "<<fitp[1]<<"     A3=  "<<fitp[2]
        <<"     A4=  "<<fitp[3]<<std::endl;

    std::cout<<"             original matrix: M_ij        "<<std::endl;
    std::cout<<"    A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2]
            <<"    A14="<<mtp[0][3]<<std::endl;
    std::cout<<"    A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2]
        <<"    A24="<<mtp[1][3]<<std::endl;
    std::cout<<"    A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2]
        <<"    A34="<<mtp[2][3]<<std::endl;
    std::cout<<"    A41="<<mtp[3][0]<<"    A42="<<mtp[3][1]<<"    A43="<<mtp[3][2]
        <<"    A44="<<mtp[3][3]<<std::endl;
    #endif

    double Org[4][4]={0},Red[4][4]={0},Final[4][4]={0};
    double Org_vec[4]={0}, Solution_vec[4]={0};
    for(int l=0; l<4;l++){
        for(int m=0; m<4; m++){
        Org[l][m]=mtp[l][m];
        }
    }

    for (int i=0; i<4; i++){
        Org_vec[i]=fitp[i];
    }

    if( MathTools::GaussJordan(mtp,4,fitp,&indxc[0],
                    &indxd[0],&ipiv[0])==false ){
        std::cerr << funcname << ": Fitting fails" << std::endl;
        return status_=false;
    }
    x0=fitp[0]; y0=fitp[2]; u0=fitp[1]; v0=fitp[3];

    #if 0
    std::cout<<"             reduced matrix        "<<std::endl;
    std::cout<<"     A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2]
        <<"     A14="<<mtp[0][3]<<std::endl;
    std::cout<<"     A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2]
        <<"     A24="<<mtp[1][3]<<std::endl;
    std::cout<<"     A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2]
        <<"     A34="<<mtp[2][3]<<std::endl;
    std::cout<<"     A41="<<mtp[3][0]<<"    A42="<<mtp[3][1]<<"    A43="<<mtp[3][2]
        <<"     A44="<<mtp[3][3]<<std::endl;
    #endif

    #if 0
    std::cout<<"             Solution from Gauss-Jordan             "<<std::endl;
    std::cout<<"     A1="<<fitp[0]<<"    A2="<<fitp[1]<<"    A3="<<fitp[2]
        <<"     A4="<<fitp[3]<<std::endl;
    #endif


    for(int l=0; l<4;l++){
        for(int m=0; m<4; m++){
        Red[l][m]=mtp[l][m];
        }
    }

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
        Solution_vec[i] +=Red[i][j]*Org_vec[j];
        }
    }

    #if 0
    std::cout<<"             Solution from explicite calculation       "<<std::endl;
    std::cout<<"    A1="<<Solution_vec[0]<<"    A2="<<Solution_vec[1]<<"    A3="<<Solution_vec[2]
        <<"    A4="<<Solution_vec[3]<<std::endl;
    #endif

    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
        for(int k=0; k<4; k++){
        Final[i][j] += Red[i][k]*Org[k][j];
        }
        }
    }
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
        if(Final[i][j]<1.0e-10) Final[i][j]=0.0;
        }
    }

    #if 0
    std::cout<< ""<<std::endl;
    std::cout<<"            final matrix        "<<std::endl;
    std::cout<<"       A11="<<std::setw(10)<<Final[0][0]<<"       A12="<<std::setw(10)<<Final[0][1]
        <<"       A13="<<std::setw(10)<<Final[0][2]<<"       A14="<<std::setw(10)<<Final[0][3]
        <<std::endl;
    std::cout<<"       A21="<<std::setw(10)<<Final[1][0]<<"       A22="<<std::setw(10)<<Final[1][1]
        <<"       A23="<<std::setw(10)<<Final[1][2]<<"       A24="<<std::setw(10)<<Final[1][3]
        <<std::endl;
    std::cout<<"       A31="<<std::setw(10)<<Final[2][0]<<"       A32="<<std::setw(10)<<Final[2][1]
        <<"       A33="<<std::setw(10)<<Final[2][2]<<"       A34="<<std::setw(10)<<Final[2][3]
        <<std::endl;
    std::cout<<"       A41="<<std::setw(10)<<Final[3][0]<<"       A42="<<std::setw(10)<<Final[3][1]
        <<"       A43="<<std::setw(10)<<Final[3][2]<<"       A44="<<std::setw(10)<<Final[3][3]
        <<std::endl;

    #endif

    double chisqr=0.0;
    for( std::size_t i=0; i<nn; ++i ){
        double res_ = res[i], zz=z[i];
        double scal=CalcX(zz)*ct[i]+CalcY(zz)*st[i];
        chisqr += res_*(w[i]-scal)*(w[i]-scal);

    #if 0
        if(1){
        //      std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
        //     std::cout<<"x coordinate = "<<(x0+u0*zz)<<std::endl;
        //     std::cout<<"y coordinate = "<<(y0+v0*zz)<<std::endl;
        std::cout<<std::setw(10)<<"layer = "<<i<<
        std::setw(10)<<"scal = "<<scal<<
        std::setw(10)<<"sdata = "<<w[i]<<std::endl;
    std::cout<<std::setw(10)<<"Res = "<<w[i]-scal<<std::endl;
    std::cout<<std::setw(10)<<"X = "<< CalcX(zz)<<" Y = "<< CalcY(zz)<<std::endl;
        std::cout<<std::setw(10)<<"chisqr = "<<chisqr<<std::endl;
        }
    #endif
    }
    chisqr /= nn-4.;
    /*  
    if(chisqr<2){
        std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
        std::cout << "chisqr = " << chisqr << " nn-4 = " << nn-4 << std::endl;
        std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
    }
    */
    chisqr_=chisqr;
    for( std::size_t i=0; i<nn; ++i ){
        DCLTrackHit *hitp = dclthits[i];
        if( hitp ){
        double zz = hitp->GetGlobalZ();
        /*  
        if(chisqr<2){
        std::cout<<std::setw(10)<<"z = "<< zz <<std::endl;
        std::cout<<std::setw(10)<<"X = "<< CalcX(zz)<<" Y = "<< CalcY(zz)<<std::endl;
        }
        */
        hitp->SetCalPosition( CalcX(zz), CalcY(zz) );
        
        }
    }

    // std::cout << "***********************************************************" << std::endl;

    return status_=true;
} // bool DCLocalTrack::AngleCorrection( void )

void DCLocalTrack::CalcHitPositions(){
    for(auto hit : dclthits){
        double z = hit->GetGlobalZ();
        hit->SetCalPosition(CalcX(z), CalcY(z));
    }
} // void DCLocalTrack::CalcHitPositions()