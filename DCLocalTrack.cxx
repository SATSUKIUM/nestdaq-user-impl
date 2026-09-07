#include "DCLocalTrack.h"
#include "DCLTrackHit.h"
#include "DCMathConstants.h"
#include "DCMathTools.h"

using namedaq::DCLocalTrack;

DCLTrackHit* DCLocalTrack::GetHit(std::size_t iHit) const {
    if(iHit < hits.size()){
        return hits[iHit];
    }
    else{
        return nullptr;
    } 
} // DCLTrackHit* nestdaq::DCLocalTrack::GetHit(std::size_t iHit) const

bool DCLocalTrack::DoFit(){
// True/False condition
// 1. Number of hits >= DCLocalMinNHits
// 2. later....
    const std::string_view funcname = "[nestdaq::DCLocalTrack::DoFit] ";

    std::size_t nHits = dclthits.size();
    if(nHits < DCLocalMinNHits){
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
    size_t 
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
    x0_=fitp[0]; y0_=fitp[2]; u0_=fitp[1]; v0_=fitp[3];

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
      //     std::cout<<"x coordinate = "<<(x0_+u0_*z_)<<std::endl;
      //     std::cout<<"y coordinate = "<<(y0_+v0_*z_)<<std::endl;
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
            dclthit->SetCalPosition(CalcX(dclthit->GetGlobalZ()), CalcY(dclthit->GetGlobalZ()), u0_, v0_);
        }
    } // for( std::size_t i=0; i<nn; ++i )

    return status_=true;
} // bool nestdaq::DCLocalTrack::DoFit()