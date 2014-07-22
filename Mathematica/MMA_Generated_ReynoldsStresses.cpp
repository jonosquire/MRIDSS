// Automatically generated temporary variables - class definition
dcmplxVecM rey_TiLap2TkxbTky, rey_TdzTiLap2Tkxb, rey_TiLap2Tky, rey_TdzTiLap2;




// Automatically generated temporary variables - class constructor
rey_TiLap2TkxbTky = dcmplxVecM(NZ_);
rey_TdzTiLap2Tkxb = dcmplxVecM(NZ_);
rey_TiLap2Tky = dcmplxVecM(NZ_);
rey_TdzTiLap2 = dcmplxVecM(NZ_);




//Assign automatically generated variables in equations 
// Vector variable definition Reynolds stress (automatic)
rey_TiLap2TkxbTky = ilap2tmp_.cast<dcmplx>().matrix()*(kxctmp_*kyctmp_);
rey_TdzTiLap2Tkxb = (ilap2tmp_.cast<dcmplx>()*kz_).matrix()*kxctmp_;
rey_TiLap2Tky = ilap2tmp_.cast<dcmplx>().matrix()*kyctmp_;
rey_TdzTiLap2 = (ilap2tmp_.cast<dcmplx>()*kz_).matrix();



////////////////////////////////////////////////////////
///  AUTOMATICALLY GENERATED REYNOLDS STRESSES      /////
///       see GenerateC++Equations.nb in MMA        /////

// Bx reynolds stress
reynolds_mat_tmp_ = C31_*rey_TdzTiLap2Tkxb.asDiagonal()-C32_*rey_TiLap2Tky.asDiagonal()-rey_TdzTiLap2Tkxb.asDiagonal()*C31_-rey_TiLap2Tky.asDiagonal()*C41_;


// By reynolds stress
reynolds_mat_tmp_ = rey_TdzTiLap2.asDiagonal()*C41_*rey_TdzTiLap2Tkxb.asDiagonal()+rey_TdzTiLap2Tkxb.asDiagonal()*C31_*rey_TiLap2TkxbTky.asDiagonal()+rey_TdzTiLap2Tkxb.asDiagonal()*C32_*rey_TdzTiLap2.asDiagonal()+rey_TiLap2TkxbTky.asDiagonal()*C32_*rey_TiLap2Tky.asDiagonal()+rey_TiLap2Tky.asDiagonal()*C41_*rey_TiLap2TkxbTky.asDiagonal()+rey_TiLap2Tky.asDiagonal()*C42_*rey_TdzTiLap2.asDiagonal()-rey_TdzTiLap2.asDiagonal()*C42_*rey_TiLap2Tky.asDiagonal()-rey_TiLap2TkxbTky.asDiagonal()*C31_*rey_TdzTiLap2Tkxb.asDiagonal();

/// AUTOMATICALLY GENERATED REYNOLDS STRESSES - END /////
////////////////////////////////////////////////////////
