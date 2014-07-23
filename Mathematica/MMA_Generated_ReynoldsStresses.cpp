// Automatically generated temporary variables - class definition
dcmplxVecM rey_TdzTiLap2TkxbTky, rey_TiLap2TkxbTkyP2, rey_TiLap2TkxbP2Tky, rey_TiLap2TkxbTky, rey_TdzTiLap2Tkxb, rey_TdzTiLap2Tky, rey_TdzP2TiLap2, rey_TiLap2Tky, rey_TdzTiLap2, rey_Tdz;




// Automatically generated temporary variables - class constructor
rey_TdzTiLap2TkxbTky = dcmplxVecM(NZ_);
rey_TiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
rey_TiLap2TkxbP2Tky = dcmplxVecM(NZ_);
rey_TiLap2TkxbTky = dcmplxVecM(NZ_);
rey_TdzTiLap2Tkxb = dcmplxVecM(NZ_);
rey_TdzTiLap2Tky = dcmplxVecM(NZ_);
rey_TdzP2TiLap2 = dcmplxVecM(NZ_);
rey_TiLap2Tky = dcmplxVecM(NZ_);
rey_TdzTiLap2 = dcmplxVecM(NZ_);
rey_Tdz = dcmplxVecM(NZ_);




//Assign automatically generated variables in equations 
// Vector variable definition Reynolds stress (automatic)
rey_TdzTiLap2TkxbTky = (ilap2tmp_.cast<dcmplx>()*kz_).matrix()*(kxctmp_*kyctmp_);
rey_TiLap2TkxbTkyP2 = ilap2tmp_.cast<dcmplx>().matrix()*(kxctmp_*pow(kyctmp_,2));
rey_TiLap2TkxbP2Tky = ilap2tmp_.cast<dcmplx>().matrix()*(kyctmp_*pow(kxctmp_,2));
rey_TiLap2TkxbTky = ilap2tmp_.cast<dcmplx>().matrix()*(kxctmp_*kyctmp_);
rey_TdzTiLap2Tkxb = (ilap2tmp_.cast<dcmplx>()*kz_).matrix()*kxctmp_;
rey_TdzTiLap2Tky = (ilap2tmp_.cast<dcmplx>()*kz_).matrix()*kyctmp_;
rey_TdzP2TiLap2 = (ilap2tmp_.cast<dcmplx>()*kz_.pow(2)).matrix();
rey_TiLap2Tky = ilap2tmp_.cast<dcmplx>().matrix()*kyctmp_;
rey_TdzTiLap2 = (ilap2tmp_.cast<dcmplx>()*kz_).matrix();
rey_Tdz = kz_.matrix();



////////////////////////////////////////////////////////
///  AUTOMATICALLY GENERATED REYNOLDS STRESSES      /////
///       see GenerateC++Equations.nb in MMA        /////

// Bx reynolds stress
reynolds_mat_tmp_ = C31_*rey_TdzTiLap2Tkxb.asDiagonal() - C32_*rey_TiLap2Tky.asDiagonal() - rey_TdzTiLap2Tkxb.asDiagonal()*C31_ - rey_TiLap2Tky.asDiagonal()*C41_;


// By reynolds stress
reynolds_mat_tmp_ = rey_TdzTiLap2.asDiagonal()*C41_*rey_TdzTiLap2Tkxb.asDiagonal() - rey_TdzTiLap2.asDiagonal()*C42_*rey_TiLap2Tky.asDiagonal() + rey_TdzTiLap2Tkxb.asDiagonal()*C31_*rey_TiLap2TkxbTky.asDiagonal() + rey_TdzTiLap2Tkxb.asDiagonal()*C32_*rey_TdzTiLap2.asDiagonal() - rey_TiLap2TkxbTky.asDiagonal()*C31_*rey_TdzTiLap2Tkxb.asDiagonal() + rey_TiLap2TkxbTky.asDiagonal()*C32_*rey_TiLap2Tky.asDiagonal() + rey_TiLap2Tky.asDiagonal()*C41_*rey_TiLap2TkxbTky.asDiagonal() + rey_TiLap2Tky.asDiagonal()*C42_*rey_TdzTiLap2.asDiagonal();


// Ux reynolds stress
reynolds_mat_tmp_ = -kxctmp_*C11_ + kxctmp_*C33_ + kyctmp_*C11_*rey_TiLap2TkxbTky.asDiagonal() + kyctmp_*C12_*rey_TdzTiLap2.asDiagonal() - kyctmp_*C33_*rey_TiLap2TkxbTky.asDiagonal() - kyctmp_*C34_*rey_TdzTiLap2.asDiagonal() + rey_Tdz.asDiagonal()*C11_*rey_TdzTiLap2Tkxb.asDiagonal() - rey_Tdz.asDiagonal()*C12_*rey_TiLap2Tky.asDiagonal() - rey_Tdz.asDiagonal()*C33_*rey_TdzTiLap2Tkxb.asDiagonal() + rey_Tdz.asDiagonal()*C34_*rey_TiLap2Tky.asDiagonal();


// Uy reynolds stress
reynolds_mat_tmp_ = rey_TdzP2TiLap2.asDiagonal()*C21_*rey_TdzTiLap2Tkxb.asDiagonal() - rey_TdzP2TiLap2.asDiagonal()*C22_*rey_TiLap2Tky.asDiagonal() - rey_TdzP2TiLap2.asDiagonal()*C43_*rey_TdzTiLap2Tkxb.asDiagonal() + rey_TdzP2TiLap2.asDiagonal()*C44_*rey_TiLap2Tky.asDiagonal() - rey_TdzTiLap2Tkxb.asDiagonal()*C21_ + rey_TdzTiLap2Tkxb.asDiagonal()*C43_ - rey_TdzTiLap2TkxbTky.asDiagonal()*C11_*rey_TdzTiLap2Tkxb.asDiagonal() + rey_TdzTiLap2TkxbTky.asDiagonal()*C12_*rey_TiLap2Tky.asDiagonal() + rey_TdzTiLap2TkxbTky.asDiagonal()*C33_*rey_TdzTiLap2Tkxb.asDiagonal() - rey_TdzTiLap2TkxbTky.asDiagonal()*C34_*rey_TiLap2Tky.asDiagonal() + rey_TdzTiLap2Tky.asDiagonal()*C21_*rey_TiLap2TkxbTky.asDiagonal() + rey_TdzTiLap2Tky.asDiagonal()*C22_*rey_TdzTiLap2.asDiagonal() - rey_TdzTiLap2Tky.asDiagonal()*C43_*rey_TiLap2TkxbTky.asDiagonal() - rey_TdzTiLap2Tky.asDiagonal()*C44_*rey_TdzTiLap2.asDiagonal() + rey_TiLap2TkxbP2Tky.asDiagonal()*C11_ - rey_TiLap2TkxbP2Tky.asDiagonal()*C33_ - rey_TiLap2TkxbTkyP2.asDiagonal()*C11_*rey_TiLap2TkxbTky.asDiagonal() - rey_TiLap2TkxbTkyP2.asDiagonal()*C12_*rey_TdzTiLap2.asDiagonal() + rey_TiLap2TkxbTkyP2.asDiagonal()*C33_*rey_TiLap2TkxbTky.asDiagonal() + rey_TiLap2TkxbTkyP2.asDiagonal()*C34_*rey_TdzTiLap2.asDiagonal();

/// AUTOMATICALLY GENERATED REYNOLDS STRESSES - END /////
////////////////////////////////////////////////////////
